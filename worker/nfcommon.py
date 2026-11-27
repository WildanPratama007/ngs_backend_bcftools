"""Bagian yang dipakai bersama oleh semua worker berbasis Nextflow (vcf, trimming, alignment, variant_calling)."""

import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional

from . import config
from .errors import REFERENCE_NOT_FOUND, WorkerError
from .s3io import S3
from .spec import Reference
from .status import StatusReporter

log = logging.getLogger("worker.nf")

# Nama proses Nextflow -> nama stage di status.json. Nama stage mengikuti PIPELINE_CONTRACT:
# validate_input ... alignment ... finalize
PROCESS_STAGE = {
    "FASTQC":          "fastqc",
    "TRIMMOMATIC_SE":  "trimming",
    "TRIMMOMATIC_PE":  "trimming",
    "BWA_ALIGN_SE":    "alignment",
    "BWA_ALIGN_PE":    "alignment",
    "MARKDUPLICATES":  "alignment",
    "COVERAGE_STATS":  "alignment",
    "GET_CHROMS":      "variant_calling",
    "CALL_VARIANTS":   "variant_calling",
    "GATHER_VCF":      "variant_calling",
}

# Progress (%) saat sebuah proses mulai. Pipeline penuh (`worker vcf`) memakai skala 0-100 untuk seluruh job.
PROGRESS_FULL = {
    "FASTQC": 15, "TRIMMOMATIC_SE": 20, "TRIMMOMATIC_PE": 20, "BWA_ALIGN_SE": 30, "BWA_ALIGN_PE": 30,
    "MARKDUPLICATES": 55, "COVERAGE_STATS": 62, "GET_CHROMS": 65, "CALL_VARIANTS": 75, "GATHER_VCF": 88,
}
# Worker per langkah: progress 0-100 untuk langkah itu saja.
PROGRESS_STEP = {
    "trimming":        {"TRIMMOMATIC_SE": 30, "TRIMMOMATIC_PE": 30},
    "alignment":       {"BWA_ALIGN_SE": 20, "BWA_ALIGN_PE": 20, "MARKDUPLICATES": 60, "COVERAGE_STATS": 80},
    "variant_calling": {"GET_CHROMS": 15, "CALL_VARIANTS": 35, "GATHER_VCF": 80},
}

_SUBMITTED_RE = re.compile(r"Submitted process >\s*(?:\w+:)*([A-Z_]+)")
_ERR_PROC_RE = re.compile(r"Error executing process >\s*'(?:\w+:)*([A-Z_]+)")
_BLOCK_END = ("Work dir:", "Tip:", "Container:", "Command output:", "Command wrapper:", "Command executed:")


class Monitor:
    """Membaca output Nextflow: pindah stage di status.json dan tangkap pesan error dari task yang gagal."""

    def __init__(self, status: StatusReporter, progress: Optional[Dict[str, int]] = None):
        self.status = status
        self.progress = progress or PROGRESS_FULL
        self.failed_process: Optional[str] = None
        self._err: List[str] = []
        self._head: List[str] = []
        self._in_err = self._captured = self._in_head = False

    def feed(self, raw: str):
        s = raw.strip()
        if not s:
            return
        log.debug("[nf] %s", s)
        m = _ERR_PROC_RE.search(s)
        if m and not self.failed_process:
            self.failed_process = m.group(1)

        if s.startswith("ERROR ~") and not self._head:
            self._in_head = True
            self._head.append(s[len("ERROR ~"):].strip())
        elif self._in_head:
            if s.startswith("Caused by") or len(self._head) >= 6:
                self._in_head = False
            else:
                self._head.append(s)

        if not self._captured:
            if s.startswith("Command error:"):
                self._in_err, self._err = True, []
            elif self._in_err:
                if s.startswith(_BLOCK_END):
                    self._in_err, self._captured = False, True
                else:
                    self._err.append(s)

        m = _SUBMITTED_RE.search(s)
        if m and m.group(1) in PROCESS_STAGE and m.group(1) in self.progress:
            self.status.stage(PROCESS_STAGE[m.group(1)], self.progress[m.group(1)])

    def message(self, rc: int) -> str:
        if self._err:
            return "\n".join(self._err[-15:])
        if self._head:
            return "\n".join(self._head)
        return f"Nextflow berakhir dengan exit code {rc}"

    def failure_details(self, rc: int, fallback_stage: str) -> dict:
        return {"failedProcess": self.failed_process,
                "stage": PROCESS_STAGE.get(self.failed_process or "", fallback_stage),
                "exitCode": rc}


def ensure_reference(reference: Reference, s3: S3, scratch: Path) -> str:
    """Kembalikan path .fa referensi. Cari di REF_BASE; kalau tidak ada, unduh dari REF_S3_PREFIX."""
    fa = config.REF_FILES[reference.value]
    needed = [fa + suf for suf in config.REF_SUFFIXES]

    def complete(base: Path) -> bool:
        return all((base / n).exists() for n in needed)

    base = Path(config.REF_BASE)
    if complete(base):
        return str(base / fa)

    if config.REF_S3_PREFIX:
        local = scratch / "ref"
        if not complete(local):
            prefix = config.REF_S3_PREFIX.rstrip("/") + "/"
            log.info("Mengunduh referensi %s dari %s", fa, prefix)
            for name in needed:
                s3.download(prefix + name, local / name)
        if complete(local):
            return str(local / fa)

    missing = [n for n in needed if not (base / n).exists()]
    raise WorkerError(REFERENCE_NOT_FOUND,
                      f"Referensi {reference.value} tidak lengkap di {config.REF_BASE} (kurang: {', '.join(missing)}) "
                      f"dan REF_S3_PREFIX {'tidak bisa dipakai' if config.REF_S3_PREFIX else 'tidak diisi'}",
                      details={"reference": reference.value, "missing": missing})


def nextflow_env(plan: dict) -> Dict[str, str]:
    env = os.environ.copy()
    if config.JAVA_HOME and os.path.isdir(config.JAVA_HOME):
        env["JAVA_HOME"] = config.JAVA_HOME
        env["PATH"] = f"{config.JAVA_HOME}/bin:{env.get('PATH', '')}"
    env.pop("JAVA_CMD", None)
    env["NXF_VER"] = config.NEXTFLOW_VERSION
    env["NXF_ANSI_LOG"] = "false"
    env.setdefault("NXF_OPTS", "-Xms128m -Xmx1g")
    env["NGS_MAX_CPUS"] = str(plan["cpus"])
    env["NGS_MAX_MEMORY"] = f"{plan['memory_mb']} MB"
    return env
