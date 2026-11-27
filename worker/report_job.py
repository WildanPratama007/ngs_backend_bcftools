"""`worker report`: FASTQ (S3) -> FastQC -> laporan HTML/ZIP + summary.json (S3)."""

import json
import logging
import os
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Optional

from . import config, proc, resources
from .errors import DISK_INSUFFICIENT, NO_OUTPUT, PIPELINE_FAILED, WorkerError
from .s3io import S3
from .spec import Layout, ReportSpec, parse_s3_uri
from .status import StatusReporter

log = logging.getLogger("worker.report")


def parse_fastqc_data(text: str) -> dict:
    """Ambil angka penting dari fastqc_data.txt."""
    modules: Dict[str, str] = {}
    basic: Dict[str, str] = {}
    qual_rows: List[tuple] = []
    adapter_max = 0.0
    current: Optional[str] = None
    for line in text.splitlines():
        if line.startswith(">>END_MODULE"):
            current = None
        elif line.startswith(">>"):
            name, _, status = line[2:].partition("\t")
            modules[name] = status.strip().upper()
            current = name
        elif line.startswith("#") or not line.strip():
            continue
        elif current == "Basic Statistics":
            k, _, v = line.partition("\t")
            basic[k] = v
        elif current == "Per sequence quality scores":
            q, _, c = line.partition("\t")
            try:
                qual_rows.append((float(q), float(c)))
            except ValueError:
                pass
        elif current == "Adapter Content":
            for cell in line.split("\t")[1:]:
                try:
                    adapter_max = max(adapter_max, float(cell))
                except ValueError:
                    pass

    total = int(basic.get("Total Sequences", "0") or 0)
    q30 = None
    denom = sum(c for _, c in qual_rows)
    if denom > 0:
        q30 = round(100.0 * sum(c for q, c in qual_rows if q >= 30) / denom, 2)
    return {
        "totalSequences": total,
        "sequenceLength": basic.get("Sequence length"),
        "pctGC": float(basic["%GC"]) if basic.get("%GC") else None,
        "pctSeqMeanQ30": q30,            # % read dengan rata-rata kualitas >= Q30
        "poorQualityFlagged": int(basic.get("Sequences flagged as poor quality", "0") or 0),
        "encoding": basic.get("Encoding"),
        "adapterMaxPct": round(adapter_max, 2),
        "modules": modules,              # nama modul FastQC -> PASS / WARN / FAIL
    }


def _read_zip(zip_path: Path) -> dict:
    with zipfile.ZipFile(zip_path) as z:
        member = next((n for n in z.namelist() if n.endswith("/fastqc_data.txt")), None)
        if not member:
            raise WorkerError(NO_OUTPUT, f"fastqc_data.txt tidak ada di {zip_path.name}")
        return parse_fastqc_data(z.read(member).decode("utf-8", errors="replace"))


def run(spec: ReportSpec, status: StatusReporter, s3: S3) -> None:
    scratch = Path(config.SCRATCH_DIR)
    job_dir = scratch / "jobs" / spec.jobId
    in_dir, out_dir = job_dir / "in", job_dir / "out"
    if job_dir.exists():
        shutil.rmtree(job_dir, ignore_errors=True)
    for d in (in_dir, out_dir):
        d.mkdir(parents=True, exist_ok=True)

    status.stage("validate_input", 2)
    targets = [("R1", spec.input.r1)] + ([("R2", spec.input.r2)] if spec.input.r2 else [])
    sizes = {rt: s3.size_of(uri) for rt, uri in targets}
    free = shutil.disk_usage(job_dir).free
    need = int(sum(sizes.values()) * 1.5)   # FASTQ + hasil FastQC (kecil)
    if not config.SKIP_DISK_CHECK and free < need:
        raise WorkerError(DISK_INSUFFICIENT,
                          f"Disk tidak cukup: bebas {free / 1e9:.1f} GB, dibutuhkan sekitar {need / 1e9:.1f} GB",
                          details={"freeBytes": free, "neededBytes": need})

    status.stage("download_input", 10)
    local: Dict[str, Path] = {}
    for rt, uri in targets:
        dest = in_dir / Path(parse_s3_uri(uri)[1]).name
        log.info("Mengunduh %s (%.2f GB)", uri, sizes[rt] / 1e9)
        s3.download(uri, dest)
        local[rt] = dest

    status.stage("fastqc", 30)
    plan = resources.plan(spec.resources.threads)
    threads = max(1, min(plan["threads"], len(local)))   # FastQC paralel per file
    env = os.environ.copy()
    if config.JAVA_HOME and os.path.isdir(config.JAVA_HOME):
        env["JAVA_HOME"] = config.JAVA_HOME
        env["PATH"] = f"{config.JAVA_HOME}/bin:{env.get('PATH', '')}"
    cmd = [config.FASTQC_PATH, "--threads", str(threads), "--outdir", str(out_dir)] + [str(p) for p in local.values()]
    log.info("Menjalankan: %s", " ".join(cmd))
    tail: List[str] = []

    def on_line(line: str):
        tail.append(line)
        del tail[:-20]

    rc = proc.run(cmd, env, str(job_dir), on_line)
    if rc != 0:
        raise WorkerError(PIPELINE_FAILED, "FastQC gagal:\n" + "\n".join(tail[-10:]),
                          details={"stage": "fastqc", "exitCode": rc})

    # summary.json dari fastqc_data.txt tiap file
    status.stage("finalize", 85)
    files = []
    for rt, path in local.items():
        base = path.name
        for suf in (".gz", ".fastq", ".fq"):   # sama dengan cara FastQC menamai hasilnya
            if base.endswith(suf):
                base = base[: -len(suf)]
        html, zp = out_dir / f"{base}_fastqc.html", out_dir / f"{base}_fastqc.zip"
        if not html.exists() or not zp.exists():
            raise WorkerError(NO_OUTPUT, f"Output FastQC untuk {path.name} tidak ditemukan")
        entry = {"readType": rt, "fileName": path.name, "sizeBytes": sizes[rt],
                 "report": {"html": f"{spec.report_dir}{html.name}", "zip": f"{spec.report_dir}{zp.name}"}}
        entry.update(_read_zip(zp))
        files.append(entry)

    summary = {
        "schemaVersion": 1,
        "jobId": spec.jobId,
        "layout": spec.input.layout.value,
        "tool": "FastQC",
        "files": files,
        "overall": {
            "modulesFail": sorted({m for f in files for m, s in f["modules"].items() if s == "FAIL"}),
            "modulesWarn": sorted({m for f in files for m, s in f["modules"].items() if s == "WARN"}),
            "pairedCountMatch": (files[0]["totalSequences"] == files[1]["totalSequences"]) if len(files) == 2 else None,
        },
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    status.stage("upload_output", 93)
    outputs = []
    uploads = [(summary_path, "summary.json")] + [(out_dir / f["report"][k].split("/")[-1], f["report"][k])
                                                  for f in files for k in ("html", "zip")]
    for path, rel in uploads:
        uri = spec.status_prefix + rel
        outputs.append({"key": rel, "uri": uri, "sizeBytes": s3.upload(path, uri)})
        log.info("Upload %s", uri)
    status.succeed(outputs, {"files": len(files), "overall": summary["overall"]})
