"""
NGS API — Pipeline Runner (BCFTools)

Dua jalur eksekusi yang terpisah:
- run_pipeline(): pipeline penuh (pipeline/main.nf), dipakai POST /api/pipeline/submit
- run_step():     satu step wizard (pipeline/main_wizard.nf -entry ...), dipakai /api/analyses/{id}/*/run

Keduanya memakai _spawn_nextflow() untuk menjalankan proses dan membaca output.
"""

import asyncio
import glob
import json
import logging
import os
import re
from typing import Optional

from .models import StageStatus, AnalysisMode, EndType, ReferenceGenome
from .config import (
    NEXTFLOW_BIN, NEXTFLOW_SCRIPT, NEXTFLOW_SCRIPT_WIZARD, NEXTFLOW_PROFILE,
    NEXTFLOW_JAVA_HOME, NEXTFLOW_VERSION,
    REF_GENOMES, RESULTS_DIR, UPLOAD_DIR, DEFAULT_THREADS, BIN_DIR,
    SAMTOOLS_BIN, BCFTOOLS_BIN, FASTQC_PATH, resolve_adapter_file, MAX_CONCURRENT_JOBS,
)
from . import job_manager

logger = logging.getLogger("ngs.pipeline")

# Semaphore dari config — bukan hardcode
_semaphore = asyncio.Semaphore(MAX_CONCURRENT_JOBS)
_queue: list[str] = []


def sanitize_sample_name(name: str) -> str:
    """Sanitasi nama sampel: hanya A-Z a-z 0-9 . - _ (sesuai Sample_Name_Sanitizer spec)."""
    sanitized = re.sub(r'[^A-Za-z0-9._\-]', '_', name)
    sanitized = re.sub(r'_+', '_', sanitized).strip('_')
    return sanitized or "sample"


def get_queue_status() -> dict:
    return {
        "max_concurrent":  MAX_CONCURRENT_JOBS,
        "active_jobs":     MAX_CONCURRENT_JOBS - _semaphore._value,
        "available_slots": _semaphore._value,
        "queued_jobs":     list(_queue),
        "queued_count":    len(_queue),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Pemetaan proses Nextflow → stage job
# ═══════════════════════════════════════════════════════════════════════════════

STAGE_OF_PROCESS = {
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

# Progress (%) saat sebuah proses mulai (Submitted). Hanya naik, tidak turun.
PROGRESS_OF_PROCESS = {
    "FASTQC":          50,
    "TRIMMOMATIC_SE":  50,
    "TRIMMOMATIC_PE":  50,
    "BWA_ALIGN_SE":    30,
    "BWA_ALIGN_PE":    30,
    "MARKDUPLICATES":  60,
    "COVERAGE_STATS":  85,
    "GET_CHROMS":      10,
    "CALL_VARIANTS":   50,
    "GATHER_VCF":      90,
}

# Stage yang harus sudah selesai begitu stage ini mulai (Nextflow tidak mencetak "completed")
PREREQ_STAGE = {"alignment": "trimming", "variant_calling": "alignment"}
ALL_STAGES   = ("fastqc", "trimming", "alignment", "variant_calling")

_SUBMITTED_RE = re.compile(r"Submitted process >\s*(?:\w+:)*([A-Z_]+)")
_ERR_PROC_RE  = re.compile(r"Error executing process >\s*'(?:\w+:)*([A-Z_]+)")


class _RunMonitor:
    """Membaca output Nextflow: update progress stage dan kumpulkan pesan error."""

    _BLOCK_END = ("Work dir:", "Tip:", "Container:", "Command output:", "Command wrapper:", "Command executed:")

    def __init__(self, job_id: str):
        self.job_id          = job_id
        self.failed_process: Optional[str] = None
        self._in_err         = False
        self._captured       = False
        self._err: list[str] = []
        self._in_head        = False
        self._head: list[str] = []

    def feed(self, raw: str):
        s = raw.strip()
        logger.debug("[NF] %s", s)

        m = _ERR_PROC_RE.search(s)
        if m and not self.failed_process:
            self.failed_process = m.group(1)

        # Error di level Nextflow (script/konfigurasi), bukan di dalam task
        if s.startswith("ERROR ~") and not self._head:
            self._in_head = True
            self._head.append(s[len("ERROR ~"):].strip())
        elif self._in_head:
            if not s or s.startswith("Caused by") or len(self._head) >= 6:
                self._in_head = False
            else:
                self._head.append(s)

        # Blok "Command error:" dari task yang gagal (stderr script, termasuk pesan "ERROR: ...")
        if not self._captured:
            if s.startswith("Command error:"):
                self._in_err = True
                self._err = []
            elif self._in_err:
                if s.startswith(self._BLOCK_END):
                    self._in_err, self._captured = False, True
                elif s:
                    self._err.append(s)

        self._update_progress(s)

    def message(self, returncode: Optional[int] = None) -> str:
        if self._err:
            return "\n".join(self._err[-15:])
        if self._head:
            return "\n".join(self._head)
        return f"Pipeline exited with error (exit code {returncode})"

    def _update_progress(self, s: str):
        m = _SUBMITTED_RE.search(s)
        if m and m.group(1) in STAGE_OF_PROCESS:
            name  = m.group(1)
            stage = STAGE_OF_PROCESS[name]
            job   = job_manager.get_job(self.job_id)
            if not job:
                return

            prereq = PREREQ_STAGE.get(stage)
            if prereq and job["stages"][prereq]["status"] == "running":
                job_manager.update_stage(self.job_id, prereq, StageStatus.COMPLETED, 100)

            cur = job["stages"][stage]
            pct = PROGRESS_OF_PROCESS[name]
            if cur["status"] not in ("completed", "failed", "cancelled") and (cur["status"] != "running" or pct > cur["progress"]):
                job_manager.update_stage(self.job_id, stage, StageStatus.RUNNING, max(pct, cur["progress"]))

        # Ringkasan coverage sementara (JSON lengkap dibaca setelah selesai)
        if "📊 COVERAGE:" in s:
            try:
                mean_cov = float(s.split("📊 COVERAGE:")[1].replace("x", "").split("|")[0].strip())
                job_manager.update_coverage(self.job_id, {"mean_depth": mean_cov})
            except (ValueError, IndexError):
                pass


def _nextflow_env() -> dict:
    """Environment subprocess Nextflow: Java 17 (bukan Java conda) dan versi Nextflow tetap."""
    env = os.environ.copy()
    if NEXTFLOW_JAVA_HOME and os.path.isdir(NEXTFLOW_JAVA_HOME):
        env["JAVA_HOME"] = NEXTFLOW_JAVA_HOME
        env["PATH"]      = f"{NEXTFLOW_JAVA_HOME}/bin:{env.get('PATH', '')}"
    env.pop("JAVA_CMD", None)

    # Profil conda_env memanggil `conda activate`; pastikan `conda` ada di PATH walau
    # uvicorn dijalankan dari shell yang environment-nya belum diaktifkan.
    extra = []
    for d in (BIN_DIR, os.path.normpath(os.path.join(BIN_DIR, "..", "..", "..", "bin"))):
        if os.path.isfile(os.path.join(d, "conda")) and d not in env.get("PATH", "").split(":"):
            extra.append(d)
    if extra:
        env["PATH"] = f"{env.get('PATH', '')}:{':'.join(extra)}"

    if NEXTFLOW_VERSION:
        env["NXF_VER"] = NEXTFLOW_VERSION
    # Output non-ANSI: satu baris "Submitted process > NAMA" per task. Parser progress bergantung padanya.
    env["NXF_ANSI_LOG"] = "false"
    return env


async def _spawn_nextflow(cmd: list, job_id: str, step: str, cwd: str):
    """Jalankan Nextflow, baca stdout+stderr, kembalikan (returncode, monitor)."""
    process = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=cwd,
        env=_nextflow_env(),
    )
    job_manager.register_nextflow_pid(job_id, process.pid, step)
    logger.info("Nextflow dimulai job=%s step=%s pid=%s", job_id, step, process.pid)

    monitor = _RunMonitor(job_id)

    async def read_stream(stream):
        async for line in stream:
            monitor.feed(line.decode("utf-8", errors="replace"))
            job_manager.update_elapsed(job_id)

    await asyncio.gather(read_stream(process.stdout), read_stream(process.stderr))
    await process.wait()
    job_manager.unregister_nextflow_pid(job_id, step)
    logger.info("Nextflow selesai job=%s step=%s rc=%s", job_id, step, process.returncode)
    return process.returncode, monitor


def _load_coverage_json(job_id: str, outdir: str):
    """Baca coverage JSON hasil step Alignment dan simpan ke job."""
    job    = job_manager.get_job(job_id)
    sample = job.get("sample_name", "") if job else ""
    files  = glob.glob(os.path.join(outdir, sample, "bwa", "*.coverage.json")) or \
             glob.glob(os.path.join(outdir, "*", "bwa", "*.coverage.json"))
    if not files:
        return
    try:
        with open(files[0]) as f:
            job_manager.update_coverage(job_id, json.load(f))
    except Exception:
        logger.exception("Gagal membaca coverage JSON %s", files[0])


def _find_upload_bed(job_id: str) -> Optional[str]:
    beds = sorted(glob.glob(os.path.join(UPLOAD_DIR, job_id, "*.bed")))
    return beds[0] if beds else None


def find_trimmed_reads(job_id: str) -> Optional[dict]:
    """
    Cari hasil Trimming untuk job. Kembalikan {"se": path} atau {"r1": path, "r2": path},
    atau None kalau belum ada.
    """
    job = job_manager.get_job(job_id)
    if not job:
        return None
    sample   = job.get("sample_name", "")
    is_pe    = job.get("end_type") == "PE"
    dirs     = [os.path.join(RESULTS_DIR, job_id, sample, "trimmed")]
    dirs    += sorted(glob.glob(os.path.join(RESULTS_DIR, job_id, "*", "trimmed")))
    for d in dirs:
        if not os.path.isdir(d):
            continue
        if is_pe:
            r1 = sorted(glob.glob(os.path.join(d, "*_R1.trimmed.fastq")))
            r2 = sorted(glob.glob(os.path.join(d, "*_R2.trimmed.fastq")))
            if r1 and r2:
                return {"r1": r1[0], "r2": r2[0]}
        else:
            se = [f for f in sorted(glob.glob(os.path.join(d, "*.trimmed.fastq")))
                  if not f.endswith(("_R1.trimmed.fastq", "_R2.trimmed.fastq"))]
            if se:
                return {"se": se[0]}
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline penuh (POST /api/pipeline/submit)
# ═══════════════════════════════════════════════════════════════════════════════

async def run_pipeline(
    job_id:        str,
    fastq_path:    str = None,
    fastq_r1:      str = None,
    fastq_r2:      str = None,
    bed_file:      str = None,
    analysis_mode: AnalysisMode = AnalysisMode.WGS,
    end_type:      EndType      = EndType.SE,
    reference:     ReferenceGenome = ReferenceGenome.GRCh38,
    sample_sex:    str = None,
    adapter_kit:   str = None,
    # Trimming
    trim_quality:  int  = 20,
    trim_min_len:  int  = 36,
    trim_adapter:  bool = True,
    trim_poly_g:   bool = False,
    # Alignment
    mark_duplicates: bool = True,
    # Variant calling
    min_qual:   int = 20,
    min_dp:     int = 10,
    min_gq:     int = 20,
    max_depth:  int = 500,
    min_bq:     int = 20,
    min_mq:     int = 20,
):
    _queue.append(job_id)
    job_manager.begin_step(job_id)
    job_manager.update_job_status(job_id, StageStatus.PENDING)

    try:
        async with _semaphore:
            if job_id in _queue:
                _queue.remove(job_id)

            if job_manager.is_cancel_requested(job_id):   # dibatalkan saat masih antre
                _mark_unfinished_cancelled(job_id)
                job_manager.update_job_status(job_id, StageStatus.CANCELLED)
                return

            await _execute_pipeline(
                job_id=job_id, fastq_path=fastq_path, fastq_r1=fastq_r1, fastq_r2=fastq_r2,
                bed_file=bed_file, analysis_mode=analysis_mode, end_type=end_type,
                reference=reference, sample_sex=sample_sex, adapter_kit=adapter_kit,
                trim_quality=trim_quality, trim_min_len=trim_min_len,
                trim_adapter=trim_adapter, trim_poly_g=trim_poly_g,
                mark_duplicates=mark_duplicates,
                min_qual=min_qual, min_dp=min_dp, min_gq=min_gq,
                max_depth=max_depth, min_bq=min_bq, min_mq=min_mq,
            )
    finally:
        if job_id in _queue:
            _queue.remove(job_id)
        job_manager.end_step(job_id)


def _mark_unfinished_cancelled(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        return
    for stage, st in job["stages"].items():
        if st["status"] in ("pending", "running"):
            job_manager.update_stage(job_id, stage, StageStatus.CANCELLED, st.get("progress", 0))


async def _execute_pipeline(
    job_id, fastq_path, fastq_r1, fastq_r2, bed_file, analysis_mode, end_type,
    reference, sample_sex, adapter_kit, trim_quality, trim_min_len, trim_adapter,
    trim_poly_g, mark_duplicates, min_qual, min_dp, min_gq, max_depth, min_bq, min_mq,
):
    job_manager.update_job_status(job_id, StageStatus.RUNNING)
    job        = job_manager.get_job(job_id)
    ref_genome = REF_GENOMES.get(reference.value, REF_GENOMES["GRCh38"])
    outdir     = os.path.join(RESULTS_DIR, job_id)
    ploidy     = "GRCh37" if reference == ReferenceGenome.GRCh37 else "GRCh38"
    work_dir   = os.path.join(outdir, ".work")
    os.makedirs(outdir, exist_ok=True)

    cmd = [
        NEXTFLOW_BIN, "-log", os.path.join(outdir, ".nextflow.log"),
        "run", NEXTFLOW_SCRIPT,
        "-profile",       NEXTFLOW_PROFILE,
        "-w",             work_dir,
        "--ref_genome",   ref_genome,
        "--outdir",       outdir,
        "--sample_name",  job["sample_name"],
        "--ploidy",       ploidy,
        "--threads",      str(DEFAULT_THREADS),
        "--bin_dir",      BIN_DIR,
        "--samtools_bin", SAMTOOLS_BIN,
        "--bcftools_bin", BCFTOOLS_BIN,
        "--fastqc_path",  FASTQC_PATH,
        "--adapter_file", resolve_adapter_file(adapter_kit, end_type.value),
        "--trim_sliding_window_quality", str(trim_quality),
        "--trim_min_len",                str(trim_min_len),
        "--trim_adapter",                str(trim_adapter).lower(),
        "--trim_poly_g",                 str(trim_poly_g).lower(),
        "--mark_duplicates",             str(mark_duplicates).lower(),
        "--min_qual",    str(min_qual),
        "--min_dp",      str(min_dp),
        "--min_gq",      str(min_gq),
        "--max_depth",   str(max_depth),
        "--min_bq",      str(min_bq),
        "--min_mq",      str(min_mq),
        "-resume",
    ]
    if sample_sex:
        cmd += ["--sample_sex", sample_sex]

    if end_type == EndType.PE and fastq_r1 and fastq_r2:
        cmd += ["--reads_r1", fastq_r1, "--reads_r2", fastq_r2]
    else:
        cmd += ["--reads", fastq_path]

    if analysis_mode in (AnalysisMode.WES, AnalysisMode.TARGETED) and bed_file:
        cmd += ["--target_bed", bed_file]

    try:
        rc, monitor = await _spawn_nextflow(cmd, job_id, "full", os.path.dirname(NEXTFLOW_SCRIPT))

        if rc == 0:
            for stage in ALL_STAGES:
                job_manager.update_stage(job_id, stage, StageStatus.COMPLETED, 100)
            _load_coverage_json(job_id, outdir)
            job_manager.update_job_status(job_id, StageStatus.COMPLETED)
        elif job_manager.is_cancel_requested(job_id):
            _mark_unfinished_cancelled(job_id)
            job_manager.update_job_status(job_id, StageStatus.CANCELLED)
        else:
            msg   = monitor.message(rc)
            stage = STAGE_OF_PROCESS.get(monitor.failed_process or "")
            if not stage:
                current = job_manager.get_job(job_id)
                running = [s for s in ALL_STAGES if current and current["stages"][s]["status"] == "running"]
                stage   = running[-1] if running else "variant_calling"
            job_manager.update_stage(job_id, stage, StageStatus.FAILED, error=msg)
            job_manager.set_nextflow_log(job_id, msg)
            job_manager.update_job_status(job_id, StageStatus.FAILED)

    except Exception as e:
        logger.exception("Pipeline penuh gagal job=%s", job_id)
        job_manager.update_job_status(job_id, StageStatus.FAILED)
        job_manager.update_stage(job_id, "variant_calling", StageStatus.FAILED, error=str(e))
        job_manager.set_nextflow_log(job_id, f"Exception: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# Wizard: satu step per eksekusi
# ═══════════════════════════════════════════════════════════════════════════════

_STEP_CONFIG = {
    "qc":              {"stage": "fastqc",          "entry": "QC"},
    "trimming":        {"stage": "trimming",        "entry": "TRIMMING"},
    "alignment":       {"stage": "alignment",       "entry": "ALIGNMENT"},
    "variant_calling": {"stage": "variant_calling", "entry": "VARIANT_CALLING"},
}


async def run_step(job_id: str, step: str, **kwargs):
    """
    Antre (batas MAX_CONCURRENT_JOBS, sama dengan pipeline penuh), lalu jalankan satu step wizard.
    Selama antre, stage berstatus "pending". Parameter step lihat _run_step_locked().
    """
    if step not in _STEP_CONFIG:
        return
    stage = _STEP_CONFIG[step]["stage"]
    tag   = f"{job_id}:{step}"

    job_manager.begin_step(job_id)
    _queue.append(tag)
    try:
        job_manager.reset_stage(job_id, stage)
        job_manager.update_stage(job_id, stage, StageStatus.PENDING, 0)
        job_manager.refresh_job_status(job_id)

        async with _semaphore:
            if tag in _queue:
                _queue.remove(tag)
            if job_manager.is_cancel_requested(job_id):   # dibatalkan saat masih antre
                job_manager.update_stage(job_id, stage, StageStatus.CANCELLED)
                return
            await _run_step_locked(job_id, step, **kwargs)
    finally:
        if tag in _queue:
            _queue.remove(tag)
        job_manager.end_step(job_id)
        job_manager.refresh_job_status(job_id)


async def _run_step_locked(
    job_id:          str,
    step:            str,         # "qc", "trimming", "alignment", "variant_calling"
    fastq_path:      str = None,
    fastq_r1:        str = None,
    fastq_r2:        str = None,
    bam_file:        str = None,
    end_type:        EndType      = EndType.SE,
    reference:       ReferenceGenome = ReferenceGenome.GRCh38,
    sample_sex:      str  = None,
    trim_quality:    int  = 20,
    trim_min_len:    int  = 36,
    trim_adapter:    bool = True,
    trim_poly_g:     bool = False,
    adapter_kit:     str  = None,   # truseq_pe, truseq_se, nextera, ... (None = ikut job / SE-PE)
    mark_duplicates: bool = True,
    min_qual:        int  = 20,
    min_dp:          int  = 10,
    min_gq:          int  = 20,
    max_depth:       int  = 500,
    min_bq:          int  = 20,
    min_mq:          int  = 20,
):
    """Jalankan satu step pipeline via Nextflow entry point (pipeline/main_wizard.nf)."""
    if step not in _STEP_CONFIG:
        return

    stage = _STEP_CONFIG[step]["stage"]
    entry = _STEP_CONFIG[step]["entry"]

    job_manager.update_stage(job_id, stage, StageStatus.RUNNING, 0)
    job_manager.refresh_job_status(job_id)

    try:
        job        = job_manager.get_job(job_id) or {}
        if step == "trimming":
            _clear_trimmed_outputs(job_id)
        ref_val    = job.get("reference", "GRCh38")
        ref_genome = REF_GENOMES.get(ref_val, REF_GENOMES["GRCh38"])
        end_type_v = job.get("end_type", "SE")
        sample     = job.get("sample_name", "")
        outdir     = os.path.join(RESULTS_DIR, job_id)
        work_dir   = os.path.join(outdir, f".work_{step}")
        ploidy     = "GRCh37" if "GRCh37" in ref_val else "GRCh38"
        kit        = adapter_kit or job.get("adapter_kit")
        os.makedirs(outdir, exist_ok=True)

        cmd = [
            NEXTFLOW_BIN, "-log", os.path.join(outdir, f".nextflow_{step}.log"),
            "run", NEXTFLOW_SCRIPT_WIZARD,
            "-profile",    NEXTFLOW_PROFILE,
            "-entry",      entry,
            "-w",          work_dir,
            "--ref_genome", ref_genome,
            "--outdir",    outdir,
            "--step_tag",  step,
            "--sample_name", sample,
            "--ploidy",    ploidy,
            "--threads",   str(DEFAULT_THREADS),
            "--bin_dir",   BIN_DIR,
            "--samtools_bin", SAMTOOLS_BIN,
            "--bcftools_bin", BCFTOOLS_BIN,
            "--fastqc_path",  FASTQC_PATH,
            "--trim_sliding_window_quality", str(trim_quality),
            "--trim_min_len",    str(trim_min_len),
            "--trim_adapter",    str(trim_adapter).lower(),
            "--trim_poly_g",     str(trim_poly_g).lower(),
            "--adapter_file",    resolve_adapter_file(kit, end_type_v),
            "--mark_duplicates", str(mark_duplicates).lower(),
            "--min_qual",        str(min_qual),
            "--min_dp",          str(min_dp),
            "--min_gq",          str(min_gq),
            "--max_depth",       str(max_depth),
            "--min_bq",          str(min_bq),
            "--min_mq",          str(min_mq),
        ]

        # BED untuk WES/Targeted: dipakai coverage (alignment) dan variant calling
        if step in ("alignment", "variant_calling") and job.get("analysis_mode") in ("WES", "Targeted"):
            bed = job.get("bed_path") or _find_upload_bed(job_id)
            if not bed:
                raise RuntimeError(f"Mode {job.get('analysis_mode')} membutuhkan BED file, tapi tidak ditemukan di upload.")
            cmd += ["--target_bed", bed]

        if step in ("qc", "trimming"):
            # Path upload disimpan di job saat Step 1; scan folder hanya fallback
            if end_type_v == "PE":
                r1 = fastq_r1 or job.get("fastq_r1")
                r2 = fastq_r2 or job.get("fastq_r2")
                if not (r1 and r2):
                    r1, r2 = _find_uploaded_pair(job_id)
                cmd += ["--reads_r1", r1, "--reads_r2", r2]
            else:
                cmd += ["--reads", fastq_path or job.get("fastq_path") or _find_uploaded_se(job_id)]

        elif step == "alignment":
            trimmed = find_trimmed_reads(job_id)
            if not trimmed:
                raise RuntimeError("Hasil Trimming tidak ditemukan. Jalankan Step 2 dulu.")
            if "se" in trimmed:
                cmd += ["--trimmed_reads", trimmed["se"]]
            else:
                cmd += ["--trimmed_r1", trimmed["r1"], "--trimmed_r2", trimmed["r2"]]

        elif step == "variant_calling":
            if not bam_file:
                raise RuntimeError("BAM file wajib untuk variant calling.")
            cmd += ["--bam_file", bam_file]
            if sample_sex:
                cmd += ["--sample_sex", sample_sex]

        rc, monitor = await _spawn_nextflow(cmd, job_id, step, os.path.dirname(NEXTFLOW_SCRIPT_WIZARD))

        if rc == 0:
            if step == "alignment":
                _load_coverage_json(job_id, outdir)
            job_manager.update_stage(job_id, stage, StageStatus.COMPLETED, 100)
        elif job_manager.is_cancel_requested(job_id):
            job_manager.update_stage(job_id, stage, StageStatus.CANCELLED)
        else:
            job_manager.update_stage(job_id, stage, StageStatus.FAILED, error=monitor.message(rc))

    except Exception as e:
        logger.exception("Step %s gagal job=%s", step, job_id)
        job_manager.update_stage(job_id, stage, StageStatus.FAILED, error=str(e))
    finally:
        job_manager.refresh_job_status(job_id)


def _clear_trimmed_outputs(job_id: str):
    """Hapus hasil Trimming lama supaya run yang gagal tidak meninggalkan file basi untuk Alignment."""
    for d in glob.glob(os.path.join(RESULTS_DIR, job_id, "*", "trimmed")):
        for f in glob.glob(os.path.join(d, "*")):
            try:
                os.remove(f)
            except OSError:
                pass


def _uploaded_fastqs(job_id: str) -> list:
    d = os.path.join(UPLOAD_DIR, job_id)
    if not os.path.isdir(d):
        return []
    return sorted(os.path.join(d, f) for f in os.listdir(d)
                  if f.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')))


def _find_uploaded_se(job_id: str) -> str:
    files = _uploaded_fastqs(job_id)
    if not files:
        raise RuntimeError("FASTQ upload tidak ditemukan. Jalankan Step 1 dulu.")
    return files[0]


def _find_uploaded_pair(job_id: str):
    files = _uploaded_fastqs(job_id)
    r1 = [f for f in files if re.search(r'[_.]R?1([_.]|$)', os.path.basename(f))]
    r2 = [f for f in files if re.search(r'[_.]R?2([_.]|$)', os.path.basename(f))]
    if not (r1 and r2):
        raise RuntimeError("Pasangan FASTQ R1/R2 tidak ditemukan di upload. Jalankan Step 1 dulu.")
    return r1[0], r2[0]
