"""
NGS API — Pipeline Runner (BCFTools)
Menjalankan Nextflow BCFTools pipeline secara async dengan:
- Level 1: Semaphore 5 job paralel (multi-sample)
- Level 2: Scatter-gather BCFTools per kromosom (via main.nf)
"""

import asyncio
import os
from .models import StageStatus, AnalysisMode, EndType, ReferenceGenome
from .config import (
    NEXTFLOW_BIN, NEXTFLOW_SCRIPT, NEXTFLOW_PROFILE,
    REF_GENOMES, RESULTS_DIR, DEFAULT_THREADS, BIN_DIR,
    SAMTOOLS_BIN, BCFTOOLS_BIN, FASTQC_PATH
)
from . import job_manager

# ── Level 1: Semaphore — maksimal 5 job berjalan bersamaan ───────────────────
MAX_CONCURRENT_JOBS = 5
_semaphore = asyncio.Semaphore(MAX_CONCURRENT_JOBS)

# ── Queue tracker ─────────────────────────────────────────────────────────────
_queue: list[str] = []   # job_id yang sedang menunggu slot


def get_queue_status() -> dict:
    """Return status antrian saat ini."""
    return {
        "max_concurrent":  MAX_CONCURRENT_JOBS,
        "active_jobs":     MAX_CONCURRENT_JOBS - _semaphore._value,
        "available_slots": _semaphore._value,
        "queued_jobs":     list(_queue),
        "queued_count":    len(_queue),
    }


async def run_pipeline(
    job_id:        str,
    fastq_path:    str = None,
    fastq_r1:      str = None,
    fastq_r2:      str = None,
    bed_file:      str = None,
    analysis_mode: AnalysisMode = AnalysisMode.WGS,
    end_type:      EndType      = EndType.SE,
    reference:     ReferenceGenome = ReferenceGenome.GRCh38,
):
    """
    Jalankan pipeline dengan semaphore.
    Jika semua slot penuh, job akan mengantri otomatis.
    """
    # Masuk antrian jika slot penuh
    _queue.append(job_id)
    job_manager.update_job_status(job_id, StageStatus.PENDING)

    async with _semaphore:
        # Keluar dari antrian — dapat slot
        if job_id in _queue:
            _queue.remove(job_id)

        await _execute_pipeline(
            job_id        = job_id,
            fastq_path    = fastq_path,
            fastq_r1      = fastq_r1,
            fastq_r2      = fastq_r2,
            bed_file      = bed_file,
            analysis_mode = analysis_mode,
            end_type      = end_type,
            reference     = reference,
        )


async def _execute_pipeline(
    job_id:        str,
    fastq_path:    str = None,
    fastq_r1:      str = None,
    fastq_r2:      str = None,
    bed_file:      str = None,
    analysis_mode: AnalysisMode = AnalysisMode.WGS,
    end_type:      EndType      = EndType.SE,
    reference:     ReferenceGenome = ReferenceGenome.GRCh38,
):
    """Eksekusi Nextflow pipeline untuk 1 job."""

    job_manager.update_job_status(job_id, StageStatus.RUNNING)
    ref_genome = REF_GENOMES.get(reference.value, REF_GENOMES["GRCh38"])
    outdir     = os.path.join(RESULTS_DIR, job_id)

    # Deteksi ploidy dari referensi
    ploidy = "GRCh37" if reference == ReferenceGenome.GRCh37 else "GRCh38"

    # Susun command Nextflow
    cmd = [
        NEXTFLOW_BIN, "run", NEXTFLOW_SCRIPT,
        "-profile",       NEXTFLOW_PROFILE,
        "--ref_genome",   ref_genome,
        "--outdir",       outdir,
        "--ploidy",       ploidy,
        "--threads",      str(DEFAULT_THREADS),
        "--bin_dir",      BIN_DIR,
        "--samtools_bin", SAMTOOLS_BIN,
        "--bcftools_bin", BCFTOOLS_BIN,
        "--fastqc_path",  FASTQC_PATH,
        "-resume",
    ]

    # Mode SE vs PE
    if end_type == EndType.PE and fastq_r1 and fastq_r2:
        cmd += ["--reads_r1", fastq_r1, "--reads_r2", fastq_r2]
    else:
        cmd += ["--reads", fastq_path]

    # WES/Targeted — tambah BED file
    if analysis_mode in (AnalysisMode.WES, AnalysisMode.TARGETED) and bed_file:
        cmd += ["--target_bed", bed_file]

    # Aktifkan scatter-gather untuk semua mode (default on)
    cmd += ["--scatter_gather", "true"]

    try:
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=os.path.dirname(NEXTFLOW_SCRIPT),
        )

        # Baca stdout dan stderr secara bersamaan
        async def read_stream(stream):
            async for line in stream:
                line_str = line.decode("utf-8", errors="replace").rstrip()
                _parse_nextflow_output(job_id, line_str)
                job_manager.update_elapsed(job_id)

        await asyncio.gather(
            read_stream(process.stdout),
            read_stream(process.stderr),
        )

        await process.wait()

        if process.returncode == 0:
            job_manager.update_job_status(job_id, StageStatus.COMPLETED)
        else:
            job_manager.update_job_status(job_id, StageStatus.FAILED)

    except Exception as e:
        job_manager.update_job_status(job_id, StageStatus.FAILED)
        job_manager.update_stage(job_id, "variant_calling", StageStatus.FAILED, error=str(e))
        # Simpan error ke job untuk debugging
        if job_id in job_manager._jobs:
            job_manager._jobs[job_id]["nextflow_log"] = f"Exception: {str(e)}"
            job_manager._save_job(job_id)


def _parse_nextflow_output(job_id: str, line: str):
    """Parse output Nextflow untuk update status tiap tahap."""

    stage_map = {
        # FastQC
        "FASTQC":           "fastqc",
        # Trimming
        "TRIMMOMATIC_SE":   "trimming",
        "TRIMMOMATIC_PE":   "trimming",
        # Alignment
        "BWA_ALIGN_SE":     "alignment",
        "BWA_ALIGN_PE":     "alignment",
        # Coverage & Chrom detection
        "COVERAGE_STATS":   "variant_calling",
        "GET_CHROMS":       "variant_calling",
        # BCFTools variant calling (scatter-gather)
        "CALL_VARIANTS":    "variant_calling",
        "GATHER_VCF":       "variant_calling",
    }

    # Progress mapping per proses BCFTools (urutan pipeline)
    bcftools_progress = {
        "COVERAGE_STATS":  10,
        "GET_CHROMS":      15,
        "CALL_VARIANTS":   50,
        "GATHER_VCF":      90,
    }

    # Format A: "Submitted process > PROCESS_NAME (sample)"  ← Nextflow 24.x DEBUG/INFO
    # Format B: "Task completed > TaskHandler[... name: PROCESS_NAME ..."
    # Format C: "[100%] process ~ PROCESS_NAME"  ← Nextflow progress bar
    # Format D: "process > PROCESS_NAME [100%]"

    for nf_name, stage in stage_map.items():
        if nf_name not in line:
            continue

        # Submitted — tandai running
        if "Submitted process" in line or "Submitted" in line:
            pct = bcftools_progress.get(nf_name, 10)
            current = job_manager.get_job(job_id)
            if current:
                cur_progress = current["stages"][stage]["progress"]
                if pct > cur_progress:
                    job_manager.update_stage(job_id, stage, StageStatus.RUNNING, pct)

        # Completed — tandai selesai per tahap
        elif "Task completed" in line or ("COMPLETED" in line and "exit: 0" in line):
            if nf_name == "FASTQC":
                job_manager.update_stage(job_id, "fastqc", StageStatus.COMPLETED, 100)
            elif nf_name in ("TRIMMOMATIC_SE", "TRIMMOMATIC_PE"):
                job_manager.update_stage(job_id, "trimming", StageStatus.COMPLETED, 100)
            elif nf_name in ("BWA_ALIGN_SE", "BWA_ALIGN_PE"):
                job_manager.update_stage(job_id, "alignment", StageStatus.COMPLETED, 100)
            elif nf_name == "GATHER_VCF":
                job_manager.update_stage(job_id, "variant_calling", StageStatus.COMPLETED, 100)
            elif nf_name in ("COVERAGE_STATS", "GET_CHROMS", "CALL_VARIANTS"):
                pct = bcftools_progress.get(nf_name, 50)
                job_manager.update_stage(job_id, "variant_calling", StageStatus.RUNNING, pct)

        # Format [100%] dari progress bar
        elif "[100%]" in line and ("✔" in line or "cached" in line.lower()):
            if nf_name == "FASTQC":
                job_manager.update_stage(job_id, "fastqc", StageStatus.COMPLETED, 100)
            elif nf_name in ("TRIMMOMATIC_SE", "TRIMMOMATIC_PE"):
                job_manager.update_stage(job_id, "trimming", StageStatus.COMPLETED, 100)
            elif nf_name in ("BWA_ALIGN_SE", "BWA_ALIGN_PE"):
                job_manager.update_stage(job_id, "alignment", StageStatus.COMPLETED, 100)
            elif nf_name == "GATHER_VCF":
                job_manager.update_stage(job_id, "variant_calling", StageStatus.COMPLETED, 100)

        # Error
        elif "failed" in line.lower() or ("exit:" in line and "exit: 0" not in line):
            job_manager.update_stage(job_id, stage, StageStatus.FAILED, error=line)

    # Deteksi pipeline selesai dari output VCF log
    if "✅ VCF Output" in line:
        job_manager.update_stage(job_id, "variant_calling", StageStatus.COMPLETED, 100)

    # Deteksi coverage dari output
    if "📊 COVERAGE:" in line:
        try:
            cov_str = line.split("📊 COVERAGE:")[1].replace("x", "").strip()
            coverage = float(cov_str)
            job_manager.update_coverage(job_id, coverage)
        except (ValueError, IndexError):
            pass
