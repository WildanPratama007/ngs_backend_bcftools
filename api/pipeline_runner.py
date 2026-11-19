"""
NGS API — Pipeline Runner (BCFTools)
Fix:
- MAX_CONCURRENT_JOBS dari config (tidak hardcode)
- Streaming upload (chunk) untuk file besar
- Sanitasi sample_name untuk mencegah path traversal & shell injection
- Pemulihan state via job_manager.recover_jobs()
"""

import asyncio
import os
import re
from .models import StageStatus, AnalysisMode, EndType, ReferenceGenome
from .config import (
    NEXTFLOW_BIN, NEXTFLOW_SCRIPT, NEXTFLOW_PROFILE,
    REF_GENOMES, RESULTS_DIR, DEFAULT_THREADS, BIN_DIR,
    SAMTOOLS_BIN, BCFTOOLS_BIN, FASTQC_PATH, MAX_CONCURRENT_JOBS
)
from . import job_manager

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
    # Trimming parameters — user-configurable
    trim_quality:  int  = 20,    # SLIDINGWINDOW quality (Q10-Q40)
    trim_min_len:  int  = 36,    # MINLEN (20-100)
    trim_adapter:  bool = True,  # pakai ILLUMINACLIP
    trim_poly_g:   bool = False, # tambah poly-G (NovaSeq)
):
    _queue.append(job_id)
    job_manager.update_job_status(job_id, StageStatus.PENDING)

    async with _semaphore:
        if job_id in _queue:
            _queue.remove(job_id)

        await _execute_pipeline(
            job_id          = job_id,
            fastq_path      = fastq_path,
            fastq_r1        = fastq_r1,
            fastq_r2        = fastq_r2,
            bed_file        = bed_file,
            analysis_mode   = analysis_mode,
            end_type        = end_type,
            reference       = reference,
            sample_sex      = sample_sex,
            trim_quality    = trim_quality,
            trim_min_len    = trim_min_len,
            trim_adapter    = trim_adapter,
            trim_poly_g     = trim_poly_g,
            mark_duplicates = mark_duplicates,
            min_qual        = min_qual,
            min_dp          = min_dp,
            min_gq          = min_gq,
            max_depth       = max_depth,
            min_bq          = min_bq,
            min_mq          = min_mq,
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
    sample_sex:    str = None,
    trim_quality:  int  = 20,
    trim_min_len:  int  = 36,
    trim_adapter:  bool = True,
    trim_poly_g:   bool = False,
    mark_duplicates: bool = True,
    # Variant calling parameters
    min_qual:   int = 20,    # Minimum QUAL score
    min_dp:     int = 10,    # Minimum depth (DP)
    min_gq:     int = 20,    # Minimum Genotype Quality (GQ)
    max_depth:  int = 500,   # Max depth per posisi (500=WGS, 2000=WES panel)
    min_bq:     int = 20,    # Minimum base quality di mpileup
    min_mq:     int = 20,    # Minimum mapping quality di mpileup
):
    job_manager.update_job_status(job_id, StageStatus.RUNNING)
    ref_genome = REF_GENOMES.get(reference.value, REF_GENOMES["GRCh38"])
    outdir     = os.path.join(RESULTS_DIR, job_id)
    ploidy     = "GRCh37" if reference == ReferenceGenome.GRCh37 else "GRCh38"
    work_dir   = os.path.join(outdir, ".work")

    cmd = [
        NEXTFLOW_BIN, "run", NEXTFLOW_SCRIPT,
        "-profile",       NEXTFLOW_PROFILE,
        "-w",             work_dir,
        "--ref_genome",   ref_genome,
        "--outdir",       outdir,
        "--ploidy",       ploidy,
        "--threads",      str(DEFAULT_THREADS),
        "--bin_dir",      BIN_DIR,
        "--samtools_bin", SAMTOOLS_BIN,
        "--bcftools_bin", BCFTOOLS_BIN,
        "--fastqc_path",  FASTQC_PATH,
        # Trimming params dari user
        "--trim_sliding_window_quality", str(trim_quality),
        "--trim_min_len",                str(trim_min_len),
        "--trim_adapter",                str(trim_adapter).lower(),
        "--trim_poly_g",                 str(trim_poly_g).lower(),
        "--mark_duplicates",             str(mark_duplicates).lower(),
        # Variant calling params dari user
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
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=os.path.dirname(NEXTFLOW_SCRIPT),
        )

        # Simpan PID untuk cancel
        job_manager.register_nextflow_pid(job_id, process.pid)

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
            # Parse coverage JSON file untuk update extended stats
            _load_coverage_json(job_id, outdir)
            job_manager.update_job_status(job_id, StageStatus.COMPLETED)
        else:
            job_manager.update_job_status(job_id, StageStatus.FAILED)

    except Exception as e:
        job_manager.update_job_status(job_id, StageStatus.FAILED)
        job_manager.update_stage(job_id, "variant_calling", StageStatus.FAILED, error=str(e))
        if job_id in job_manager._jobs:
            job_manager._jobs[job_id]["nextflow_log"] = f"Exception: {str(e)}"
            job_manager._save_job(job_id)


def _load_coverage_json(job_id: str, outdir: str):
    """Baca coverage JSON file setelah pipeline selesai dan update job."""
    import json, glob
    pattern = os.path.join(outdir, "**", "*.coverage.json")
    files = glob.glob(pattern, recursive=True)
    if not files:
        return
    try:
        with open(files[0]) as f:
            cov_data = json.load(f)
        job_manager.update_coverage(job_id, cov_data)
    except Exception:
        pass


def _parse_nextflow_output(job_id: str, line: str):
    """Parse output Nextflow untuk update status tiap tahap."""

    stage_map = {
        "FASTQC":           "fastqc",
        "TRIMMOMATIC_SE":   "trimming",
        "TRIMMOMATIC_PE":   "trimming",
        "BWA_ALIGN_SE":     "alignment",
        "BWA_ALIGN_PE":     "alignment",
        "MARKDUPLICATES":   "alignment",
        "COVERAGE_STATS":   "variant_calling",
        "GET_CHROMS":       "variant_calling",
        "CALL_VARIANTS":    "variant_calling",
        "GATHER_VCF":       "variant_calling",
    }

    bcftools_progress = {
        "MARKDUPLICATES":  5,
        "COVERAGE_STATS":  10,
        "GET_CHROMS":      15,
        "CALL_VARIANTS":   50,
        "GATHER_VCF":      90,
    }

    for nf_name, stage in stage_map.items():
        if nf_name not in line:
            continue

        if "Submitted process" in line or "Submitted" in line:
            pct = bcftools_progress.get(nf_name, 10)
            current = job_manager.get_job(job_id)
            if current:
                cur_progress = current["stages"][stage]["progress"]
                if pct > cur_progress:
                    job_manager.update_stage(job_id, stage, StageStatus.RUNNING, pct)

        elif "Task completed" in line or ("COMPLETED" in line and "exit: 0" in line):
            if nf_name == "FASTQC":
                job_manager.update_stage(job_id, "fastqc", StageStatus.COMPLETED, 100)
            elif nf_name in ("TRIMMOMATIC_SE", "TRIMMOMATIC_PE"):
                job_manager.update_stage(job_id, "trimming", StageStatus.COMPLETED, 100)
            elif nf_name == "MARKDUPLICATES":
                job_manager.update_stage(job_id, "alignment", StageStatus.COMPLETED, 100)
            elif nf_name in ("BWA_ALIGN_SE", "BWA_ALIGN_PE"):
                job_manager.update_stage(job_id, "alignment", StageStatus.RUNNING, 50)
            elif nf_name == "GATHER_VCF":
                job_manager.update_stage(job_id, "variant_calling", StageStatus.COMPLETED, 100)
            elif nf_name in ("COVERAGE_STATS", "GET_CHROMS", "CALL_VARIANTS"):
                pct = bcftools_progress.get(nf_name, 50)
                job_manager.update_stage(job_id, "variant_calling", StageStatus.RUNNING, pct)

        elif "[100%]" in line and ("✔" in line or "cached" in line.lower()):
            if nf_name == "FASTQC":
                job_manager.update_stage(job_id, "fastqc", StageStatus.COMPLETED, 100)
            elif nf_name in ("TRIMMOMATIC_SE", "TRIMMOMATIC_PE"):
                job_manager.update_stage(job_id, "trimming", StageStatus.COMPLETED, 100)
            elif nf_name == "MARKDUPLICATES":
                job_manager.update_stage(job_id, "alignment", StageStatus.COMPLETED, 100)
            elif nf_name == "GATHER_VCF":
                job_manager.update_stage(job_id, "variant_calling", StageStatus.COMPLETED, 100)

        elif "failed" in line.lower() or ("exit:" in line and "exit: 0" not in line):
            job_manager.update_stage(job_id, stage, StageStatus.FAILED, error=line)

    if "✅ VCF Output" in line:
        job_manager.update_stage(job_id, "variant_calling", StageStatus.COMPLETED, 100)

    # Deteksi coverage dari output — sekarang format JSON
    if "📊 COVERAGE:" in line:
        try:
            cov_str = line.split("📊 COVERAGE:")[1].replace("x", "").split("|")[0].strip()
            mean_cov = float(cov_str)
            # Simpan sebagai object minimal kalau JSON file belum di-parse
            # Job manager akan diupdate lagi saat JSON file dibaca
            job_manager.update_coverage(job_id, {"mean_depth": mean_cov})
        except (ValueError, IndexError):
            pass
