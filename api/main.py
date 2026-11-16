"""
NGS E2E Secondary Analysis — REST API (BCFTools)
Fix:
- Streaming upload (chunk) untuk file besar — tidak load ke memori
- Sanitasi file.filename sebelum os.path.join (fix path traversal)
- Fix sample name extraction untuk format Illumina (S1_S1_L001_R1_001.fastq.gz)
- Fix download trimming PE (cari _R1.trimmed.fastq)
- Fix download VCF (.vcf.gz fallback)
- Recovery job saat startup
- Error response format sesuai API Contract
"""

import os
import re
import asyncio
from typing import Optional
from datetime import datetime

from fastapi import FastAPI, UploadFile, File, Form, HTTPException, Request
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
import aiofiles

from .models import (
    AnalysisMode, EndType, ReferenceGenome, VariantType,
    SubmitResponse, PipelineStatus, StageStatus
)
from .config import UPLOAD_DIR, RESULTS_DIR
from . import job_manager
from . import pipeline_runner
from .pipeline_runner import sanitize_sample_name

# ── App Init ──────────────────────────────────────────────────────────────────
app = FastAPI(
    title       = "NGS Secondary Analysis API — BCFTools",
    description = "REST API untuk pipeline NGS: FastQC → Trimming → BWA MEM → MarkDup → BCFTools mpileup/call/norm/filter",
    version     = "1.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins     = ["*"],
    allow_credentials = True,
    allow_methods     = ["*"],
    allow_headers     = ["*"],
)


@app.on_event("startup")
async def startup():
    """Recover jobs saat server start — mark job running sebagai failed."""
    job_manager.recover_jobs()


# ── Helper ────────────────────────────────────────────────────────────────────

def _safe_filename(filename: str) -> str:
    """
    Sanitasi filename untuk mencegah path traversal.
    Hanya ambil basename, hapus karakter berbahaya.
    """
    # Ambil basename saja (cegah ../../../etc/passwd)
    name = os.path.basename(filename)
    # Hapus karakter selain alfanumerik, titik, strip, underscore
    name = re.sub(r'[^A-Za-z0-9._\-]', '_', name)
    return name or "upload"


async def save_upload(file: UploadFile, dest_dir: str) -> str:
    """
    Simpan file upload ke disk dengan streaming (chunk 1MB).
    Fix: tidak load seluruh file ke memori.
    Fix: sanitasi filename sebelum os.path.join.
    """
    os.makedirs(dest_dir, exist_ok=True)
    safe_name = _safe_filename(file.filename or "upload")
    dest      = os.path.join(dest_dir, safe_name)
    chunk_size = 1024 * 1024  # 1 MB per chunk
    async with aiofiles.open(dest, "wb") as f:
        while True:
            chunk = await file.read(chunk_size)
            if not chunk:
                break
            await f.write(chunk)
    return dest


def _extract_sample_name(filename: str, end_type: EndType) -> str:
    """
    Fix: ekstrak sample name dari format Illumina maupun format sederhana.
    Illumina: S1_S1_L001_R1_001.fastq.gz → S1_S1_L001
    Sederhana: NA12878_R1.fastq → NA12878
    """
    name = os.path.basename(filename)
    # Hapus ekstensi
    name = re.sub(r'\.(fastq\.gz|fq\.gz|fastq|fq)$', '', name, flags=re.IGNORECASE)
    if end_type == EndType.PE:
        # Format Illumina: _R1_001, _R1_002, _1_001, dll
        name = re.sub(r'[_\.]R[12][_\.]?\d*$', '', name)
        name = re.sub(r'[_\.]?[12][_\.]?\d*$', '', name)
        # Format sederhana: _R1, _1
        name = re.sub(r'[_\.]R[12]$', '', name)
        name = re.sub(r'[_\.][12]$', '', name)
    return sanitize_sample_name(name) or "sample"


def get_stage_output(job_id: str, stage: str) -> Optional[str]:
    """Cari file output untuk tahap tertentu."""
    job = job_manager.get_job(job_id)
    if not job:
        return None
    sample = job["sample_name"]
    base   = os.path.join(RESULTS_DIR, job_id, sample)

    paths = {
        "fastqc": [
            os.path.join(base, "fastqc", f"{sample}_fastqc.html"),
            os.path.join(base, "fastqc", f"{sample}_fastqc.zip"),
        ],
        "trimming": [
            # SE
            os.path.join(base, "trimmed", f"{sample}.trimmed.fastq"),
            # PE R1
            os.path.join(base, "trimmed", f"{sample}_R1.trimmed.fastq"),
        ],
        "alignment": [
            os.path.join(base, "bwa", f"{sample}.markdup.bam"),
            os.path.join(base, "bwa", f"{sample}.aligned.sorted.bam"),
        ],
        "variant_calling": [
            # bgzipped (output baru)
            os.path.join(base, "vcf", f"{sample}.filtered.vcf.gz"),
            # plain vcf (fallback)
            os.path.join(base, "vcf", f"{sample}.filtered.vcf"),
        ],
    }

    for path in paths.get(stage, []):
        if os.path.exists(path):
            return path
    return None


def _validate_fastq_filename(filename: str) -> bool:
    """
    Fix #23: validasi ekstensi filename sesuai API Contract.
    Hanya .fastq.gz atau .fq.gz (huruf kecil), max 255 karakter.
    """
    if not filename or len(filename) > 255:
        return False
    return bool(re.search(r'\.(fastq\.gz|fq\.gz)$', filename))


def _error_response(code: str, message: str, http_status: int, details: dict = None):
    """Format error sesuai API Contract section 4."""
    body = {"error": {"code": code, "message": message}}
    if details:
        body["error"]["details"] = details
    raise HTTPException(status_code=http_status, detail=body)


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Submit Pipeline
# ═══════════════════════════════════════════════════════════════════════════════

@app.post("/api/pipeline/submit", response_model=SubmitResponse, tags=["Pipeline"])
async def submit_pipeline(
    analysis_mode: AnalysisMode    = Form(AnalysisMode.WGS),
    end_type:      EndType         = Form(EndType.SE),
    reference:     ReferenceGenome = Form(ReferenceGenome.GRCh38),
    variant_type:  VariantType     = Form(VariantType.GERMLINE),
    sample_sex:    Optional[str]   = Form(None),   # Fix #19: "male" atau "female"
    fastq:         Optional[UploadFile] = File(None),
    fastq_r1:      Optional[UploadFile] = File(None),
    fastq_r2:      Optional[UploadFile] = File(None),
    bed_file:      Optional[UploadFile] = File(None),
):
    if variant_type == VariantType.SOMATIC:
        _error_response("VALIDATION_FAILED", "Somatic variant calling belum didukung.", 422)

    if end_type == EndType.PE and (not fastq_r1 or not fastq_r2):
        _error_response("VALIDATION_FAILED", "PE mode membutuhkan fastq_r1 dan fastq_r2.", 400)
    if end_type == EndType.SE and not fastq:
        _error_response("VALIDATION_FAILED", "SE mode membutuhkan fastq.", 400)
    if analysis_mode in (AnalysisMode.WES, AnalysisMode.TARGETED) and not bed_file:
        _error_response("VALIDATION_FAILED", f"{analysis_mode.value} membutuhkan bed_file.", 400)

    # Fix #23: validasi ekstensi filename
    files_to_validate = []
    if end_type == EndType.PE:
        files_to_validate = [fastq_r1, fastq_r2]
    else:
        files_to_validate = [fastq]

    invalid_files = []
    for f in files_to_validate:
        if f and not _validate_fastq_filename(f.filename or ""):
            invalid_files.append(f.filename)
    if invalid_files:
        _error_response(
            "VALIDATION_FAILED",
            "Ekstensi file tidak valid. Hanya .fastq.gz atau .fq.gz yang didukung.",
            400,
            {"files": invalid_files}
        )

    # Fix #19: validasi nilai sex
    if sample_sex and sample_sex.lower() not in ("male", "female"):
        _error_response("VALIDATION_FAILED", "sample_sex harus 'male' atau 'female'.", 400)
    sex = sample_sex.lower() if sample_sex else None

    # Ekstrak dan sanitasi sample name
    source_filename = fastq_r1.filename if end_type == EndType.PE else fastq.filename
    sample_name = _extract_sample_name(source_filename or "sample", end_type)

    job_id  = job_manager.create_job(sample_name, analysis_mode, end_type, reference, variant_type)
    job_dir = os.path.join(UPLOAD_DIR, job_id)

    fastq_path = fastq_r1_path = fastq_r2_path = bed_path = None

    if end_type == EndType.PE:
        fastq_r1_path = await save_upload(fastq_r1, job_dir)
        fastq_r2_path = await save_upload(fastq_r2, job_dir)
    else:
        fastq_path = await save_upload(fastq, job_dir)

    if bed_file:
        bed_path = await save_upload(bed_file, job_dir)

    asyncio.create_task(
        pipeline_runner.run_pipeline(
            job_id        = job_id,
            fastq_path    = fastq_path,
            fastq_r1      = fastq_r1_path,
            fastq_r2      = fastq_r2_path,
            bed_file      = bed_path,
            analysis_mode = analysis_mode,
            end_type      = end_type,
            reference     = reference,
            sample_sex    = sex,            # Fix #19
        )
    )

    return SubmitResponse(
        job_id  = job_id,
        message = f"Pipeline submitted untuk sampel '{sample_name}'",
        status  = "running",
    )


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Status & Jobs
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/pipeline/{job_id}/status", tags=["Pipeline"])
async def get_pipeline_status(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        _error_response("NOT_FOUND", f"Job '{job_id}' tidak ditemukan.", 404)
    job_manager.update_elapsed(job_id)
    return job_manager.get_job(job_id)


@app.get("/api/pipeline/jobs", tags=["Pipeline"])
async def list_all_jobs():
    return job_manager.list_jobs()


@app.get("/api/pipeline/queue", tags=["Pipeline"])
async def get_queue():
    return pipeline_runner.get_queue_status()


@app.delete("/api/pipeline/{job_id}", tags=["Pipeline"])
async def delete_job(job_id: str):
    """Hapus job: stop Nextflow + hapus folder uploads/results."""
    success = job_manager.delete_job(job_id)
    if not success:
        _error_response("NOT_FOUND", f"Job '{job_id}' tidak ditemukan.", 404)
    return {"message": f"Job '{job_id}' berhasil dihapus."}


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Download per Tahap
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/pipeline/{job_id}/fastqc/download", tags=["FastQC"])
async def download_fastqc(job_id: str):
    path = get_stage_output(job_id, "fastqc")
    if not path:
        _error_response("NOT_FOUND", "FastQC output belum tersedia.", 404)
    return FileResponse(path, filename=os.path.basename(path), media_type="text/html")


@app.get("/api/pipeline/{job_id}/fastqc/status", tags=["FastQC"])
async def fastqc_status(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        _error_response("NOT_FOUND", "Job tidak ditemukan.", 404)
    return job["stages"]["fastqc"]


@app.get("/api/pipeline/{job_id}/trimming/download", tags=["Trimming"])
async def download_trimming(job_id: str):
    path = get_stage_output(job_id, "trimming")
    if not path:
        _error_response("NOT_FOUND", "Trimming output belum tersedia.", 404)
    return FileResponse(path, filename=os.path.basename(path), media_type="application/octet-stream")


@app.get("/api/pipeline/{job_id}/trimming/status", tags=["Trimming"])
async def trimming_status(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        _error_response("NOT_FOUND", "Job tidak ditemukan.", 404)
    return job["stages"]["trimming"]


@app.get("/api/pipeline/{job_id}/alignment/download", tags=["BWA Alignment"])
async def download_alignment(job_id: str):
    path = get_stage_output(job_id, "alignment")
    if not path:
        _error_response("NOT_FOUND", "Alignment output belum tersedia.", 404)
    return FileResponse(path, filename=os.path.basename(path), media_type="application/octet-stream")


@app.get("/api/pipeline/{job_id}/alignment/status", tags=["BWA Alignment"])
async def alignment_status(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        _error_response("NOT_FOUND", "Job tidak ditemukan.", 404)
    return job["stages"]["alignment"]


@app.get("/api/pipeline/{job_id}/variant-calling/download", tags=["Variant Calling"])
async def download_vcf(job_id: str):
    path = get_stage_output(job_id, "variant_calling")
    if not path:
        _error_response("NOT_FOUND", "VCF output belum tersedia.", 404)
    # Tentukan media type berdasarkan ekstensi
    media_type = "application/octet-stream" if path.endswith(".gz") else "text/plain"
    return FileResponse(path, filename=os.path.basename(path), media_type=media_type)


@app.get("/api/pipeline/{job_id}/variant-calling/status", tags=["Variant Calling"])
async def variant_calling_status(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        _error_response("NOT_FOUND", "Job tidak ditemukan.", 404)
    return job["stages"]["variant_calling"]


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Health & Info
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/health", tags=["System"])
async def health_check():
    return {
        "status":   "ok",
        "version":  "1.1.0",
        "backend":  "BCFTools",
        "pipeline": "FastQC → Trimmomatic → BWA MEM → MarkDup → BCFTools mpileup/call/norm/filter",
        "time":     datetime.utcnow().isoformat() + "Z",
    }


@app.get("/api/references", tags=["System"])
async def list_references():
    from .config import REF_GENOMES
    result = {}
    for key, path in REF_GENOMES.items():
        result[key] = {
            "path":      path,
            "available": os.path.exists(path),
        }
    return result


@app.get("/", tags=["System"])
async def root():
    return {
        "name":    "NGS Secondary Analysis API — BCFTools",
        "version": "1.1.0",
        "backend": "BCFTools mpileup/call/norm/filter",
        "docs":    "/docs",
        "health":  "/api/health",
    }
