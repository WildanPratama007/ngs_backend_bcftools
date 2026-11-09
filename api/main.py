"""
NGS E2E Secondary Analysis — REST API (BCFTools)
FastAPI Backend untuk pipeline NGS: FastQC → Trimming → BWA MEM → BCFTools Variant Calling
"""

import os
import asyncio
import uuid
from typing import Optional, List
from datetime import datetime

from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import aiofiles

from .models import (
    AnalysisMode, EndType, ReferenceGenome, VariantType,
    SubmitResponse, PipelineStatus, StageStatus
)
from .config import UPLOAD_DIR, RESULTS_DIR
from . import job_manager
from . import pipeline_runner

# ── App Init ──────────────────────────────────────────────────────────────────
app = FastAPI(
    title       = "NGS Secondary Analysis API — BCFTools",
    description = "REST API untuk pipeline NGS: FastQC → Trimming → BWA MEM → BCFTools mpileup/call/filter",
    version     = "1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins     = ["*"],
    allow_credentials = True,
    allow_methods     = ["*"],
    allow_headers     = ["*"],
)

# ── Helper ────────────────────────────────────────────────────────────────────
async def save_upload(file: UploadFile, dest_dir: str) -> str:
    """Simpan file upload ke disk, return path."""
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, file.filename)
    async with aiofiles.open(dest, "wb") as f:
        content = await file.read()
        await f.write(content)
    return dest


def get_stage_output(job_id: str, stage: str) -> Optional[str]:
    """Cari file output untuk tahap tertentu."""
    job     = job_manager.get_job(job_id)
    if not job:
        return None
    sample  = job["sample_name"]
    base    = os.path.join(RESULTS_DIR, job_id, sample)

    paths = {
        "fastqc":          [
            os.path.join(base, "fastqc", f"{sample}_fastqc.html"),
            os.path.join(base, "fastqc", f"{sample}_fastqc.zip"),
        ],
        "trimming":        [os.path.join(base, "trimmed", f"{sample}.trimmed.fastq")],
        "alignment":       [os.path.join(base, "bwa", f"{sample}.aligned.sorted.bam")],
        "variant_calling": [os.path.join(base, "vcf", f"{sample}.filtered.vcf")],
    }

    for path in paths.get(stage, []):
        if os.path.exists(path):
            return path
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Submit Pipeline
# ═══════════════════════════════════════════════════════════════════════════════

@app.post("/api/pipeline/submit", response_model=SubmitResponse, tags=["Pipeline"])
async def submit_pipeline(
    analysis_mode: AnalysisMode = Form(AnalysisMode.WGS),
    end_type:      EndType      = Form(EndType.SE),
    reference:     ReferenceGenome = Form(ReferenceGenome.GRCh38),
    variant_type:  VariantType  = Form(VariantType.GERMLINE),
    fastq:         Optional[UploadFile] = File(None),
    fastq_r1:      Optional[UploadFile] = File(None),
    fastq_r2:      Optional[UploadFile] = File(None),
    bed_file:      Optional[UploadFile] = File(None),
):
    # Validasi somatic — belum didukung
    if variant_type == VariantType.SOMATIC:
        raise HTTPException(
            status_code=422,
            detail="Somatic variant calling belum didukung. Coming soon."
        )

    # Validasi input file
    if end_type == EndType.PE and (not fastq_r1 or not fastq_r2):
        raise HTTPException(status_code=400, detail="PE mode membutuhkan fastq_r1 dan fastq_r2.")
    if end_type == EndType.SE and not fastq:
        raise HTTPException(status_code=400, detail="SE mode membutuhkan fastq.")
    if analysis_mode in (AnalysisMode.WES, AnalysisMode.TARGETED) and not bed_file:
        raise HTTPException(status_code=400, detail=f"{analysis_mode} membutuhkan bed_file.")

    # Tentukan sample name
    if end_type == EndType.PE:
        sample_name = fastq_r1.filename.replace("_R1.fastq.gz","").replace("_1.fastq.gz","").replace(".fastq","")
    else:
        sample_name = fastq.filename.replace(".fastq.gz","").replace(".fastq","")

    # Buat job
    job_id   = job_manager.create_job(sample_name, analysis_mode, end_type, reference, variant_type)
    job_dir  = os.path.join(UPLOAD_DIR, job_id)

    # Simpan file upload
    fastq_path = fastq_r1_path = fastq_r2_path = bed_path = None

    if end_type == EndType.PE:
        fastq_r1_path = await save_upload(fastq_r1, job_dir)
        fastq_r2_path = await save_upload(fastq_r2, job_dir)
    else:
        fastq_path = await save_upload(fastq, job_dir)

    if bed_file:
        bed_path = await save_upload(bed_file, job_dir)

    # Jalankan pipeline dengan asyncio.create_task — tidak blocking
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
        )
    )

    return SubmitResponse(
        job_id  = job_id,
        message = f"Pipeline submitted untuk sampel '{sample_name}'",
        status  = "running",
    )


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Status
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/pipeline/{job_id}/status", tags=["Pipeline"])
async def get_pipeline_status(job_id: str):
    """Cek status pipeline dan semua tahap."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' tidak ditemukan.")
    job_manager.update_elapsed(job_id)
    return job_manager.get_job(job_id)


@app.get("/api/pipeline/jobs", tags=["Pipeline"])
async def list_all_jobs():
    """List semua job yang pernah disubmit."""
    return job_manager.list_jobs()


@app.get("/api/pipeline/queue", tags=["Pipeline"])
async def get_queue():
    """Status antrian pipeline saat ini."""
    return pipeline_runner.get_queue_status()


@app.delete("/api/pipeline/{job_id}", tags=["Pipeline"])
async def delete_job(job_id: str):
    """Hapus job dari sistem."""
    success = job_manager.delete_job(job_id)
    if not success:
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' tidak ditemukan.")
    return {"message": f"Job '{job_id}' berhasil dihapus."}


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Download per Tahap
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/pipeline/{job_id}/fastqc/download", tags=["FastQC"])
async def download_fastqc(job_id: str):
    """Download FastQC HTML report."""
    path = get_stage_output(job_id, "fastqc")
    if not path:
        raise HTTPException(status_code=404, detail="FastQC output belum tersedia.")
    return FileResponse(path, filename=os.path.basename(path), media_type="text/html")


@app.get("/api/pipeline/{job_id}/fastqc/status", tags=["FastQC"])
async def fastqc_status(job_id: str):
    """Status tahap FastQC."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job tidak ditemukan.")
    return job["stages"]["fastqc"]


@app.get("/api/pipeline/{job_id}/trimming/download", tags=["Trimming"])
async def download_trimming(job_id: str):
    """Download Trimmed FASTQ."""
    path = get_stage_output(job_id, "trimming")
    if not path:
        raise HTTPException(status_code=404, detail="Trimming output belum tersedia.")
    return FileResponse(path, filename=os.path.basename(path), media_type="application/octet-stream")


@app.get("/api/pipeline/{job_id}/trimming/status", tags=["Trimming"])
async def trimming_status(job_id: str):
    """Status tahap Trimming."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job tidak ditemukan.")
    return job["stages"]["trimming"]


@app.get("/api/pipeline/{job_id}/alignment/download", tags=["BWA Alignment"])
async def download_alignment(job_id: str):
    """Download Sorted BAM."""
    path = get_stage_output(job_id, "alignment")
    if not path:
        raise HTTPException(status_code=404, detail="Alignment output belum tersedia.")
    return FileResponse(path, filename=os.path.basename(path), media_type="application/octet-stream")


@app.get("/api/pipeline/{job_id}/alignment/status", tags=["BWA Alignment"])
async def alignment_status(job_id: str):
    """Status tahap Alignment."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job tidak ditemukan.")
    return job["stages"]["alignment"]


@app.get("/api/pipeline/{job_id}/variant-calling/download", tags=["Variant Calling"])
async def download_vcf(job_id: str):
    """Download Filtered VCF."""
    path = get_stage_output(job_id, "variant_calling")
    if not path:
        raise HTTPException(status_code=404, detail="VCF output belum tersedia.")
    return FileResponse(path, filename=os.path.basename(path), media_type="text/plain")


@app.get("/api/pipeline/{job_id}/variant-calling/status", tags=["Variant Calling"])
async def variant_calling_status(job_id: str):
    """Status tahap Variant Calling."""
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job tidak ditemukan.")
    return job["stages"]["variant_calling"]


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Health Check & Info
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/health", tags=["System"])
async def health_check():
    return {
        "status":   "ok",
        "version":  "1.0.0",
        "backend":  "BCFTools",
        "pipeline": "FastQC → Trimmomatic → BWA MEM → BCFTools mpileup/call/filter",
        "time":     datetime.now().isoformat(),
    }


@app.get("/api/references", tags=["System"])
async def list_references():
    """List referensi genom yang tersedia."""
    from .config import REF_GENOMES
    result = {}
    for key, path in REF_GENOMES.items():
        # Convert WSL path ke Windows path untuk cek keberadaan file
        win_path = path.replace("/mnt/c/", "C:/").replace("/", "/")
        result[key] = {
            "path":      path,
            "available": os.path.exists(path) or os.path.exists(win_path),
        }
    return result


@app.get("/", tags=["System"])
async def root():
    return {
        "name":     "NGS Secondary Analysis API — BCFTools",
        "version":  "1.0.0",
        "backend":  "BCFTools mpileup/call/filter",
        "docs":     "/docs",
        "health":   "/api/health",
    }
