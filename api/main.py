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

import glob
import logging
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
    AnalysisMode, EndType, ReferenceGenome, VariantType, TrimPreset,
    SubmitResponse, PipelineStatus, StageStatus
)
from .config import UPLOAD_DIR, RESULTS_DIR, ADAPTER_KIT_MAP, CORS_ORIGINS, PIPELINE_MODE
from . import job_manager
from . import pipeline_runner
from . import preflight
from .contract import router as contract_router, cancel_run
from .errors import error_response, install_handlers
from .pipeline_runner import sanitize_sample_name

logging.basicConfig(level=os.getenv("LOG_LEVEL", "INFO"),
                    format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger("ngs.api")

# ── App Init ──────────────────────────────────────────────────────────────────
app = FastAPI(
    title       = "NGS Secondary Analysis API — BCFTools",
    description = "REST API untuk pipeline NGS: FastQC → Trimming → BWA MEM → MarkDup → BCFTools mpileup/call/norm/filter",
    version     = "1.1.0",
)

install_handlers(app)
app.include_router(contract_router)

# Origin dari env CORS_ORIGINS. Credentials hanya diaktifkan untuk daftar origin eksplisit
# (browser menolak kombinasi wildcard "*" dengan credentials).
_wildcard = "*" in CORS_ORIGINS
app.add_middleware(
    CORSMiddleware,
    allow_origins     = ["*"] if _wildcard else CORS_ORIGINS,
    allow_credentials = not _wildcard,
    allow_methods     = ["*"],
    allow_headers     = ["*"],
)


@app.on_event("startup")
async def startup():
    """Recover jobs saat server start dan cek tool pipeline."""
    job_manager.recover_jobs()
    preflight.run_preflight()
    if _wildcard:
        logger.warning("CORS_ORIGINS='*' — isi domain frontend di production.")


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
        # Hapus penanda read hanya sekali: _R1, _R1_001, .R1, _1, _1_001
        # (regex lama dijalankan berurutan dan ikut memotong digit di akhir nama sampel)
        name = re.sub(r'[_.]R?1(?:_\d+)?$', '', name)
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

    # FastQC menamai output dari nama file input (PE: dua laporan, R1 dan R2) — ambil yang ada
    if stage == "fastqc":
        reports = sorted(glob.glob(os.path.join(base, "fastqc", "*_fastqc.html")))
        if reports:
            return reports[0]
    return None


def _validate_fastq_filename(filename: str) -> bool:
    """
    Validasi ekstensi filename sesuai API Contract.
    Production: .fastq.gz atau .fq.gz
    Testing: juga terima .fastq dan .fq (uncompressed)
    """
    if not filename or len(filename) > 255:
        return False
    return bool(re.search(r'\.(fastq\.gz|fq\.gz|fastq|fq)$', filename, re.IGNORECASE))


_error_response = error_response   # envelope kontrak: lihat api/errors.py


def _require_pipeline_ready():
    """Tolak pekerjaan baru kalau pipeline dimatikan atau tool wajib tidak lengkap (kontrak: 503)."""
    if PIPELINE_MODE == "disabled":
        _error_response("PIPELINE_DISABLED", "Pipeline dinonaktifkan di server ini.", 503)
    missing = preflight.missing_required()
    if missing:
        logger.error("Pipeline belum siap: %s", "; ".join(missing))
        _error_response("PIPELINE_NOT_CONFIGURED", "Tool pipeline belum lengkap di server.", 503)


def _require_valid_id(analysis_id: str):
    if not job_manager.is_valid_id(analysis_id):
        _error_response("VALIDATION_FAILED",
                        "analysis_id hanya boleh huruf, angka, '_', '-', '.' (maks 64, diawali huruf/angka).", 400,
                        {"fields": [{"path": "analysis_id", "reason": "format tidak valid"}]})


def _validate_trim_params(quality: int, min_len: int):
    if not (10 <= quality <= 40):
        _error_response("VALIDATION_FAILED", "trim_quality harus antara 10-40.", 400)
    if not (20 <= min_len <= 100):
        _error_response("VALIDATION_FAILED", "trim_min_len harus antara 20-100.", 400)


def _validate_variant_params(min_qual: int, min_dp: int, min_gq: int,
                             max_depth: int, min_bq: int, min_mq: int):
    if not (0 <= min_qual <= 999):
        _error_response("VALIDATION_FAILED", "min_qual harus antara 0-999.", 400)
    if not (1 <= min_dp <= 500):
        _error_response("VALIDATION_FAILED", "min_dp harus antara 1-500.", 400)
    if not (0 <= min_gq <= 99):
        _error_response("VALIDATION_FAILED", "min_gq harus antara 0-99.", 400)
    if not (100 <= max_depth <= 10000):
        _error_response("VALIDATION_FAILED", "max_depth harus antara 100-10000.", 400)
    if not (0 <= min_bq <= 40):
        _error_response("VALIDATION_FAILED", "min_bq harus antara 0-40.", 400)
    if not (0 <= min_mq <= 60):
        _error_response("VALIDATION_FAILED", "min_mq harus antara 0-60.", 400)


def _validate_adapter_kit(adapter_kit: Optional[str]):
    if adapter_kit and adapter_kit not in ADAPTER_KIT_MAP:
        _error_response("VALIDATION_FAILED", f"adapter_kit tidak dikenal: '{adapter_kit}'.", 400,
                        {"allowed": sorted(ADAPTER_KIT_MAP)})


def _validate_sample_sex(sample_sex: Optional[str]) -> Optional[str]:
    if sample_sex and sample_sex.lower() not in ("male", "female"):
        _error_response("VALIDATION_FAILED", "sample_sex harus 'male' atau 'female'.", 400)
    return sample_sex.lower() if sample_sex else None


def _get_job_or_404(analysis_id: str) -> dict:
    job = job_manager.get_job(analysis_id)
    if not job:
        _error_response("NOT_FOUND", f"Analysis '{analysis_id}' tidak ditemukan. Jalankan Step 1 dulu.", 404)
    return job


def _ensure_not_running(job: dict, stage: str):
    if job["stages"][stage]["status"] == "running":
        _error_response("STATUS_CONFLICT", f"Step '{stage}' sedang berjalan.", 409)


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Submit Pipeline
# ═══════════════════════════════════════════════════════════════════════════════

@app.post("/api/pipeline/submit", response_model=SubmitResponse, tags=["Pipeline"])
async def submit_pipeline(
    analysis_mode: AnalysisMode    = Form(AnalysisMode.WGS),
    end_type:      EndType         = Form(EndType.SE),
    reference:     ReferenceGenome = Form(ReferenceGenome.GRCh38),
    variant_type:  VariantType     = Form(VariantType.GERMLINE),
    sample_sex:    Optional[str]   = Form(None),
    # ── Trimming parameters — user-configurable ──────────────────────────────
    # Preset shortcut: "research" = Q20/36bp, "clinical" = Q30/36bp
    trim_preset:   Optional[TrimPreset] = Form(None),
    trim_quality:  Optional[int]  = Form(None),   # Q10-Q40, default 20
    trim_min_len:  Optional[int]  = Form(None),   # 20-100bp, default 36
    trim_adapter:  Optional[bool] = Form(None),   # default True
    trim_poly_g:   Optional[bool] = Form(None),   # default False (NovaSeq only)
    adapter_kit:   Optional[str]  = Form(None),   # truseq_pe, truseq_se, nextera, ... (default: ikut SE/PE)
    # ── Alignment parameters ─────────────────────────────────────────────────
    mark_duplicates: Optional[bool] = Form(None),
    # ── Variant calling parameters ───────────────────────────────────────────
    min_qual:   Optional[int]  = Form(None),   # default 20
    min_dp:     Optional[int]  = Form(None),   # default 10
    min_gq:     Optional[int]  = Form(None),   # default 20
    max_depth:  Optional[int]  = Form(None),   # default 500 (WGS), 2000 (WES panel)
    min_bq:     Optional[int]  = Form(None),   # default 20 (minimum base quality)
    min_mq:     Optional[int]  = Form(None),   # default 20 (minimum mapping quality)
    fastq:         Optional[UploadFile] = File(None),
    fastq_r1:      Optional[UploadFile] = File(None),
    fastq_r2:      Optional[UploadFile] = File(None),
    bed_file:      Optional[UploadFile] = File(None),
):
    _require_pipeline_ready()
    if variant_type == VariantType.SOMATIC:
        _error_response("VALIDATION_FAILED", "Somatic variant calling belum didukung.", 400)

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
    sex = _validate_sample_sex(sample_sex)

    # ── Resolve trimming parameters ──────────────────────────────────────────
    # Preset menimpa nilai individual bila keduanya diberikan
    if trim_preset == TrimPreset.CLINICAL:
        final_trim_quality = 30
        final_trim_min_len = 36
        final_trim_adapter = True
        final_trim_poly_g  = False
    elif trim_preset == TrimPreset.RESEARCH:
        final_trim_quality = 20
        final_trim_min_len = 36
        final_trim_adapter = True
        final_trim_poly_g  = False
    else:
        # Pakai nilai individual, fallback ke default clinical bila tidak diisi
        final_trim_quality = trim_quality if trim_quality is not None else 20
        final_trim_min_len = trim_min_len if trim_min_len is not None else 36
        final_trim_adapter = trim_adapter if trim_adapter is not None else True
        final_trim_poly_g  = trim_poly_g  if trim_poly_g  is not None else False

    # Validasi range
    _validate_trim_params(final_trim_quality, final_trim_min_len)
    _validate_adapter_kit(adapter_kit)

    # Resolve mark_duplicates — default True
    final_mark_duplicates = mark_duplicates if mark_duplicates is not None else True

    # Resolve variant calling params — default clinical grade
    final_min_qual  = min_qual  if min_qual  is not None else 20
    final_min_dp    = min_dp    if min_dp    is not None else 10
    final_min_gq    = min_gq    if min_gq    is not None else 20
    final_max_depth = max_depth if max_depth is not None else 500
    final_min_bq    = min_bq    if min_bq    is not None else 20
    final_min_mq    = min_mq    if min_mq    is not None else 20

    _validate_variant_params(final_min_qual, final_min_dp, final_min_gq,
                             final_max_depth, final_min_bq, final_min_mq)

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
            sample_sex    = sex,
            adapter_kit   = adapter_kit,
            trim_quality  = final_trim_quality,
            trim_min_len    = final_trim_min_len,
            trim_adapter    = final_trim_adapter,
            trim_poly_g     = final_trim_poly_g,
            mark_duplicates = final_mark_duplicates,
            min_qual        = final_min_qual,
            min_dp          = final_min_dp,
            min_gq          = final_min_gq,
            max_depth       = final_max_depth,
            min_bq          = final_min_bq,
            min_mq          = final_min_mq,
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
# WIZARD ENDPOINTS — Step-by-step (sesuai UI wizard)
# Setiap step dijalankan secara terpisah dengan button "Run" di UI
# ═══════════════════════════════════════════════════════════════════════════════

@app.post("/api/analyses/{analysis_id}/qc/run", tags=["Wizard"])
async def run_qc(
    analysis_id: str,
    fastq:         Optional[UploadFile] = File(None),
    fastq_r1:      Optional[UploadFile] = File(None),
    fastq_r2:      Optional[UploadFile] = File(None),
    bed_file:      Optional[UploadFile] = File(None),
    end_type:      EndType         = Form(EndType.SE),
    reference:     ReferenceGenome = Form(ReferenceGenome.GRCh38),
    analysis_mode: AnalysisMode    = Form(AnalysisMode.WGS),
    adapter_kit:   Optional[str]   = Form(None),  # truseq_pe, nextera, truseq_se, dll
):
    """
    Step 1 — Upload FASTQ & jalankan FastQC.
    Upload selesai sebelum response dikembalikan; FastQC lanjut di background.
    Trimming (Step 2) bisa langsung dijalankan tanpa menunggu FastQC selesai.
    """
    _require_valid_id(analysis_id)
    _require_pipeline_ready()
    if end_type == EndType.PE and (not fastq_r1 or not fastq_r2):
        _error_response("VALIDATION_FAILED", "PE mode membutuhkan fastq_r1 dan fastq_r2.", 400)
    if end_type == EndType.SE and not fastq:
        _error_response("VALIDATION_FAILED", "SE mode membutuhkan fastq.", 400)

    files_to_validate = [fastq_r1, fastq_r2] if end_type == EndType.PE else [fastq]
    for f in [x for x in files_to_validate if x]:
        if not _validate_fastq_filename(f.filename or ""):
            _error_response("VALIDATION_FAILED", f"File {f.filename} harus berekstensi .fastq.gz atau .fq.gz", 400)

    needs_bed = analysis_mode in (AnalysisMode.WES, AnalysisMode.TARGETED)
    if needs_bed and not bed_file:
        _error_response("VALIDATION_FAILED", f"{analysis_mode.value} membutuhkan bed_file.", 400)
    _validate_adapter_kit(adapter_kit)

    source_filename = fastq_r1.filename if end_type == EndType.PE else fastq.filename
    sample_name = _extract_sample_name(source_filename or "sample", end_type)

    # Gunakan analysis_id dari frontend (atau buat baru)
    job_id  = analysis_id
    job_dir = os.path.join(UPLOAD_DIR, job_id)

    existing = job_manager.get_job(job_id)
    if existing:
        for stage in ("fastqc", "trimming", "alignment", "variant_calling"):
            _ensure_not_running(existing, stage)

    # Simpan FASTQ
    fastq_path = fastq_r1_path = fastq_r2_path = None
    if end_type == EndType.PE:
        fastq_r1_path = await save_upload(fastq_r1, job_dir)
        fastq_r2_path = await save_upload(fastq_r2, job_dir)
    else:
        fastq_path = await save_upload(fastq, job_dir)

    bed_path = await save_upload(bed_file, job_dir) if needs_bed else None

    # Buat job baru, atau perbarui metadata kalau analysis_id dipakai ulang
    if not existing:
        job_manager.create_job_with_id(job_id, sample_name, analysis_mode, end_type, reference, VariantType.GERMLINE)
    if existing:  # input baru → hasil step lama tidak berlaku lagi
        for stage in ("fastqc", "trimming", "alignment", "variant_calling"):
            job_manager.reset_stage(job_id, stage)
    job_manager.update_job_fields(
        job_id,
        coverage      = None,
        sample_name   = sample_name,
        analysis_mode = analysis_mode.value,
        end_type      = end_type.value,
        reference     = reference.value,
        adapter_kit   = adapter_kit,
        fastq_path    = fastq_path,
        fastq_r1      = fastq_r1_path,
        fastq_r2      = fastq_r2_path,
        bed_path      = bed_path,
    )

    # Jalankan hanya QC entry point
    asyncio.create_task(
        pipeline_runner.run_step(
            job_id      = job_id,
            step        = "qc",
            fastq_path  = fastq_path,
            fastq_r1    = fastq_r1_path,
            fastq_r2    = fastq_r2_path,
            end_type    = end_type,
            reference   = reference,
        )
    )

    return {"analysis_id": job_id, "step": "qc", "status": "running",
            "message": f"FastQC dimulai untuk sampel '{sample_name}'"}


@app.post("/api/analyses/{analysis_id}/trimming/run", tags=["Wizard"])
async def run_trimming(
    analysis_id: str,
    trim_preset:  Optional[TrimPreset] = Form(None),
    trim_quality: Optional[int]  = Form(None),
    trim_min_len: Optional[int]  = Form(None),
    trim_adapter: Optional[bool] = Form(None),
    trim_poly_g:  Optional[bool] = Form(None),
    adapter_kit:  Optional[str]  = Form(None),   # truseq_pe, nextera, dll (default: kit dari Step 1)
):
    """
    Step 2 — Trimming dengan parameter yang bisa diatur user.
    Hanya butuh upload FASTQ dari Step 1 (selesai saat /qc/run mengembalikan 200).
    Tidak perlu menunggu FastQC selesai; QC dan Trimming boleh jalan bersamaan.
    """
    _require_pipeline_ready()
    job = _get_job_or_404(analysis_id)
    _ensure_not_running(job, "trimming")

    # Cek step 1 sudah selesai (FASTQ sudah ada)
    if not os.path.exists(os.path.join(UPLOAD_DIR, analysis_id)):
        _error_response("STATUS_CONFLICT", "Upload belum selesai. Jalankan Step 1 (QC) dulu.", 409)

    # Resolve trim params
    if trim_preset == TrimPreset.CLINICAL:
        q, ml = 30, 36
    elif trim_preset == TrimPreset.RESEARCH:
        q, ml = 20, 36
    else:
        q  = trim_quality if trim_quality is not None else 20
        ml = trim_min_len if trim_min_len is not None else 36
    _validate_trim_params(q, ml)
    _validate_adapter_kit(adapter_kit)

    # Trimming baru → hasil Alignment dan Variant Calling sebelumnya tidak berlaku lagi
    for downstream in ("alignment", "variant_calling"):
        if job["stages"][downstream]["status"] == "running":
            _error_response("STATUS_CONFLICT", f"Step '{downstream}' sedang berjalan.", 409)
        job_manager.reset_stage(analysis_id, downstream)

    asyncio.create_task(
        pipeline_runner.run_step(
            job_id       = analysis_id,
            step         = "trimming",
            trim_quality = q,
            trim_min_len = ml,
            trim_adapter = trim_adapter if trim_adapter is not None else True,
            trim_poly_g  = trim_poly_g  if trim_poly_g  is not None else False,
            adapter_kit  = adapter_kit,
        )
    )

    return {"analysis_id": analysis_id, "step": "trimming", "status": "running"}


@app.post("/api/analyses/{analysis_id}/alignment/run", tags=["Wizard"])
async def run_alignment(
    analysis_id:     str,
    mark_duplicates: Optional[bool] = Form(None),
):
    """
    Step 3 — Alignment (BWA MEM + MarkDuplicates + Coverage). SE dan PE.
    Membutuhkan Step 2 (Trimming) sudah selesai.
    """
    _require_pipeline_ready()
    job = _get_job_or_404(analysis_id)
    _ensure_not_running(job, "alignment")

    trimming_status = job["stages"]["trimming"]["status"]
    if trimming_status != "completed":
        _error_response("STATUS_CONFLICT", "Trimming belum selesai. Jalankan Step 2 dulu.", 409)
    if not pipeline_runner.find_trimmed_reads(analysis_id):
        _error_response("NOT_FOUND", "File hasil Trimming tidak ditemukan di disk. Jalankan Step 2 ulang.", 404)

    # Alignment baru → Variant Calling sebelumnya tidak berlaku lagi
    _ensure_not_running(job, "variant_calling")
    job_manager.reset_stage(analysis_id, "variant_calling")

    asyncio.create_task(
        pipeline_runner.run_step(
            job_id          = analysis_id,
            step            = "alignment",
            mark_duplicates = mark_duplicates if mark_duplicates is not None else True,
        )
    )

    return {"analysis_id": analysis_id, "step": "alignment", "status": "running"}


@app.post("/api/analyses/{analysis_id}/variant-calling/run", tags=["Wizard"])
async def run_variant_calling(
    analysis_id: str,
    min_qual:    Optional[int] = Form(None),
    min_dp:      Optional[int] = Form(None),
    min_gq:      Optional[int] = Form(None),
    max_depth:   Optional[int] = Form(None),
    min_bq:      Optional[int] = Form(None),
    min_mq:      Optional[int] = Form(None),
    sample_sex:  Optional[str] = Form(None),   # male/female → ploidy chrX/chrY
):
    """
    Step 4 — Variant Calling (BCFTools mpileup → call → norm → filter).
    Membutuhkan Step 3 (Alignment) sudah selesai.
    """
    _require_pipeline_ready()
    job = _get_job_or_404(analysis_id)
    _ensure_not_running(job, "variant_calling")

    if job["stages"]["alignment"]["status"] != "completed":
        _error_response("STATUS_CONFLICT", "Alignment belum selesai. Jalankan Step 3 dulu.", 409)

    sex = _validate_sample_sex(sample_sex)
    p = {
        "min_qual":  min_qual  if min_qual  is not None else 20,
        "min_dp":    min_dp    if min_dp    is not None else 10,
        "min_gq":    min_gq    if min_gq    is not None else 20,
        "max_depth": max_depth if max_depth is not None else 500,
        "min_bq":    min_bq    if min_bq    is not None else 20,
        "min_mq":    min_mq    if min_mq    is not None else 20,
    }
    _validate_variant_params(**p)

    # Cari BAM file dari step sebelumnya
    sample = job["sample_name"]
    bam_path = os.path.join(RESULTS_DIR, analysis_id, sample, "bwa", f"{sample}.markdup.bam")
    if not os.path.exists(bam_path):
        _error_response("NOT_FOUND", "BAM file dari step Alignment tidak ditemukan.", 404)

    asyncio.create_task(
        pipeline_runner.run_step(
            job_id    = analysis_id,
            step      = "variant_calling",
            bam_file  = bam_path,
            sample_sex = sex,
            **p,
        )
    )

    return {"analysis_id": analysis_id, "step": "variant_calling", "status": "running"}


@app.post("/api/analyses/{analysis_id}/cancel", tags=["Wizard"])
async def cancel_analysis(analysis_id: str):
    """
    Batalkan step yang sedang antre/berjalan. Step antre langsung CANCELLED; step berjalan
    dihentikan dan menjadi CANCELLED setelah Nextflow berhenti (jobStatus CANCEL_REQUESTED).
    """
    job = _get_job_or_404(analysis_id)
    return {"analysis_id": analysis_id, "jobStatus": cancel_run(job)}


@app.get("/api/analyses/{analysis_id}/status", tags=["Wizard"])
async def get_analysis_status(analysis_id: str):
    """Status semua step untuk wizard UI."""
    job = job_manager.get_job(analysis_id)
    if not job:
        _error_response("NOT_FOUND", f"Analysis '{analysis_id}' tidak ditemukan.", 404)
    job_manager.update_elapsed(analysis_id)
    return job_manager.get_job(analysis_id)


# ═══════════════════════════════════════════════════════════════════════════════
# ENDPOINT: Health & Info
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/health", tags=["System"])
async def health_check(refresh: bool = False):
    """`?refresh=true` menjalankan ulang pemeriksaan tool (mis. setelah install bcftools baru)."""
    result = preflight.run_preflight() if refresh else preflight.last_result()
    return {
        "status":        "ok" if result.get("ok") else "degraded",
        "version":       "1.1.0",
        "backend":       "BCFTools",
        "pipelineMode":  PIPELINE_MODE,
        "pipeline":      "FastQC → Trimmomatic → BWA MEM → MarkDup → BCFTools mpileup/call/norm/filter",
        "missing":       preflight.missing_required(),
        "checks":        result,
        "time":          datetime.utcnow().isoformat() + "Z",
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
