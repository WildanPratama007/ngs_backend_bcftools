"""
NGS API — Job Manager
Mengelola state setiap pipeline job (in-memory + file JSON)
"""

import json
import os
import uuid
from datetime import datetime
from typing import Optional, Dict
from .models import PipelineStatus, StageStatus, AnalysisMode, EndType, ReferenceGenome, VariantType
from .config import JOBS_DIR


# In-memory store
_jobs: Dict[str, dict] = {}


def create_job(
    sample_name:   str,
    analysis_mode: AnalysisMode,
    end_type:      EndType,
    reference:     ReferenceGenome,
    variant_type:  VariantType,
) -> str:
    """Buat job baru dan return job_id."""
    job_id = str(uuid.uuid4())[:8] + "-" + sample_name[:16].replace(" ", "_")

    job = PipelineStatus(
        job_id        = job_id,
        sample_name   = sample_name,
        status        = StageStatus.PENDING,
        analysis_mode = analysis_mode,
        end_type      = end_type,
        reference     = reference,
        variant_type  = variant_type,
    )

    _jobs[job_id] = job.model_dump(mode="json")
    _save_job(job_id)
    return job_id


def get_job(job_id: str) -> Optional[dict]:
    """Ambil job by ID."""
    if job_id in _jobs:
        return _jobs[job_id]
    # Coba load dari file
    path = os.path.join(JOBS_DIR, f"{job_id}.json")
    if os.path.exists(path):
        with open(path) as f:
            job = json.load(f)
            _jobs[job_id] = job
            return job
    return None


def list_jobs() -> list:
    """List semua job."""
    # Load dari file jika belum di memory
    for fname in os.listdir(JOBS_DIR):
        if fname.endswith(".json"):
            jid = fname.replace(".json", "")
            if jid not in _jobs:
                path = os.path.join(JOBS_DIR, fname)
                with open(path) as f:
                    _jobs[jid] = json.load(f)
    return list(_jobs.values())


def update_job_status(job_id: str, status: StageStatus):
    """Update status keseluruhan job."""
    if job_id in _jobs:
        _jobs[job_id]["status"] = status.value
        _save_job(job_id)


def update_stage(job_id: str, stage: str, status: StageStatus, progress: int = 0, error: str = None):
    """Update status tahap tertentu."""
    if job_id not in _jobs:
        return
    now = datetime.now().isoformat()
    _jobs[job_id]["stages"][stage]["status"]   = status.value
    _jobs[job_id]["stages"][stage]["progress"] = progress
    if status == StageStatus.RUNNING and not _jobs[job_id]["stages"][stage].get("started_at"):
        _jobs[job_id]["stages"][stage]["started_at"] = now
    if status in (StageStatus.COMPLETED, StageStatus.FAILED):
        _jobs[job_id]["stages"][stage]["finished_at"] = now
    if error:
        _jobs[job_id]["stages"][stage]["error"] = error
    _save_job(job_id)


def update_elapsed(job_id: str):
    """Hitung elapsed time sejak job dibuat. Stop jika sudah completed/failed."""
    if job_id not in _jobs:
        return
    # Jangan update jika sudah selesai
    final_statuses = ("completed", "failed")
    if _jobs[job_id].get("status") in final_statuses:
        return
    created = datetime.fromisoformat(_jobs[job_id]["created_at"])
    delta   = datetime.now() - created
    hours, remainder = divmod(int(delta.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)
    _jobs[job_id]["elapsed"] = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    _save_job(job_id)


def update_coverage(job_id: str, coverage: float):
    """Update mean coverage setelah alignment selesai."""
    if job_id not in _jobs:
        return
    _jobs[job_id]["coverage"] = round(coverage, 2)
    _save_job(job_id)


def delete_job(job_id: str) -> bool:
    """Hapus job dari memory dan disk."""
    if job_id in _jobs:
        del _jobs[job_id]
    path = os.path.join(JOBS_DIR, f"{job_id}.json")
    if os.path.exists(path):
        os.remove(path)
        return True
    return False


def _save_job(job_id: str):
    """Simpan job ke file JSON."""
    path = os.path.join(JOBS_DIR, f"{job_id}.json")
    with open(path, "w") as f:
        json.dump(_jobs[job_id], f, indent=2, default=str)
