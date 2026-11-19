"""
NGS API — Job Manager
Fix:
- Pemulihan state saat restart (load jobs dari file JSON)
- Hapus job: stop Nextflow + hapus folder uploads/results
- State running saat restart di-mark failed
"""

import json
import os
import shutil
import signal
import uuid
from datetime import datetime
from typing import Optional, Dict
from .models import PipelineStatus, StageStatus, AnalysisMode, EndType, ReferenceGenome, VariantType
from .config import JOBS_DIR, UPLOAD_DIR, RESULTS_DIR

# In-memory store
_jobs: Dict[str, dict] = {}

# Track PID Nextflow per job untuk cancel
_nextflow_pids: Dict[str, int] = {}


def register_nextflow_pid(job_id: str, pid: int):
    """Simpan PID proses Nextflow untuk job ini."""
    _nextflow_pids[job_id] = pid


def recover_jobs():
    """
    Panggil saat startup — load semua job dari disk.
    Job yang berstatus 'running' saat server mati di-mark failed.
    """
    if not os.path.exists(JOBS_DIR):
        return
    for fname in os.listdir(JOBS_DIR):
        if not fname.endswith(".json"):
            continue
        path = os.path.join(JOBS_DIR, fname)
        try:
            with open(path) as f:
                job = json.load(f)
            jid = job.get("job_id", fname.replace(".json", ""))
            # Job running saat server mati → mark failed
            if job.get("status") == "running":
                job["status"] = "failed"
                job["nextflow_log"] = "Server restart — job interrupted"
                for stage in job.get("stages", {}).values():
                    if stage.get("status") == "running":
                        stage["status"] = "failed"
                        stage["finished_at"] = datetime.now().isoformat()
                with open(path, "w") as f:
                    json.dump(job, f, indent=2, default=str)
            _jobs[jid] = job
        except Exception:
            pass


def create_job(
    sample_name:   str,
    analysis_mode: AnalysisMode,
    end_type:      EndType,
    reference:     ReferenceGenome,
    variant_type:  VariantType,
) -> str:
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
    if job_id in _jobs:
        return _jobs[job_id]
    path = os.path.join(JOBS_DIR, f"{job_id}.json")
    if os.path.exists(path):
        with open(path) as f:
            job = json.load(f)
            _jobs[job_id] = job
            return job
    return None


def list_jobs() -> list:
    for fname in os.listdir(JOBS_DIR):
        if fname.endswith(".json"):
            jid = fname.replace(".json", "")
            if jid not in _jobs:
                try:
                    with open(os.path.join(JOBS_DIR, fname)) as f:
                        _jobs[jid] = json.load(f)
                except Exception:
                    pass
    return list(_jobs.values())


def update_job_status(job_id: str, status: StageStatus):
    if job_id in _jobs:
        _jobs[job_id]["status"] = status.value
        _save_job(job_id)


def update_stage(job_id: str, stage: str, status: StageStatus, progress: int = 0, error: str = None):
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
    if job_id not in _jobs:
        return
    if _jobs[job_id].get("status") in ("completed", "failed"):
        return
    created = datetime.fromisoformat(_jobs[job_id]["created_at"])
    delta   = datetime.now() - created
    hours, remainder = divmod(int(delta.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)
    _jobs[job_id]["elapsed"] = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    _save_job(job_id)


def update_coverage(job_id: str, coverage: dict):
    """Update extended coverage stats setelah alignment selesai."""
    if job_id not in _jobs:
        return
    _jobs[job_id]["coverage"] = coverage
    _save_job(job_id)


def delete_job(job_id: str) -> bool:
    """
    Hapus job: stop Nextflow process + hapus folder uploads & results + hapus JSON.
    """
    # Stop Nextflow process bila masih jalan
    pid = _nextflow_pids.pop(job_id, None)
    if pid:
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass

    # Hapus dari memory
    if job_id in _jobs:
        del _jobs[job_id]

    # Hapus JSON
    json_path = os.path.join(JOBS_DIR, f"{job_id}.json")
    found = os.path.exists(json_path)
    if found:
        os.remove(json_path)

    # Hapus folder uploads
    upload_path = os.path.join(UPLOAD_DIR, job_id)
    if os.path.exists(upload_path):
        shutil.rmtree(upload_path, ignore_errors=True)

    # Hapus folder results
    result_path = os.path.join(RESULTS_DIR, job_id)
    if os.path.exists(result_path):
        shutil.rmtree(result_path, ignore_errors=True)

    return found


def _save_job(job_id: str):
    path = os.path.join(JOBS_DIR, f"{job_id}.json")
    with open(path, "w") as f:
        json.dump(_jobs[job_id], f, indent=2, default=str)
