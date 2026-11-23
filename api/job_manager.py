"""
NGS API — Job Manager
Fix:
- Pemulihan state saat restart (load jobs dari file JSON)
- Hapus job: stop Nextflow + hapus folder uploads/results
- State running saat restart di-mark failed
"""

import json
import os
import re
import shutil
import signal
import uuid
from datetime import datetime
from typing import Optional, Dict
from .models import PipelineStatus, StageStatus, AnalysisMode, EndType, ReferenceGenome, VariantType
from .config import JOBS_DIR, UPLOAD_DIR, RESULTS_DIR

# In-memory store
_jobs: Dict[str, dict] = {}

# Track PID Nextflow per job per step (QC dan Trimming bisa jalan bersamaan)
_nextflow_pids: Dict[str, Dict[str, int]] = {}

# Step yang sedang antre/berjalan per job, dan job yang diminta batal
_active_steps: Dict[str, int] = {}
_cancel_requested: set = set()

# ID dipakai sebagai nama folder di uploads/ dan results/, jadi harus aman sebagai satu komponen path:
# diawali huruf/angka, tanpa "/", spasi, atau "..".
_ID_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$')


def is_valid_id(job_id) -> bool:
    return isinstance(job_id, str) and bool(_ID_RE.match(job_id)) and ".." not in job_id


def begin_step(job_id: str):
    _active_steps[job_id] = _active_steps.get(job_id, 0) + 1


def end_step(job_id: str):
    """Dipanggil saat sebuah step selesai/dibatalkan. Flag batal dibersihkan setelah semua step berhenti."""
    n = _active_steps.get(job_id, 0) - 1
    if n > 0:
        _active_steps[job_id] = n
        return
    _active_steps.pop(job_id, None)
    _cancel_requested.discard(job_id)
    if job_id in _jobs and _jobs[job_id].get("cancel_requested"):
        _jobs[job_id]["cancel_requested"] = False
        _save_job(job_id)


def is_cancel_requested(job_id: str) -> bool:
    return job_id in _cancel_requested


def request_cancel(job_id: str) -> str:
    """
    Batalkan semua step job ini. Step yang masih antre langsung CANCELLED; step yang berjalan
    dihentikan (SIGTERM ke Nextflow) dan berakhir CANCELLED setelah proses selesai.
    Mengembalikan "CANCELLED" (hanya yang antre) atau "CANCEL_REQUESTED" (ada proses berjalan).
    """
    _cancel_requested.add(job_id)
    pids = _nextflow_pids.get(job_id, {})
    job = _jobs.get(job_id)
    if job:
        for st in job["stages"].values():
            if st["status"] == StageStatus.PENDING.value:
                st["status"] = StageStatus.CANCELLED.value
                st["finished_at"] = datetime.now().isoformat()
        job["cancel_requested"] = bool(pids)
        _save_job(job_id)
    for pid in list(pids.values()):
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    if _active_steps.get(job_id, 0) == 0:   # tidak ada yang antre/berjalan → tidak perlu flag
        _cancel_requested.discard(job_id)
    return "CANCEL_REQUESTED" if pids else "CANCELLED"


def register_nextflow_pid(job_id: str, pid: int, step: str = "full"):
    """Simpan PID proses Nextflow untuk job + step ini."""
    _nextflow_pids.setdefault(job_id, {})[step] = pid


def unregister_nextflow_pid(job_id: str, step: str = "full"):
    _nextflow_pids.get(job_id, {}).pop(step, None)


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
            # Stage yang running/pending saat server mati tidak mungkin dilanjutkan → failed.
            # (job-level "pending" sengaja tidak disentuh: di wizard artinya "sebagian step selesai")
            changed = False
            for stage in job.get("stages", {}).values():
                if stage.get("status") in ("running", "pending"):
                    stage["status"] = "failed"
                    stage["error"] = stage.get("error") or "Server restart — step interrupted"
                    stage["finished_at"] = datetime.now().isoformat()
                    changed = True
            if job.get("cancel_requested"):
                job["cancel_requested"] = False
                changed = True
            if job.get("status") == "running" or changed:
                job["status"] = "failed"
                job["nextflow_log"] = job.get("nextflow_log") or "Server restart — job interrupted"
                changed = True
            if changed:
                with open(path, "w") as f:
                    json.dump(job, f, indent=2, default=str)
            _jobs[jid] = job
        except Exception:
            pass


def create_job_with_id(
    job_id:        str,
    sample_name:   str,
    analysis_mode: AnalysisMode,
    end_type:      EndType,
    reference:     ReferenceGenome,
    variant_type:  VariantType,
) -> str:
    """Buat job dengan ID yang sudah ditentukan (dari frontend analysis_id)."""
    if not is_valid_id(job_id):
        raise ValueError(f"ID tidak valid: {job_id!r}")
    job = PipelineStatus(
        job_id        = job_id,
        sample_name   = sample_name,
        status        = StageStatus.NOT_STARTED,
        analysis_mode = analysis_mode,
        end_type      = end_type,
        reference     = reference,
        variant_type  = variant_type,
    )
    _jobs[job_id] = job.model_dump(mode="json")
    # Set semua stage ke not_started
    for stage in _jobs[job_id]["stages"].values():
        stage["status"] = StageStatus.NOT_STARTED.value
    _save_job(job_id)
    return job_id


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
    if not is_valid_id(job_id):   # jangan pernah membentuk path dari ID mentah
        return None
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
    if status in (StageStatus.COMPLETED, StageStatus.FAILED, StageStatus.CANCELLED):
        _jobs[job_id]["stages"][stage]["finished_at"] = now
    if error:
        _jobs[job_id]["stages"][stage]["error"] = error
    _save_job(job_id)


def update_job_fields(job_id: str, **fields):
    """Simpan metadata tambahan di job (mis. adapter_kit, path upload)."""
    if job_id not in _jobs:
        return
    _jobs[job_id].update(fields)
    _save_job(job_id)


def reset_stage(job_id: str, stage: str):
    """Bersihkan state stage sebelum dijalankan ulang (retry)."""
    if job_id not in _jobs:
        return
    _jobs[job_id]["stages"][stage].update(
        {"status": StageStatus.NOT_STARTED.value, "progress": 0,
         "started_at": None, "finished_at": None, "error": None}
    )
    _save_job(job_id)


def set_nextflow_log(job_id: str, message: str):
    if job_id in _jobs:
        _jobs[job_id]["nextflow_log"] = message
        _save_job(job_id)


def refresh_job_status(job_id: str):
    """Turunkan status level-job dari status stage (dipakai oleh wizard)."""
    if job_id not in _jobs:
        return
    statuses = [s["status"] for s in _jobs[job_id]["stages"].values()]
    if "running" in statuses:
        new = StageStatus.RUNNING
    elif "failed" in statuses:
        new = StageStatus.FAILED
    elif "cancelled" in statuses:
        new = StageStatus.CANCELLED
    elif all(s == "completed" for s in statuses):
        new = StageStatus.COMPLETED
    elif any(s == "completed" for s in statuses):
        new = StageStatus.PENDING
    else:
        new = StageStatus.NOT_STARTED
    update_job_status(job_id, new)


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
    # ID yang bukan job dikenal dan bukan ID aman tidak boleh dipakai membentuk path (rmtree di bawah)
    if job_id not in _jobs and not is_valid_id(job_id):
        return False

    # Stop Nextflow process bila masih jalan
    for pid in _nextflow_pids.pop(job_id, {}).values():
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
