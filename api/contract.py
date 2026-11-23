"""
NGS API — Endpoint sesuai API_Contract Secondary Analysis (bagian yang tidak butuh AWS).

Prefix: /api/secondary-analysis
  GET  /analyses/{analysisId}/status   → getStatus   (Status_View, bagian 3.7)
  POST /status/batch                   → getStatusBatch (bagian 3.8)
  POST /analyses/{analysisId}/cancel   → cancel (bagian 3.6)

Yang sengaja tidak ada di sini karena bergantung pada AWS/layanan lain (lihat README):
registerFiles/startUpload/signParts/completeUpload/cancelUpload (S3 multipart), submit/retry
(memakai FASTQ dari S3 + Job_Message), autentikasi Cognito, role, institusi, langganan, kuota.
"""

from datetime import datetime, timezone
from typing import List, Optional

from fastapi import APIRouter
from pydantic import BaseModel, ConfigDict

from . import job_manager
from .config import PIPELINE_MODE
from .errors import error_response

router = APIRouter(prefix="/api/secondary-analysis", tags=["Contract"])

STAGES = ("fastqc", "trimming", "alignment", "variant_calling")
TERMINAL = ("SUCCEEDED", "FAILED", "CANCELLED")

# Pemetaan JobStatus → AnalysisStatus (kontrak bagian 2.5)
_ANALYSIS_STATUS = {
    "QUEUED": "PROCESSING", "RUNNING": "PROCESSING", "CANCEL_REQUESTED": "PROCESSING",
    "SUCCEEDED": "QC_REVIEW", "FAILED": "FAILED", "CANCELLED": "CANCELLED",
}


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")   # field tak dikenal → VALIDATION_FAILED


class StatusBatchRequest(_Strict):
    analysisIds: List[str]


class CancelRequest(_Strict):
    jobId: str


def _iso(ts: Optional[str]) -> Optional[str]:
    """Waktu job disimpan sebagai waktu lokal tanpa zona; kontrak meminta ISO 8601 UTC dengan milidetik."""
    if not ts:
        return None
    dt = datetime.fromisoformat(ts)
    dt = dt.astimezone(timezone.utc)   # naive → dianggap waktu lokal server
    return dt.strftime("%Y-%m-%dT%H:%M:%S.") + f"{dt.microsecond // 1000:03d}Z"


def derive_job_status(job: dict) -> Optional[str]:
    """
    JobStatus kontrak dari status stage. None = belum pernah dijalankan (run: null).
    Di mode wizard, jeda antar step (sebagian stage selesai, tidak ada yang jalan) dilaporkan RUNNING
    karena Run belum selesai.
    """
    vals = [job["stages"][s]["status"] for s in STAGES]
    if all(v == "not_started" for v in vals):
        return None
    active = [v for v in vals if v in ("running", "pending")]
    if job.get("cancel_requested") and active:
        return "CANCEL_REQUESTED"
    if "running" in vals:
        return "RUNNING"
    if "pending" in vals:
        return "QUEUED"
    if "failed" in vals:
        return "FAILED"
    if "cancelled" in vals:
        return "CANCELLED"
    if all(v == "completed" for v in vals):
        return "SUCCEEDED"
    return "RUNNING"


def build_status_view(job: dict) -> dict:
    job_status = derive_job_status(job)
    if job_status is None:
        return {"analysisStatus": "READY", "pipelineMode": PIPELINE_MODE, "run": None}

    stages   = [(s, job["stages"][s]) for s in STAGES]
    progress = int(sum(st.get("progress") or 0 for _, st in stages) / len(stages))
    current  = next((n for n, st in reversed(stages) if st["status"] == "running"), None) \
            or next((n for n, st in stages if st["status"] == "pending"), None) \
            or next((n for n, st in stages if st["status"] in ("failed", "cancelled")), None) \
            or next((n for n, st in reversed(stages) if st["status"] != "not_started"), STAGES[0])

    started  = sorted(st["started_at"] for _, st in stages if st.get("started_at"))
    finished = sorted(st["finished_at"] for _, st in stages if st.get("finished_at"))
    terminal = job_status in TERMINAL

    error_msg = None
    if job_status == "FAILED":
        failed = next((st for _, st in stages if st["status"] == "failed"), None)
        error_msg = (failed.get("error") or "Pipeline gagal")[:1000] if failed else "Pipeline gagal"

    return {
        "analysisStatus": _ANALYSIS_STATUS.get(job_status, "PROCESSING"),
        "pipelineMode":   PIPELINE_MODE,
        "run": {
            "jobId":        job["job_id"],
            "jobStatus":    job_status,
            "stage":        current,
            "progress":     100 if job_status == "SUCCEEDED" else min(progress, 99 if not terminal else 100),
            "message":      f"{current}: {progress}%",
            "errorMessage": error_msg,
            "submittedAt":  _iso(job.get("created_at")),
            "startedAt":    _iso(started[0]) if started else None,
            "finishedAt":   _iso(finished[-1]) if terminal and finished else None,
            "mode":         PIPELINE_MODE if PIPELINE_MODE in ("mock", "aws") else "mock",
        },
    }


def _not_found():
    # Isi respons sama untuk id tidak ada, id tidak valid, atau milik pihak lain (kontrak bagian 3)
    error_response("NOT_FOUND", "Analisis tidak ditemukan.", 404)


def cancel_run(job: dict) -> str:
    """Batalkan Run job. Mengembalikan jobStatus hasil ("CANCELLED" atau "CANCEL_REQUESTED")."""
    status = derive_job_status(job)
    if status is None or status in TERMINAL or status == "CANCEL_REQUESTED":
        error_response("STATUS_CONFLICT", "Run tidak bisa dibatalkan pada status ini.", 409,
                       {"currentStatus": build_status_view(job)["analysisStatus"]})
    return job_manager.request_cancel(job["job_id"])


@router.get("/analyses/{analysis_id}/status")
async def get_status(analysis_id: str):
    job = job_manager.get_job(analysis_id)
    if not job:
        _not_found()
    return build_status_view(job)


@router.post("/status/batch")
async def get_status_batch(req: StatusBatchRequest):
    ids = req.analysisIds
    if not (1 <= len(ids) <= 25):
        error_response("VALIDATION_FAILED", "analysisIds harus 1 sampai 25 id.", 400,
                       {"fields": [{"path": "analysisIds", "reason": "jumlah harus 1 sampai 25"}]})
    if len(set(ids)) != len(ids):
        error_response("VALIDATION_FAILED", "analysisIds tidak boleh duplikat.", 400,
                       {"fields": [{"path": "analysisIds", "reason": "id duplikat"}]})

    items = []
    for aid in ids:
        job = job_manager.get_job(aid)
        if job:
            items.append({"analysisId": aid, "found": True, "status": build_status_view(job)})
        else:
            items.append({"analysisId": aid, "found": False})
    return {"items": items}


@router.post("/analyses/{analysis_id}/cancel")
async def cancel(analysis_id: str, req: CancelRequest):
    job = job_manager.get_job(analysis_id)
    if not job:
        _not_found()
    if req.jobId != job["job_id"]:
        error_response("STATUS_CONFLICT", "jobId bukan Run terakhir.", 409,
                       {"currentStatus": build_status_view(job)["analysisStatus"]})
    return {"jobStatus": cancel_run(job)}
