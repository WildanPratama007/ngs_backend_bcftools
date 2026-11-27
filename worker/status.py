"""
status.json: ditulis ke <output.prefix>status.json saat mulai, saat tahap berganti, dan di akhir
(termasuk saat gagal), supaya EventBridge/Lambda/DynamoDB bisa membacanya tanpa memparsing log.
"""

import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

from .s3io import S3

log = logging.getLogger("worker.status")

SCHEMA_VERSION = 1


def _now() -> str:
    dt = datetime.now(timezone.utc)
    return dt.strftime("%Y-%m-%dT%H:%M:%S.") + f"{dt.microsecond // 1000:03d}Z"


class StatusReporter:
    def __init__(self, task: str, job_id: str, prefix_uri: Optional[str], s3: Optional[S3], local_dir: Path):
        self.s3 = s3
        self.uri = f"{prefix_uri}status.json" if prefix_uri else None
        self.local = local_dir / "status.json"
        self._last_push = 0.0
        self.doc = {
            "schemaVersion": SCHEMA_VERSION,
            "task": task,
            "jobId": job_id,
            "status": "RUNNING",
            "stage": "validate_input",
            "progress": 0,
            "errorCode": None,
            "errorMessage": None,
            "details": None,
            "startedAt": _now(),
            "updatedAt": _now(),
            "finishedAt": None,
            "outputs": [],
            "metrics": None,
        }

    def _write(self, force: bool = False):
        self.doc["updatedAt"] = _now()
        try:
            self.local.parent.mkdir(parents=True, exist_ok=True)
            self.local.write_text(json.dumps(self.doc, indent=2))
        except OSError:
            pass
        # Dibatasi (maks 1 tulis / 5 detik) kecuali dipaksa; kegagalan menulis status tidak boleh mematikan job
        if self.s3 and self.uri and (force or time.time() - self._last_push >= 5):
            try:
                self.s3.put_json(self.uri, self.doc)
                self._last_push = time.time()
            except Exception as e:  # noqa: BLE001
                log.warning("status.json gagal ditulis ke S3: %s", e)

    def start(self):
        log.info("status RUNNING stage=%s", self.doc["stage"])
        self._write(force=True)

    def stage(self, name: str, progress: Optional[int] = None):
        changed = name != self.doc["stage"]
        self.doc["stage"] = name
        if progress is not None and progress >= self.doc["progress"]:
            self.doc["progress"] = min(progress, 99)
        if changed:
            log.info("stage=%s progress=%s", name, self.doc["progress"])
        self._write(force=changed)

    def succeed(self, outputs: List[dict], metrics: Optional[dict] = None):
        self.doc.update(status="SUCCEEDED", stage="finalize", progress=100, outputs=outputs,
                        metrics=metrics, finishedAt=_now())
        self._write(force=True)
        log.info("status SUCCEEDED (%d file output)", len(outputs))

    def fail(self, code: str, message: str, details: Optional[dict] = None):
        self.doc.update(status="FAILED", errorCode=code, errorMessage=message[:4000],
                        details=details or None, finishedAt=_now())
        self._write(force=True)
        log.error("status FAILED code=%s: %s", code, message)
