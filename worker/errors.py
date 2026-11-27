"""Kode error worker. Kode ini muncul di status.json (errorCode) dan menentukan exit code proses."""

from typing import Optional

# Exit code proses
EXIT_OK = 0
EXIT_FAILED = 1        # pipeline / I/O gagal
EXIT_INVALID_SPEC = 2  # spesifikasi job tidak valid (tidak ada yang dijalankan)
EXIT_TERMINATED = 143  # dihentikan (SIGTERM, mis. ECS stop atau Spot)

INVALID_SPEC = "INVALID_SPEC"
INPUT_NOT_FOUND = "INPUT_NOT_FOUND"
DOWNLOAD_FAILED = "DOWNLOAD_FAILED"
DISK_INSUFFICIENT = "DISK_INSUFFICIENT"
REFERENCE_NOT_FOUND = "REFERENCE_NOT_FOUND"
PIPELINE_FAILED = "PIPELINE_FAILED"
NO_OUTPUT = "NO_OUTPUT"
UPLOAD_FAILED = "UPLOAD_FAILED"
TERMINATED = "TERMINATED"
INTERNAL_ERROR = "INTERNAL_ERROR"


class WorkerError(Exception):
    def __init__(self, code: str, message: str, exit_code: int = EXIT_FAILED, details: Optional[dict] = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.exit_code = exit_code
        self.details = details or {}
