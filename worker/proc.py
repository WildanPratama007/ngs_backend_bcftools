"""Menjalankan proses anak (Nextflow/FastQC) dengan penanganan SIGTERM dari ECS."""

import os
import signal
import subprocess
from typing import Callable, Optional

from .errors import EXIT_TERMINATED, TERMINATED, WorkerError

_current: Optional[subprocess.Popen] = None


def _on_term(signum, frame):
    # ECS mengirim SIGTERM saat task dihentikan (stop, deploy, atau interupsi Spot), lalu SIGKILL setelah stopTimeout.
    if _current is not None and _current.poll() is None:
        try:
            os.killpg(os.getpgid(_current.pid), signal.SIGTERM)
        except ProcessLookupError:
            pass
    raise WorkerError(TERMINATED, "Task dihentikan (SIGTERM) sebelum selesai", EXIT_TERMINATED)


def install_signal_handlers():
    signal.signal(signal.SIGTERM, _on_term)
    signal.signal(signal.SIGINT, _on_term)


def run(cmd: list, env: dict, cwd: str, on_line: Callable[[str], None]) -> int:
    """Jalankan perintah, kirim tiap baris output (stdout+stderr) ke on_line, kembalikan return code."""
    global _current
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace",
                            cwd=cwd, env=env, start_new_session=True)
    _current = proc
    try:
        for line in proc.stdout:
            on_line(line.rstrip("\n"))
        return proc.wait()
    finally:
        if proc.poll() is None:   # keluar karena exception: jangan tinggalkan proses yatim
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass
            proc.wait()
        _current = None
