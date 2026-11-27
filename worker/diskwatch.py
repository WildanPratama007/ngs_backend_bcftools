"""Mengukur puncak pemakaian disk scratch selama job (selisih ruang bebas filesystem), untuk menyetel DISK_FACTOR."""

import shutil
import threading
from pathlib import Path


class DiskPeak:
    def __init__(self, path: Path, interval: float = 10.0):
        self.path = path
        self.interval = interval
        self.start_free = shutil.disk_usage(path).free
        self.min_free = self.start_free
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True)

    def sample(self):
        self._sample()

    def _sample(self):
        try:
            self.min_free = min(self.min_free, shutil.disk_usage(self.path).free)
        except OSError:
            pass

    def _loop(self):
        while not self._stop.wait(self.interval):
            self._sample()

    def start(self):
        self._thread.start()
        return self

    def stop(self) -> dict:
        self._stop.set()
        self._sample()
        return {"freeAtStartBytes": self.start_free, "peakUsedBytes": max(0, self.start_free - self.min_free)}
