"""
NGS API — Pemeriksaan lingkungan saat startup.

Versi tool yang salah (mis. bcftools < 1.21) baru ketahuan saat step gagal di tengah jalan.
Di sini semuanya dicek di awal, dan hasilnya dipakai /api/health serta guard submit.
"""

import logging
import os
import re
import subprocess

from .config import (
    NEXTFLOW_BIN, NEXTFLOW_JAVA_HOME, FASTQC_PATH, BIN_DIR, SAMTOOLS_BIN, BCFTOOLS_BIN,
    REF_GENOMES, ADAPTER_KIT_MAP,
)

logger = logging.getLogger("ngs.preflight")

MIN_BCFTOOLS = (1, 21)   # concat --rm-dups butuh -a dan versi baru
MIN_SAMTOOLS = (1, 13)

_last: dict = {"ok": False, "checks": {}}


def _version(cmd: list) -> str:
    try:
        # samtools --version memuat karakter non-UTF8 (©) di baris berikutnya
        out = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=15)
        return (out.stdout or out.stderr).strip().splitlines()[0] if (out.stdout or out.stderr) else ""
    except Exception as e:  # tool tidak ada / tidak bisa dijalankan
        return f"ERROR: {e}"


def _parse(text: str):
    m = re.search(r"(\d+)\.(\d+)", text)
    return (int(m.group(1)), int(m.group(2))) if m else None


def _check_min_version(name: str, binary: str, minimum: tuple) -> dict:
    if not os.path.isfile(binary):
        return {"ok": False, "detail": f"{binary} tidak ditemukan"}
    line = _version([binary, "--version"])
    ver  = _parse(line)
    if not ver:
        return {"ok": False, "detail": f"versi {name} tidak terbaca: {line}"}
    need = ".".join(map(str, minimum))
    return {"ok": ver >= minimum, "detail": f"{line} (minimal {need})"}


def _check_file(path: str, what: str) -> dict:
    ok = bool(path) and os.path.exists(path)
    return {"ok": ok, "detail": path if ok else f"{what} tidak ditemukan: {path}"}


def run_preflight() -> dict:
    required = {
        "bcftools":    _check_min_version("bcftools", BCFTOOLS_BIN, MIN_BCFTOOLS),
        "samtools":    _check_min_version("samtools", SAMTOOLS_BIN, MIN_SAMTOOLS),
        "nextflow":    _check_file(NEXTFLOW_BIN, "nextflow"),
        "java":        _check_file(os.path.join(NEXTFLOW_JAVA_HOME or "", "bin", "java"), "java 17"),
        "fastqc":      _check_file(FASTQC_PATH, "fastqc"),
        "trimmomatic": _check_file(os.path.join(BIN_DIR, "trimmomatic"), "trimmomatic"),
        "bwa":         _check_file(os.path.join(BIN_DIR, "bwa"), "bwa"),
    }
    optional = {
        f"reference:{name}": _check_file(path, f"referensi {name}") for name, path in REF_GENOMES.items()
    }
    for kit, path in ADAPTER_KIT_MAP.items():
        optional[f"adapter:{kit}"] = _check_file(path, f"adapter {kit}")

    result = {
        "ok":       all(c["ok"] for c in required.values()),
        "required": required,
        "optional": optional,
    }
    _last.clear()
    _last.update(result)

    for name, c in {**required, **optional}.items():
        if not c["ok"]:
            level = logger.error if name in required else logger.warning
            level("Preflight %s: %s", name, c["detail"])
    return result


def last_result() -> dict:
    return _last


def missing_required() -> list:
    return [f"{n}: {c['detail']}" for n, c in _last.get("required", {}).items() if not c["ok"]]
