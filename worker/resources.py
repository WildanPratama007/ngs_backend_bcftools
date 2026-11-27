"""Deteksi CPU dan memori yang benar-benar tersedia untuk container (batas cgroup, bukan host)."""

import json
import os
import urllib.request
from pathlib import Path


def _read(path: str):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def _ecs_task_limits() -> dict:
    """
    Batas CPU (vCPU) dan memori (MiB) task dari endpoint metadata ECS. Di Fargate batas task tidak terlihat dari
    cgroup container (yang terbaca adalah memori mesin induk), jadi sumber ini yang dipercaya bila ada.
    """
    base = os.environ.get("ECS_CONTAINER_METADATA_URI_V4")
    if not base:
        return {}
    try:
        with urllib.request.urlopen(base + "/task", timeout=3) as resp:
            limits = json.load(resp).get("Limits") or {}
        return {k: v for k, v in limits.items() if isinstance(v, (int, float)) and v > 0}
    except Exception:  # noqa: BLE001 — metadata tidak tersedia: pakai cgroup
        return {}


def cpu_limit() -> int:
    meta_cpu = _ecs_task_limits().get("CPU")
    if meta_cpu:
        return max(1, int(meta_cpu))
    # cgroup v2: "max 100000" atau "400000 100000"
    cpu_max = _read("/sys/fs/cgroup/cpu.max")
    if cpu_max:
        quota, _, period = cpu_max.partition(" ")
        if quota != "max" and period.isdigit() and int(period) > 0:
            return max(1, int(int(quota) / int(period)))
    # cgroup v1
    quota, period = _read("/sys/fs/cgroup/cpu/cpu.cfs_quota_us"), _read("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    if quota and period and quota.lstrip("-").isdigit() and int(quota) > 0 and period.isdigit():
        return max(1, int(int(quota) / int(period)))
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def memory_limit_bytes() -> int:
    candidates = []
    meta_mem = _ecs_task_limits().get("Memory")   # MiB
    if meta_mem:
        candidates.append(int(meta_mem) << 20)
    v2 = _read("/sys/fs/cgroup/memory.max")
    if v2 and v2 != "max" and v2.isdigit():
        candidates.append(int(v2))
    v1 = _read("/sys/fs/cgroup/memory/memory.limit_in_bytes")
    if v1 and v1.isdigit():
        candidates.append(int(v1))
    meminfo = _read("/proc/meminfo")
    if meminfo:
        for line in meminfo.splitlines():
            if line.startswith("MemTotal:"):
                candidates.append(int(line.split()[1]) * 1024)
                break
    sane = [c for c in candidates if c < (1 << 50)]   # limit "tanpa batas" cgroup v1 bernilai sangat besar
    return min(sane) if sane else 4 << 30


def plan(threads_override=None) -> dict:
    """
    Rencana sumber daya untuk Nextflow:
      cpus     : CPU container
      memory   : memori yang boleh dipakai proses (sisakan untuk OS, JVM Nextflow, dan I/O)
      threads  : thread per alat (BWA/samtools/Trimmomatic), maks 16
      sort_mem : memori per thread untuk `samtools sort`; threads x sort_mem <= ~40% memori
    """
    cpus = int(os.environ.get("NGS_MAX_CPUS") or cpu_limit())
    mem_total = memory_limit_bytes()
    headroom = max(1 << 30, int(mem_total * 0.08))
    mem_usable = max(512 << 20, mem_total - headroom)
    if os.environ.get("NGS_MAX_MEMORY_MB"):
        mem_usable = int(os.environ["NGS_MAX_MEMORY_MB"]) << 20

    threads = int(threads_override or os.environ.get("THREADS") or min(cpus, 16))
    threads = max(1, min(threads, cpus))
    sort_mb = int(mem_usable * 0.4 / threads / (1 << 20))
    sort_mb = max(256, min(sort_mb, 4096))
    return {
        "cpus": cpus,
        "memory_mb": mem_usable >> 20,
        "memory_total_mb": mem_total >> 20,
        "threads": threads,
        "sort_mem": f"{sort_mb}M",
    }
