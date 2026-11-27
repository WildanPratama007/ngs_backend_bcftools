"""`worker vcf`: FASTQ (S3) -> pipeline penuh (pipeline/main.nf) -> VCF, BAM, QC (S3)."""

import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from . import config, proc, resources
from .diskwatch import DiskPeak
from .nfcommon import Monitor, ensure_reference, nextflow_env
from .errors import (DISK_INSUFFICIENT, NO_OUTPUT, PIPELINE_FAILED, REFERENCE_NOT_FOUND, WorkerError)
from .s3io import S3
from .spec import JobSpec, Layout, Reference, parse_s3_uri
from .status import StatusReporter

log = logging.getLogger("worker.vcf")
def _build_command(spec: JobSpec, plan: dict, ref_fa: str, local: Dict[str, str], out_dir: Path, work_dir: Path) -> list:
    q, minlen = spec.trimming.effective()
    vc = spec.variantCalling
    cmd = [
        config.NEXTFLOW_BIN, "-log", str(out_dir / ".nextflow.log"),
        "run", str(Path(config.PIPELINE_DIR) / "main.nf"),
        "-profile", config.NEXTFLOW_PROFILE,
        "-w", str(work_dir),
        "--outdir", str(out_dir),
        "--ref_genome", ref_fa,
        "--sample_name", spec.sample,
        "--ploidy", "GRCh37" if spec.reference == Reference.GRCh37 else "GRCh38",
        "--threads", str(plan["threads"]),
        "--sort_mem", plan["sort_mem"],
        "--bin_dir", config.BIN_DIR,
        "--samtools_bin", config.SAMTOOLS_BIN,
        "--bcftools_bin", config.BCFTOOLS_BIN,
        "--fastqc_path", config.FASTQC_PATH,
        "--adapter_file", config.resolve_adapter(spec.trimming.adapterKit.value if spec.trimming.adapterKit else None,
                                                 spec.input.layout.value),
        "--trim_sliding_window_quality", str(q),
        "--trim_min_len", str(minlen),
        "--trim_adapter", str(spec.trimming.adapter).lower(),
        "--trim_poly_g", str(spec.trimming.polyG).lower(),
        "--mark_duplicates", str(spec.alignment.markDuplicates).lower(),
        "--min_qual", str(vc.minQual), "--min_dp", str(vc.minDp), "--min_gq", str(vc.minGq),
        "--max_depth", str(vc.maxDepth), "--min_bq", str(vc.minBq), "--min_mq", str(vc.minMq),
        # hemat disk: hard link (satu filesystem) dan tidak mem-publish FASTQ trimmed / BAM sebelum markdup
        "--publish_mode", "link", "--publish_intermediate", "false",
        "--step_tag", "pipeline",
    ]
    if spec.input.layout == Layout.PAIRED:
        cmd += ["--reads_r1", local["r1"], "--reads_r2", local["r2"]]
    else:
        cmd += ["--reads", local["r1"]]
    if spec.sampleSex:
        cmd += ["--sample_sex", spec.sampleSex.value]
    if "bed" in local:
        cmd += ["--target_bed", local["bed"]]
    return cmd


def _variant_counts(vcf: Path) -> Optional[dict]:
    try:
        out = subprocess.run([config.BCFTOOLS_BIN, "query", "-f", "%FILTER\n", str(vcf)],
                             capture_output=True, text=True, errors="replace", timeout=900, check=True).stdout
    except (subprocess.SubprocessError, OSError):
        return None
    filters = out.splitlines()
    return {"total": len(filters), "pass": sum(1 for f in filters if f in ("PASS", "."))}


def _collect_outputs(spec: JobSpec, out_dir: Path) -> List[tuple]:
    """Daftar (path lokal, key relatif terhadap prefix output)."""
    s = spec.sample
    base = out_dir / s
    vcf = base / "vcf" / f"{s}.filtered.vcf.gz"
    if not vcf.exists() or not (base / "vcf" / f"{s}.filtered.vcf.gz.tbi").exists():
        raise WorkerError(NO_OUTPUT, f"Pipeline selesai tetapi VCF tidak ditemukan: {vcf.name}")

    files = [(vcf, f"vcf/{vcf.name}"), (base / "vcf" / f"{vcf.name}.tbi", f"vcf/{vcf.name}.tbi")]
    if spec.output.includeBam:
        bam = base / "bwa" / f"{s}.markdup.bam"
        if not bam.exists() or not (base / "bwa" / f"{s}.markdup.bam.bai").exists():
            raise WorkerError(NO_OUTPUT, f"Pipeline selesai tetapi BAM tidak ditemukan: {bam.name}")
        files += [(bam, f"bam/{bam.name}"), (base / "bwa" / f"{bam.name}.bai", f"bam/{bam.name}.bai")]
    for name in (f"{s}.coverage.json", f"{s}.markdup.metrics.txt"):
        if (base / "bwa" / name).exists():
            files.append((base / "bwa" / name, f"qc/{name}"))
    for p in sorted((base / "fastqc").glob("*")) if (base / "fastqc").exists() else []:
        files.append((p, f"qc/{p.name}"))
    if (out_dir / ".nextflow.log").exists():
        files.append((out_dir / ".nextflow.log", "logs/nextflow.log"))
    return files


# ── job ───────────────────────────────────────────────────────────────────────

def run(spec: JobSpec, status: StatusReporter, s3: S3) -> None:
    scratch = Path(config.SCRATCH_DIR)
    job_dir = scratch / "jobs" / spec.jobId
    in_dir, out_dir, work_dir = job_dir / "in", job_dir / "out", job_dir / "work"
    if job_dir.exists():
        shutil.rmtree(job_dir, ignore_errors=True)
    for d in (in_dir, out_dir, work_dir):
        d.mkdir(parents=True, exist_ok=True)

    # 1. periksa input dan ruang disk sebelum membuang waktu
    status.stage("validate_input", 2)
    uris = {"r1": spec.input.r1}
    if spec.input.r2:
        uris["r2"] = spec.input.r2
    if spec.bed:
        uris["bed"] = spec.bed
    sizes = {k: s3.size_of(u) for k, u in uris.items()}
    fastq_bytes = sizes["r1"] + sizes.get("r2", 0)
    free = shutil.disk_usage(job_dir).free
    need = config.fastq_disk_need("vcf", [(uris[k], sizes[k]) for k in ("r1", "r2") if k in uris])
    eff_factor = need / fastq_bytes if fastq_bytes else 0.0
    log.info("Input %.2f GB, disk bebas %.2f GB, perkiraan kebutuhan %.2f GB (faktor %.1f)",
             fastq_bytes / 1e9, free / 1e9, need / 1e9, eff_factor)
    if not config.SKIP_DISK_CHECK and free < need:
        raise WorkerError(DISK_INSUFFICIENT,
                          f"Disk tidak cukup: bebas {free / 1e9:.1f} GB, perkiraan kebutuhan {need / 1e9:.1f} GB "
                          f"({eff_factor:g} x ukuran FASTQ {fastq_bytes / 1e9:.1f} GB). "
                          f"Perbesar volume task atau atur DISK_FACTOR.",
                          details={"freeBytes": free, "neededBytes": need, "fastqBytes": fastq_bytes})

    plan = resources.plan(spec.resources.threads)
    log.info("Sumber daya: %s", plan)
    disk = DiskPeak(job_dir).start()
    ref_fa = ensure_reference(spec.reference, s3, scratch)

    # 2. unduh input
    status.stage("download_input", 5)
    local: Dict[str, str] = {}
    for key, uri in uris.items():
        dest = in_dir / Path(parse_s3_uri(uri)[1]).name
        log.info("Mengunduh %s (%.2f GB)", uri, sizes[key] / 1e9)
        s3.download(uri, dest)
        local[key] = str(dest)

    # 3. pipeline
    status.stage("fastqc", 10)
    monitor = Monitor(status)
    cmd = _build_command(spec, plan, ref_fa, local, out_dir, work_dir)
    log.info("Menjalankan: %s", " ".join(cmd))
    rc = proc.run(cmd, nextflow_env(plan), config.PIPELINE_DIR, monitor.feed)
    if rc != 0:
        raise WorkerError(PIPELINE_FAILED, monitor.message(rc),
                          details=monitor.failure_details(rc, status.doc["stage"]))

    # 4. kumpulkan dan unggah hasil (work dir tidak dibutuhkan lagi; hard link menjaga file hasil tetap ada)
    status.stage("upload_output", 93)
    files = _collect_outputs(spec, out_dir)
    disk.sample()   # titik terpenuh biasanya tepat sebelum work dir dibuang
    shutil.rmtree(work_dir, ignore_errors=True)
    outputs = []
    for path, rel in files:
        uri = spec.output.prefix + rel
        outputs.append({"key": rel, "uri": uri, "sizeBytes": s3.upload(path, uri)})
        log.info("Upload %s", uri)

    # 5. ringkasan (tidak menggagalkan job bila gagal dibaca)
    metrics = {}
    cov = out_dir / spec.sample / "bwa" / f"{spec.sample}.coverage.json"
    if cov.exists():
        try:
            import json
            metrics["coverage"] = json.loads(cov.read_text())
        except ValueError:
            pass
    counts = _variant_counts(out_dir / spec.sample / "vcf" / f"{spec.sample}.filtered.vcf.gz")
    if counts:
        metrics["variants"] = counts
    metrics["resources"] = plan
    peak = disk.stop()
    peak["fastqBytes"] = fastq_bytes
    peak["peakOverFastq"] = round(peak["peakUsedBytes"] / fastq_bytes, 2) if fastq_bytes else None   # bandingkan dengan DISK_FACTOR
    metrics["scratch"] = peak
    log.info("Puncak disk scratch %.2f GB = %.2f x ukuran FASTQ (DISK_FACTOR=%g)",
             peak["peakUsedBytes"] / 1e9, peak["peakOverFastq"] or 0, config.DISK_FACTOR)
    status.succeed(outputs, metrics)
