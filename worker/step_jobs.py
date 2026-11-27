"""
Worker per langkah wizard (opsi B: dipicu platform lewat ECS RunTask):

    worker trimming         FASTQ mentah (S3)         -> <analysisPrefix>trimmed/
    worker alignment        <analysisPrefix>trimmed/  -> <analysisPrefix>alignment/
    worker variant_calling  <analysisPrefix>alignment/-> <analysisPrefix>variants/

(Step 1 = `worker report`, hasil ke <analysisPrefix>qc/.) Tiap langkah menulis status.json di foldernya sendiri.
Menjalankan sebuah langkah menghapus hasil lamanya dan hasil langkah di hilirnya, supaya tidak ada berkas basi
yang terbaca langkah berikutnya.
"""

import json
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from . import config, proc, resources
from .diskwatch import DiskPeak
from .errors import DISK_INSUFFICIENT, INPUT_NOT_FOUND, NO_OUTPUT, PIPELINE_FAILED, WorkerError
from .nfcommon import PROGRESS_STEP, Monitor, ensure_reference, nextflow_env
from .s3io import S3
from .spec import (AlignmentStepSpec, Layout, Reference, TrimmingStepSpec, VariantStepSpec, parse_s3_uri)
from .status import StatusReporter

log = logging.getLogger("worker.step")

_SE_STATS = re.compile(r"Input Reads:\s*(\d+)\s+Surviving:\s*(\d+)\s*\(([\d.]+)%\)\s+Dropped:\s*(\d+)\s*\(([\d.]+)%\)")
_PE_STATS = re.compile(r"Input Read Pairs:\s*(\d+)\s+Both Surviving:\s*(\d+)\s*\(([\d.]+)%\)\s+Forward Only Surviving:\s*(\d+)"
                       r"\s*\(([\d.]+)%\)\s+Reverse Only Surviving:\s*(\d+)\s*\(([\d.]+)%\)\s+Dropped:\s*(\d+)\s*\(([\d.]+)%\)")


# ── nama berkas yang menjadi kontrak antar langkah ────────────────────────────

def trimmed_names(sample: str, layout: Layout) -> List[str]:
    if layout == Layout.PAIRED:
        return [f"{sample}_R1.trimmed.fastq.gz", f"{sample}_R2.trimmed.fastq.gz"]
    return [f"{sample}.trimmed.fastq.gz"]


def bam_names(sample: str) -> List[str]:
    return [f"{sample}.markdup.bam", f"{sample}.markdup.bam.bai"]


# ── persiapan (dipanggil sebelum status.json pertama ditulis) ─────────────────

def prepare(task: str, spec, s3: S3) -> None:
    """Periksa input ada, lalu hapus hasil lama langkah ini dan langkah hilirnya."""
    if task == "trimming":
        for uri in [spec.input.r1] + ([spec.input.r2] if spec.input.r2 else []):
            s3.size_of(uri)
    elif task == "alignment":
        tprefix = spec.step_prefix("trimming")
        missing = [n for n in trimmed_names(spec.sample, spec.layout) if not s3.exists(tprefix + n)]
        if missing:
            raise WorkerError(INPUT_NOT_FOUND, f"Hasil trimming tidak ditemukan di {tprefix} ({', '.join(missing)}). "
                                               f"Jalankan step Trimming dulu, dengan layout dan sampleName yang sama.",
                              details={"missing": missing, "prefix": tprefix})
    elif task == "variant_calling":
        aprefix = spec.step_prefix("alignment")
        missing = [n for n in bam_names(spec.sample) if not s3.exists(aprefix + n)]
        if missing:
            raise WorkerError(INPUT_NOT_FOUND, f"BAM tidak ditemukan di {aprefix} ({', '.join(missing)}). "
                                               f"Jalankan step Alignment dulu, dengan sampleName yang sama.",
                              details={"missing": missing, "prefix": aprefix})

    for prefix in [spec.status_prefix] + spec.downstream_prefixes:
        n = s3.delete_prefix(prefix)
        if n:
            log.info("Menghapus %d objek lama di %s", n, prefix)


# ── utilitas bersama ──────────────────────────────────────────────────────────

def _dirs(job_id: str) -> Tuple[Path, Path, Path, Path]:
    job_dir = Path(config.SCRATCH_DIR) / "jobs" / job_id
    if job_dir.exists():
        shutil.rmtree(job_dir, ignore_errors=True)
    in_dir, out_dir, work_dir = job_dir / "in", job_dir / "out", job_dir / "work"
    for d in (in_dir, out_dir, work_dir):
        d.mkdir(parents=True, exist_ok=True)
    return job_dir, in_dir, out_dir, work_dir


def _check_disk(step: str, job_dir: Path, input_bytes: int, need: Optional[int] = None):
    free = shutil.disk_usage(job_dir).free
    if need is None:
        need = int(input_bytes * config.disk_factor(step))
    factor = need / input_bytes if input_bytes else 0.0
    log.info("Input %.2f GB, disk bebas %.2f GB, perkiraan kebutuhan %.2f GB (faktor %.1f)",
             input_bytes / 1e9, free / 1e9, need / 1e9, factor)
    if not config.SKIP_DISK_CHECK and free < need:
        raise WorkerError(DISK_INSUFFICIENT,
                          f"Disk tidak cukup: bebas {free / 1e9:.1f} GB, perkiraan kebutuhan {need / 1e9:.1f} GB "
                          f"({factor:g} x ukuran input {input_bytes / 1e9:.1f} GB). Perbesar volume task atau atur DISK_FACTOR.",
                          details={"freeBytes": free, "neededBytes": need, "inputBytes": input_bytes})


def _base_cmd(entry: str, spec, plan: dict, out_dir: Path, work_dir: Path) -> list:
    return [
        config.NEXTFLOW_BIN, "-log", str(out_dir / ".nextflow.log"),
        "run", str(Path(config.PIPELINE_DIR) / "main_wizard.nf"),
        "-profile", config.NEXTFLOW_PROFILE,
        "-entry", entry,
        "-w", str(work_dir),
        "--outdir", str(out_dir),
        "--sample_name", spec.sample,
        "--threads", str(plan["threads"]),
        "--sort_mem", plan["sort_mem"],
        "--bin_dir", config.BIN_DIR,
        "--samtools_bin", config.SAMTOOLS_BIN,
        "--bcftools_bin", config.BCFTOOLS_BIN,
        "--fastqc_path", config.FASTQC_PATH,
        "--publish_mode", "link",
        "--step_tag", spec.STEP,
    ]


def _run_nextflow(cmd: list, plan: dict, status: StatusReporter, step: str):
    monitor = Monitor(status, PROGRESS_STEP[step])
    log.info("Menjalankan: %s", " ".join(cmd))
    rc = proc.run(cmd, nextflow_env(plan), config.PIPELINE_DIR, monitor.feed)
    if rc != 0:
        raise WorkerError(PIPELINE_FAILED, monitor.message(rc), details=monitor.failure_details(rc, step))


def _upload(s3: S3, files: List[Tuple[Path, str]], spec, status: StatusReporter) -> List[dict]:
    """Unggah (path lokal, nama relatif terhadap folder langkah). `key` di hasil relatif terhadap analysisPrefix."""
    status.stage("upload_output", 93)
    outputs = []
    for path, rel in files:
        uri = spec.status_prefix + rel
        size = s3.upload(path, uri)
        outputs.append({"key": uri[len(spec.analysisPrefix):], "uri": uri, "sizeBytes": size})
        log.info("Upload %s", uri)
    return outputs


def _finish(status: StatusReporter, outputs: List[dict], metrics: dict, plan: dict, disk: DiskPeak, input_bytes: int):
    metrics["resources"] = plan
    peak = disk.stop()
    peak["inputBytes"] = input_bytes
    peak["peakOverInput"] = round(peak["peakUsedBytes"] / input_bytes, 2) if input_bytes else None
    metrics["scratch"] = peak
    status.succeed(outputs, metrics)


def parse_trimmomatic_log(text: str) -> Optional[dict]:
    m = _PE_STATS.search(text)
    if m:
        g = m.groups()
        return {"inputPairs": int(g[0]), "bothSurviving": int(g[1]), "bothSurvivingPct": float(g[2]),
                "forwardOnly": int(g[3]), "reverseOnly": int(g[5]), "dropped": int(g[7]), "droppedPct": float(g[8])}
    m = _SE_STATS.search(text)
    if m:
        g = m.groups()
        return {"inputReads": int(g[0]), "surviving": int(g[1]), "survivingPct": float(g[2]),
                "dropped": int(g[3]), "droppedPct": float(g[4])}
    return None


# ── Step 2: trimming ──────────────────────────────────────────────────────────

def run_trimming(spec: TrimmingStepSpec, status: StatusReporter, s3: S3) -> None:
    job_dir, in_dir, out_dir, work_dir = _dirs(spec.jobId)
    status.stage("validate_input", 2)
    uris = {"r1": spec.input.r1}
    if spec.input.r2:
        uris["r2"] = spec.input.r2
    sizes = {k: s3.size_of(u) for k, u in uris.items()}
    input_bytes = sum(sizes.values())
    _check_disk("trimming", job_dir, input_bytes,
                need=config.fastq_disk_need("trimming", [(uris[k], sizes[k]) for k in uris]))
    plan = resources.plan(spec.resources.threads)
    disk = DiskPeak(job_dir).start()

    status.stage("download_input", 5)
    local: Dict[str, str] = {}
    for k, uri in uris.items():
        dest = in_dir / Path(parse_s3_uri(uri)[1]).name
        log.info("Mengunduh %s (%.2f GB)", uri, sizes[k] / 1e9)
        s3.download(uri, dest)
        local[k] = str(dest)

    q, minlen = spec.trimming.effective()
    cmd = _base_cmd("TRIMMING", spec, plan, out_dir, work_dir) + [
        "--adapter_file", config.resolve_adapter(spec.trimming.adapterKit.value if spec.trimming.adapterKit else None,
                                                 spec.layout.value),
        "--trim_sliding_window_quality", str(q), "--trim_min_len", str(minlen),
        "--trim_adapter", str(spec.trimming.adapter).lower(), "--trim_poly_g", str(spec.trimming.polyG).lower(),
        "--trimmed_ext", ".fastq.gz",
    ]
    cmd += (["--reads_r1", local["r1"], "--reads_r2", local["r2"]] if spec.layout == Layout.PAIRED else ["--reads", local["r1"]])
    status.stage("trimming", 10)
    _run_nextflow(cmd, plan, status, "trimming")

    base = out_dir / spec.sample / "trimmed"
    files: List[Tuple[Path, str]] = []
    for name in trimmed_names(spec.sample, spec.layout):
        p = base / name
        if not p.exists() or p.stat().st_size == 0:
            raise WorkerError(NO_OUTPUT, f"Hasil trimming tidak ditemukan: {name}")
        files.append((p, name))

    stats = None
    for lg in sorted(work_dir.rglob("trimmomatic.log")):
        stats = parse_trimmomatic_log(lg.read_text(errors="replace"))
        if stats:
            break
    metrics = {"sample": spec.sample, "layout": spec.layout.value, "trimming": stats}
    stats_path = out_dir / "stats.json"
    stats_path.write_text(json.dumps({"sample": spec.sample, "layout": spec.layout.value, "trimming": stats,
                                      "parameters": {"quality": q, "minLen": minlen, "adapter": spec.trimming.adapter,
                                                     "polyG": spec.trimming.polyG}}, indent=2))
    files.append((stats_path, "stats.json"))
    if (out_dir / ".nextflow.log").exists():
        files.append((out_dir / ".nextflow.log", "nextflow.log"))

    outputs = _upload(s3, files, spec, status)
    shutil.rmtree(work_dir, ignore_errors=True)
    _finish(status, outputs, metrics, plan, disk, input_bytes)


# ── Step 3: alignment ─────────────────────────────────────────────────────────

def _flagstat(bam: Path, threads: int) -> Optional[dict]:
    try:
        out = subprocess.run([config.SAMTOOLS_BIN, "flagstat", "-@", str(threads), str(bam)],
                             capture_output=True, text=True, errors="replace", timeout=900, check=True).stdout
    except (subprocess.SubprocessError, OSError):
        return None

    m_total = re.search(r"^(\d+) \+ \d+ in total", out, re.M)
    m_mapped = re.search(r"^(\d+) \+ \d+ mapped \(([\d.]+)%", out, re.M)
    m_dups = re.search(r"^(\d+) \+ \d+ duplicates", out, re.M)
    m_pp = re.search(r"^(\d+) \+ \d+ properly paired \(([\d.]+)%", out, re.M)
    return {"totalReads": int(m_total.group(1)) if m_total else None,
            "mapped": int(m_mapped.group(1)) if m_mapped else None,
            "mappedPct": float(m_mapped.group(2)) if m_mapped else None,
            "duplicates": int(m_dups.group(1)) if m_dups else None,
            "properlyPairedPct": float(m_pp.group(2)) if m_pp else None}


def run_alignment(spec: AlignmentStepSpec, status: StatusReporter, s3: S3) -> None:
    job_dir, in_dir, out_dir, work_dir = _dirs(spec.jobId)
    status.stage("validate_input", 2)
    tprefix = spec.step_prefix("trimming")
    uris = {n: tprefix + n for n in trimmed_names(spec.sample, spec.layout)}
    if spec.bed:
        uris["__bed__"] = spec.bed
    sizes = {k: s3.size_of(u) for k, u in uris.items()}
    input_bytes = sum(v for k, v in sizes.items() if k != "__bed__")
    _check_disk("alignment", job_dir, input_bytes)
    plan = resources.plan(spec.resources.threads)
    disk = DiskPeak(job_dir).start()
    ref_fa = ensure_reference(spec.reference, s3, Path(config.SCRATCH_DIR))

    status.stage("download_input", 5)
    local: Dict[str, str] = {}
    for k, uri in uris.items():
        dest = in_dir / Path(parse_s3_uri(uri)[1]).name
        log.info("Mengunduh %s (%.2f GB)", uri, sizes[k] / 1e9)
        s3.download(uri, dest)
        local[k] = str(dest)

    names = trimmed_names(spec.sample, spec.layout)
    cmd = _base_cmd("ALIGNMENT", spec, plan, out_dir, work_dir) + [
        "--ref_genome", ref_fa, "--mark_duplicates", str(spec.alignment.markDuplicates).lower(),
        "--publish_intermediate", "false",
    ]
    cmd += (["--trimmed_r1", local[names[0]], "--trimmed_r2", local[names[1]]] if spec.layout == Layout.PAIRED
            else ["--trimmed_reads", local[names[0]]])
    if "__bed__" in local:
        cmd += ["--target_bed", local["__bed__"]]
    status.stage("alignment", 10)
    _run_nextflow(cmd, plan, status, "alignment")

    base = out_dir / spec.sample / "bwa"
    bam = base / f"{spec.sample}.markdup.bam"
    if not bam.exists() or not (base / f"{bam.name}.bai").exists():
        raise WorkerError(NO_OUTPUT, f"BAM hasil alignment tidak ditemukan: {bam.name}")
    files = [(bam, bam.name), (base / f"{bam.name}.bai", f"{bam.name}.bai")]
    for name in (f"{spec.sample}.coverage.json", f"{spec.sample}.markdup.metrics.txt"):
        if (base / name).exists():
            files.append((base / name, name))
    if (out_dir / ".nextflow.log").exists():
        files.append((out_dir / ".nextflow.log", "nextflow.log"))

    metrics = {"sample": spec.sample, "layout": spec.layout.value, "flagstat": _flagstat(bam, plan["threads"])}
    cov = base / f"{spec.sample}.coverage.json"
    if cov.exists():
        try:
            metrics["coverage"] = json.loads(cov.read_text())
        except ValueError:
            pass

    outputs = _upload(s3, files, spec, status)
    shutil.rmtree(work_dir, ignore_errors=True)
    _finish(status, outputs, metrics, plan, disk, input_bytes)


# ── Step 4: variant calling ───────────────────────────────────────────────────

def run_variant_calling(spec: VariantStepSpec, status: StatusReporter, s3: S3) -> None:
    job_dir, in_dir, out_dir, work_dir = _dirs(spec.jobId)
    status.stage("validate_input", 2)
    aprefix = spec.step_prefix("alignment")
    uris = {n: aprefix + n for n in bam_names(spec.sample)}
    if spec.bed:
        uris["__bed__"] = spec.bed
    sizes = {k: s3.size_of(u) for k, u in uris.items()}
    input_bytes = sum(v for k, v in sizes.items() if k != "__bed__")
    _check_disk("variant_calling", job_dir, input_bytes)
    plan = resources.plan(spec.resources.threads)
    disk = DiskPeak(job_dir).start()
    ref_fa = ensure_reference(spec.reference, s3, Path(config.SCRATCH_DIR))

    status.stage("download_input", 5)
    local: Dict[str, str] = {}
    for k, uri in uris.items():
        dest = in_dir / Path(parse_s3_uri(uri)[1]).name
        log.info("Mengunduh %s (%.2f GB)", uri, sizes[k] / 1e9)
        s3.download(uri, dest)
        local[k] = str(dest)

    vc = spec.variantCalling
    cmd = _base_cmd("VARIANT_CALLING", spec, plan, out_dir, work_dir) + [
        "--ref_genome", ref_fa, "--ploidy", "GRCh37" if spec.reference == Reference.GRCh37 else "GRCh38",
        "--bam_file", local[f"{spec.sample}.markdup.bam"],
        "--min_qual", str(vc.minQual), "--min_dp", str(vc.minDp), "--min_gq", str(vc.minGq),
        "--max_depth", str(vc.maxDepth), "--min_bq", str(vc.minBq), "--min_mq", str(vc.minMq),
    ]
    if spec.sampleSex:
        cmd += ["--sample_sex", spec.sampleSex.value]
    if "__bed__" in local:
        cmd += ["--target_bed", local["__bed__"]]
    status.stage("variant_calling", 10)
    _run_nextflow(cmd, plan, status, "variant_calling")

    base = out_dir / spec.sample / "vcf"
    vcf = base / f"{spec.sample}.filtered.vcf.gz"
    if not vcf.exists() or not (base / f"{vcf.name}.tbi").exists():
        raise WorkerError(NO_OUTPUT, f"VCF hasil variant calling tidak ditemukan: {vcf.name}")
    files = [(vcf, vcf.name), (base / f"{vcf.name}.tbi", f"{vcf.name}.tbi")]
    if (out_dir / ".nextflow.log").exists():
        files.append((out_dir / ".nextflow.log", "nextflow.log"))

    metrics = {"sample": spec.sample}
    try:
        filt = subprocess.run([config.BCFTOOLS_BIN, "query", "-f", "%FILTER\n", str(vcf)], capture_output=True, text=True,
                              errors="replace", timeout=900, check=True).stdout.splitlines()
        metrics["variants"] = {"total": len(filt), "pass": sum(1 for f in filt if f in ("PASS", "."))}
    except (subprocess.SubprocessError, OSError):
        pass

    outputs = _upload(s3, files, spec, status)
    shutil.rmtree(work_dir, ignore_errors=True)
    _finish(status, outputs, metrics, plan, disk, input_bytes)
