"""
NGS Backend — Konfigurasi dari environment variables (.env)
"""

import os
from pathlib import Path
from dotenv import load_dotenv

# Load .env dari root backend
load_dotenv(Path(__file__).parent.parent / ".env")

# ── Tools ──────────────────────────────────────────────────────────────────────
NEXTFLOW_BIN  = os.getenv("NEXTFLOW_BIN",  "/home/mahasiswa/miniconda3/envs/ngs_platform/bin/nextflow")
FASTQC_PATH   = os.getenv("FASTQC_PATH",   "/mnt/c/Users/Mahasiswa/Fastt/FastQC/fastqc")
BIN_DIR       = os.getenv("BIN_DIR",       "/home/mahasiswa/miniconda3/envs/ngs_platform/bin")
SAMTOOLS_BIN  = os.getenv("SAMTOOLS_BIN",  "/usr/bin/samtools")
BCFTOOLS_BIN  = os.getenv("BCFTOOLS_BIN",  "/usr/bin/bcftools")

# ── Referensi Genom ────────────────────────────────────────────────────────────
REF_BASE = os.getenv("REF_BASE", "/mnt/c/Users/Mahasiswa/Fastt/ref_genome")

REF_GENOMES = {
    "GRCh38_chr22": f"{REF_BASE}/chr22.fa",
    "GRCh38":       f"{REF_BASE}/hg38.fa",
    "GRCh37":       f"{REF_BASE}/hg19.fa",
}

# ── Nextflow Pipeline ──────────────────────────────────────────────────────────
NEXTFLOW_SCRIPT  = str(Path(__file__).parent.parent / "pipeline" / "main.nf")
NEXTFLOW_CONFIG  = str(Path(__file__).parent.parent / "pipeline" / "nextflow.config")
NEXTFLOW_PROFILE = os.getenv("NEXTFLOW_PROFILE", "conda_env")

# ── Direktori Data ─────────────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).parent.parent
UPLOAD_DIR  = os.getenv("UPLOAD_DIR",  str(BASE_DIR / "data" / "uploads"))
RESULTS_DIR = os.getenv("RESULTS_DIR", str(BASE_DIR / "data" / "results"))
JOBS_DIR    = os.getenv("JOBS_DIR",    str(BASE_DIR / "data" / "jobs"))

for d in [UPLOAD_DIR, RESULTS_DIR, JOBS_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Pipeline ───────────────────────────────────────────────────────────────────
DEFAULT_TRIMMER_OPTS  = "SLIDINGWINDOW:4:20 MINLEN:25 LEADING:3 TRAILING:3"
DEFAULT_THREADS       = int(os.getenv("DEFAULT_THREADS", "4"))
MAX_CONCURRENT_JOBS   = int(os.getenv("MAX_CONCURRENT_JOBS", "5"))
