"""
NGS Backend — Konfigurasi dari environment variables (.env)
Fix: MAX_CONCURRENT_JOBS dibaca dari env, bukan hardcode
"""

import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent.parent / ".env")

# ── Tools ──────────────────────────────────────────────────────────────────────
NEXTFLOW_BIN  = os.getenv("NEXTFLOW_BIN",  "/home/mahasiswa/miniconda3/envs/ngs_platform/bin/nextflow")
FASTQC_PATH   = os.getenv("FASTQC_PATH",   "/mnt/c/Users/Mahasiswa/Fastt/FastQC/fastqc")
BIN_DIR       = os.getenv("BIN_DIR",       "/home/mahasiswa/miniconda3/envs/ngs_platform/bin")
SAMTOOLS_BIN  = os.getenv("SAMTOOLS_BIN",  "/usr/bin/samtools")
BCFTOOLS_BIN  = os.getenv("BCFTOOLS_BIN",  "/usr/local/bin/bcftools")
ADAPTER_FILE     = os.getenv("ADAPTER_FILE",    "/home/mahasiswa/miniconda3/envs/ngs_platform/share/trimmomatic/adapters/TruSeq3-PE-2.fa")
ADAPTER_FILE_SE  = os.getenv("ADAPTER_FILE_SE", "/home/mahasiswa/miniconda3/envs/ngs_platform/share/trimmomatic/adapters/TruSeq3-SE.fa")
ADAPTER_NEXTERA  = os.getenv("ADAPTER_NEXTERA", "/home/mahasiswa/miniconda3/envs/ngs_platform/share/trimmomatic/adapters/NexteraPE-PE.fa")

# Java/Nextflow untuk subprocess pipeline (Nextflow 24.x butuh Java 17+, bukan Java dari conda)
NEXTFLOW_JAVA_HOME = os.getenv("NEXTFLOW_JAVA_HOME", "/usr/lib/jvm/java-17-openjdk-amd64")
NEXTFLOW_VERSION   = os.getenv("NEXTFLOW_VERSION",   "24.10.4")

# Map adapter_kit → file path
ADAPTER_KIT_MAP = {
    "truseq_pe":       ADAPTER_FILE,
    "truseq_se":       ADAPTER_FILE_SE,
    "truseq_pcr_free": ADAPTER_FILE,
    "nextera":         ADAPTER_NEXTERA,
    "nextera_dna_flex":ADAPTER_NEXTERA,
}


def resolve_adapter_file(adapter_kit, end_type_value: str) -> str:
    """
    Pilih file adapter. Kit eksplisit menang; kalau tidak ada, pakai TruSeq sesuai SE/PE.
    Catatan: TruSeq3-SE.fa hanya berisi adapter SE, jadi kit PE tidak dipakai untuk data SE
    kecuali user memilihnya secara eksplisit.
    """
    if adapter_kit and adapter_kit in ADAPTER_KIT_MAP:
        return ADAPTER_KIT_MAP[adapter_kit]
    return ADAPTER_FILE_SE if end_type_value == "SE" else ADAPTER_FILE

# ── Referensi Genom ────────────────────────────────────────────────────────────
REF_BASE = os.getenv("REF_BASE", "/mnt/c/Users/Mahasiswa/Fastt/ref_genome")

REF_GENOMES = {
    "GRCh38_chr22": f"{REF_BASE}/chr22.fa",
    "GRCh38":       f"{REF_BASE}/hg38.fa",
    "GRCh37":       f"{REF_BASE}/hg19.fa",
}

# ── Nextflow Pipeline ──────────────────────────────────────────────────────────
NEXTFLOW_SCRIPT       = str(Path(__file__).parent.parent / "pipeline" / "main.nf")
NEXTFLOW_SCRIPT_WIZARD = str(Path(__file__).parent.parent / "pipeline" / "main_wizard.nf")
NEXTFLOW_CONFIG  = str(Path(__file__).parent.parent / "pipeline" / "nextflow.config")
NEXTFLOW_PROFILE = os.getenv("NEXTFLOW_PROFILE", "conda_env")

# ── Server ─────────────────────────────────────────────────────────────────────
# Origin yang boleh memanggil API, pisahkan dengan koma. Produksi: isi domain frontend.
CORS_ORIGINS = [o.strip() for o in os.getenv("CORS_ORIGINS", "*").split(",") if o.strip()]

# Mode yang dilaporkan ke frontend (enum kontrak: disabled | mock | aws).
# "disabled" menolak submit dengan PIPELINE_DISABLED. Set "aws" di deployment AWS.
PIPELINE_MODE = os.getenv("PIPELINE_MODE", "mock")

# ── Direktori Data ─────────────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).parent.parent
UPLOAD_DIR  = os.getenv("UPLOAD_DIR",  str(BASE_DIR / "data" / "uploads"))
RESULTS_DIR = os.getenv("RESULTS_DIR", str(BASE_DIR / "data" / "results"))
JOBS_DIR    = os.getenv("JOBS_DIR",    str(BASE_DIR / "data" / "jobs"))

for d in [UPLOAD_DIR, RESULTS_DIR, JOBS_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Pipeline ───────────────────────────────────────────────────────────────────
DEFAULT_TRIMMER_OPTS = "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25 LEADING:3 TRAILING:3"
DEFAULT_THREADS      = int(os.getenv("DEFAULT_THREADS", "4"))
# Fix: baca dari env, bukan hardcode di pipeline_runner
MAX_CONCURRENT_JOBS  = int(os.getenv("MAX_CONCURRENT_JOBS", "5"))
