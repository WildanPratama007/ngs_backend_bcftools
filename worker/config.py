"""
Konfigurasi worker dari environment. Nilai bawaan sama dengan ENV di Dockerfile,
jadi di dalam image tidak ada yang wajib diatur kecuali akses S3 (IAM task role).
"""

import os
from pathlib import Path


def env(name: str, default: str = "") -> str:
    return os.environ.get(name, default)


# ── S3 ────────────────────────────────────────────────────────────────────────
AWS_REGION       = env("AWS_REGION") or env("AWS_DEFAULT_REGION") or "ap-southeast-3"   # Jakarta
S3_ENDPOINT_URL  = env("S3_ENDPOINT_URL")   # kosong = AWS. Isi untuk MinIO/LocalStack saat uji lokal
S3_SSE           = env("S3_SSE")            # "AES256" atau "aws:kms" (opsional; bucket default encryption tetap berlaku)
S3_KMS_KEY_ID    = env("S3_KMS_KEY_ID")

# ── Direktori kerja ───────────────────────────────────────────────────────────
# Disk lokal sementara (EBS/ephemeral). Semua input, work dir Nextflow, dan output ditaruh di sini.
SCRATCH_DIR      = env("SCRATCH_DIR", "/scratch")
DISK_FACTOR      = float(env("DISK_FACTOR", "8"))   # `worker vcf`: perkiraan ruang kerja = faktor x total ukuran FASTQ (.gz)
# Worker per langkah: faktor x total ukuran input langkah itu. Perkiraan awal, belum diukur pada data besar.
# Env DISK_FACTOR (bila diisi) menimpa semuanya.
_STEP_DISK_FACTOR = {"trimming": 3.0, "alignment": 6.0, "variant_calling": 3.0}


def disk_factor(step: str) -> float:
    return DISK_FACTOR if env("DISK_FACTOR") else _STEP_DISK_FACTOR.get(step, DISK_FACTOR)


# Faktor di atas dihitung terhadap FASTQ .gz. FASTQ teks biasa kira-kira sekian kali lebih besar daripada .gz
# (perkiraan, belum diukur), jadi kebutuhan disk dihitung per file: input .gz = faktor x ukuran;
# input teks = ukuran (file itu sendiri) + (faktor - 1) x ukuran / rasio, supaya tidak terlalu pesimis.
FASTQ_PLAIN_RATIO = max(1.0, float(env("FASTQ_PLAIN_RATIO", "4")))


def fastq_disk_need(step: str, files) -> int:
    """Perkiraan ruang kerja (byte) untuk daftar (nama/URI file FASTQ, ukuran byte)."""
    f, r = disk_factor(step), FASTQ_PLAIN_RATIO
    total = 0.0
    for name, size in files:
        total += size * f if name.lower().endswith(".gz") else size * (1 + (f - 1) / r)
    return int(total)


SKIP_DISK_CHECK  = env("SKIP_DISK_CHECK").lower() in ("1", "true", "yes")

# ── Tool & pipeline ───────────────────────────────────────────────────────────
PIPELINE_DIR     = env("PIPELINE_DIR", str(Path(__file__).resolve().parent.parent / "pipeline"))
NEXTFLOW_BIN     = env("NEXTFLOW_BIN", "nextflow")
NEXTFLOW_PROFILE = env("NEXTFLOW_PROFILE", "local_tools")
JAVA_HOME        = env("NEXTFLOW_JAVA_HOME", "/opt/java")
NEXTFLOW_VERSION = env("NEXTFLOW_VERSION", "24.10.4")
FASTQC_PATH      = env("FASTQC_PATH", "/usr/bin/fastqc")
BIN_DIR          = env("BIN_DIR", "/usr/local/bin")
SAMTOOLS_BIN     = env("SAMTOOLS_BIN", "/usr/bin/samtools")
BCFTOOLS_BIN     = env("BCFTOOLS_BIN", "/usr/local/bin/bcftools")

# ── Referensi genom ───────────────────────────────────────────────────────────
REF_BASE         = env("REF_BASE", "/data/ref_genome")
# Bila file referensi tidak ada di REF_BASE, unduh dari sini (mis. s3://bucket/ref_genome/).
# Nama file di S3 sama dengan lokal: hg38.fa, hg38.fa.fai, hg38.fa.bwt, ...
REF_S3_PREFIX    = env("REF_S3_PREFIX")
REF_FILES        = {"GRCh38_chr22": "chr22.fa", "GRCh38": "hg38.fa", "GRCh37": "hg19.fa"}
REF_SUFFIXES     = ["", ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]

# ── Adapter Trimmomatic ───────────────────────────────────────────────────────
_ADAPTER_PE      = env("ADAPTER_FILE",     "/usr/share/trimmomatic/TruSeq3-PE-2.fa")
_ADAPTER_SE      = env("ADAPTER_FILE_SE",  "/usr/share/trimmomatic/TruSeq3-SE.fa")
_ADAPTER_NEXTERA = env("ADAPTER_NEXTERA",  "/usr/share/trimmomatic/NexteraPE-PE.fa")
ADAPTER_KIT_MAP  = {
    "truseq_pe":        _ADAPTER_PE,
    "truseq_se":        _ADAPTER_SE,
    "truseq_pcr_free":  _ADAPTER_PE,
    "nextera":          _ADAPTER_NEXTERA,
    "nextera_dna_flex": _ADAPTER_NEXTERA,
}


def resolve_adapter(adapter_kit, layout: str) -> str:
    if adapter_kit:
        return ADAPTER_KIT_MAP[adapter_kit]
    return _ADAPTER_SE if layout == "SINGLE" else _ADAPTER_PE
