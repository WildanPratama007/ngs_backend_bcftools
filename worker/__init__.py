"""
Worker ECS: task yang berjalan sekali lalu selesai (bukan server API).

    python -m worker vcf     # FASTQ (S3) -> VCF + BAM + QC (S3)
    python -m worker report  # FASTQ (S3) -> laporan FastQC + summary.json (S3)

Spesifikasi job dikirim sebagai JSON (env JOB_SPEC / JOB_SPEC_S3_URI / --spec-file).
Dokumentasi lengkap: docs/ECS_JOB_SPEC.md
"""

__version__ = "1.0.0"
