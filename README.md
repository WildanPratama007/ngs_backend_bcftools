# 🧬 NGS Secondary Analysis — Backend BCFTools

Backend E2E pipeline NGS berbasis **BCFTools** dengan REST API.  
Pipeline: `FastQC → Trimmomatic → BWA MEM → BCFTools mpileup → bcftools call → bcftools filter`

Port default: **8001**

---

## 📁 Struktur

```
ngs_backend_bcftools/
├── api/                    ← FastAPI REST API
│   ├── main.py             ← Endpoints
│   ├── models.py           ← Pydantic models
│   ├── job_manager.py      ← Job state management
│   ├── pipeline_runner.py  ← Async Nextflow runner (semaphore 5 job)
│   └── config.py           ← Path & konfigurasi dari .env
│
├── pipeline/               ← Nextflow pipeline (BCFTools)
│   ├── main.nf             ← Pipeline DSL2
│   └── nextflow.config     ← Konfigurasi & resource allocation
│
├── scripts/
│   └── setup.sh            ← Setup environment & dependencies
│
├── data/                   ← Auto-dibuat, tidak masuk git
│   ├── uploads/            ← FASTQ input dari user
│   ├── results/            ← Output pipeline (BAM, VCF)
│   └── jobs/               ← Job state JSON
│
├── requirements.txt
├── .env.example
└── README.md
```

---

## 🚀 Quick Start

```bash
# 1. Setup environment
bash scripts/setup.sh

# 2. Copy dan edit konfigurasi
cp .env.example .env

# 3. Jalankan API
conda activate ngs_platform
uvicorn api.main:app --host 0.0.0.0 --port 8001 --reload
```

Dokumentasi API: **http://localhost:8001/docs**

---

## 📥 Persiapan Referensi Genom

BCFTools **tidak membutuhkan** file `.dict` (berbeda dengan GATK). Cukup `.fa`, `.fai`, dan BWA index.

### Download Referensi

**GRCh38 chr22 (untuk testing — cepat):**
```bash
# Download chr22 ~51MB
curl -L -o ref_genome/chr22.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip ref_genome/chr22.fa.gz

# Build BWA index (~3 menit)
bwa index ref_genome/chr22.fa
samtools faidx ref_genome/chr22.fa
# Tidak perlu .dict untuk BCFTools
```

**GRCh38 Full hg38 (untuk production):**
```bash
# Download hg38 ~938MB compressed
curl -L -C - -o ref_genome/hg38.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ref_genome/hg38.fa.gz

# Build index (~2 jam)
bwa index ref_genome/hg38.fa
samtools faidx ref_genome/hg38.fa
```

**GRCh37 Full hg19 (untuk data lama):**
```bash
curl -L -C - -o ref_genome/hg19.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip ref_genome/hg19.fa.gz
bwa index ref_genome/hg19.fa
samtools faidx ref_genome/hg19.fa
```

### File Index yang Dibutuhkan per Referensi

| File | Dibuat oleh | Dibutuhkan untuk |
|---|---|---|
| `*.fa` | Download | Semua tools |
| `*.fa.fai` | `samtools faidx` | BCFTools mpileup, BWA |
| `*.fa.amb` | `bwa index` | BWA MEM |
| `*.fa.ann` | `bwa index` | BWA MEM |
| `*.fa.bwt` | `bwa index` | BWA MEM |
| `*.fa.pac` | `bwa index` | BWA MEM |
| `*.fa.sa` | `bwa index` | BWA MEM |

> ✅ Tidak perlu `.dict` — BCFTools tidak membutuhkannya (berbeda dengan GATK)

---

## 🔄 Alur Pipeline & Output per Tahap

```
Input: FASTQ (SE atau PE)
```

### Tahap 1 — FastQC (Quality Control)
```
Input:  sample.fastq (SE) atau sample_R1/R2.fastq (PE)
Output: results/{job_id}/{sample}/fastqc/
        ├── {sample}_fastqc.html  ← Laporan QC visual
        └── {sample}_fastqc.zip   ← Data QC untuk MultiQC
```

### Tahap 2 — Trimmomatic (Adapter Trimming)
```
Input:  sample.fastq
Output: results/{job_id}/{sample}/trimmed/
        ├── {sample}.trimmed.fastq         ← SE: reads bersih
        ├── {sample}_R1.trimmed.fastq      ← PE: R1 paired
        ├── {sample}_R1.unpaired.fastq     ← PE: R1 yang pasangannya dibuang
        ├── {sample}_R2.trimmed.fastq      ← PE: R2 paired
        └── {sample}_R2.unpaired.fastq     ← PE: R2 yang pasangannya dibuang
```

### Tahap 3 — BWA MEM + Samtools (Alignment)
```
Input:  {sample}.trimmed.fastq + ref_genome.fa
Output: results/{job_id}/{sample}/bwa/
        ├── {sample}.aligned.sorted.bam     ← BAM sorted
        ├── {sample}.aligned.sorted.bam.bai ← Index BAM
        └── {sample}.coverage.txt           ← Mean coverage (📊 COVERAGE: XX.Xx)
```

### Tahap 4 — BCFTools mpileup (Scatter per Kromosom)
```
Input:  {sample}.aligned.sorted.bam + ref_genome.fa
Output: work/ (internal Nextflow, scatter per kromosom)
        └── {sample}.{chrom}.bcf      ← Pileup per kromosom
        └── {sample}.{chrom}.bcf.csi  ← Index BCF
```

### Tahap 5 — BCFTools concat + call + filter (Gather)
```
Input:  semua {sample}.{chrom}.bcf
Output: results/{job_id}/{sample}/vcf/
        └── {sample}.filtered.vcf  ← VCF final yang sudah difilter
```

Pipeline dalam proses GATHER_VCF:
```
bcftools concat → bcftools call -m -v → bcftools filter -e 'QUAL<20 || DP<5'
```

### Output Utama untuk Tertiary Analysis
```
✅ {sample}.filtered.vcf  ← File yang perlu dikirim ke tertiary analysis
```

---

## 🔧 Tools, Framework & Library

### Pipeline Tools

| Tool | Versi | Fungsi | Instalasi |
|---|---|---|---|
| **Nextflow** | 24.10.6 | Workflow manager | conda bioconda |
| **FastQC** | 0.11.x | Quality control reads | Lokal binary |
| **Trimmomatic** | 0.39 | Adapter & quality trimming | conda bioconda |
| **BWA MEM** | 0.7.19 | Short read alignment ke GRCh38/GRCh37 | conda bioconda |
| **Samtools** | 1.13 | BAM sort, index, coverage stats | apt (`/usr/bin/samtools`) |
| **BCFTools mpileup** | 1.13 | Pileup per kromosom dari BAM | apt (`/usr/bin/bcftools`) |
| **BCFTools call** | 1.13 | Variant calling dari pileup (`-m -v`) | apt (`/usr/bin/bcftools`) |
| **BCFTools filter** | 1.13 | Filter berdasarkan QUAL dan DP | apt (`/usr/bin/bcftools`) |
| **BCFTools concat** | 1.13 | Gabung BCF per kromosom (scatter-gather) | apt (`/usr/bin/bcftools`) |

> ⚠️ Samtools dan BCFTools menggunakan versi dari **apt** (`/usr/bin/`) karena versi conda mengalami konflik `libcrypto.so.1.0.0` di Ubuntu versi baru.

### API Framework & Library

| Library | Versi | Fungsi |
|---|---|---|
| **FastAPI** | 0.111.0 | REST API framework |
| **Uvicorn** | 0.29.0 | ASGI server |
| **Pydantic** | 2.7.1 | Data validation & models |
| **python-multipart** | 0.0.9 | File upload handling |
| **aiofiles** | 23.2.1 | Async file I/O |
| **python-dotenv** | 1.0.1 | Environment variable loading |

### Parameter Pipeline — Nilai Default

| Parameter | Nilai Default | Alasan |
|---|---|---|
| `trimmer_opts` | `SLIDINGWINDOW:4:20 MINLEN:25 LEADING:3 TRAILING:3` | Balance kualitas dan jumlah reads tersisa |
| `threads` | `4` | Aman untuk 6-core lokal, scalable ke 16 di AWS |
| `samtools sort -m` | `4G` per thread | Aman untuk RAM 16GB+ |
| `ploidy` | `GRCh38` | Standard human diploid (bcftools call --ploidy) |
| `min_qual` | `20` | Filter QUAL<20 (probabilitas error <1%) |
| `min_dp` | `5` | Filter DP<5 (minimum read depth per posisi) |
| `-m` (bcftools call) | multiallelic model | Lebih akurat dari consensus model untuk human |
| `-v` (bcftools call) | variant sites only | Hanya output posisi yang ada variant |

### Referensi Genom yang Didukung

| Referensi | URL Download | Index | Keterangan |
|---|---|---|---|
| GRCh38 chr22 | [UCSC hg38 chr22](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz) | ✅ `.fai` + BWA | Default testing |
| GRCh38 Full | [UCSC hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) | ✅ | Production |
| GRCh37 Full | [UCSC hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) | ✅ | Data lama |

---

## 🌐 API Endpoints

| Method | Endpoint | Deskripsi |
|---|---|---|
| `POST` | `/api/pipeline/submit` | Submit pipeline (SE/PE, WGS/WES/Targeted) |
| `GET` | `/api/pipeline/{id}/status` | Cek status semua tahap |
| `GET` | `/api/pipeline/jobs` | List semua job |
| `GET` | `/api/pipeline/queue` | Status antrian (max 5 job paralel) |
| `GET` | `/api/pipeline/{id}/fastqc/download` | Download FastQC HTML |
| `GET` | `/api/pipeline/{id}/trimming/download` | Download Trimmed FASTQ |
| `GET` | `/api/pipeline/{id}/alignment/download` | Download BAM |
| `GET` | `/api/pipeline/{id}/variant-calling/download` | Download VCF |
| `GET` | `/api/references` | List referensi genom tersedia |
| `GET` | `/api/health` | Health check |
| `DELETE` | `/api/pipeline/{id}` | Hapus job |

---

## ⚙️ Status Fitur

| Fitur | Status | Keterangan |
|---|---|---|
| Single-End (SE) | ✅ Didukung | Berjalan end-to-end |
| Paired-End (PE) | ✅ Didukung | PE dengan samtools fixmate |
| WGS | ✅ Didukung | Scatter-gather per kromosom |
| WES (dengan BED) | ✅ Ada di kode | `-R bed_file` di mpileup |
| Targeted (dengan BED) | ✅ Ada di kode | `-R bed_file` di mpileup |
| GRCh38 | ✅ | hg38.fa + index tersedia |
| GRCh37 | ✅ | hg19.fa + index tersedia |
| Multi-job paralel | ✅ | Semaphore 5 job bersamaan |
| Scatter-gather BCFTools | ✅ | mpileup per kromosom, concat + call + filter |
| Coverage stats otomatis | ✅ | `📊 COVERAGE: XX.Xx` di log |

---

## 🖥️ Analisis Komputasi

### Estimasi Waktu Pipeline

| Konfigurasi | WES 100MB | WES 15GB (30x) | WGS 174GB |
|---|---|---|---|
| **Lokal saat ini** (`/mnt/c/`, 4 thread) | ~8-10 menit | ~2-3 jam | ~10-15 jam |
| Lokal optimal (`/home/`, 6 thread) | ~4-6 menit | ~1-1.5 jam | ~5-8 jam |
| AWS c7i.4xlarge EBS gp3 (16 thread) | **~1-2 menit** | **~15-25 menit** | **~1.5-2 jam** |

> BCFTools umumnya **2-3x lebih cepat** dari GATK untuk data coverage rendah-sedang karena tidak ada MarkDuplicates, HaplotypeCaller yang berat, atau GVCF intermediate.

### BCFTools vs GATK — Perbandingan Singkat

| Aspek | BCFTools | GATK HaplotypeCaller |
|---|---|---|
| Kecepatan | ✅ Lebih cepat (~2-3x) | ❌ Lebih lambat |
| Akurasi SNP | ✅ Baik untuk coverage >20x | ✅ Lebih akurat untuk coverage rendah |
| Akurasi INDEL | ⚠️ Kurang akurat untuk INDEL panjang | ✅ Lebih akurat |
| Kebutuhan resource | ✅ RAM lebih hemat | ❌ RAM lebih besar |
| File `.dict` | ✅ Tidak perlu | ❌ Wajib |
| BQSR | ✅ Tidak perlu | ⚠️ Disarankan untuk klinis |
| Cocok untuk | WGS coverage >20x, screening | Clinical grade, INDEL analysis |

### Storage I/O

| Kondisi | Throughput | Dampak |
|---|---|---|
| `/mnt/c/` (NTFS via WSL) | ~100-200 MB/s | BCFTools baca/tulis lambat |
| Linux native `/home/` | ~400-500 MB/s | ~2x lebih cepat |
| AWS EBS gp3 | ~1,000 MB/s | ~5x lebih cepat |

```bash
# Pindah data/ ke Linux FS untuk performa lebih baik
# Edit .env:
UPLOAD_DIR=/home/mahasiswa/ngs_data/uploads
RESULTS_DIR=/home/mahasiswa/ngs_data/results
```

---

## 📌 Catatan Teknis

- `samtools` dan `bcftools` pakai `/usr/bin/` (apt) — bukan conda, untuk menghindari konflik `libcrypto.so.1.0.0`
- Pipeline di-test pada Nextflow 24.10.6, BCFTools 1.13, WSL2 Ubuntu
- BCFTools tidak memerlukan MarkDuplicates sebagai prerequisite — langsung mpileup dari BAM sorted
- Filter default `QUAL<20 || DP<5` sudah cukup untuk screening; untuk klinis pertimbangkan `QUAL<30 || DP<10`
- Coverage dihitung via `samtools depth` → mean depth seluruh genome/region
