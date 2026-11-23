# 🧬 NGS Secondary Analysis — Backend BCFTools

Backend E2E pipeline NGS berbasis **BCFTools** dengan REST API.  
Pipeline: `FastQC → Trimmomatic → BWA MEM → MarkDuplicates → BCFTools mpileup/call/norm/filter`

Port default: **8001** | Versi: **1.2.0**

---

## 📁 Struktur

```
ngs_backend_bcftools/
├── api/
│   ├── main.py             ← Endpoints + validasi parameter
│   ├── models.py           ← Pydantic models + TrimPreset + CoverageStats
│   ├── job_manager.py      ← Job state + recovery saat restart
│   ├── pipeline_runner.py  ← Async Nextflow runner
│   └── config.py           ← Path & konfigurasi dari .env
│
├── pipeline/
│   ├── main.nf             ← Pipeline DSL2 (8 proses)
│   └── nextflow.config     ← Default parameter pipeline
│
├── scripts/
│   └── setup.sh
│
├── data/
│   ├── uploads/
│   ├── results/
│   └── jobs/
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

Dokumentasi API interaktif: **http://localhost:8001/docs**

---

## 📥 Persiapan Referensi Genom

BCFTools **tidak membutuhkan** file `.dict`. Cukup `.fa`, `.fai`, dan BWA index.

**GRCh38 chr22 (testing):**
```bash
curl -L -o ref_genome/chr22.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip ref_genome/chr22.fa.gz
bwa index ref_genome/chr22.fa
samtools faidx ref_genome/chr22.fa
```

**GRCh38 Full hg38 (production):**
```bash
curl -L -C - -o ref_genome/hg38.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ref_genome/hg38.fa.gz
bwa index ref_genome/hg38.fa
samtools faidx ref_genome/hg38.fa
```

| Referensi | URL | Keterangan |
|---|---|---|
| GRCh38 chr22 | [UCSC chr22](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz) | Default testing |
| GRCh38 Full | [UCSC hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) | Production |
| GRCh37 Full | [UCSC hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) | Data lama |

---

## 🔄 Alur Pipeline & Output per Tahap

```
Input: FASTQ (SE atau PE)
  │
  ▼ Step 1: FastQC — Quality Control
  │   → {sample}/fastqc/{sample}_fastqc.html
  │   → {sample}/fastqc/{sample}_fastqc.zip
  │
  ▼ Step 2: Trimmomatic — Adapter & Quality Trimming
  │   Parameter: trim_quality, trim_min_len, trim_adapter, trim_poly_g
  │   → {sample}/trimmed/{sample}_R1.trimmed.fastq  (PE)
  │   → {sample}/trimmed/{sample}_R2.trimmed.fastq  (PE)
  │   → {sample}/trimmed/{sample}.trimmed.fastq     (SE)
  │
  ▼ Step 3: BWA MEM — Alignment ke Referensi
  │   Read Group: ID={sample}.L001, SM={sample}, PL=ILLUMINA, LB={sample}
  │   → {sample}/bwa/{sample}.aligned.sorted.bam
  │
  ▼ Step 4: samtools markdup — Mark PCR Duplicates
  │   Parameter: mark_duplicates (on/off)
  │   → {sample}/bwa/{sample}.markdup.bam
  │   → {sample}/bwa/{sample}.markdup.metrics.txt
  │
  ▼ Step 5: Coverage Stats — Metrik Coverage
  │   → {sample}/bwa/{sample}.coverage.json
  │   (mean_depth, median_depth, pct_bases_10/20/30x, warning)
  │
  ▼ Step 6: Deteksi Kromosom dari BAM header
  │   (otomatis, hanya kromosom utama chr1-chr22, chrX, chrY, chrM)
  │
  ▼ Step 7: BCFTools mpileup (scatter per kromosom)
  │   Parameter: max_depth, min_bq, min_mq
  │   → work/.work/{sample}.{chrom}.bcf  (internal)
  │
  ▼ Step 8: BCFTools concat → call → norm → filter (gather)
      Parameter: min_qual, min_dp, min_gq, sample_sex
      Filter: QUAL>=min_qual && FORMAT/DP>=min_dp && FORMAT/GQ>=min_gq
      → {sample}/vcf/{sample}.filtered.vcf.gz  ✅ OUTPUT UTAMA
      → {sample}/vcf/{sample}.filtered.vcf.gz.tbi
```

---

## 🔧 Tools

| Tool | Versi | Fungsi | Instalasi |
|---|---|---|---|
| **Nextflow** | 24.10.6 | Workflow manager | conda |
| **FastQC** | 0.11.x | Quality control reads | local binary |
| **Trimmomatic** | 0.39 | Adapter & quality trimming | conda |
| **BWA MEM** | 0.7.19 | Short read alignment | conda |
| **Samtools** | 1.13 | Sort, index, markdup, coverage | apt `/usr/bin/` |
| **BCFTools** | ≥ 1.21 | mpileup, call, norm, filter, concat | build dari source → `/usr/local/bin/bcftools` |

> BCFTools ≥ 1.21 wajib: pipeline memakai `concat -a --rm-dups all` (opsi `--rm-dups` baru ada di versi baru dan hanya valid bersama `-a`). Versi apt 1.13 tidak cukup.

---

## ⚙️ Parameter Lengkap — User-Configurable

Semua parameter di bawah bisa diatur saat submit job via API atau Swagger UI (`/docs`).

---

### Step 2 — Trimming

| Parameter API | Type | Default | Range | Deskripsi |
|---|---|---|---|---|
| `trim_preset` | enum | — | `research`, `clinical` | Shortcut: research=Q20/36bp, clinical=Q30/36bp |
| `trim_quality` | int | 20 | 10–40 | Phred quality threshold SLIDINGWINDOW |
| `trim_min_len` | int | 36 | 20–100 | Minimum read length setelah trimming (bp) |
| `trim_adapter` | bool | true | true/false | Aktifkan ILLUMINACLIP adapter removal |
| `trim_poly_g` | bool | false | true/false | Tambah poly-G ke adapter list (NovaSeq/NextSeq) |

**Panduan memilih `trim_quality`:**

| Kondisi | Nilai | Alasan |
|---|---|---|
| Coverage rendah (<15x) | 15 | Pertahankan sebanyak mungkin reads |
| WGS standard (30x) | 20 | Balance (default) |
| WES panel | 25 | Kualitas lebih penting dari jumlah reads |
| Clinical grade | 30 | False positive minimal |

**Adapter file yang tersedia:**
```
TruSeq3-PE-2.fa    ← Default (Illumina TruSeq PE, HiSeq/NovaSeq)
TruSeq3-PE.fa      ← TruSeq PE versi lama
TruSeq2-PE.fa      ← Data lama (HiSeq 2000)
NexteraPE-PE.fa    ← Nextera XT kit
TruSeq3-SE.fa      ← Single-End
```
Adapter dipilih lewat parameter `adapter_kit` (wizard: dikirim di Step 1 dan dipakai di Step 2; `/api/pipeline/submit`: parameter form).
Kalau kosong, dipakai TruSeq SE untuk data SE dan TruSeq PE untuk data PE.

| `adapter_kit` | File adapter (override lewat `.env`) |
|---|---|
| `truseq_pe`, `truseq_pcr_free` | `ADAPTER_FILE` (TruSeq3-PE-2.fa) |
| `truseq_se` | `ADAPTER_FILE_SE` (TruSeq3-SE.fa) |
| `nextera`, `nextera_dna_flex` | `ADAPTER_NEXTERA` (NexteraPE-PE.fa) |

`trim_poly_g=true` menambah satu langkah `ILLUMINACLIP` dengan sekuens poly-G yang dibuat otomatis di tiap task (untuk NovaSeq/NextSeq).

Kalau Trimmomatic menyisakan 0 reads, error yang dikembalikan di `stages.trimming.error` menyebutkan apakah `trim_quality` terlalu tinggi.

---

### Step 4 — Alignment & MarkDuplicates

| Parameter API | Type | Default | Deskripsi |
|---|---|---|---|
| `mark_duplicates` | bool | true | Tandai PCR duplicate. `false` untuk PCR-free library atau amplicon sequencing |
| `sample_sex` | string | null | `male` atau `female` — ploidy chrX/chrY. Dipakai saat variant calling (wizard: dikirim di Step 4) |

**Kapan `mark_duplicates=false`:**

| Library type | Nilai | Alasan |
|---|---|---|
| Standard PCR (TruSeq) | `true` | Ada PCR duplicate dari library prep |
| PCR-free library | `false` | Tidak ada amplifikasi → tidak ada duplicate |
| Amplicon sequencing | `false` | Semua reads memang dari amplicon yang sama |

**Kenapa `sample_sex` penting:**
- Laki-laki: chrX/chrY hemizigot (1 copy) → ploidy 1
- Perempuan: chrX diploid (2 copy), chrY tidak ada → ploidy 2/0
- Tanpa ini, variant chrX laki-laki dilaporkan sebagai heterozygous padahal seharusnya hemizigot → penting untuk diagnosis X-linked disease (hemofilia, buta warna, dll)

---

### Step 7 — BCFTools mpileup

| Parameter API | Type | Default | Range | Deskripsi |
|---|---|---|---|---|
| `max_depth` | int | 500 | 100–10000 | Max read depth per posisi. 500=WGS, 2000=WES panel dalam |
| `min_bq` | int | 20 | 0–40 | Minimum base quality untuk dihitung di pileup |
| `min_mq` | int | 20 | 0–60 | Minimum mapping quality read |

**Panduan `max_depth`:**

| Use case | Nilai | Alasan |
|---|---|---|
| WGS standard (30x) | 500 | Coverage tidak akan melebihi ini |
| WES panel (200x) | 1000 | Cukup untuk panel standar |
| Deep panel (1000x) | 2000 | Sesuai UI SenusaBio default |

**Panduan `min_bq` dan `min_mq`:**

| Mode | min_bq | min_mq | Dampak |
|---|---|---|---|
| Permissive | 13 | 0 | Semua reads dihitung, false positive lebih tinggi |
| Standard | 20 | 20 | Balance (default) |
| Strict | 25 | 30 | False positive minimal, coverage turun sedikit |

---

### Step 8 — Variant Calling & Filter

| Parameter API | Type | Default | Range | Deskripsi |
|---|---|---|---|---|
| `min_qual` | int | 20 | 0–999 | Minimum QUAL score variant |
| `min_dp` | int | 10 | 1–500 | Minimum read depth (FORMAT/DP) |
| `min_gq` | int | 20 | 0–99 | Minimum Genotype Quality (FORMAT/GQ) |

**Filter expression yang dibangun otomatis:**
```
QUAL<{min_qual} || FORMAT/DP<{min_dp} || FORMAT/GQ<{min_gq}
```
Contoh dengan default: `QUAL<20 || FORMAT/DP<10 || FORMAT/GQ<20`

**Panduan filter per use case:**

| Use case | min_qual | min_dp | min_gq |
|---|---|---|---|
| Screening populasi | 10 | 3 | 10 |
| WGS research | 20 | 10 | 20 |
| WES clinical | 30 | 10 | 20 |
| Deep panel | 20 | 20 | 30 |

> Variant yang tidak lolos filter **tidak dihapus** — hanya diberi tag `LOWQUAL`. Downstream tools tetap bisa akses semua variant.

---

## 📊 Coverage Report (Output Baru)

Coverage tidak lagi hanya angka mean. Response API sekarang berisi object lengkap:

```json
{
  "coverage": {
    "mean_depth":     22.4,
    "median_depth":   21.0,
    "pct_bases_10x":  96.2,
    "pct_bases_20x":  85.3,
    "pct_bases_30x":  61.7,
    "total_bases":    51000000,
    "covered_bases":  48900000,
    "warning":        null
  }
}
```

**Warning otomatis:**
- `< 10x`: "Coverage sangat rendah — hasil variant calling tidak reliable"
- `< 20x`: "Coverage rendah — beberapa variant mungkin terlewat"

---

## 🌐 API Endpoints

### Submit Pipeline

```
POST /api/pipeline/submit
```

**Contoh submit lengkap:**
```bash
curl -X POST http://localhost:8001/api/pipeline/submit \
  -F "analysis_mode=WES" \
  -F "end_type=PE" \
  -F "reference=GRCh38" \
  -F "sample_sex=male" \
  # Trimming
  -F "trim_preset=clinical" \
  # Alignment
  -F "mark_duplicates=true" \
  # Variant calling
  -F "min_qual=20" \
  -F "min_dp=10" \
  -F "min_gq=20" \
  -F "max_depth=2000" \
  -F "min_bq=20" \
  -F "min_mq=20" \
  # Files
  -F "fastq_r1=@sample_R1.fastq.gz" \
  -F "fastq_r2=@sample_R2.fastq.gz" \
  -F "bed_file=@targets.bed"
```

> File harus berekstensi `.fastq.gz` atau `.fq.gz`

### Wizard — Satu Step per Request

Pipeline yang sama bisa dijalankan step demi step (dipakai `test_wizard.html`). Nextflow script-nya terpisah dari pipeline penuh: `pipeline/main_wizard.nf` (entry `QC`, `TRIMMING`, `ALIGNMENT`, `VARIANT_CALLING`), sedangkan `pipeline/main.nf` hanya untuk `/api/pipeline/submit`.

| Step | Endpoint | Syarat |
|---|---|---|
| 1 Upload + FastQC | `POST /api/analyses/{id}/qc/run` | — (form: `end_type`, `reference`, `analysis_mode`, `adapter_kit`, FASTQ, `bed_file` untuk WES/Targeted) |
| 2 Trimming | `POST /api/analyses/{id}/trimming/run` | Upload selesai. **Tidak menunggu FastQC** |
| 3 Alignment | `POST /api/analyses/{id}/alignment/run` | Trimming `completed` (SE dan PE) |
| 4 Variant Calling | `POST /api/analyses/{id}/variant-calling/run` | Alignment `completed`; param `sample_sex`, `min_qual`, `min_dp`, `min_gq`, `max_depth`, `min_bq`, `min_mq` |
| Status | `GET /api/analyses/{id}/status` | Status + progress + error tiap stage |

- Step 1 dan 2 independen dan boleh berjalan bersamaan. Hanya Step 2 → 3 → 4 yang berurutan.
- Menjalankan ulang sebuah step menghapus hasil step di hilirnya (mis. Trimming ulang me-reset Alignment dan Variant Calling).
- Step yang sedang `running` menolak request ganda dengan `409 STATUS_CONFLICT`.
- Untuk WES/Targeted, BED yang di-upload di Step 1 dipakai oleh coverage (Step 3) dan variant calling (Step 4).
- Pesan error dari task Nextflow yang gagal (blok `Command error:`) disimpan di `stages.<stage>.error`.
- `analysis_id`: huruf/angka/`_`/`-`/`.`, maksimal 64 karakter, diawali huruf atau angka (dipakai sebagai nama folder). ID lain ditolak dengan `VALIDATION_FAILED`.

### Endpoint sesuai API Contract (`/api/secondary-analysis`)

Bagian kontrak yang tidak butuh AWS (`api/contract.py`):

| Operasi kontrak | Endpoint | Catatan |
|---|---|---|
| getStatus | `GET /api/secondary-analysis/analyses/{id}/status` | `Status_View`: `analysisStatus`, `pipelineMode`, `run` (jobId, jobStatus, stage, progress, message, errorMessage, submittedAt, startedAt, finishedAt, mode). Waktu ISO 8601 UTC |
| getStatusBatch | `POST /api/secondary-analysis/status/batch` | 1–25 id unik; id tak ditemukan → `found: false` |
| cancel | `POST /api/secondary-analysis/analyses/{id}/cancel` body `{"jobId": "..."}` | Run antre → `CANCELLED`; Run berjalan → `CANCEL_REQUESTED` (Nextflow dihentikan, stage jadi `cancelled`) |

- Pembatalan juga tersedia tanpa body untuk wizard: `POST /api/analyses/{id}/cancel`.
- Semua error memakai envelope `{"error": {"code", "message", "details"}}`. Validasi request → `VALIDATION_FAILED` (400), tool pipeline tidak lengkap → `PIPELINE_NOT_CONFIGURED` (503), `PIPELINE_MODE=disabled` → `PIPELINE_DISABLED` (503).
- `jobId` = id analisis. `pipelineMode`/`mode` dari env `PIPELINE_MODE` (enum kontrak tidak punya nilai "local"; set `aws` di deployment AWS).
- Di mode wizard, jeda antar step (sebagian step selesai, tidak ada yang berjalan) dilaporkan `RUNNING` karena Run belum selesai.
- Step yang menunggu slot (`MAX_CONCURRENT_JOBS`) berstatus `pending` / `QUEUED`. Batas ini berlaku untuk pipeline penuh dan step wizard.

**Belum ada karena bergantung pada AWS/layanan lain:** registerFiles, startUpload, signParts, completeUpload, cancelUpload (S3 multipart), submit dan retry (FASTQ dari S3 + Job_Message SQS), autentikasi Cognito, role baca/tulis, pemisahan institusi, langganan, kuota, serta tabel transisi `DRAFT/UPLOADING/READY`.

### Endpoints Lain

| Method | Endpoint | Deskripsi |
|---|---|---|
| `GET` | `/api/pipeline/{id}/status` | Status + coverage report |
| `GET` | `/api/pipeline/jobs` | List semua job |
| `GET` | `/api/pipeline/queue` | Status antrian |
| `GET` | `/api/pipeline/{id}/fastqc/download` | Download FastQC HTML |
| `GET` | `/api/pipeline/{id}/trimming/download` | Download Trimmed FASTQ |
| `GET` | `/api/pipeline/{id}/alignment/download` | Download BAM (post-markdup) |
| `GET` | `/api/pipeline/{id}/variant-calling/download` | Download VCF (.vcf.gz) |
| `GET` | `/api/references` | List referensi tersedia |
| `GET` | `/api/health` | Health check + hasil pemeriksaan tool (bcftools ≥ 1.21, samtools, Java 17, Nextflow, FastQC, Trimmomatic, BWA). `?refresh=true` mengulang pemeriksaan. `status: degraded` bila ada tool wajib yang kurang |
| `DELETE` | `/api/pipeline/{id}` | Hapus job + file |

---

## ⚙️ Status Fitur

| Fitur | Status | Keterangan |
|---|---|---|
| Single-End (SE) | ✅ | End-to-end |
| Paired-End (PE) | ✅ | Dengan samtools fixmate |
| WGS | ✅ | Scatter-gather per kromosom |
| WES (BED) | ✅ | `-r chrom -T bed` per task, no duplikat |
| Targeted (BED) | ✅ | Sama dengan WES |
| GRCh38 | ✅ | hg38.fa tersedia |
| GRCh37 | ✅ | hg19.fa tersedia |
| Trimming user-configurable | ✅ | Preset + 4 parameter individual |
| MarkDuplicates toggle | ✅ | On/off via `mark_duplicates` |
| Read Group per lane | ✅ | Format `{sample}.L001` |
| Extended coverage stats | ✅ | mean, median, pct_bases_10/20/30x, warning |
| Ploidy chrX/chrY per sex | ✅ | Via `sample_sex` |
| Variant calling user-configurable | ✅ | 6 parameter: min_qual/dp/gq, max_depth, min_bq/mq |
| bcftools norm | ✅ | Left-align + split multiallelic |
| Output .vcf.gz + .tbi | ✅ | Siap untuk annotation |
| Work dir per job | ✅ | Tidak bentrok antar job paralel |
| Streaming upload | ✅ | Chunk 1MB, tidak OOM |
| Recovery saat restart | ✅ | Job running → failed |
| Validasi ekstensi file | ✅ | `.fastq.gz`/`.fq.gz` only |
| Docker profile | ✅ | bcftools v1.13 |
| Somatic calling | ❌ | Belum didukung |
| Multi-lane | ❌ | Planned |
| Autentikasi | ❌ | Planned (AWS Cognito) |

---

## ⚠️ Keterbatasan — Multi-lane

Pipeline ini hanya menerima **satu pasang FASTQ per job**.

Mesin Illumina sering menghasilkan beberapa file per sampel (multi-lane):
```
sample_L001_R1.fastq.gz + sample_L001_R2.fastq.gz  ← Lane 1
sample_L002_R1.fastq.gz + sample_L002_R2.fastq.gz  ← Lane 2
```

**Jika sampel kamu multi-lane, gabungkan dulu di lab sebelum submit:**
```bash
cat sample_L001_R1.fastq.gz sample_L002_R1.fastq.gz > sample_R1.fastq.gz
cat sample_L001_R2.fastq.gz sample_L002_R2.fastq.gz > sample_R2.fastq.gz
```

> Dukungan multi-lane direncanakan di versi berikutnya.

---

## 🖥️ Estimasi Waktu Pipeline

| Konfigurasi | WES 100MB | WES 15GB (30x) | WGS 174GB |
|---|---|---|---|
| Lokal `/mnt/c/` 4 thread | ~10-15 menit | ~2-3 jam | ~12-15 jam |
| Lokal `/home/` 6 thread | ~5-8 menit | ~1-1.5 jam | ~6-8 jam |
| AWS c7i.4xlarge EBS gp3 | **~1-2 menit** | **~15-25 menit** | **~1.5-2 jam** |

---

## 📌 Catatan Teknis

- `samtools` pakai `/usr/bin/samtools` (apt) — bukan conda (versi conda gagal di `libcrypto`); `bcftools` pakai `/usr/local/bin/bcftools` (≥ 1.21)
- Nextflow dijalankan dengan Java 17 (`NEXTFLOW_JAVA_HOME`) dan `NXF_VER` (`NEXTFLOW_VERSION`, default 24.10.4)
- Pipeline di-test pada Nextflow 24.10.4, BCFTools 1.21, WSL2 Ubuntu 22.04
- Coverage WES/Targeted dihitung dari `samtools depth -b {bed}` (persentase ≥10/20/30x terukur); WGS memakai `samtools coverage` (persentase diperkirakan dari mean depth)
- `--sample_name` dikirim API ke Nextflow supaya nama folder hasil sama dengan yang dicari backend
- Output VCF berformat `.vcf.gz` + `.tbi` (bgzipped + tabix indexed)
- Variant ditandai `LOWQUAL` (tidak dihapus) bila tidak lolos filter
- Work dir per job: `data/results/{job_id}/.work/` — aman untuk paralel
- `sample_sex` penting untuk diagnosis X-linked disease
- `max_depth=0` di bcftools bukan unlimited — pipeline menggunakan 8000 sebagai unlimited
