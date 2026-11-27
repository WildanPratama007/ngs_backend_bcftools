# Parameter Job untuk Container ECS (FASTQ → VCF dan FASTQ → Report)

Satu image (`ngs-backend-bcftools`) berisi dua task yang **berjalan sekali lalu selesai**. Tidak ada server API di mode ini.
Input dibaca dari S3, hasil ditulis ke S3, dan status akhir ditulis ke `status.json`.

| Task | Command container | Fungsi | Compute |
|---|---|---|---|
| FASTQ → VCF | `python -m worker vcf` | FastQC → Trimmomatic → BWA MEM → MarkDup → coverage → BCFTools (VCF + BAM + QC) | Besar: 8–16+ vCPU, ≥ 32 GB, disk besar |
| FASTQ → Report | `python -m worker report` | FastQC (HTML/ZIP) + `summary.json` | Kecil: 2 vCPU, 4 GB, disk ≈ 1,5 × FASTQ |

Image tanpa `command` menjalankan API (`uvicorn`). Untuk worker, **selalu isi `command`** di task definition atau di override `RunTask`.
Docker `HEALTHCHECK` di image hanya untuk mode API; ECS memakai `healthCheck` di task definition, dan untuk `docker run` lokal pakai `--no-healthcheck`.

---

## 1. Cara mengirim parameter

Spesifikasi job adalah **satu JSON**. Tiga cara memberikannya (urutan prioritas):

| Cara | Kapan |
|---|---|
| `--spec-file path.json` | Uji lokal |
| Env `JOB_SPEC` = isi JSON | Spesifikasi kecil. Override `RunTask` dibatasi sekitar 8 KB (cek dokumentasi ECS) |
| Env `JOB_SPEC_S3_URI` = `s3://bucket/key.json` | **Disarankan.** Tulis JSON ke S3 lebih dulu, kirim hanya URI-nya. Tanpa masalah escape JSON-dalam-JSON dan batas ukuran |

Contoh `RunTask` (region Jakarta):

```bash
aws ecs run-task --region ap-southeast-3 \
  --cluster ngs --launch-type EC2 \
  --task-definition ngs-worker-vcf \
  --overrides '{"containerOverrides":[{"name":"worker","environment":[
      {"name":"JOB_SPEC_S3_URI","value":"s3://ngs-data/specs/ana_001.json"}]}]}'
```

Task `report` identik, dengan task definition `ngs-worker-report` (command `python -m worker report`).

Periksa spesifikasi tanpa menjalankan apa pun (tidak butuh S3):

```bash
python -m worker vcf --validate-only --spec-file spec.json     # exit 0 = valid, 2 = tidak valid
```

---

## 2. Spesifikasi `worker vcf`

Field yang tidak dikenal **ditolak** (exit 2). Rentang sama dengan API (`/api/pipeline/submit`).

```json
{
  "jobId": "ana_001",
  "sampleName": "S1",
  "input":  { "layout": "PAIRED",
              "r1": "s3://ngs-data/in/S1_L001_R1_001.fastq.gz",
              "r2": "s3://ngs-data/in/S1_L001_R2_001.fastq.gz" },
  "output": { "prefix": "s3://ngs-data/runs/ana_001/", "includeBam": true },
  "reference": "GRCh38",
  "analysisMode": "WES",
  "bed": "s3://ngs-data/refs/targets.bed",
  "sampleSex": "male",
  "trimming":       { "preset": null, "quality": 30, "minLen": 36, "adapter": true, "adapterKit": "truseq_pe", "polyG": false },
  "alignment":      { "markDuplicates": true },
  "variantCalling": { "minQual": 20, "minDp": 10, "minGq": 20, "maxDepth": 500, "minBq": 20, "minMq": 20 },
  "resources":      { "threads": null }
}
```

| Field | Tipe | Wajib | Default | Aturan |
|---|---|---|---|---|
| `jobId` | string | ya | — | `A-Za-z0-9_.-`, maks 64, diawali huruf/angka, tanpa `..` |
| `sampleName` | string | tidak | dari nama file R1 | Aturan sama dengan `jobId`. Dipakai sebagai read group `SM` dan nama file hasil |
| `input.layout` | `PAIRED` \| `SINGLE` | ya | — | |
| `input.r1` | URI S3 | ya | — | Berakhiran `.fastq.gz`, `.fq.gz`, `.fastq`, atau `.fq` (huruf kecil). FASTQ teks biasa boleh; perkiraan disk memakai asumsi teks ±4× lebih besar dari `.gz` (`FASTQ_PLAIN_RATIO`, belum diukur) |
| `input.r2` | URI S3 | wajib bila `PAIRED` | — | Tidak boleh ada bila `SINGLE`; tidak boleh sama dengan `r1` |
| `output.prefix` | URI S3 | ya | — | `s3://bucket/prefix`. Harus ada prefix (tidak boleh akar bucket). Garis miring akhir ditambahkan otomatis |
| `output.includeBam` | bool | tidak | `true` | Upload `markdup.bam` + `.bai` |
| `reference` | `GRCh38` \| `GRCh37` \| `GRCh38_chr22` | ya | — | `GRCh38_chr22` hanya untuk uji |
| `analysisMode` | `WGS` \| `WES` \| `Targeted` | tidak | `WGS` | |
| `bed` | URI S3 | wajib bila `WES`/`Targeted` | — | Berakhiran `.bed`; ditolak bila `WGS` |
| `sampleSex` | `male` \| `female` | tidak | `null` | Ploidy chrX/chrY (`bcftools call -S`). Tanpa nilai ini semua kromosom diploid |
| `trimming.preset` | `research` \| `clinical` | tidak | `null` | Menimpa `quality`/`minLen`: research = Q20, clinical = Q30 (minLen 36) |
| `trimming.quality` | int | tidak | 20 | 10–40 (SLIDINGWINDOW) |
| `trimming.minLen` | int | tidak | 36 | 20–100 |
| `trimming.adapter` | bool | tidak | `true` | ILLUMINACLIP |
| `trimming.adapterKit` | enum | tidak | SE→`truseq_se`, PE→`truseq_pe` | `truseq_pe`, `truseq_se`, `truseq_pcr_free`, `nextera`, `nextera_dna_flex` |
| `trimming.polyG` | bool | tidak | `false` | Untuk NovaSeq/NextSeq |
| `alignment.markDuplicates` | bool | tidak | `true` | `false` untuk PCR-free / amplicon |
| `variantCalling.minQual` | int | tidak | 20 | 0–999 |
| `variantCalling.minDp` | int | tidak | 10 | 1–500 |
| `variantCalling.minGq` | int | tidak | 20 | 0–99 |
| `variantCalling.maxDepth` | int | tidak | 500 | 100–10000 |
| `variantCalling.minBq` | int | tidak | 20 | 0–40 |
| `variantCalling.minMq` | int | tidak | 20 | 0–60 |
| `resources.threads` | int | tidak | otomatis | 1–64. Otomatis = jumlah CPU container (maks 16) |

Catatan:
- **Satu pasang FASTQ per job.** Multi-lane belum didukung; gabungkan lane sebelum diunggah.
- Batas CPU dan memori Nextflow dibaca otomatis dari container (cgroup), tidak perlu dikirim.
- Variant ditandai `LOWQUAL` di kolom FILTER bila tidak lolos `minQual/minDp/minGq`; tidak dihapus.

## 3. Spesifikasi `worker report`

```json
{
  "jobId": "rep_001",
  "input":  { "layout": "PAIRED", "r1": "s3://ngs-data/in/S1_R1.fastq.gz", "r2": "s3://ngs-data/in/S1_R2.fastq.gz" },
  "output": { "prefix": "s3://ngs-data/reports/rep_001/" },
  "resources": { "threads": null }
}
```

Field sama dengan di atas (`jobId`, `input`, `output.prefix`, `resources`). Tidak ada parameter lain. Cocok untuk WES maupun WGS, karena hanya membaca FASTQ.

---

## 4. Hasil di S3

### `worker vcf` → `<output.prefix>`

| Key | Isi |
|---|---|
| `vcf/<sample>.filtered.vcf.gz` + `.tbi` | **Hasil utama** (bgzip + tabix, sudah norm dan difilter) |
| `bam/<sample>.markdup.bam` + `.bai` | Bila `includeBam: true` |
| `qc/<sample>.coverage.json` | mean/median depth, %≥10x/20x/30x, covered bases, peringatan |
| `qc/<sample>.markdup.metrics.txt` | Metrik duplikat |
| `qc/<file>_fastqc.html` / `.zip` | FastQC tiap FASTQ |
| `logs/nextflow.log` | Log Nextflow |
| `status.json` | Status akhir (lihat bagian 5) |

Coverage WES/Targeted dihitung di region BED (terukur). WGS memakai `samtools coverage`; persentase ≥10x/20x/30x diperkirakan dari mean depth.

### `worker report` → `<output.prefix>`

| Key | Isi |
|---|---|
| `summary.json` | Ringkasan semua file (bagian 6) |
| `qc/<file>_fastqc.html` / `.zip` | Laporan FastQC per FASTQ |
| `status.json` | Status akhir |

> Di `worker report` nama laporan mengikuti nama file: `qc/<nama>_fastqc.html`. Bila file FASTQ bernama sama di dua prefix berbeda, hasilnya tetap terpisah karena `output.prefix` berbeda.

## 5. `status.json`

Ditulis saat mulai, tiap pergantian tahap (maks 1 kali per 5 detik untuk update progres), dan **selalu** di akhir, termasuk saat gagal. Baca file ini (atau event ECS task stopped) untuk memperbarui DynamoDB/UI.

```json
{
  "schemaVersion": 1, "task": "vcf", "jobId": "ana_001",
  "status": "SUCCEEDED",
  "stage": "finalize", "progress": 100,
  "errorCode": null, "errorMessage": null, "details": null,
  "startedAt": "2026-11-25T08:00:00.000Z", "updatedAt": "...", "finishedAt": "...",
  "outputs": [ { "key": "vcf/S1.filtered.vcf.gz", "uri": "s3://ngs-data/runs/ana_001/vcf/S1.filtered.vcf.gz", "sizeBytes": 1578 } ],
  "metrics": {
    "coverage":  { "mean_depth": 98.2, "pct_bases_20x": 97.1, "...": "..." },
    "variants":  { "total": 52311, "pass": 49870 },
    "resources": { "cpus": 16, "memory_mb": 28000, "threads": 16, "sort_mem": "700M" },
    "scratch":   { "peakUsedBytes": 123456789012, "fastqBytes": 15000000000, "peakOverFastq": 8.2 }
  }
}
```

- `status`: `RUNNING` | `SUCCEEDED` | `FAILED`.
- `stage` (urut): `validate_input` → `download_input` → `fastqc` → `trimming` → `alignment` → `variant_calling` → `upload_output` → `finalize`.
- `progress`: 0–100, hanya naik. 100 hanya saat `SUCCEEDED`.
- `metrics.scratch.peakOverFastq` adalah **puncak pemakaian disk dibagi ukuran FASTQ**. Bandingkan dengan `DISK_FACTOR` (bawaan 8) setelah run nyata pertama, lalu setel (bagian 8).

| `errorCode` | Exit code | Arti |
|---|---|---|
| `INVALID_SPEC` | 2 | Spesifikasi tidak valid. `details.fields` memuat `[{path, reason}]`. Tidak ada yang dijalankan |
| `INPUT_NOT_FOUND` | 1 | FASTQ/BED tidak ada, atau tidak ada izin baca |
| `DOWNLOAD_FAILED` | 1 | Unduhan gagal atau ukuran tidak cocok |
| `DISK_INSUFFICIENT` | 1 | Ruang disk < `DISK_FACTOR` × ukuran FASTQ. Dicek sebelum mengunduh |
| `REFERENCE_NOT_FOUND` | 1 | Referensi tidak lengkap di `REF_BASE` dan `REF_S3_PREFIX` |
| `PIPELINE_FAILED` | 1 | Proses Nextflow/FastQC gagal. `errorMessage` memuat pesan task; `details.failedProcess`, `details.stage` |
| `NO_OUTPUT` | 1 | Pipeline selesai tetapi file hasil tidak ada |
| `UPLOAD_FAILED` | 1 | Unggahan gagal atau ukuran tidak cocok |
| `TERMINATED` | 143 | Task menerima SIGTERM (stop, deploy, atau interupsi Spot) |
| `INTERNAL_ERROR` | 1 | Kesalahan tak terduga (lihat log CloudWatch) |

Bila task dimatikan paksa (SIGKILL setelah `stopTimeout`, OOM kill), `status.json` terakhir bisa tetap `RUNNING`. Pemanggil harus memakai event ECS "task stopped" sebagai sumber kebenaran kedua.

## 6. `summary.json` (report)

```json
{
  "schemaVersion": 1, "jobId": "rep_001", "layout": "PAIRED", "tool": "FastQC",
  "files": [ {
    "readType": "R1", "fileName": "S1_R1.fastq.gz", "sizeBytes": 1593311,
    "totalSequences": 25000, "sequenceLength": "150", "pctGC": 47.0,
    "pctSeqMeanQ30": 0.0, "poorQualityFlagged": 0, "encoding": "Sanger / Illumina 1.9",
    "adapterMaxPct": 0.0,
    "modules": { "Basic Statistics": "PASS", "Per base sequence quality": "FAIL", "...": "..." },
    "report": { "html": "qc/S1_R1_fastqc.html", "zip": "qc/S1_R1_fastqc.zip" }
  } ],
  "overall": { "modulesFail": ["Per base sequence quality"], "modulesWarn": [], "pairedCountMatch": true }
}
```

- `pctSeqMeanQ30`: persen **read** yang rata-rata kualitasnya ≥ Q30 (dari modul "Per sequence quality scores"). Ini bukan persen basis ≥ Q30.
- `adapterMaxPct`: nilai tertinggi di modul Adapter Content (semua adapter dan posisi).
- `pairedCountMatch`: jumlah read R1 dan R2 sama (`null` untuk SINGLE).
- Nama field adalah rancangan worker ini, bukan dari spec QC_Report eksternal. Sesuaikan bila ada schema resmi.

---

## 7. Variabel environment container

| Env | Default | Fungsi |
|---|---|---|
| `JOB_SPEC`, `JOB_SPEC_S3_URI` | — | Spesifikasi job (bagian 1) |
| `AWS_REGION` | `ap-southeast-3` | Region S3 (Jakarta). Kredensial dari **IAM task role**, jangan taruh kunci di env |
| `SCRATCH_DIR` | `/scratch` | Disk kerja sementara: input, work Nextflow, output. Harus bisa ditulis user container |
| `DISK_FACTOR` | `8` | Perkiraan kebutuhan disk = faktor × ukuran FASTQ. Job gagal cepat (`DISK_INSUFFICIENT`) bila ruang bebas lebih kecil |
| `SKIP_DISK_CHECK` | — | `true` menonaktifkan pengecekan di atas |
| `REF_BASE` | `/data/ref_genome` | Folder referensi lokal (`hg38.fa` + `.fai` + index BWA, dst.) |
| `REF_S3_PREFIX` | — | Bila referensi tidak ada di `REF_BASE`, unduh dari sini (mis. `s3://ngs-data/ref_genome/`). Nama file S3 sama dengan lokal |
| `THREADS` | otomatis | Thread per alat (maks 16) bila `resources.threads` kosong |
| `NGS_MAX_CPUS`, `NGS_MAX_MEMORY_MB` | otomatis | Paksa batas yang dilihat Nextflow (biasanya tidak perlu) |
| `S3_SSE`, `S3_KMS_KEY_ID` | — | `AES256` atau `aws:kms` untuk objek hasil |
| `S3_ENDPOINT_URL` | — | Hanya untuk uji lokal (MinIO/moto/LocalStack). **Kosongkan di AWS** |
| `LOG_LEVEL` | `INFO` | Log ke stdout → CloudWatch (`awslogs`) |

Referensi di S3 (`REF_S3_PREFIX`) harus berisi, untuk GRCh38: `hg38.fa`, `hg38.fa.fai`, `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` (GRCh37: `hg19.fa*`). Mengunduh hg38 + index (± 8–9 GB) di tiap task menambah waktu awal; alternatifnya menyimpannya di snapshot EBS/AMI dan mengisi `REF_BASE`.

## 8. Task definition ECS (panduan)

- **IAM task role** (minimum): `s3:GetObject` pada prefix input, `bed`, `specs`, dan `ref_genome`; `s3:PutObject` pada `output.prefix`. Tambahkan `s3:ListBucket` bila ingin pesan "tidak ditemukan" yang jelas (tanpa itu S3 membalas 403). Bila memakai SSE-KMS, beri izin `kms:Decrypt`/`kms:GenerateDataKey`.
- **Jaringan:** pakai S3 *gateway VPC endpoint* agar transfer tidak lewat NAT Gateway (biaya data processing). Bucket dan ECR di `ap-southeast-3`.
- **Disk:** `SCRATCH_DIR` butuh ruang besar (perkiraan awal 8 × ukuran FASTQ; ukur dengan `metrics.scratch` pada run nyata pertama). Pilihan: Fargate ephemeral storage (maks 200 GiB), atau volume EBS yang di-mount ke `/scratch` (ECS-managed EBS per task atau volume host). **Volume yang di-mount biasanya milik root**; container berjalan sebagai uid 1000 dan akan berhenti dengan pesan jelas bila `/scratch` tidak bisa ditulis. Atur `user` task ke `0` atau sesuaikan izin volume.
- **Ukuran awal (perkiraan, belum diukur):** report 2 vCPU/4 GB. VCF WES 16 vCPU/32 GB. VCF WGS ≥ 32 vCPU/64 GB dan disk berukuran TB.
- `stopTimeout`: 120 detik (worker menulis `status.json` dan menghentikan Nextflow saat SIGTERM; tanpa ini ECS memakai 30 detik).
- Tanpa Spot untuk percobaan pertama. Karena tiap job mandiri dan hasilnya idempoten (ditulis ulang ke prefix yang sama), retry setelah interupsi Spot aman, tetapi pekerjaan diulang dari awal.
- Log: `awslogs` driver ke CloudWatch.

## 9. Uji lokal tanpa AWS

```bash
# validasi spesifikasi
docker run --rm -i ngs-backend-bcftools:1.3.0 python -m worker vcf --validate-only --spec-file /dev/stdin < spec.json

# server S3 palsu (moto) di host
pip install "moto[server]" && moto_server -H 0.0.0.0 -p 5000

# jalankan worker terhadap server itu (host.docker.internal di Docker Desktop)
docker run --rm --no-healthcheck \
  -e S3_ENDPOINT_URL=http://host.docker.internal:5000 \
  -e AWS_ACCESS_KEY_ID=test -e AWS_SECRET_ACCESS_KEY=test \
  -e REF_S3_PREFIX=s3://ngs-data/ref_genome/ \
  -e JOB_SPEC_S3_URI=s3://ngs-data/specs/ana_001.json \
  ngs-backend-bcftools:1.3.0 python -m worker vcf
```

Yang sudah diuji dengan cara ini (data kecil, referensi chr22): VCF SE Targeted dengan BED dan referensi dari S3, VCF PE, report PE, serta jalur gagal (input hilang, spesifikasi tidak valid, trimming 0 read, disk tidak cukup, SIGTERM saat berjalan: exit 143, `status.json` `TERMINATED`).

## 10. Batasan yang diketahui

- Belum diuji pada S3 asli, IAM role, EBS per task, atau ECS di Jakarta. Semua itu perlu dicoba di akun AWS.
- Belum diuji pada data skala WES/WGS sungguhan. `DISK_FACTOR` dan ukuran compute di atas adalah perkiraan.
- FASTQ trimmed masih ditulis tanpa kompresi di work dir (sementara). Kompresinya akan menurunkan kebutuhan disk.
- Satu pasang FASTQ per job; belum multi-lane.
- Versi alat di image (apt Ubuntu 22.04): BWA 0.7.17, samtools 1.13, Trimmomatic 0.39, FastQC 0.11.9, bcftools 1.21 (build dari source).

---

# Bagian B: Worker per Langkah (wizard) — SE dan PE

Untuk wizard yang pengguna jalankan **satu langkah demi satu** (tinjau hasil, ubah parameter, ulangi), tiap langkah adalah task ECS sendiri.
Penyerahan antar langkah lewat S3, dan platform hanya perlu tahu `analysisPrefix` serta nama sampel. Yang didukung: `SINGLE` dan `PAIRED` (tanpa multi-lane).

| Step wizard | Command | Membaca | Menulis (di bawah `analysisPrefix`) |
|---|---|---|---|
| ① Upload & QC | `python -m worker report` | FASTQ mentah (`input.r1/r2`) | `qc/` |
| ② Trimming | `python -m worker trimming` | FASTQ mentah (`input.r1/r2`) | `trimmed/` |
| ③ Alignment | `python -m worker alignment` | `trimmed/` | `alignment/` |
| ④ Variant Calling | `python -m worker variant_calling` | `alignment/` | `variants/` |
| ⑤ Results | (tidak ada task) | — | Unduh dari S3 dengan URL bertanda tangan |

Step ① dan ② tidak saling bergantung dan boleh berjalan bersamaan. Step ③ butuh ② selesai, dan ④ butuh ③ selesai.
`analysisPrefix` berbentuk `s3://bucket/analyses/<id>/`.

## Susunan S3 per analisis

```
analyses/<id>/
  input/                    ← FASTQ + BED (lokasi unggahan platform; spesifikasi menyebut URI-nya langsung)
  qc/                       ← status.json, summary.json, <file>_fastqc.html, <file>_fastqc.zip
  trimmed/                  ← status.json, stats.json, nextflow.log,
                              SE: <sample>.trimmed.fastq.gz
                              PE: <sample>_R1.trimmed.fastq.gz, <sample>_R2.trimmed.fastq.gz
  alignment/                ← status.json, nextflow.log, <sample>.markdup.bam, .bai, <sample>.coverage.json, <sample>.markdup.metrics.txt
  variants/                 ← status.json, nextflow.log, <sample>.filtered.vcf.gz, .tbi
```

Nama berkas di atas adalah **kontrak antar langkah**: langkah berikutnya mencarinya persis dengan nama itu. Karena itu `sampleName` dan `layout` harus **sama di semua langkah** satu analisis
(pakai `metrics.sample` dari `status.json` trimming bila nama diturunkan otomatis).

## Spesifikasi tiap langkah

Semua langkah memakai `jobId` (ID analisis, aturan sama seperti di atas), `analysisPrefix`, dan `resources.threads` (opsional).

**① report** (sama seperti Bagian 3, tetapi hasil ke `qc/`):
```json
{ "jobId": "ana_001", "analysisPrefix": "s3://ngs-data/analyses/ana_001/",
  "input": { "layout": "PAIRED", "r1": "s3://ngs-data/analyses/ana_001/input/S1_R1.fastq.gz", "r2": "s3://ngs-data/analyses/ana_001/input/S1_R2.fastq.gz" } }
```
Untuk lokasi hasil bebas, pakai `output.prefix` sebagai pengganti `analysisPrefix` (isi salah satu).

**② trimming**
```json
{ "jobId": "ana_001", "analysisPrefix": "s3://ngs-data/analyses/ana_001/", "sampleName": "S1",
  "input":    { "layout": "PAIRED", "r1": "s3://.../S1_R1.fastq.gz", "r2": "s3://.../S1_R2.fastq.gz" },
  "trimming": { "preset": null, "quality": 30, "minLen": 36, "adapter": true, "adapterKit": "truseq_pe", "polyG": false } }
```
Field `trimming` sama seperti di `worker vcf`. `sampleName` boleh dikosongkan (diturunkan dari nama file R1). Hasil trimming selalu ter-gzip.

**③ alignment**
```json
{ "jobId": "ana_001", "analysisPrefix": "s3://ngs-data/analyses/ana_001/", "sampleName": "S1",
  "layout": "PAIRED", "reference": "GRCh38", "analysisMode": "WES", "bed": "s3://.../targets.bed",
  "alignment": { "markDuplicates": true } }
```
`sampleName` **wajib**. `bed` wajib untuk `WES`/`Targeted` dan dipakai untuk menghitung coverage di region target.

**④ variant_calling**
```json
{ "jobId": "ana_001", "analysisPrefix": "s3://ngs-data/analyses/ana_001/", "sampleName": "S1",
  "reference": "GRCh38", "analysisMode": "WES", "bed": "s3://.../targets.bed", "sampleSex": "male",
  "variantCalling": { "minQual": 20, "minDp": 10, "minGq": 20, "maxDepth": 500, "minBq": 20, "minMq": 20 } }
```
`sampleName` **wajib**. `reference`, `analysisMode`, dan `bed` harus sama dengan langkah ③.

## Status per langkah

Tiap langkah menulis `status.json` di foldernya (`qc/`, `trimmed/`, `alignment/`, `variants/`), dengan skema yang sama seperti Bagian 5 (`progress` 0–100 untuk langkah itu saja).
Pemantauan per langkah: baca `<analysisPrefix><folder>/status.json`, dan gunakan event ECS *task stopped* sebagai sumber kebenaran kedua.

`metrics` per langkah:

| Langkah | `metrics` |
|---|---|
| trimming | `sample`, `layout`, `trimming` (SE: `inputReads`, `surviving`, `survivingPct`, `dropped`; PE: `inputPairs`, `bothSurviving`, `bothSurvivingPct`, `forwardOnly`, `reverseOnly`, `dropped`) |
| alignment | `sample`, `flagstat` (`totalReads`, `mapped`, `mappedPct`, `duplicates`, `properlyPairedPct`), `coverage` |
| variant_calling | `sample`, `variants` (`total`, `pass`) |
| semua | `resources`, `scratch` (puncak pemakaian disk) |

## Perilaku menjalankan ulang

Saat sebuah langkah dimulai, worker **menghapus hasil lama langkah itu dan semua langkah di hilirnya** sebelum bekerja:

| Langkah dijalankan | Folder yang dibersihkan |
|---|---|
| trimming | `trimmed/`, `alignment/`, `variants/` |
| alignment | `alignment/`, `variants/` |
| variant_calling | `variants/` |

Dengan begitu run yang gagal tidak meninggalkan berkas basi yang terbaca langkah berikutnya. Konsekuensinya, jika langkah ② gagal pada percobaan ulang, hasil ③ dan ④ sebelumnya ikut hilang.
Folder `qc/` tidak disentuh oleh langkah lain. Task role perlu izin `s3:DeleteObject` pada `analyses/*` (sudah ada di `aws/task-role-policy.json`).

Jika input langkah belum ada (mis. alignment sebelum trimming), task berakhir `FAILED` dengan `errorCode: INPUT_NOT_FOUND` dan pesan yang menyebut langkah yang harus dijalankan lebih dulu; tidak ada yang dihapus.

## Ukuran task awal per langkah (perkiraan, belum diukur pada data besar)

| Langkah | Task definition | Fargate (uji chr22) | Catatan untuk WES/WGS |
|---|---|---|---|
| report | `ngs-worker-report` | 2 vCPU / 4 GB | Cukup untuk WES dan WGS (disk ≈ 1,5 × FASTQ) |
| trimming | `ngs-worker-trimming` | 4 vCPU / 8 GB | Disk ≈ 3 × FASTQ |
| alignment | `ngs-worker-alignment` | 4 vCPU / 16 GB | Langkah terberat: 16 vCPU / 32 GB, disk ≈ 6 × FASTQ trimmed (EC2) |
| variant calling | `ngs-worker-variant-calling` | 4 vCPU / 16 GB | 8–16 vCPU / 32 GB, disk ≈ 3 × BAM |

Faktor disk bawaan: trimming 3, alignment 6, variant calling 3 (`DISK_FACTOR` di env menimpa semuanya). Ukur dengan `metrics.scratch.peakOverInput` pada run nyata.
