# Panduan AWS: dari Nol sampai Task Berjalan di ECS (Jakarta)

Panduan ini membuat bucket S3, repository ECR, role IAM, ECS cluster, dan dua task definition,
lalu menjalankan **uji pertama (FastQC → VCF chr22)** di Fargate. Setelah itu ada catatan untuk WES/WGS.

> **Status:** perintah di sini **belum pernah dijalankan di akun AWS** (saya tidak punya akses). Image dan worker sudah
> diuji lokal dengan S3 palsu. Bila ada perintah yang gagal, baca bagian *Troubleshooting* di akhir dan kirim pesan errornya.
>
> Semua perintah memakai **PowerShell** di Windows dan **AWS CLI v2**. Jalankan dari folder `ngs_backend_bcftools`.
> Biaya uji pertama kecil (menit Fargate + beberapa MB di S3), tetapi **hapus resource setelah selesai** (bagian 11).

| Langkah | Isi |
|---|---|
| 0 | Persiapan akun, CLI, variabel |
| 1 | Bucket S3 + unggah referensi dan data uji |
| 2 | ECR + push image |
| 3 | Role IAM |
| 4 | Log group dan jaringan |
| 5 | ECS cluster + task definition |
| 6–8 | Uji: report, lalu VCF, lalu baca hasil |
| 9–10 | Naik ke WES/WGS, dan pemantauan |
| 11 | Bersih-bersih |

---

## 0. Persiapan

1. **Aktifkan region Jakarta.** `ap-southeast-3` adalah region *opt-in*: konsol AWS → pojok kanan atas nama akun → *Account* → *AWS Regions* → aktifkan **Asia Pacific (Jakarta)**. Butuh beberapa menit. Sebelum ini, semua perintah di region itu gagal dengan error token tidak valid.
2. **Pasang AWS CLI v2** (`aws --version`) dan pastikan Docker Desktop berjalan.
3. **Login CLI.** Pilih salah satu:
   - Punya IAM Identity Center (SSO): `aws configure sso`, lalu `aws sso login`.
   - Hanya punya IAM user: konsol → IAM → Users → *Security credentials* → *Create access key*, lalu `aws configure` (isi Access Key ID, Secret, region `ap-southeast-3`, output `json`). **Jangan tempel kunci itu di chat atau commit.** Nonaktifkan setelah selesai. Jangan membuat access key untuk akun *root*.
   - Identitas ini perlu izin membuat IAM role, S3, ECR, ECS, EC2 (VPC), dan CloudWatch Logs. Untuk percobaan pertama, termudah memakai administrator, atau minta admin akun yang menyiapkannya.
4. **Set variabel** (dipakai di semua langkah berikutnya; ulangi bila membuka terminal baru):

```powershell
$REGION  = "ap-southeast-3"
$ACCOUNT = aws sts get-caller-identity --query Account --output text      # harus menampilkan 12 digit
$BUCKET  = "ngs-data-$ACCOUNT"                                            # nama bucket harus unik sedunia
$REPO    = "ngs-backend-bcftools"
$ECR     = "$ACCOUNT.dkr.ecr.$REGION.amazonaws.com"
$IMAGE   = "${ECR}/${REPO}:1.3.0"

# Mengisi placeholder __ACCOUNT__, __BUCKET__, __ECR__, __REGION__ di file template
function Render($src, $dst) {
  (Get-Content $src -Raw).Replace("__ACCOUNT__",$ACCOUNT).Replace("__BUCKET__",$BUCKET).Replace("__ECR__",$ECR).Replace("__REGION__",$REGION) |
    Set-Content $dst -Encoding ascii
}
$TMP = $env:TEMP
```

Jika `aws sts get-caller-identity` gagal, perbaiki login dulu sebelum lanjut.

---

## 1. Bucket S3, referensi, dan data uji

```powershell
# Bucket (di luar us-east-1 wajib LocationConstraint). Bucket baru otomatis terenkripsi (SSE-S3)
# dan Block Public Access aktif; cek di konsol S3 > bucket > Permissions.
aws s3api create-bucket --bucket $BUCKET --region $REGION --create-bucket-configuration LocationConstraint=$REGION

# Hapus sisa multipart upload yang tidak selesai setelah 7 hari (menghindari biaya tersembunyi)
aws s3api put-bucket-lifecycle-configuration --bucket $BUCKET --lifecycle-configuration file://aws/s3-lifecycle.json

# Referensi (untuk uji: chr22; untuk produksi: hg38 beserta semua index BWA-nya)
aws s3 cp ..\ref_genome\ s3://$BUCKET/ref_genome/ --recursive --exclude "*" --include "chr22.fa*" --region $REGION

# Data uji: satu FASTQ single-end
aws s3 cp ..\datas\test_samples\latihan_S1_L001_R1_001.fastq.gz s3://$BUCKET/in/ --region $REGION

aws s3 ls s3://$BUCKET/ --recursive --human-readable --region $REGION
```

Untuk produksi tambahkan referensi hg38 (± 8–9 GB): `--include "hg38.fa*"`. File referensi wajib lengkap:
`hg38.fa`, `.fai`, `.amb`, `.ann`, `.bwt`, `.pac`, `.sa`. (`.dict` tidak dibutuhkan.)

> Path `..\ref_genome` mengasumsikan struktur folder proyek Anda (`Fastt\ref_genome`). Sesuaikan jika berbeda.

---

## 2. ECR dan push image

Image `ngs-backend-bcftools:1.3.0` harus ada di Docker Anda (build lokal): `docker build -t ngs-backend-bcftools:1.3.0 .`

```powershell
aws ecr create-repository --repository-name $REPO --region $REGION --image-scanning-configuration scanOnPush=true

aws ecr get-login-password --region $REGION | docker login --username AWS --password-stdin $ECR

docker tag ngs-backend-bcftools:1.3.0 $IMAGE
docker push $IMAGE          # ± 350 MB terkompresi di ECR; push berikutnya hanya mengunggah layer yang berubah

aws ecr describe-images --repository-name $REPO --region $REGION --query "imageDetails[].imageTags"
```

> **`docker login` gagal di Docker Desktop Windows** dengan pesan `error saving credentials ... The stub received bad data`
> (token ECR terlalu besar untuk credential store Windows), lalu `docker push` berkata `authorization token has expired`.
> Solusi: pakai folder konfigurasi Docker sementara yang berisi token langsung, tanpa credential store, lalu hapus setelah push:
>
> ```powershell
> $cfg = "$env:TEMP\dockercfg"; New-Item -ItemType Directory -Force $cfg | Out-Null
> $tok = aws ecr get-login-password --region $REGION
> $b64 = [Convert]::ToBase64String([Text.Encoding]::ASCII.GetBytes("AWS:$tok"))
> [IO.File]::WriteAllText("$cfg\config.json", '{"auths":{"' + $ECR + '":{"auth":"' + $b64 + '"}}}', (New-Object Text.UTF8Encoding $false))
> Remove-Variable tok,b64
> docker --config $cfg push $IMAGE
> Remove-Item $cfg -Recurse -Force
> ```
>
> Jangan menempelkan token ke chat atau menyimpannya di repo.

**Setelah push: daftarkan revisi task definition baru.** Task definition menunjuk tag image, jadi tag baru tidak otomatis dipakai.
Ambil revisi yang aktif, ganti `image`, lalu daftarkan sebagai revisi baru (role, CPU/RAM, log, dan `stopTimeout` ikut tersalin):

```powershell
foreach ($fam in "ngs-worker-report","ngs-worker-trimming","ngs-worker-alignment","ngs-worker-variant-calling","ngs-worker-vcf") {
  $d = (aws ecs describe-task-definition --task-definition $fam --region $REGION --query taskDefinition --output json) | ConvertFrom-Json
  foreach ($k in "taskDefinitionArn","revision","status","requiresAttributes","compatibilities","registeredAt","registeredBy","deregisteredAt") { $d.PSObject.Properties.Remove($k) }
  $d.containerDefinitions[0].image = $IMAGE
  [IO.File]::WriteAllText("$TMP\$fam.json", ($d | ConvertTo-Json -Depth 20), (New-Object Text.UTF8Encoding $false))
  aws ecs register-task-definition --cli-input-json "file://$TMP\$fam.json" --region $REGION --query "taskDefinition.taskDefinitionArn" --output text
}
```

`run-task` tanpa nomor revisi memakai revisi terbaru. Revisi dan tag lama sebaiknya dipertahankan satu versi sebagai cadangan rollback, lalu dihapus setelah versi baru terbukti jalan.

Image dibuild untuk **linux/amd64**. Jika Anda membuild di mesin ARM (mis. Mac M1) gunakan `docker build --platform linux/amd64 ...`.

---

## 3. Role IAM

Ada dua role untuk task. Task role memberi **kode worker** izin ke S3 (inilah pengganti access key). Execution role dipakai ECS untuk menarik image dan menulis log.

```powershell
# 3a. Task role + izin S3 (hanya ke bucket ini)
aws iam create-role --role-name ngs-worker-task-role --assume-role-policy-document file://aws/trust-ecs-tasks.json
Render aws\task-role-policy.json "$TMP\task-role-policy.json"
aws iam put-role-policy --role-name ngs-worker-task-role --policy-name ngs-worker-s3 --policy-document file://$TMP\task-role-policy.json

# 3b. Execution role
aws iam create-role --role-name ngs-ecs-execution-role --assume-role-policy-document file://aws/trust-ecs-tasks.json
aws iam attach-role-policy --role-name ngs-ecs-execution-role --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
```

Bila memakai SSE-KMS dengan key sendiri, tambahkan `kms:Decrypt` dan `kms:GenerateDataKey` pada key itu ke task role.

**Izin untuk yang menjalankan task** (Anda, Lambda, atau Step Functions) selain `ecs:RunTask`, `ecs:DescribeTasks`, `ecs:StopTask`:
`iam:PassRole` untuk kedua role di atas dengan kondisi `iam:PassedToService = ecs-tasks.amazonaws.com`. Jika Anda administrator, ini sudah tercakup.

---

## 4. Log group dan jaringan

```powershell
aws logs create-log-group --log-group-name /ecs/ngs-worker --region $REGION
aws logs put-retention-policy --log-group-name /ecs/ngs-worker --retention-in-days 30 --region $REGION

# Jaringan: untuk uji pertama pakai VPC bawaan (subnet publik). Untuk produksi pakai subnet privat.
$VPC     = aws ec2 describe-vpcs --filters Name=isDefault,Values=true --query "Vpcs[0].VpcId" --output text --region $REGION
$SUBNET  = (aws ec2 describe-subnets --filters Name=vpc-id,Values=$VPC --query "Subnets[0].SubnetId" --output text --region $REGION)
$SG      = aws ec2 create-security-group --group-name ngs-worker-sg --description "ngs worker egress only" --vpc-id $VPC --query GroupId --output text --region $REGION
$RTB     = aws ec2 describe-route-tables --filters Name=vpc-id,Values=$VPC --query "RouteTables[0].RouteTableId" --output text --region $REGION

# S3 gateway endpoint: transfer ke S3 tidak lewat internet/NAT (gratis)
aws ec2 create-vpc-endpoint --vpc-id $VPC --vpc-endpoint-type Gateway --service-name "com.amazonaws.$REGION.s3" --route-table-ids $RTB --region $REGION

"VPC=$VPC SUBNET=$SUBNET SG=$SG"
```

Security group baru tidak punya aturan *inbound* dan mengizinkan semua *outbound*. Itu yang dibutuhkan task (tidak ada server yang menerima koneksi).
Jika `$VPC` bernilai `None`, akun Anda tidak punya VPC bawaan di region ini; buat VPC lewat konsol (*VPC → Create VPC → VPC and more*).

---

## 5. ECS cluster dan task definition

```powershell
aws ecs create-cluster --cluster-name ngs --region $REGION

Render aws\task-def-report.json "$TMP\td-report.json"
Render aws\task-def-vcf.json    "$TMP\td-vcf.json"
aws ecs register-task-definition --cli-input-json file://$TMP\td-report.json --region $REGION
aws ecs register-task-definition --cli-input-json file://$TMP\td-vcf.json    --region $REGION
```

Isi template (`aws/task-def-*.json`):

| | `ngs-worker-report` | `ngs-worker-vcf` (uji) |
|---|---|---|
| Compute | Fargate 2 vCPU / 4 GB | Fargate 4 vCPU / 16 GB |
| Disk ephemeral | 50 GiB | 50 GiB (maks Fargate 200 GiB) |
| Command | `python -m worker report` | `python -m worker vcf` |
| `stopTimeout` | 120 detik | 120 detik |

Ukuran VCF di atas hanya untuk **uji chr22**. Untuk WES/WGS lihat bagian 9.

---

## 6. Uji 1: report (FastQC)

```powershell
Render aws\spec-report.example.json "$TMP\rep_001.json"
aws s3 cp "$TMP\rep_001.json" s3://$BUCKET/specs/rep_001.json --region $REGION

$ov = @{ containerOverrides = @(@{ name = "worker"; environment = @(@{ name = "JOB_SPEC_S3_URI"; value = "s3://$BUCKET/specs/rep_001.json" }) }) } | ConvertTo-Json -Depth 6 -Compress
Set-Content "$TMP\ov-report.json" $ov -Encoding ascii

$TASK = aws ecs run-task --cluster ngs --launch-type FARGATE --task-definition ngs-worker-report --count 1 `
  --network-configuration "awsvpcConfiguration={subnets=[$SUBNET],securityGroups=[$SG],assignPublicIp=ENABLED}" `
  --overrides file://$TMP\ov-report.json --region $REGION --query "tasks[0].taskArn" --output text
$TASK
```

`assignPublicIp=ENABLED` dibutuhkan agar task di subnet publik bisa menarik image dari ECR. (Di subnet privat dengan VPC endpoint tidak perlu.)

Pantau:

```powershell
aws ecs describe-tasks --cluster ngs --tasks $TASK --region $REGION --query "tasks[0].[lastStatus,stopCode,stoppedReason,containers[0].exitCode]" --output text
aws logs tail /ecs/ngs-worker --follow --region $REGION          # Ctrl+C untuk berhenti
```

Status melalui `PROVISIONING → PENDING → RUNNING → STOPPED`. Exit code `0` = sukses. Lalu baca hasilnya:

```powershell
aws s3 ls s3://$BUCKET/reports/rep_001/ --region $REGION
aws s3 cp s3://$BUCKET/reports/rep_001/status.json - --region $REGION
aws s3 cp s3://$BUCKET/reports/rep_001/summary.json - --region $REGION
```

`status.json` harus berisi `"status": "SUCCEEDED"`. Bila `FAILED`, lihat `errorCode` dan `errorMessage` (tabel kode di `docs/ECS_JOB_SPEC.md`).

---

## 7. Uji 2: FASTQ → VCF (chr22)

```powershell
Render aws\spec-vcf-chr22.example.json "$TMP\vcf_001.json"
aws s3 cp "$TMP\vcf_001.json" s3://$BUCKET/specs/vcf_001.json --region $REGION

$ov = @{ containerOverrides = @(@{ name = "worker"; environment = @(@{ name = "JOB_SPEC_S3_URI"; value = "s3://$BUCKET/specs/vcf_001.json" }) }) } | ConvertTo-Json -Depth 6 -Compress
Set-Content "$TMP\ov-vcf.json" $ov -Encoding ascii

$TASK = aws ecs run-task --cluster ngs --launch-type FARGATE --task-definition ngs-worker-vcf --count 1 `
  --network-configuration "awsvpcConfiguration={subnets=[$SUBNET],securityGroups=[$SG],assignPublicIp=ENABLED}" `
  --overrides file://$TMP\ov-vcf.json --region $REGION --query "tasks[0].taskArn" --output text

aws logs tail /ecs/ngs-worker --follow --region $REGION
```

Spesifikasi ini memakai `trimming.quality = 15` karena data uji sintetis; kualitas 20 ke atas membuang semua read pada file ini.
Perkiraan waktu: beberapa menit (termasuk mengunduh referensi chr22 ± 140 MB dari S3).

---

## 8. Membaca hasil

```powershell
aws s3 ls s3://$BUCKET/runs/vcf_001/ --recursive --human-readable --region $REGION
aws s3 cp s3://$BUCKET/runs/vcf_001/status.json - --region $REGION
aws s3 cp s3://$BUCKET/runs/vcf_001/vcf/ .\hasil_vcf\ --recursive --region $REGION
```

Yang diharapkan: `vcf/*.filtered.vcf.gz` + `.tbi`, `bam/*.markdup.bam` + `.bai`, `qc/` (coverage, FastQC), `logs/nextflow.log`, dan `status.json`
dengan `metrics.variants` (sekitar 7 varian untuk data uji ini) dan `metrics.scratch.peakOverFastq`.

---

## 9. Naik ke WES / WGS

1. **Referensi:** unggah `hg38.fa*` ke `s3://$BUCKET/ref_genome/` (bagian 1). Spec memakai `"reference": "GRCh38"`.
2. **Ganti ke EC2 untuk disk besar.** Fargate maksimal 200 GiB ephemeral; WES dan WGS biasanya butuh lebih (perkiraan awal 8 × ukuran FASTQ, ukur dengan `metrics.scratch`). Cara termudah lewat konsol: *ECS → Clusters → Create cluster → Infrastructure: Amazon EC2 instances*, pilih instance (WES: `c7i.4xlarge` 16 vCPU/32 GB, sesuai `AWS_COST_ANALYSIS.md`), ukuran root volume EBS (WES 300–500 GB, gp3), dan biarkan konsol membuat Auto Scaling Group serta instance role.
3. **Task definition EC2:** salin template, ubah `requiresCompatibilities` ke `["EC2"]`, hapus `ephemeralStorage` dan `runtimePlatform`, tambahkan `cpu`/`memory` sesuai instance (sisakan 4–6 GB RAM untuk OS), dan bila perlu `"user": "0"` bila `/scratch` dipasang dari volume host milik root. Jalankan dengan `--launch-type EC2` (atau capacity provider cluster), tanpa `assignPublicIp` bila di subnet privat.
4. **Jalankan satu sampel nyata dulu**, baca `metrics.scratch.peakOverFastq` di `status.json`, lalu setel `DISK_FACTOR` dan ukuran volume.
5. Spot hanya untuk task (bukan server), setelah satu run nyata sukses; interupsi berarti job diulang dari awal.

---

## 10. Pemantauan dan otomatisasi (opsional)

- **Event task selesai:** buat aturan EventBridge untuk `ECS Task State Change` dengan `lastStatus = STOPPED` di cluster `ngs`, target Lambda/SQS. Gunakan ini bersama `status.json` (event tetap datang walaupun container dimatikan paksa dan `status.json` tertinggal `RUNNING`).
- **Anggaran:** AWS Budgets dengan peringatan, mis. USD 20/bulan, supaya task yang lupa dihentikan tidak membengkak.
- **Pemicu otomatis** (S3 upload → RunTask) butuh Lambda/Step Functions; itu di luar dua container ini.

---

## 11. Bersih-bersih

```powershell
aws ecs delete-cluster --cluster ngs --region $REGION                              # setelah semua task STOPPED dan task definition dideregister
aws ecs deregister-task-definition --task-definition ngs-worker-report:1 --region $REGION
aws ecs deregister-task-definition --task-definition ngs-worker-vcf:1    --region $REGION
aws ecr delete-repository --repository-name $REPO --force --region $REGION
aws logs delete-log-group --log-group-name /ecs/ngs-worker --region $REGION
aws s3 rb s3://$BUCKET --force --region $REGION                                    # menghapus bucket dan SEMUA isinya
aws iam delete-role-policy --role-name ngs-worker-task-role --policy-name ngs-worker-s3
aws iam delete-role --role-name ngs-worker-task-role
aws iam detach-role-policy --role-name ngs-ecs-execution-role --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
aws iam delete-role --role-name ngs-ecs-execution-role
aws ec2 delete-security-group --group-id $SG --region $REGION
```

VPC endpoint S3 gratis; hapus lewat konsol bila tidak dipakai. Jika memakai access key, **nonaktifkan dan hapus** di IAM.

---

## Troubleshooting

| Gejala | Penyebab umum | Perbaikan |
|---|---|---|
| `The security token included in the request is invalid` / `InvalidClientTokenId` di Jakarta | Region belum diaktifkan, atau kredensial salah | Langkah 0.1; ulangi `aws sts get-caller-identity` |
| `AccessDenied` saat `run-task` / `iam:PassRole` | Identitas Anda tidak boleh meneruskan role | Beri `iam:PassRole` untuk kedua role |
| Task `STOPPED` dengan `CannotPullContainerError` | Image tidak ada di ECR, atau task tidak punya jalur internet/endpoint | Cek `$IMAGE` sama dengan di task definition; pakai `assignPublicIp=ENABLED` di subnet publik, atau buat VPC endpoint ECR |
| `ResourceInitializationError` | Execution role belum ada/izin kurang, atau tidak bisa menjangkau ECR/CloudWatch | Langkah 3b; cek jaringan |
| `status.json`: `INPUT_NOT_FOUND` | Key salah, atau task role tidak boleh membaca | Cek `aws s3 ls`; cek policy task role dan nama bucket di spec |
| `UPLOAD_FAILED` | Task role tidak punya `PutObject`/`GetObject` di prefix output | Langkah 3a (policy sudah mencakup `runs/*` dan `reports/*`) |
| `DISK_INSUFFICIENT` | Ephemeral storage terlalu kecil | Naikkan `ephemeralStorage` (maks 200 GiB) atau pindah ke EC2; atau setel `DISK_FACTOR` |
| `PIPELINE_FAILED` + trimming 0 read | `trimming.quality` terlalu tinggi untuk data | Turunkan (data uji ini: 15) |
| Exit code 137, tanpa `status.json` baru | Kehabisan memori (OOM) | Naikkan memori task |
| `Invalid JSON` dari AWS CLI untuk `file://` | File ber-BOM | Gunakan `Set-Content -Encoding ascii` seperti di fungsi `Render` |

---

## Tambahan: task definition per langkah (wizard)

Selain `ngs-worker-report` dan `ngs-worker-vcf`, ada tiga task definition untuk wizard per langkah. Cara mendaftarkannya sama seperti langkah 5 (gunakan `Render`):

```powershell
foreach ($f in "task-def-step-trimming","task-def-step-alignment","task-def-step-variant") {
  Render "aws\$f.json" "$TMP\$f.json"
  aws ecs register-task-definition --cli-input-json file://$TMP\$f.json --region $REGION
}
```

| Template | Family | Command | Ukuran awal (Fargate, uji chr22) |
|---|---|---|---|
| `task-def-step-trimming.json` | `ngs-worker-trimming` | `python -m worker trimming` | 4 vCPU / 8 GB |
| `task-def-step-alignment.json` | `ngs-worker-alignment` | `python -m worker alignment` | 4 vCPU / 16 GB |
| `task-def-step-variant.json` | `ngs-worker-variant-calling` | `python -m worker variant_calling` | 4 vCPU / 16 GB |

Task role harus memuat izin untuk folder `analyses/*` (baca, tulis, **hapus**). Template `aws/task-role-policy.json` sudah memuatnya; terapkan ulang dengan `aws iam put-role-policy`
(langkah 3a) bila role dibuat sebelum folder `analyses/` ada. Spesifikasi tiap langkah, susunan S3, dan perilaku jalan-ulang ada di `docs/ECS_JOB_SPEC.md` (Bagian B).
