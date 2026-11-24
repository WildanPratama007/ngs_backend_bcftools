#!/bin/bash
# ── NGS Backend Setup Script (BCFTools) ───────────────────────────────────────
# Jalankan dari folder ngs_backend_bcftools:  bash scripts/setup.sh
# Target: Ubuntu 22.04 / WSL2 dengan conda. Butuh sudo untuk apt dan make install.

set -e

BCFTOOLS_MIN="1.21"          # concat --rm-dups butuh -a dan versi baru (apt jammy hanya 1.13)
BCFTOOLS_VERSION="1.21"      # versi yang dibuild bila yang terpasang terlalu lama

echo "🧬 NGS Backend Setup (BCFTools)"
echo "==============================="

# 1. Cek conda
if ! command -v conda &> /dev/null; then
    echo "❌ Conda tidak ditemukan. Install Miniconda dulu."
    exit 1
fi

# 2. Buat conda env
echo "📦 Membuat conda environment ngs_platform..."
conda create -n ngs_platform python=3.11 -y 2>/dev/null || echo "  ✅ Env sudah ada"

# 3. Tools bioinformatika dari conda (bcftools dan samtools sengaja TIDAK dari conda)
echo "🔧 Menginstall Nextflow, BWA, Trimmomatic, FastQC..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ngs_platform

conda install -c conda-forge mamba -y -q
mamba install -c bioconda -c conda-forge \
    nextflow bwa trimmomatic fastqc -y -q

# 4. Paket sistem: Java 17 (dibutuhkan Nextflow 24.x), samtools, dan dependensi build htslib
echo "🔧 Menginstall paket sistem via apt..."
sudo apt-get update -q
sudo apt-get install -y -q \
    openjdk-17-jre-headless samtools \
    build-essential wget bzip2 zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev

# 5. BCFTools >= 1.21 dari source (versi apt 1.13 membuat pipeline gagal)
have_bcftools_ok() {
    command -v /usr/local/bin/bcftools &> /dev/null || return 1
    local v
    v=$(/usr/local/bin/bcftools --version | head -1 | awk '{print $2}')
    [ "$(printf '%s\n%s\n' "$BCFTOOLS_MIN" "$v" | sort -V | head -1)" = "$BCFTOOLS_MIN" ]
}

if have_bcftools_ok; then
    echo "  ✅ bcftools $(/usr/local/bin/bcftools --version | head -1 | awk '{print $2}') sudah terpasang"
else
    echo "🔨 Build bcftools ${BCFTOOLS_VERSION} dari source..."
    BUILD_DIR=$(mktemp -d)
    wget -q -O "$BUILD_DIR/bcftools.tar.bz2" \
        "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"
    tar -xjf "$BUILD_DIR/bcftools.tar.bz2" -C "$BUILD_DIR"
    (cd "$BUILD_DIR/bcftools-${BCFTOOLS_VERSION}" && ./configure -q && make -j"$(nproc)" -s && sudo make install -s)
    rm -rf "$BUILD_DIR"
fi

# 6. Python dependencies
echo "🐍 Menginstall Python dependencies..."
pip install -r requirements.txt -q

# 7. Buat folder data
echo "📁 Membuat folder data..."
mkdir -p data/uploads data/results data/jobs

# 8. Copy .env
if [ ! -f .env ]; then
    cp .env.example .env
    echo "📝 File .env dibuat dari .env.example — edit path sesuai mesin Anda"
fi

# 9. Verifikasi
echo ""
echo "✅ Verifikasi instalasi:"
echo "  Nextflow:    $(nextflow -version 2>&1 | grep -m1 version || echo 'tidak terbaca')"
echo "  Java:        $(/usr/lib/jvm/java-17-openjdk-amd64/bin/java -version 2>&1 | head -1)"
echo "  BWA:         $(bwa 2>&1 | grep Version | head -1)"
echo "  Trimmomatic: $(trimmomatic -version 2>&1 | tail -1)"
echo "  FastQC:      $(fastqc --version 2>&1 | head -1)"
echo "  Samtools:    $(/usr/bin/samtools --version | head -1)"
echo "  BCFTools:    $(/usr/local/bin/bcftools --version | head -1)"

echo ""
echo "ℹ️  Isi di .env sesuai mesin Anda:"
echo "   FASTQC_PATH=$(command -v fastqc)"
echo "   BIN_DIR=$(dirname "$(command -v bwa)")"
echo "   REF_BASE=<folder referensi genom>   (lihat README: Persiapan Referensi Genom)"
echo ""
echo "🚀 Setup selesai! Jalankan API dengan:"
echo "   conda activate ngs_platform"
echo "   uvicorn api.main:app --host 0.0.0.0 --port 8001 --reload"
echo "   Lalu cek: curl http://localhost:8001/api/health   (status harus \"ok\")"
