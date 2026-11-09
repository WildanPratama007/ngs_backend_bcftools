#!/bin/bash
# ── NGS Backend Setup Script ──────────────────────────────────────────────────
# Jalankan: bash scripts/setup.sh

set -e

echo "🧬 NGS Backend Setup"
echo "===================="

# 1. Cek conda
if ! command -v conda &> /dev/null; then
    echo "❌ Conda tidak ditemukan. Install Miniconda dulu."
    exit 1
fi

# 2. Buat conda env
echo "📦 Membuat conda environment ngs_platform..."
conda create -n ngs_platform python=3.11 -y 2>/dev/null || echo "  ✅ Env sudah ada"

# 3. Aktivasi dan install tools
echo "🔧 Menginstall bioinformatika tools..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ngs_platform

conda install -c conda-forge mamba -y -q
mamba install -c bioconda -c conda-forge \
    nextflow bwa trimmomatic gatk4 -y -q

# 4. Install samtools & bcftools via apt
echo "🔧 Menginstall samtools & bcftools via apt..."
sudo apt-get update -q
sudo apt-get install -y -q samtools bcftools

# 5. Install Python dependencies
echo "🐍 Menginstall Python dependencies..."
pip install -r requirements.txt -q

# 6. Buat folder data
echo "📁 Membuat folder data..."
mkdir -p data/uploads data/results data/jobs

# 7. Copy .env
if [ ! -f .env ]; then
    cp .env.example .env
    echo "📝 File .env dibuat dari .env.example — edit sesuai kebutuhan"
fi

# 8. Verifikasi
echo ""
echo "✅ Verifikasi instalasi:"
echo "  Nextflow: $(nextflow -version 2>&1 | head -1)"
echo "  GATK:     $(gatk --version 2>&1 | grep GATK)"
echo "  BWA:      $(bwa 2>&1 | grep Version | head -1)"
echo "  Samtools: $(/usr/bin/samtools --version | head -1)"

echo ""
echo "🚀 Setup selesai! Jalankan API dengan:"
echo "   uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload"
