# syntax=docker/dockerfile:1
#
# NGS Backend (BCFTools) — API FastAPI + pipeline Nextflow dalam satu image.
# Nextflow menjalankan semua proses secara lokal di dalam container (profil `local_tools`),
# jadi image ini sudah memuat semua tool. Referensi genom dan data job dipasang sebagai volume.
#
#   docker build -t ngs-backend-bcftools:1.1.1 .
#   docker compose up -d        # lihat docker-compose.yml
#
ARG UBUNTU=ubuntu:22.04

# ── Stage 1: build bcftools >= 1.21 (apt jammy hanya 1.13, tidak punya concat --rm-dups yang dibutuhkan) ──
FROM ${UBUNTU} AS bcftools-build
ARG BCFTOOLS_VERSION=1.21
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential ca-certificates wget bzip2 \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
 && rm -rf /var/lib/apt/lists/*
RUN wget -q -O /tmp/bcftools.tar.bz2 \
      "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
 && tar -xjf /tmp/bcftools.tar.bz2 -C /tmp \
 && cd "/tmp/bcftools-${BCFTOOLS_VERSION}" \
 && ./configure --prefix=/usr/local \
 && make -j"$(nproc)" \
 && make install

# ── Stage 2: image runtime ──────────────────────────────────────────────────────
FROM ${UBUNTU}
ARG NEXTFLOW_VERSION=24.10.4
ENV DEBIAN_FRONTEND=noninteractive

# Java 17 dipasang lebih dulu supaya fastqc/trimmomatic tidak menarik Java 11.
# procps: Nextflow memakai `ps` untuk report/trace. tini: PID 1 yang meneruskan SIGTERM dan membersihkan zombie.
RUN apt-get update \
 && apt-get install -y --no-install-recommends openjdk-17-jre-headless \
 && apt-get install -y --no-install-recommends \
      bwa samtools trimmomatic fastqc \
      python3 python3-venv python3-pip \
      tini curl ca-certificates procps \
      libcurl4 libssl3 libbz2-1.0 liblzma5 zlib1g \
 && rm -rf /var/lib/apt/lists/* \
 && ln -s "/usr/lib/jvm/java-17-openjdk-$(dpkg --print-architecture)" /opt/java \
 && /opt/java/bin/java -version

COPY --from=bcftools-build /usr/local/bin/bcftools /usr/local/bin/bcftools

# Pipeline memanggil ${BIN_DIR}/trimmomatic dan ${BIN_DIR}/bwa. Paket Debian hanya menyediakan
# TrimmomaticSE/PE, jadi dibuat pembungkus `trimmomatic` (sub-perintah SE/PE sama dengan versi conda).
RUN printf '#!/bin/sh\nexec /opt/java/bin/java -jar /usr/share/java/trimmomatic.jar "$@"\n' > /usr/local/bin/trimmomatic \
 && chmod 755 /usr/local/bin/trimmomatic \
 && ln -s /usr/bin/bwa /usr/local/bin/bwa

# Launcher Nextflow dari rilis yang sama. Framework-nya diunduh ke NXF_HOME saat build (lihat bawah),
# jadi container tidak perlu internet saat runtime.
RUN curl -fsSL -o /usr/local/bin/nextflow \
      "https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow" \
 && chmod 755 /usr/local/bin/nextflow

# Dependensi Python di virtualenv (layer ini ter-cache selama requirements.txt tidak berubah)
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt \
 && rm /tmp/requirements.txt

# User non-root. UID 1000 mengikuti user default WSL/Linux sehingga bind mount mudah dibaca/ditulis.
RUN useradd -u 1000 -m -s /bin/bash ngs \
 && mkdir -p /app /data/uploads /data/results /data/jobs /data/ref_genome /opt/nextflow /scratch \
 && chown -R ngs:ngs /app /data /opt/nextflow /scratch

ENV PATH="/opt/venv/bin:${PATH}" \
    PYTHONUNBUFFERED=1 \
    NXF_HOME=/opt/nextflow \
    NXF_DISABLE_CHECK_LATEST=true \
    NEXTFLOW_BIN=/usr/local/bin/nextflow \
    NEXTFLOW_JAVA_HOME=/opt/java \
    NEXTFLOW_VERSION=${NEXTFLOW_VERSION} \
    NEXTFLOW_PROFILE=local_tools \
    FASTQC_PATH=/usr/bin/fastqc \
    BIN_DIR=/usr/local/bin \
    SAMTOOLS_BIN=/usr/bin/samtools \
    BCFTOOLS_BIN=/usr/local/bin/bcftools \
    ADAPTER_FILE=/usr/share/trimmomatic/TruSeq3-PE-2.fa \
    ADAPTER_FILE_SE=/usr/share/trimmomatic/TruSeq3-SE.fa \
    ADAPTER_NEXTERA=/usr/share/trimmomatic/NexteraPE-PE.fa \
    REF_BASE=/data/ref_genome \
    UPLOAD_DIR=/data/uploads \
    RESULTS_DIR=/data/results \
    JOBS_DIR=/data/jobs \
    SCRATCH_DIR=/scratch \
    AWS_REGION=ap-southeast-3

USER ngs
WORKDIR /app

# Muat framework Nextflow sekali saat build (ke NXF_HOME) agar start pertama tidak menulis/mengunduh
RUN JAVA_HOME=/opt/java NXF_VER="${NEXTFLOW_VERSION}" nextflow -version

COPY --chown=ngs:ngs api      /app/api
COPY --chown=ngs:ngs pipeline /app/pipeline
COPY --chown=ngs:ngs worker   /app/worker

VOLUME ["/data"]
EXPOSE 8001

HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
  CMD python3 -c "import json,sys,urllib.request as u; d=json.load(u.urlopen('http://127.0.0.1:8001/api/health', timeout=5)); sys.exit(0 if d.get('status')=='ok' else 1)"

# Mode API (bawaan). Mode worker ECS: ganti command, mis. ["python", "-m", "worker", "vcf"]  (docs/ECS_JOB_SPEC.md)
# SATU worker saja: status job disimpan di memori proses. Jangan tambahkan --workers.
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8001"]
