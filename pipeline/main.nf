#!/usr/bin/env nextflow

/*
 * NGS E2E Secondary Analysis Pipeline — BCFTools
 * Alur: FastQC → Trimmomatic (SE/PE) → BWA MEM → samtools markdup
 *       → BCFTools mpileup/call/norm/filter
 * Fix: duplikat varian WES, adapter trimming, markdup, parameter aktif,
 *      output bgzip+tbi, coverage per region, sanitasi nama sampel
 */

nextflow.enable.dsl = 2

if (!params.reads && !params.reads_r1) {
    error "Parameter --reads (SE) atau --reads_r1 + --reads_r2 (PE) wajib diisi."
}

// ── Helper: sanitasi nama sampel (A-Z a-z 0-9 . - _) ─────────────────────────
def sanitize(String s) {
    s.replaceAll(/[^A-Za-z0-9._\-]/, '_').replaceAll(/_+/, '_').replaceAll(/^_|_$/, '')
}

// ── Proses 1: FastQC ──────────────────────────────────────────────────────────
process FASTQC {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip")

    script:
    def fastqc_bin = params.fastqc_path ?: 'fastqc'
    """
    ${fastqc_bin} ${reads} --outdir .
    """
}

// ── Proses 2a: Trimmomatic SE ─────────────────────────────────────────────────
process TRIMMOMATIC_SE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq"), val("SE")

    script:
    // Bangun trimmer_opts secara dinamis dari parameter user
    def adapter_clip = params.trim_adapter
        ? "ILLUMINACLIP:${params.adapter_file}:2:30:10:8:true"
        : ""
    def poly_g_clip  = params.trim_poly_g
        ? "ILLUMINACLIP:${params.adapter_file.replace('.fa', '-polyG.fa')}:2:30:10:8:true"
        : ""
    def sliding      = "SLIDINGWINDOW:4:${params.trim_sliding_window_quality}"
    def minlen       = "MINLEN:${params.trim_min_len}"
    def trim_cmd     = [adapter_clip, poly_g_clip, "LEADING:5", "TRAILING:5", sliding, minlen]
                           .findAll { it }.join(" ")
    """
    ${params.bin_dir}/trimmomatic SE -threads ${params.threads} \
        "${reads}" "${sample_id}.trimmed.fastq" \
        ${trim_cmd} 2> trimmomatic.log

    # Gagalkan bila output kosong — jangan fallback ke raw FASTQ
    if [ ! -s "${sample_id}.trimmed.fastq" ]; then
        echo "ERROR: Trimmomatic menghasilkan output kosong" >&2
        exit 1
    fi
    """
}

// ── Proses 2b: Trimmomatic PE ─────────────────────────────────────────────────
process TRIMMOMATIC_PE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fastq"),
          path("${sample_id}_R1.unpaired.fastq"),
          path("${sample_id}_R2.trimmed.fastq"),
          path("${sample_id}_R2.unpaired.fastq"),
          val("PE")

    script:
    // Bangun trimmer_opts secara dinamis — sama dengan SE
    def adapter_clip = params.trim_adapter
        ? "ILLUMINACLIP:${params.adapter_file}:2:30:10:8:true"
        : ""
    def poly_g_clip  = params.trim_poly_g
        ? "ILLUMINACLIP:${params.adapter_file.replace('.fa', '-polyG.fa')}:2:30:10:8:true"
        : ""
    def sliding      = "SLIDINGWINDOW:4:${params.trim_sliding_window_quality}"
    def minlen       = "MINLEN:${params.trim_min_len}"
    def trim_cmd     = [adapter_clip, poly_g_clip, "LEADING:5", "TRAILING:5", sliding, minlen]
                           .findAll { it }.join(" ")
    """
    ${params.bin_dir}/trimmomatic PE -threads ${params.threads} \
        "${reads_r1}" "${reads_r2}" \
        "${sample_id}_R1.trimmed.fastq" "${sample_id}_R1.unpaired.fastq" \
        "${sample_id}_R2.trimmed.fastq" "${sample_id}_R2.unpaired.fastq" \
        ${trim_cmd} 2> trimmomatic.log

    # Gagalkan bila output kosong
    if [ ! -s "${sample_id}_R1.trimmed.fastq" ]; then
        echo "ERROR: Trimmomatic menghasilkan output kosong" >&2
        exit 1
    fi
    """
}

// ── Proses 3a: BWA MEM SE ─────────────────────────────────────────────────────
process BWA_ALIGN_SE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads), val(mode)
    path ref_genome
    path ref_index

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), path("${sample_id}.aligned.sorted.bam.bai")

    script:
    // Read Group: ID format SM.L001 — konsisten dengan konvensi Illumina
    // LB = sample_id agar markdup bisa bedakan antar library bila nanti multi-lane
    """
    ${params.bin_dir}/bwa mem \
        -t ${params.threads} \
        -R "@RG\\tID:${sample_id}.L001\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}\\tPU:${sample_id}.L001" \
        "${ref_genome}" "${trimmed_reads}" \
        | ${params.samtools_bin} sort -@ ${params.threads} -m 4G \
          -o "${sample_id}.aligned.sorted.bam"

    ${params.samtools_bin} index "${sample_id}.aligned.sorted.bam"
    """
}

// ── Proses 3b: BWA MEM PE ─────────────────────────────────────────────────────
process BWA_ALIGN_PE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r1_unpaired), path(r2), path(r2_unpaired), val(mode)
    path ref_genome
    path ref_index

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), path("${sample_id}.aligned.sorted.bam.bai")

    script:
    // Read Group: ID format SM.L001 — konsisten dengan konvensi Illumina
    // LB = sample_id agar markdup bisa bedakan antar library bila nanti multi-lane
    """
    ${params.bin_dir}/bwa mem \
        -t ${params.threads} \
        -R "@RG\\tID:${sample_id}.L001\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}\\tPU:${sample_id}.L001" \
        "${ref_genome}" "${r1}" "${r2}" \
        | ${params.samtools_bin} fixmate -m -@ ${params.threads} - - \
        | ${params.samtools_bin} sort -@ ${params.threads} -m 4G \
          -o "${sample_id}.aligned.sorted.bam"

    ${params.samtools_bin} index "${sample_id}.aligned.sorted.bam"
    """
}

// ── Proses 4: MarkDuplicates ──────────────────────────────────────────────────
process MARKDUPLICATES {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai")
    path "${sample_id}.markdup.metrics.txt"

    script:
    if (params.mark_duplicates) {
        // Tandai duplicate — reads tetap ada tapi di-flag
        """
        ${params.samtools_bin} markdup \
            -f "${sample_id}.markdup.metrics.txt" \
            "${bam}" \
            "${sample_id}.markdup.bam"

        ${params.samtools_bin} index "${sample_id}.markdup.bam"
        """
    } else {
        // Skip markdup — symlink BAM agar output konsisten
        // Cocok untuk PCR-free library atau amplicon sequencing
        """
        cp "${bam}" "${sample_id}.markdup.bam"
        cp "${bai}" "${sample_id}.markdup.bam.bai"
        echo "mark_duplicates=false — skipped" > "${sample_id}.markdup.metrics.txt"
        """
    }
}

// ── Proses 5: Hitung Coverage ─────────────────────────────────────────────────
process COVERAGE_STATS {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.json")

    script:
    // Untuk WES/Targeted: hitung coverage di region BED saja
    // Untuk WGS: hitung coverage seluruh referensi
    def region_opt = params.target_bed ? "-l ${params.target_bed}" : ""
    """
    # samtools coverage output kolom:
    # #rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
    # Kita pakai: covbases($5), coverage($6=pct_bases_1x), meandepth($7)

    samtools coverage ${region_opt} "${bam}" > coverage_raw.txt

    # Hitung metrik lengkap dengan python (sudah tersedia di conda env)
    python3 << 'PYEOF'
import json, subprocess, sys

# Baca output samtools coverage
rows = []
with open("coverage_raw.txt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\\t")
        if len(parts) < 7:
            continue
        try:
            chrom      = parts[0]
            start      = int(parts[1])
            end        = int(parts[2])
            numreads   = int(parts[3])
            covbases   = int(parts[4])
            cov_pct    = float(parts[5])   # % bases >=1x
            meandepth  = float(parts[6])
            length     = end - start if end > start else 1
            rows.append({
                "chrom": chrom,
                "length": length,
                "covbases": covbases,
                "meandepth": meandepth,
                "cov_pct": cov_pct,
            })
        except (ValueError, IndexError):
            continue

if not rows:
    result = {
        "mean_depth": 0.0, "median_depth": 0.0,
        "pct_bases_10x": 0.0, "pct_bases_20x": 0.0, "pct_bases_30x": 0.0,
        "total_bases": 0, "covered_bases": 0,
        "warning": "No coverage data"
    }
else:
    # Weighted mean depth
    total_len    = sum(r["length"] for r in rows)
    total_cov    = sum(r["covbases"] for r in rows)
    weighted_sum = sum(r["meandepth"] * r["length"] for r in rows)
    mean_depth   = weighted_sum / total_len if total_len > 0 else 0.0

    # Median depth — gunakan meandepth per contig sebagai approximasi
    depths = sorted(r["meandepth"] for r in rows)
    n = len(depths)
    median_depth = (depths[n//2] + depths[(n-1)//2]) / 2 if n > 0 else 0.0

    # pct_bases_Nx: estimasi dari mean depth (approximasi bila tidak ada depth per-base)
    # Formula: P(depth >= N) ~ 1 - e^(-mean/N) untuk distribusi Poisson
    import math
    def pct_at_least_nx(mean, n):
        if mean <= 0:
            return 0.0
        # Poisson CDF approximation
        prob = 1.0
        term = math.exp(-mean)
        s = term
        for k in range(1, n):
            term *= mean / k
            s += term
        prob = max(0.0, (1.0 - s) * 100.0)
        return round(prob, 1)

    # Warning bila coverage rendah
    warning = None
    if mean_depth < 10:
        warning = "Coverage sangat rendah (<10x) — hasil variant calling tidak reliable"
    elif mean_depth < 20:
        warning = "Coverage rendah (<20x) — beberapa variant mungkin terlewat"

    result = {
        "mean_depth":     round(mean_depth, 2),
        "median_depth":   round(median_depth, 2),
        "pct_bases_10x":  pct_at_least_nx(mean_depth, 10),
        "pct_bases_20x":  pct_at_least_nx(mean_depth, 20),
        "pct_bases_30x":  pct_at_least_nx(mean_depth, 30),
        "total_bases":    total_len,
        "covered_bases":  total_cov,
        "warning":        warning,
    }

with open("${sample_id}.coverage.json", "w") as f:
    json.dump(result, f, indent=2)

print("📊 COVERAGE: " + str(result["mean_depth"]) + "x")
print("   pct_bases_20x: " + str(result["pct_bases_20x"]) + "%")
if result.get("warning"):
    print("⚠️  " + result["warning"])
PYEOF
    """
}

// ── Proses 6: Deteksi Kromosom dari BAM ───────────────────────────────────────
process GET_CHROMS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path(bam), path(bai), path("chroms.txt")

    script:
    """
    # Hanya ambil kromosom utama — exclude alt, random, chrUn, decoy, HLA
    ${params.samtools_bin} view -H "${bam}" \
        | grep "^@SQ" \
        | awk '{print \$2}' \
        | sed 's/SN://' \
        | grep -E "^(chr[0-9]+|chr[XYM]|[0-9]+|[XY]|MT)\$" \
        > chroms.txt

    if [ ! -s chroms.txt ]; then
        echo "ERROR: Tidak ada kromosom utama ditemukan di BAM header" >&2
        exit 1
    fi
    """
}

// ── Proses 7: BCFTools mpileup per Kromosom (Scatter) ─────────────────────────
process CALL_VARIANTS {
    tag "${sample_id}:${chrom}"

    input:
    tuple val(sample_id), path(bam), path(bai), val(chrom)
    path ref_genome
    path ref_index

    output:
    tuple val(sample_id), val(chrom), path("${sample_id}.${chrom}.bcf"), path("${sample_id}.${chrom}.bcf.csi")

    script:
    // Fix WES/Targeted: gunakan -r {chrom} bersama BED agar tidak duplikat varian
    def region_flag = params.target_bed ? "-r ${chrom} -T ${params.target_bed}" : "-r ${chrom}"
    // Fix bug: max_depth=0 di bcftools bukan unlimited, gunakan 8000 sebagai unlimited
    def max_depth   = (params.max_depth as int) > 0 ? params.max_depth : 8000
    def min_bq      = params.min_bq ?: 20
    def min_mq      = params.min_mq ?: 20
    """
    ${params.bcftools_bin} mpileup \
        -O b \
        -o "${sample_id}.${chrom}.bcf" \
        -f "${ref_genome}" \
        ${region_flag} \
        --max-depth ${max_depth} \
        --min-BQ ${min_bq} \
        --min-MQ ${min_mq} \
        -a FORMAT/AD,FORMAT/DP \
        "${bam}"

    ${params.bcftools_bin} index "${sample_id}.${chrom}.bcf"
    """
}

// ── Proses 8: Gather + Call + Norm + Filter VCF ───────────────────────────────
process GATHER_VCF {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bcf_files), path(bcf_indexes)
    path ref_genome

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), path("${sample_id}.filtered.vcf.gz.tbi")

    script:
    def min_qual = params.min_qual ?: 20
    def min_dp   = params.min_dp   ?: 10
    def min_gq   = params.min_gq   ?: 20
    // Fix #19: gunakan samples.txt untuk ploidy chrX/chrY berdasarkan jenis kelamin
    def ploidy_flag = params.sample_sex
        ? "--ploidy-file samples.txt"
        : "--ploidy ${params.ploidy}"
    """
    # Fix #19: buat samples.txt untuk bcftools call agar ploidy chrX/chrY benar
    if [ "${params.sample_sex}" = "male" ]; then
        echo -e "X\\t1\\t${sample_id}\\tM\\t1" > samples.txt
        echo -e "Y\\t1\\t${sample_id}\\tM\\t1" >> samples.txt
    elif [ "${params.sample_sex}" = "female" ]; then
        echo -e "X\\t1\\t${sample_id}\\tF\\t2" > samples.txt
        echo -e "Y\\t1\\t${sample_id}\\tF\\t0" >> samples.txt
    fi

    # Sort BCF list by chromosome order
    BCF_LIST=\$(ls *.bcf | grep -v "\\.csi" | sort -V | tr '\\n' ' ')

    # Concat tanpa --allow-overlaps untuk hindari duplikat pada WES/Targeted
    ${params.bcftools_bin} concat \
        --rm-dups all \
        -O b \
        -o "${sample_id}.merged.bcf" \
        \${BCF_LIST}

    # Call variants dengan AD dan GQ annotation
    # Fix #19: pakai samples.txt bila ada, default ploidy bila tidak
    if [ -f samples.txt ]; then
        ${params.bcftools_bin} call \
            --ploidy-file samples.txt \
            -m -v \
            -f GQ \
            -o "${sample_id}.called.vcf.gz" \
            -O z \
            "${sample_id}.merged.bcf"
    else
        ${params.bcftools_bin} call \
            --ploidy ${params.ploidy} \
            -m -v \
            -f GQ \
            -o "${sample_id}.called.vcf.gz" \
            -O z \
            "${sample_id}.merged.bcf"
    fi

    # Normalize: left-align dan split multiallelic
    ${params.bcftools_bin} norm \
        -f "${ref_genome}" \
        -m -any \
        -O z \
        -o "${sample_id}.norm.vcf.gz" \
        "${sample_id}.called.vcf.gz"

    ${params.bcftools_bin} index -t "${sample_id}.norm.vcf.gz"

    # Filter berdasarkan params (tidak hardcode)
    # Filter expression sesuai UI: QUAL>=min_qual && FORMAT/DP>=min_dp && FORMAT/GQ>=min_gq
    ${params.bcftools_bin} filter \
        -s LOWQUAL \
        -e "QUAL<${min_qual} || FORMAT/DP<${min_dp} || FORMAT/GQ<${min_gq}" \
        -O z \
        -o "${sample_id}.filtered.vcf.gz" \
        "${sample_id}.norm.vcf.gz"

    ${params.bcftools_bin} index -t "${sample_id}.filtered.vcf.gz"
    """
}

// ── Workflow Utama ─────────────────────────────────────────────────────────────
workflow {

    ref_ch       = Channel.fromPath(params.ref_genome)
    ref_index_ch = Channel.fromPath("${params.ref_genome}.*").collect()

    if (params.reads_r1 && params.reads_r2) {
        log.info "Mode: Paired-End (PE)"
        pe_reads_ch = Channel.fromPath(params.reads_r1).map { r1 ->
            def raw_sid = r1.baseName.replaceAll(/_R1.*$|_1$/, '')
            def sid     = raw_sid.replaceAll(/[^A-Za-z0-9._\-]/, '_').replaceAll(/_+/, '_').replaceAll(/^_|_$/, '')
            tuple(sid, r1, file(params.reads_r2))
        }
        fastqc_out  = FASTQC(pe_reads_ch.map { sid, r1, r2 -> tuple(sid, [r1, r2]) })
        trimmed_out = TRIMMOMATIC_PE(pe_reads_ch)
        aligned_out = BWA_ALIGN_PE(trimmed_out, ref_ch.first(), ref_index_ch)
    } else {
        log.info "Mode: Single-End (SE)"
        se_reads_ch = Channel.fromPath(params.reads).map { f ->
            def raw_sid = f.baseName.replaceAll(/\.fastq$/, '')
            def sid     = raw_sid.replaceAll(/[^A-Za-z0-9._\-]/, '_').replaceAll(/_+/, '_').replaceAll(/^_|_$/, '')
            tuple(sid, f)
        }
        fastqc_out  = FASTQC(se_reads_ch)
        trimmed_out = TRIMMOMATIC_SE(se_reads_ch)
        aligned_out = BWA_ALIGN_SE(trimmed_out, ref_ch.first(), ref_index_ch)
    }

    // MarkDuplicates
    markdup_out = MARKDUPLICATES(aligned_out)
    bam_out     = markdup_out[0]

    // Hitung coverage (per region bila WES/Targeted)
    coverage_out = COVERAGE_STATS(bam_out)
    coverage_out.map { sid, cov_file ->
        def cov = new groovy.json.JsonSlurper().parse(cov_file)
        "📊 COVERAGE: ${cov.mean_depth}x | pct_bases_20x: ${cov.pct_bases_20x}%${cov.warning ? ' | ⚠️ ' + cov.warning : ''}"
    }.view()

    // Scatter-gather BCFTools per kromosom
    bam_with_chroms = GET_CHROMS(bam_out)

    scatter_ch = bam_with_chroms.flatMap { sample_id, bam, bai, chroms_file ->
        chroms_file.readLines()
            .findAll { it.trim() != "" }
            .collect { chrom -> tuple(sample_id, bam, bai, chrom.trim()) }
    }

    per_chrom_vcf = CALL_VARIANTS(scatter_ch, ref_ch.first(), ref_index_ch)

    gathered_ch = per_chrom_vcf
        .groupTuple(by: 0)
        .map { sid, chrom_list, bcf_files, bcf_indexes -> tuple(sid, bcf_files, bcf_indexes) }

    vcf_out = GATHER_VCF(gathered_ch, ref_ch.first())

    vcf_out.map { sid, vcf, tbi ->
        "✅ VCF Output: ${params.outdir}/${sid}/vcf/${vcf.name}"
    }.view()
}
