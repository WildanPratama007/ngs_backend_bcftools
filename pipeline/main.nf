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
    """
    ${params.bin_dir}/trimmomatic SE -threads ${params.threads} \
        "${reads}" "${sample_id}.trimmed.fastq" \
        ${params.trimmer_opts} 2> trimmomatic.log

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
    """
    ${params.bin_dir}/trimmomatic PE -threads ${params.threads} \
        "${reads_r1}" "${reads_r2}" \
        "${sample_id}_R1.trimmed.fastq" "${sample_id}_R1.unpaired.fastq" \
        "${sample_id}_R2.trimmed.fastq" "${sample_id}_R2.unpaired.fastq" \
        ${params.trimmer_opts} 2> trimmomatic.log

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
    """
    ${params.bin_dir}/bwa mem \
        -t ${params.threads} \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1" \
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
    """
    ${params.bin_dir}/bwa mem \
        -t ${params.threads} \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1" \
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
    """
    ${params.samtools_bin} markdup \
        -f "${sample_id}.markdup.metrics.txt" \
        "${bam}" \
        "${sample_id}.markdup.bam"

    ${params.samtools_bin} index "${sample_id}.markdup.bam"
    """
}

// ── Proses 5: Hitung Coverage ─────────────────────────────────────────────────
process COVERAGE_STATS {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.txt")

    script:
    // Untuk WES/Targeted: hitung coverage di region BED saja
    // Untuk WGS: hitung coverage seluruh referensi
    def region_opt = params.target_bed ? "-l ${params.target_bed}" : ""
    """
    ${params.samtools_bin} coverage ${region_opt} "${bam}" \
        | awk 'NR>1 && \$3>0 {sum+=\$7*\$3; len+=\$3} END{printf "%.2f\\n", (len>0 ? sum/len : 0)}' \
        > "${sample_id}.coverage.txt"
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
    // BED di-filter per kromosom via -r sehingga setiap task hanya proses 1 kromosom
    def region_flag = params.target_bed ? "-r ${chrom} -T ${params.target_bed}" : "-r ${chrom}"
    def max_depth   = params.max_depth ?: 0  // 0 = unlimited
    """
    ${params.bcftools_bin} mpileup \
        -O b \
        -o "${sample_id}.${chrom}.bcf" \
        -f "${ref_genome}" \
        ${region_flag} \
        --max-depth ${max_depth} \
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
    def min_qual    = params.min_qual ?: 20
    def min_dp      = params.min_dp   ?: 5
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
    ${params.bcftools_bin} filter \
        -s LOWQUAL \
        -e "QUAL<${min_qual} || FORMAT/DP<${min_dp}" \
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
        "📊 COVERAGE: ${cov_file.text.trim()}x"
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
