#!/usr/bin/env nextflow

/*
 * NGS E2E Secondary Analysis Pipeline — BCFTools
 * Alur: FastQC → Trimmomatic (SE/PE) → BWA MEM → BCFTools mpileup/call/filter
 * Support: SE/PE, WGS/WES/Targeted, GRCh38/GRCh37
 */

nextflow.enable.dsl = 2

if (!params.reads && !params.reads_r1) {
    error "Parameter --reads (SE) atau --reads_r1 + --reads_r2 (PE) wajib diisi."
}

def CHROMOSOMES_HG38 = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]
def CHROMOSOMES_HG19 = (1..22).collect { "${it}" } + ["X", "Y"]

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
        ${reads} ${sample_id}.trimmed.fastq \
        ${params.trimmer_opts} 2> trimmomatic.log

    if [ ! -s ${sample_id}.trimmed.fastq ]; then
        cp ${reads} ${sample_id}.trimmed.fastq
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
        ${reads_r1} ${reads_r2} \
        ${sample_id}_R1.trimmed.fastq ${sample_id}_R1.unpaired.fastq \
        ${sample_id}_R2.trimmed.fastq ${sample_id}_R2.unpaired.fastq \
        ${params.trimmer_opts} 2> trimmomatic.log

    if [ ! -s ${sample_id}_R1.trimmed.fastq ]; then
        cp ${reads_r1} ${sample_id}_R1.trimmed.fastq
        cp ${reads_r2} ${sample_id}_R2.trimmed.fastq
        touch ${sample_id}_R1.unpaired.fastq ${sample_id}_R2.unpaired.fastq
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
        ${ref_genome} ${trimmed_reads} \
        | ${params.samtools_bin} sort -@ ${params.threads} -m 4G \
          -o ${sample_id}.aligned.sorted.bam

    ${params.samtools_bin} index ${sample_id}.aligned.sorted.bam
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
        ${ref_genome} ${r1} ${r2} \
        | ${params.samtools_bin} fixmate -m -@ ${params.threads} - - \
        | ${params.samtools_bin} sort -@ ${params.threads} -m 4G \
          -o ${sample_id}.aligned.sorted.bam

    ${params.samtools_bin} index ${sample_id}.aligned.sorted.bam
    """
}

// ── Proses 4: Hitung Coverage ─────────────────────────────────────────────────
process COVERAGE_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.txt")

    script:
    """
    # samtools coverage menghitung meandepth terhadap SELURUH panjang referensi
    # (termasuk posisi dengan depth=0), sehingga hasilnya akurat sebagai genome coverage
    ${params.samtools_bin} coverage ${bam} \
        | awk 'NR>1 && \$3>0 {sum+=\$7*\$3; len+=\$3} END{printf "%.2f\\n", (len>0 ? sum/len : 0)}' \
        > ${sample_id}.coverage.txt
    """
}

// ── Proses 5: Deteksi Kromosom dari BAM ───────────────────────────────────────
process GET_CHROMS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path(bam), path(bai), path("chroms.txt")

    script:
    """
    ${params.samtools_bin} view -H ${bam} \
        | grep "^@SQ" \
        | awk '{print \$2}' \
        | sed 's/SN://' \
        > chroms.txt
    """
}

// ── Proses 6: BCFTools per Kromosom (Scatter) ─────────────────────────────────
process CALL_VARIANTS {
    tag "${sample_id}:${chrom}"

    input:
    tuple val(sample_id), path(bam), path(bai), val(chrom)
    path ref_genome
    path ref_index

    output:
    tuple val(sample_id), val(chrom), path("${sample_id}.${chrom}.bcf"), path("${sample_id}.${chrom}.bcf.csi")

    script:
    def region_flag = params.target_bed ? "-R ${params.target_bed}" : "-r ${chrom}"
    """
    ${params.bcftools_bin} mpileup \
        -O b -o ${sample_id}.${chrom}.bcf \
        -f ${ref_genome} \
        ${region_flag} \
        ${bam}

    ${params.bcftools_bin} index ${sample_id}.${chrom}.bcf
    """
}

// ── Proses 7: Gather + Call + Filter VCF ─────────────────────────────────────
process GATHER_VCF {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bcf_files), path(bcf_indexes)
    path ref_genome

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf")

    script:
    """
    BCF_LIST=\$(ls *.bcf | grep -v ".csi" | sort -V | tr '\\n' ' ')

    ${params.bcftools_bin} concat \
        --allow-overlaps -O b \
        -o ${sample_id}.merged.bcf \
        \${BCF_LIST}

    ${params.bcftools_bin} call \
        --ploidy ${params.ploidy} -m -v \
        -o ${sample_id}.vcf \
        ${sample_id}.merged.bcf

    ${params.bcftools_bin} filter \
        -s LOWQUAL -e 'QUAL<20 || DP<5' \
        ${sample_id}.vcf \
        > ${sample_id}.filtered.vcf
    """
}

// ── Workflow Utama ─────────────────────────────────────────────────────────────
workflow {

    ref_ch       = Channel.fromPath(params.ref_genome)
    ref_index_ch = Channel.fromPath("${params.ref_genome}.*").collect()

    if (params.reads_r1 && params.reads_r2) {
        log.info "Mode: Paired-End (PE)"
        pe_reads_ch = Channel.fromPath(params.reads_r1).map { r1 ->
            def sid = r1.baseName.replaceAll(/_R1.*$|_1$/, '')
            tuple(sid, r1, file(params.reads_r2))
        }
        fastqc_out  = FASTQC(pe_reads_ch.map { sid, r1, r2 -> tuple(sid, [r1, r2]) })
        trimmed_out = TRIMMOMATIC_PE(pe_reads_ch)
        aligned_out = BWA_ALIGN_PE(trimmed_out, ref_ch.first(), ref_index_ch)
    } else {
        log.info "Mode: Single-End (SE)"
        se_reads_ch = Channel.fromPath(params.reads).map { f ->
            tuple(f.baseName.replaceAll(/\.fastq$/, ''), f)
        }
        fastqc_out  = FASTQC(se_reads_ch)
        trimmed_out = TRIMMOMATIC_SE(se_reads_ch)
        aligned_out = BWA_ALIGN_SE(trimmed_out, ref_ch.first(), ref_index_ch)
    }

    // Hitung coverage
    coverage_out = COVERAGE_STATS(aligned_out)
    coverage_out.map { sid, cov_file ->
        "📊 COVERAGE: ${cov_file.text.trim()}x"
    }.view()

    // Scatter-gather BCFTools per kromosom
    bam_with_chroms = GET_CHROMS(aligned_out)

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

    vcf_out.map { sid, vcf ->
        "✅ VCF Output: ${params.outdir}/${sid}/vcf/${vcf.name}"
    }.view()
}
