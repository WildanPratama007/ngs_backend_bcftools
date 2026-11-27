#!/usr/bin/env nextflow

/*
 * NGS Pipeline Wizard — Per-Step Execution
 * File ini dipakai oleh wizard endpoints (/api/analyses/id/step/run)
 * Setiap step dijalankan terpisah via -entry flag
 *
 * Entry points:
 *   -entry QC              -> Step 1: FastQC saja
 *   -entry TRIMMING        -> Step 2: Trimmomatic saja
 *   -entry ALIGNMENT       -> Step 3: BWA MEM + MarkDup + Coverage
 *   -entry VARIANT_CALLING -> Step 4: mpileup + call + norm + filter
 */

nextflow.enable.dsl = 2

// ── Helper: sanitasi nama sampel ──────────────────────────────────────────────
def sanitizeName(String s) {
    s.replaceAll(/[^A-Za-z0-9._\-]/, '_').replaceAll(/_+/, '_').replaceAll(/^_|_$/, '')
}

// ══════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ══════════════════════════════════════════════════════════════════════════════

process FASTQC {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/fastqc", mode: params.publish_mode

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

// ── Helper: nama sampel ───────────────────────────────────────────────────────
// --sample_name dari API adalah sumber kebenaran (sama dengan folder hasil di backend).
// Kalau tidak diberikan, nama diturunkan dari nama file.
def sampleIdSE(f) {
    params.sample_name
        ? sanitizeName(params.sample_name.toString())
        : sanitizeName(f.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, ''))
}

def sampleIdPE(f) {
    params.sample_name
        ? sanitizeName(params.sample_name.toString())
        : sanitizeName(f.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '').replaceAll(/[_.]R?1(_\d+)?$/, ''))
}

// Pesan saran saat Trimmomatic menyisakan 0 reads
def trimHint(int q) {
    q >= 25
        ? "SLIDINGWINDOW quality terlalu tinggi (Q${q}). Turunkan ke Q15-Q20 (preset Research atau Custom)."
        : "Semua reads gagal memenuhi kualitas minimum (Q${q}). Turunkan SLIDINGWINDOW ke Q10-Q15 atau periksa FastQC."
}

// Isi file poly-G untuk ILLUMINACLIP (NovaSeq/NextSeq)
def polyGSetup() {
    'echo ">polyG" > polyG.fa; echo "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" >> polyG.fa'
}

process TRIMMOMATIC_SE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/trimmed", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed${params.trimmed_ext}"), val("SE")

    script:
    if (params.trim_adapter && !params.adapter_file) {
        error "trim_adapter=true tetapi --adapter_file tidak diisi."
    }
    def q            = params.trim_sliding_window_quality as int
    def adapter_clip = params.trim_adapter ? "ILLUMINACLIP:${params.adapter_file}:2:30:10:8:true" : ""
    def poly_g_clip  = params.trim_poly_g  ? "ILLUMINACLIP:polyG.fa:2:30:10" : ""
    def poly_g_setup = params.trim_poly_g  ? polyGSetup() : ""
    def trim_cmd     = [adapter_clip, poly_g_clip, "LEADING:5", "TRAILING:5",
                        "SLIDINGWINDOW:4:${q}", "MINLEN:${params.trim_min_len}"]
                           .findAll { it }.join(" ")
    def hint         = trimHint(q)
    """
    ${poly_g_setup}
    ${params.bin_dir}/trimmomatic SE -threads ${params.threads} \
        "${reads}" "${sample_id}.trimmed${params.trimmed_ext}" \
        ${trim_cmd} 2> trimmomatic.log || { cat trimmomatic.log >&2; exit 1; }

    nonempty() { case "\$1" in *.gz) [ -n "\$(gzip -dc "\$1" | head -c 1)" ];; *) [ -s "\$1" ];; esac; }
    if ! nonempty "${sample_id}.trimmed${params.trimmed_ext}"; then
        INPUT=\$(grep -oP 'Input Reads: \\K\\d+' trimmomatic.log || echo 0)
        echo "ERROR: Trimmomatic menghasilkan 0 reads dari \${INPUT} reads input." >&2
        echo "ERROR: ${hint}" >&2
        exit 1
    fi
    """
}

process TRIMMOMATIC_PE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/trimmed", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed${params.trimmed_ext}"),
          path("${sample_id}_R1.unpaired.fastq"),
          path("${sample_id}_R2.trimmed${params.trimmed_ext}"),
          path("${sample_id}_R2.unpaired.fastq"),
          val("PE")

    script:
    if (params.trim_adapter && !params.adapter_file) {
        error "trim_adapter=true tetapi --adapter_file tidak diisi."
    }
    def q            = params.trim_sliding_window_quality as int
    def adapter_clip = params.trim_adapter ? "ILLUMINACLIP:${params.adapter_file}:2:30:10:8:true" : ""
    def poly_g_clip  = params.trim_poly_g  ? "ILLUMINACLIP:polyG.fa:2:30:10" : ""
    def poly_g_setup = params.trim_poly_g  ? polyGSetup() : ""
    def trim_cmd     = [adapter_clip, poly_g_clip, "LEADING:5", "TRAILING:5",
                        "SLIDINGWINDOW:4:${q}", "MINLEN:${params.trim_min_len}"]
                           .findAll { it }.join(" ")
    def hint         = trimHint(q)
    """
    ${poly_g_setup}
    ${params.bin_dir}/trimmomatic PE -threads ${params.threads} \
        "${reads_r1}" "${reads_r2}" \
        "${sample_id}_R1.trimmed${params.trimmed_ext}" "${sample_id}_R1.unpaired.fastq" \
        "${sample_id}_R2.trimmed${params.trimmed_ext}" "${sample_id}_R2.unpaired.fastq" \
        ${trim_cmd} 2> trimmomatic.log || { cat trimmomatic.log >&2; exit 1; }

    nonempty() { case "\$1" in *.gz) [ -n "\$(gzip -dc "\$1" | head -c 1)" ];; *) [ -s "\$1" ];; esac; }
    if ! nonempty "${sample_id}_R1.trimmed${params.trimmed_ext}"; then
        INPUT=\$(grep -oP 'Input Read Pairs: \\K\\d+' trimmomatic.log || echo 0)
        echo "ERROR: Trimmomatic (PE) menghasilkan 0 read pairs dari \${INPUT} pairs input." >&2
        echo "ERROR: ${hint}" >&2
        exit 1
    fi
    """
}

process BWA_ALIGN_SE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: params.publish_mode, enabled: params.publish_intermediate

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
        -R "@RG\\tID:${sample_id}.L001\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}\\tPU:${sample_id}.L001" \
        "${ref_genome}" "${trimmed_reads}" \
        | ${params.samtools_bin} sort -@ ${params.threads} -m ${params.sort_mem} \
          -o "${sample_id}.aligned.sorted.bam"
    ${params.samtools_bin} index "${sample_id}.aligned.sorted.bam"
    """
}

// Wizard: hanya pasangan trimmed (R1/R2); file unpaired tidak dipakai di alignment
process BWA_ALIGN_PE {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: params.publish_mode, enabled: params.publish_intermediate

    input:
    tuple val(sample_id), path(r1), path(r2), val(mode)
    path ref_genome
    path ref_index

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), path("${sample_id}.aligned.sorted.bam.bai")

    script:
    """
    ${params.bin_dir}/bwa mem \
        -t ${params.threads} \
        -R "@RG\\tID:${sample_id}.L001\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}\\tPU:${sample_id}.L001" \
        "${ref_genome}" "${r1}" "${r2}" \
        | ${params.samtools_bin} fixmate -m -@ ${params.threads} - - \
        | ${params.samtools_bin} sort -@ ${params.threads} -m ${params.sort_mem} \
          -o "${sample_id}.aligned.sorted.bam"
    ${params.samtools_bin} index "${sample_id}.aligned.sorted.bam"
    """
}

process MARKDUPLICATES {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: params.publish_mode

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai")
    path "${sample_id}.markdup.metrics.txt"

    script:
    if (params.mark_duplicates) {
        """
        ${params.samtools_bin} markdup \
            -f "${sample_id}.markdup.metrics.txt" \
            "${bam}" "${sample_id}.markdup.bam"
        ${params.samtools_bin} index "${sample_id}.markdup.bam"
        """
    } else {
        """
        cp "${bam}" "${sample_id}.markdup.bam"
        cp "${bai}" "${sample_id}.markdup.bam.bai"
        echo "mark_duplicates=false" > "${sample_id}.markdup.metrics.txt"
        """
    }
}

// Coverage di region target (WES/Targeted): hitung dari depth per-basis dalam BED.
// `samtools coverage -l` adalah min-read-length, BUKAN BED, jadi tidak bisa dipakai di sini.
def targetCoverageScript(String samtools, String bed, String bam, String sample_id) {
    """
    set -o pipefail
    cat > depth_stats.py << 'PYEOF'
import json, sys
n = s = cov = c10 = c20 = c30 = 0
hist = {}
with sys.stdin as f:
    for line in f:
        d = int(line.rsplit(None, 1)[-1])
        n += 1; s += d
        if d > 0:  cov += 1
        if d >= 10: c10 += 1
        if d >= 20: c20 += 1
        if d >= 30: c30 += 1
        k = d if d < 2000 else 2000
        hist[k] = hist.get(k, 0) + 1

if n == 0:
    result = {"mean_depth": 0.0, "median_depth": 0.0, "pct_bases_10x": 0.0, "pct_bases_20x": 0.0,
              "pct_bases_30x": 0.0, "total_bases": 0, "covered_bases": 0,
              "warning": "Tidak ada posisi di region target (periksa BED dan referensi)"}
else:
    mean = s / n
    half, acc, median = n / 2.0, 0, 0
    for k in sorted(hist):
        acc += hist[k]
        if acc >= half:
            median = k
            break
    warning = None
    if mean < 10:   warning = "Coverage sangat rendah (<10x) — hasil tidak reliable"
    elif mean < 20: warning = "Coverage rendah (<20x) — beberapa variant mungkin terlewat"
    result = {"mean_depth": round(mean, 2), "median_depth": float(median),
              "pct_bases_10x": round(100.0 * c10 / n, 1), "pct_bases_20x": round(100.0 * c20 / n, 1),
              "pct_bases_30x": round(100.0 * c30 / n, 1),
              "total_bases": n, "covered_bases": cov, "warning": warning}

with open("${sample_id}.coverage.json", "w") as f:
    json.dump(result, f, indent=2)

print("📊 COVERAGE: " + str(result["mean_depth"]) + "x | pct_bases_20x: " + str(result["pct_bases_20x"]) + "%")
if result.get("warning"): print("⚠️  " + result["warning"])
PYEOF
    ${samtools} depth -a -b "${bed}" "${bam}" | python3 depth_stats.py
    """
}

process COVERAGE_STATS {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/bwa", mode: params.publish_mode

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.json")

    script:
    if (params.target_bed) {
        targetCoverageScript(params.samtools_bin, params.target_bed.toString(), bam.toString(), sample_id)
    } else {
    """
    ${params.samtools_bin} coverage "${bam}" > coverage_raw.txt

    python3 << 'PYEOF'
import json, math

rows = []
with open("coverage_raw.txt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\\t")
        if len(p) < 7: continue
        try:
            rows.append({"length": int(p[2])-int(p[1]), "covbases": int(p[4]), "meandepth": float(p[6])})
        except: continue

if not rows:
    result = {"mean_depth":0.0,"median_depth":0.0,"pct_bases_10x":0.0,"pct_bases_20x":0.0,"pct_bases_30x":0.0,"total_bases":0,"covered_bases":0,"warning":"No coverage data"}
else:
    total_len  = sum(r["length"] for r in rows)
    total_cov  = sum(r["covbases"] for r in rows)
    mean_depth = sum(r["meandepth"]*r["length"] for r in rows) / total_len if total_len > 0 else 0
    depths     = sorted(r["meandepth"] for r in rows)
    n          = len(depths)
    median_depth = (depths[n//2] + depths[(n-1)//2]) / 2 if n > 0 else 0

    def pct_nx(mean, n):
        if mean <= 0: return 0.0
        term = math.exp(-mean); s = term
        for k in range(1, n):
            term *= mean/k; s += term
        return round(max(0.0, (1.0-s)*100), 1)

    warning = None
    if mean_depth < 10: warning = "Coverage sangat rendah (<10x) — hasil tidak reliable"
    elif mean_depth < 20: warning = "Coverage rendah (<20x) — beberapa variant mungkin terlewat"

    result = {
        "mean_depth": round(mean_depth,2), "median_depth": round(median_depth,2),
        "pct_bases_10x": pct_nx(mean_depth,10), "pct_bases_20x": pct_nx(mean_depth,20),
        "pct_bases_30x": pct_nx(mean_depth,30),
        "total_bases": total_len, "covered_bases": total_cov, "warning": warning
    }

with open("${sample_id}.coverage.json","w") as f:
    json.dump(result, f, indent=2)

print("📊 COVERAGE: " + str(result["mean_depth"]) + "x | pct_bases_20x: " + str(result["pct_bases_20x"]) + "%")
if result.get("warning"): print("⚠️  " + result["warning"])
PYEOF
    """
    }
}

process GET_CHROMS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path(bam), path(bai), path("chroms.txt")

    script:
    """
    ${params.samtools_bin} view -H "${bam}" \
        | grep "^@SQ" | awk '{print \$2}' | sed 's/SN://' \
        | grep -E "^(chr[0-9]+|chr[XYM]|[0-9]+|[XY]|MT)\$" > chroms.txt

    if [ ! -s chroms.txt ]; then
        echo "ERROR: Tidak ada kromosom utama ditemukan di BAM header." >&2
        echo "" >&2
        echo "Kemungkinan penyebab:" >&2
        echo "  ⚠️  Referensi yang dipilih tidak sesuai dengan data FASTQ" >&2
        echo "     - Data chr22 → gunakan referensi 'GRCh38 chr22 (testing)'" >&2
        echo "     - Data full WGS → gunakan 'GRCh38' atau 'GRCh37'" >&2
        echo "  ⚠️  Coverage terlalu rendah — tidak ada reads yang ter-align" >&2
        exit 1
    fi
    """
}

process CALL_VARIANTS {
    tag "${sample_id}:${chrom}"

    input:
    tuple val(sample_id), path(bam), path(bai), val(chrom)
    path ref_genome
    path ref_index

    output:
    tuple val(sample_id), val(chrom), path("${sample_id}.${chrom}.bcf"), path("${sample_id}.${chrom}.bcf.csi")

    script:
    def region_flag = params.target_bed ? "-r ${chrom} -T ${params.target_bed}" : "-r ${chrom}"
    def max_depth   = (params.max_depth as int) > 0 ? params.max_depth : 8000
    """
    ${params.bcftools_bin} mpileup -O b \
        -o "${sample_id}.${chrom}.bcf" \
        -f "${ref_genome}" \
        ${region_flag} \
        --max-depth ${max_depth} \
        --min-BQ ${params.min_bq} \
        --min-MQ ${params.min_mq} \
        -a FORMAT/AD,FORMAT/DP "${bam}"
    ${params.bcftools_bin} index "${sample_id}.${chrom}.bcf"
    """
}

process GATHER_VCF {
    tag        "$sample_id"
    publishDir "${params.outdir}/${sample_id}/vcf", mode: params.publish_mode

    input:
    tuple val(sample_id), path(bcf_files), path(bcf_indexes)
    path ref_genome

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), path("${sample_id}.filtered.vcf.gz.tbi")

    script:
    def min_qual = params.min_qual ?: 20
    def min_dp   = params.min_dp   ?: 10
    def min_gq   = params.min_gq   ?: 20
    // sample_sex -> file "NAMA SEX" untuk bcftools call -S (ploidy chrX/chrY)
    def sex_code  = params.sample_sex ? (params.sample_sex.toString().toLowerCase() == 'male' ? 'M' : 'F') : null
    def sex_setup = sex_code ? "echo \"${sample_id} ${sex_code}\" > samples.txt" : ""
    def sex_flag  = sex_code ? "-S samples.txt" : ""
    """
    ${sex_setup}
    BCF_LIST=\$(ls *.bcf | sort -V | tr '\\n' ' ')
    ${params.bcftools_bin} concat -a --rm-dups all -O b -o "${sample_id}.merged.bcf" \${BCF_LIST}

    ${params.bcftools_bin} call --ploidy ${params.ploidy} ${sex_flag} -m -v -f GQ -O z -o "${sample_id}.called.vcf.gz" "${sample_id}.merged.bcf"

    ${params.bcftools_bin} norm -f "${ref_genome}" -m -any -O z -o "${sample_id}.norm.vcf.gz" "${sample_id}.called.vcf.gz"
    ${params.bcftools_bin} index -t "${sample_id}.norm.vcf.gz"

    ${params.bcftools_bin} filter -s LOWQUAL \
        -e "QUAL<${min_qual} || FORMAT/DP<${min_dp} || FORMAT/GQ<${min_gq}" \
        -O z -o "${sample_id}.filtered.vcf.gz" "${sample_id}.norm.vcf.gz"
    ${params.bcftools_bin} index -t "${sample_id}.filtered.vcf.gz"
    """
}

// ══════════════════════════════════════════════════════════════════════════════
// ENTRY POINT 1: QC
// ══════════════════════════════════════════════════════════════════════════════
workflow QC {
    if (!params.reads && !params.reads_r1) error "Parameter --reads atau --reads_r1/r2 wajib."

    if (params.reads_r1 && params.reads_r2) {
        reads_ch = Channel.fromPath(params.reads_r1).map { r1 ->
            tuple(sampleIdPE(r1), [r1, file(params.reads_r2)])
        }
    } else {
        reads_ch = Channel.fromPath(params.reads).map { f ->
            tuple(sampleIdSE(f), f)
        }
    }

    fastqc_out = FASTQC(reads_ch)
    fastqc_out.map { sid, html, zip -> "✅ FastQC: ${params.outdir}/${sid}/fastqc/${html.name}" }.view()
}

// ══════════════════════════════════════════════════════════════════════════════
// ENTRY POINT 2: TRIMMING
// ══════════════════════════════════════════════════════════════════════════════
workflow TRIMMING {
    if (!params.reads && !params.reads_r1) error "Parameter --reads atau --reads_r1/r2 wajib."

    if (params.reads_r1 && params.reads_r2) {
        pe_ch = Channel.fromPath(params.reads_r1).map { r1 ->
            tuple(sampleIdPE(r1), r1, file(params.reads_r2))
        }
        trimmed = TRIMMOMATIC_PE(pe_ch)
        trimmed.map { sid, r1, r1u, r2, r2u, mode -> "✅ Trimming: ${params.outdir}/${sid}/trimmed/" }.view()
    } else {
        se_ch = Channel.fromPath(params.reads).map { f ->
            tuple(sampleIdSE(f), f)
        }
        trimmed = TRIMMOMATIC_SE(se_ch)
        trimmed.map { sid, fq, mode -> "✅ Trimming: ${params.outdir}/${sid}/trimmed/" }.view()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// ENTRY POINT 3: ALIGNMENT
// Input WAJIB dari step TRIMMING:
//   SE: --trimmed_reads
//   PE: --trimmed_r1 dan --trimmed_r2
// ══════════════════════════════════════════════════════════════════════════════
workflow ALIGNMENT {
    def is_pe = params.trimmed_r1 && params.trimmed_r2
    if (!is_pe && !params.trimmed_reads) {
        error "ALIGNMENT: --trimmed_reads (SE) atau --trimmed_r1 + --trimmed_r2 (PE) wajib diisi. Jalankan step TRIMMING dulu."
    }

    ref_ch       = Channel.fromPath(params.ref_genome)
    ref_index_ch = Channel.fromPath("${params.ref_genome}.*").collect()

    if (is_pe) {
        log.info "ALIGNMENT (PE): ${params.trimmed_r1} + ${params.trimmed_r2}"
        trimmed_ch = Channel.fromPath(params.trimmed_r1).map { r1 ->
            def sid = sanitizeName(r1.name.replaceAll(/_R1\.trimmed\.fastq(\.gz)?$/, ''))
            tuple(sid, r1, file(params.trimmed_r2), "PE")
        }
        aligned_out = BWA_ALIGN_PE(trimmed_ch, ref_ch.first(), ref_index_ch)
    } else {
        log.info "ALIGNMENT (SE): ${params.trimmed_reads}"
        trimmed_ch = Channel.fromPath(params.trimmed_reads).map { f ->
            def sid = sanitizeName(f.name.replaceAll(/\.trimmed\.fastq(\.gz)?$/, ''))
            tuple(sid, f, "SE")
        }
        aligned_out = BWA_ALIGN_SE(trimmed_ch, ref_ch.first(), ref_index_ch)
    }

    markdup_out  = MARKDUPLICATES(aligned_out)
    bam_out      = markdup_out[0]
    coverage_out = COVERAGE_STATS(bam_out)

    coverage_out.map { sid, cov_file ->
        def cov = new groovy.json.JsonSlurper().parse(cov_file)
        "📊 COVERAGE: ${cov.mean_depth}x | pct_bases_20x: ${cov.pct_bases_20x}%${cov.warning ? ' | ⚠️ ' + cov.warning : ''}"
    }.view()

    bam_out.map { sid, bam, bai -> "✅ Alignment: ${params.outdir}/${sid}/bwa/${bam.name}" }.view()
}

// ══════════════════════════════════════════════════════════════════════════════
// ENTRY POINT 4: VARIANT_CALLING
// ══════════════════════════════════════════════════════════════════════════════
workflow VARIANT_CALLING {
    if (!params.bam_file) error "Parameter --bam_file wajib untuk variant calling."

    ref_ch       = Channel.fromPath(params.ref_genome)
    ref_index_ch = Channel.fromPath("${params.ref_genome}.*").collect()

    bam_ch = Channel.fromPath(params.bam_file).map { bam ->
        def sid = bam.name.replaceAll(/\.markdup\.bam$/, '')
        tuple(sid, bam, file("${bam}.bai"))
    }

    bam_with_chroms = GET_CHROMS(bam_ch)

    scatter_ch = bam_with_chroms.flatMap { sid, bam, bai, chroms_file ->
        chroms_file.readLines().findAll { it.trim() }.collect { chrom ->
            tuple(sid, bam, bai, chrom.trim())
        }
    }

    per_chrom_vcf = CALL_VARIANTS(scatter_ch, ref_ch.first(), ref_index_ch)

    gathered_ch = per_chrom_vcf
        .groupTuple(by: 0)
        .map { sid, cl, bcfs, csis -> tuple(sid, bcfs, csis) }

    vcf_out = GATHER_VCF(gathered_ch, ref_ch.first())
    vcf_out.map { sid, vcf, tbi -> "✅ VCF: ${params.outdir}/${sid}/vcf/${vcf.name}" }.view()
}
