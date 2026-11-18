process BCFTOOLS_CALL {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.bcftools_stats.txt"),                    emit: stats

    script:
    """
    bcftools mpileup \\
        --fasta-ref ${fasta} \\
        --max-depth 100000 \\
        --min-BQ ${params.min_base_quality} \\
        --annotate FORMAT/AD,FORMAT/DP \\
        --output-type u ${bam} \\
    | bcftools call --ploidy 1 --multiallelic-caller --variants-only --output-type u \\
    | bcftools norm --fasta-ref ${fasta} --output-type u \\
    | bcftools +fill-tags --output-type u -- -t FORMAT/VAF \\
    | bcftools filter \\
        --include 'QUAL >= ${params.min_variant_qual} && INFO/DP >= ${params.min_depth}' \\
        --output-type z --output ${meta.id}.vcf.gz
    tabix -p vcf ${meta.id}.vcf.gz
    bcftools stats ${meta.id}.vcf.gz > ${meta.id}.bcftools_stats.txt
    """
}

process BCFTOOLS_CONSENSUS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/consensus", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(tbi), path(depth)
    path fasta
    path primer_bed

    output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.low_coverage.bed"), emit: mask

    script:
    """
    # Positions below the minimum depth become a BED mask (0-based, half-open)
    zcat ${depth} | awk -v min=${params.min_depth} 'BEGIN{OFS="\\t"}
        \$3 < min { if (chrom == \$1 && \$2 == end + 1) { end = \$2 } else { if (chrom != "") print chrom, start - 1, end; chrom = \$1; start = \$2; end = \$2 } }
        END { if (chrom != "") print chrom, start - 1, end }' > ${meta.id}.low_coverage.bed

    # Majority variants only: VAF at or above the consensus threshold
    bcftools view --include 'FORMAT/VAF[0:0] >= ${params.consensus_min_af}' --output-type z --output majority.vcf.gz ${vcf}
    tabix -p vcf majority.vcf.gz

    # For amplicon data, indels at or within 10 bp of a primer site are a common artifact of amplicon
    # junctions; they are reported in the variant table with a review flag but not applied to the consensus
    if [ "${primer_bed.name}" != "NO_FILE" ]; then
        awk 'BEGIN{OFS="\\t"} { s = \$2 - 10; if (s < 0) s = 0; print \$1, s, \$3 + 10 }' ${primer_bed} > primer_padded.bed
        bcftools view --types snps,mnps --output-type z --output snvs.vcf.gz majority.vcf.gz
        bcftools view --types indels --targets-file ^primer_padded.bed --output-type z --output indels.vcf.gz majority.vcf.gz
        tabix -p vcf snvs.vcf.gz
        tabix -p vcf indels.vcf.gz
        bcftools concat --allow-overlaps snvs.vcf.gz indels.vcf.gz --output-type z --output applied.vcf.gz
    else
        cp majority.vcf.gz applied.vcf.gz
    fi
    tabix -f -p vcf applied.vcf.gz

    # Variants inside low-coverage (masked) regions cannot be applied
    if [ -s ${meta.id}.low_coverage.bed ]; then
        bcftools view --targets-file ^${meta.id}.low_coverage.bed --output-type z --output pass.vcf.gz applied.vcf.gz
    else
        cp applied.vcf.gz pass.vcf.gz
    fi
    tabix -f -p vcf pass.vcf.gz

    bcftools consensus \\
        --fasta-ref ${fasta} \\
        --mask ${meta.id}.low_coverage.bed \\
        --mark-del - \\
        pass.vcf.gz \\
    | awk -v id="${meta.id}" '/^>/ { print ">" id; next } { print }' > ${meta.id}.consensus.fasta
    """
}
