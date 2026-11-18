process ANNOTATE_VARIANTS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(tbi)
    path gff
    path primer_bed

    output:
    tuple val(meta), path("${meta.id}.variants.tsv"), emit: table

    script:
    def primers = primer_bed.name != 'NO_FILE' ? "--primer-bed ${primer_bed}" : ''
    """
    annotate_variants.py --vcf ${vcf} --gff ${gff} --sample ${meta.id} ${primers} \\
        --consensus-min-af ${params.consensus_min_af} > ${meta.id}.variants.tsv
    """
}
