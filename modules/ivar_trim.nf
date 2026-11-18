process IVAR_TRIM {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.log'

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed

    output:
    tuple val(meta), path("${meta.id}.primertrim.bam"), emit: bam
    tuple val(meta), path("${meta.id}.ivar_trim.log"),  emit: log

    script:
    """
    ivar trim -e -i ${bam} -b ${primer_bed} -q ${params.min_base_quality} -m ${params.min_read_length} \\
        -p ${meta.id}.primertrim > ${meta.id}.ivar_trim.log
    """
}
