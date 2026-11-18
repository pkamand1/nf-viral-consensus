process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(alignment)
    val suffix

    output:
    tuple val(meta), path("${meta.id}${suffix}.bam"), path("${meta.id}${suffix}.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}${suffix}.bam ${alignment}
    samtools index ${meta.id}${suffix}.bam
    """
}
