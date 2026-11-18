process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat
    tuple val(meta), path("${meta.id}.coverage.tsv"), emit: coverage
    tuple val(meta), path("${meta.id}.depth.tsv.gz"), emit: depth

    script:
    """
    samtools flagstat ${bam} > ${meta.id}.flagstat.txt
    samtools coverage ${bam} > ${meta.id}.coverage.tsv
    samtools depth -a -J -Q 0 -q 0 ${bam} | gzip > ${meta.id}.depth.tsv.gz
    """
}
