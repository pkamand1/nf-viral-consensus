process BOWTIE2_BUILD {
    tag "${fasta.simpleName}"
    label 'process_medium'

    input:
    path fasta

    output:
    path "bowtie2_index", emit: index

    script:
    """
    mkdir bowtie2_index
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_index/reference
    """
}

process BOWTIE2_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.log'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.sam"),         emit: sam
    tuple val(meta), path("${meta.id}.bowtie2.log"), emit: log

    script:
    """
    bowtie2 \\
        -x ${index}/reference \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        --very-sensitive --no-unal \\
        --rg-id ${meta.id} --rg SM:${meta.id} \\
        --threads ${task.cpus} \\
        -S ${meta.id}.sam \\
        2> ${meta.id}.bowtie2.log
    """
}
