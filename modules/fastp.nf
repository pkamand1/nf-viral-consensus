process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/qc/fastp", mode: params.publish_dir_mode, pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_{1,2}.trim.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),          emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),          emit: html

    script:
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_1.trim.fastq.gz --out2 ${meta.id}_2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred ${params.min_base_quality} \\
        --cut_tail --cut_tail_mean_quality ${params.min_base_quality} \\
        --length_required ${params.min_read_length} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --report_title "${meta.id}"
    """
}
