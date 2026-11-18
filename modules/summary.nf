process SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: params.publish_dir_mode

    input:
    path coverage
    path depth
    path consensus
    path variants

    output:
    path "sample_summary.tsv", emit: samples
    path "all_variants.tsv",   emit: variants

    script:
    """
    summarize_samples.py --coverage ${coverage} --depth ${depth} --consensus ${consensus} \\
        --variants ${variants} --min-depth ${params.min_depth} \\
        --min-genome-fraction ${params.min_genome_fraction}
    """
}

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path inputs
    path config

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_report_data",  emit: data

    script:
    """
    multiqc --force --config ${config} --filename multiqc_report.html .
    """
}
