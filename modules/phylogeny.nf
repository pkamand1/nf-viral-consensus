process MAFFT {
    label 'process_medium'
    publishDir "${params.outdir}/phylogeny", mode: params.publish_dir_mode

    input:
    path consensus
    path reference

    output:
    path "alignment.fasta",      emit: alignment
    path "included_samples.txt", emit: included

    script:
    """
    # Keep consensus sequences that pass the genome-completeness threshold (fraction of N bases)
    for f in ${consensus}; do
        awk -v max_n=${1 - params.min_genome_fraction} '
            /^>/ { h = \$0; next } { s = s \$0 }
            END { n = gsub(/[Nn]/, "", s); if (n / (length(s) + n) <= max_n) { print h; print s } }' "\$f"
    done > consensus_all.fasta
    grep '^>' consensus_all.fasta | sed 's/^>//' > included_samples.txt

    # Align each consensus to the reference, keeping reference coordinates
    mafft --6merpair --thread ${task.cpus} --keeplength --addfragments consensus_all.fasta ${reference} > alignment.fasta
    """
}

process IQTREE {
    label 'process_medium'
    publishDir "${params.outdir}/phylogeny", mode: params.publish_dir_mode

    input:
    path alignment

    output:
    path "tree.treefile", emit: tree
    path "tree.iqtree",   emit: report

    script:
    """
    n_seq=\$(grep -c '^>' ${alignment})
    if [ "\$n_seq" -lt 3 ]; then
        echo "Fewer than three sequences passed QC; no tree was built." > tree.iqtree
        echo "();" > tree.treefile
        exit 0
    fi
    bootstrap=""
    if [ "\$n_seq" -ge 4 ] && [ ${params.tree_bootstrap} -gt 0 ]; then
        bootstrap="-B ${params.tree_bootstrap}"
    fi
    iqtree3 -s ${alignment} -st DNA -m ${params.tree_model} \$bootstrap -T ${task.cpus} --prefix tree --seed 12345 -redo
    """
}
