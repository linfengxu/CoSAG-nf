process EXTRACT_COSAG_CONTIGS {
    tag "extract_cosag_contigs"
    label 'process_low'

    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path(cluster_tnf_json)   // cluster_data_tnf.json

    output:
    path("selected_fastas/*.fasta"), emit: fastas
    path("selected_fastas/selected_contigs.tsv"), emit: manifest
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p selected_fastas

    extract_cosag_contigs.py \\
        --json               ${cluster_tnf_json} \\
        --outdir             selected_fastas \\
        --min_completeness   ${params.min_completeness} \\
        --max_contamination  ${params.max_contamination} \\
        --mode               copy \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
