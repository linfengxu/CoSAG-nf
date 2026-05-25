process FILTER_HIGH_CONTAM {
    tag "filter_high_contam"
    label 'process_low'

    //container 'python:3.11-slim'
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path(updated_json)

    output:
    path("filtered_clusters/*.json"), emit: filtered_clusters
    path("filter_contam.log"),        emit: log
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    filter_high_contam.py \\
        --updated_json      ${updated_json} \\
        --out_dir           filtered_clusters \\
        --out_log           filter_contam.log \\
        --max_contamination ${params.max_contamination} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
