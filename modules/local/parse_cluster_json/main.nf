process PARSE_CLUSTER_JSON {
    tag "parse_cluster_json"
    label 'process_low'

    //container 'python:3.11-slim'
    container 'quay.io/xulf2022/python3.8_bio:v1'
    input:
    path(cluster_json)

    output:
    path("clusters.tsv"), emit: tsv
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    parse_cluster_json.py \\
        --cluster_json ${cluster_json} \\
        --out_tsv      clusters.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}