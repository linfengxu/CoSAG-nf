process FILTER_SAGS {
    tag "filter_sags"
    label 'process_low'

    //container 'python:3.11-slim'   // 替换为你的 Python Dockerfile 镜像
    container 'quay.io/xulf2022/python3.8_bio:v1'
    input:
    path(quality_report)   // CheckM2 quality_report.tsv

    output:
    path("sags_pass.tsv"),     emit: pass_tsv
    path("sags_fail.tsv"),     emit: fail_tsv, optional: true
    path("filter_stats.json"), emit: stats
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    filter_sags.py \\
        --input              ${quality_report} \\
        --output_pass        sags_pass.tsv \\
        --output_fail        sags_fail.tsv \\
        --stats_json         filter_stats.json \\
        --max_contamination  ${params.max_contamination} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}