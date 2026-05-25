process FORMAT_REPORT {
    tag "format_report"
    label 'process_low'

    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path(final_json)

    output:
    path("cosag_report.html"), emit: report
    path("versions.yml"),      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    format_report.py \\
        ${final_json} \\
        -o cosag_report.html \\
        -t "CoSAG Analysis Report" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
