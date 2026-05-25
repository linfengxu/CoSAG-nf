process SOURMASH_QUALITY_FILTER {
    tag "sourmash_quality_filter"
    label 'process_low'

    //container 'quay.io/biocontainers/python:3.11'
    //container '/cpfs01/projects-HDD/cfff-86962b7a8e68_HDD/public/singularity_sif/python3.8_v1_Bio_updated.sif'
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path(similarity_matrix_csv)

    output:
    path("filtered_similarity_matrix.csv"), emit: filtered_matrix
    path("connectivity_report.txt"),        emit: report
    path("low_quality_sags.txt"),           emit: low_quality,  optional: true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sourmash_quality_filter.py \\
        --matrix_csv       ${similarity_matrix_csv} \\
        --out_matrix       filtered_similarity_matrix.csv \\
        --out_report       connectivity_report.txt \\
        --out_low_quality  low_quality_sags.txt \\
        --min_similarity   ${params.sourmash_min_similarity} \\
        --min_connections  ${params.sourmash_min_connections} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}