process UPDATE_GTDBTK_RESULTS {
    tag "update_gtdbtk_${stage}"
    label 'process_low'

    //container 'python:3.11-slim'
    container 'quay.io/xulf2022/python3.8_bio:v1'
    //container '/cpfs01/projects-HDD/cfff-86962b7a8e68_HDD/public/singularity_sif/python3.8_v1_Bio_updated.sif'
    input:
    path(main_json)
    path(bac_summary)    // gtdbtk.bac120.summary.tsv  (optional: true)
    path(ar_summary)     // gtdbtk.ar53.summary.tsv    (optional: true)
    val(stage)           // 'sag' | 'coassembly'

    output:
    path("cluster_data_gtdbtk.json"), emit: updated_json
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def bac_arg  = !bac_summary.name.startsWith('NO_') ? "--bac_summary ${bac_summary}" : ''
    def ar_arg   = !ar_summary.name.startsWith('NO_')  ? "--ar_summary  ${ar_summary}"  : ''
    """
    update_gtdbtk_results.py \\
        --main_json  ${main_json} \\
        ${bac_arg} \\
        ${ar_arg} \\
        --out_json   cluster_data_gtdbtk.json \\
        --stage      ${stage} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
