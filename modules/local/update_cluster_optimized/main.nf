process UPDATE_CLUSTER_OPTIMIZED {
    tag "update_cluster_optimized"
    label 'process_low'

    //container 'python:3.11-slim'
   // container '/cpfs01/projects-HDD/cfff-86962b7a8e68_HDD/public/singularity_sif/python3.8_v1_Bio_updated.sif'
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path(main_json)         // cluster_data_updated.json
    path(opt_jsons)         // cluster_*_tnf_optimized.json collected

    output:
    path("cluster_data_tnf.json"), emit: updated_json
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    update_cluster_optimized.py \\
        --main_json  ${main_json} \\
        --opt_jsons  ${opt_jsons} \\
        --out_json   cluster_data_tnf.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
