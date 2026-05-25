process UPDATE_CLUSTER_JSON {
    tag "update_cluster_json"
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path cluster_json
    path checkm2_tsv
    path contig_list

    output:
    path "updated_clusters.json",        emit: updated_json
    path "per_cluster_json/*.json",      emit: filtered_clusters
    path "versions.yml",                 emit: versions

    script:
    """
    update_cluster_json.py \\
        --json          ${cluster_json} \\
        --checkm2       ${checkm2_tsv} \\
        --contig-list   ${contig_list} \\
        --out           updated_clusters.json \\
        --per-cluster-dir per_cluster_json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
