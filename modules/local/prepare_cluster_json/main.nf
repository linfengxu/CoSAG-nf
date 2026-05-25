process PREPARE_CLUSTER_JSON {
    tag "prepare_cluster_json"
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path cluster_tsv
    path samplesheet
    path sag_contigs_mapping
    path checkm2_tsv

    output:
    path "round1_clusters.json", emit: json
    path "versions.yml",         emit: versions

    script:
    """
    prepare_cluster_json.py \\
        --samples   ${samplesheet} \\
        --clusters  ${cluster_tsv} \\
        --contigs   ${sag_contigs_mapping} \\
        --checkm2   ${checkm2_tsv} \\
        --outdir    . \\
        --out       round1_clusters.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
