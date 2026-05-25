
process PREPARE_ROUND_JSON {
    tag "prepare_round${round}_json"
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path prev_json      // 上一轮 updated_clusters.json
    path clusters_tsv   // 新一轮聚类结果 TSV
    path round1_json    // Round 1 JSON（Round 3+ 传入，Round 2 传 NO_FILE）
    val  round

    output:
    path "round${round}_clusters.json",  emit: json
    path "cluster_*/merged_R1.fastq.gz", emit: merged_r1
    path "cluster_*/merged_R2.fastq.gz", emit: merged_r2
    path "versions.yml",                 emit: versions

    script:
    def r1_arg = round1_json.name != 'NO_FILE' ? "--round1-json ${round1_json}" : ""
    """
    prepare_round_json.py \\
        --prev-json   ${prev_json} \\
        --clusters    ${clusters_tsv} \\
        --round       ${round} \\
        --outdir      . \\
        --out         round${round}_clusters.json \\
        --threads     ${task.cpus} \\
        ${r1_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
