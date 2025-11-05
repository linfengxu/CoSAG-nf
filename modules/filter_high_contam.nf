process FILTER_HIGH_CONTAM {
    tag "filter_high_contam"
    publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/prepared_data/high_contam", mode: 'copy'
    container "${params.python.container}"

    
    input:
    path cluster_json
    
    output:
    path "output_dir/cluster_*.json", emit: filtered_clusters
    path "filter_log.txt", emit: log
    
    script:
    """
    # 创建输出目录
    mkdir -p output_dir
    
    # 运行 Python 脚本，传递 JSON 文件和输出目录路径
    python ${projectDir}/bin/filter_cluster.py ${cluster_json} output_dir > filter_log.txt 2>&1
    """
}