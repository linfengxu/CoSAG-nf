#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PREPARE_CLUSTER_JSON {
    tag "prepare_cluster_json"
    //publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/prepared_data", mode: 'copy'
    container "${params.python.container}"
    
    input:
    path cluster_tsv
    path samplesheet
    path sag_contigs_mapping  // SAG ID到contigs文件的映射
    
    output:
    path "cluster_data.json", emit: cluster_json
    path "prepare.log", emit: log
    
    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import json
    import sys
    from pathlib import Path
    import os
    
    # 读取聚类结果
    print("Reading cluster data from ${cluster_tsv}")
    cluster_df = pd.read_csv("${cluster_tsv}", sep='\\t')
    
    # 读取样本表
    print("Reading sample sheet from ${samplesheet}")
    sample_df = pd.read_csv("${samplesheet}", sep='\\t')
    
    # 创建样本ID到reads文件的映射
    sample_to_reads = {}
    for _, row in sample_df.iterrows():
        sample_to_reads[row['sampleID']] = {
            'read1': row['forwardReads'],
            'read2': row['reverseReads']
        }
    
    # 读取SAG ID到contigs文件的映射
    sag_to_contigs = {}
    
    print(f"Reading SAG contigs mapping from ${sag_contigs_mapping}")
    
    with open('${sag_contigs_mapping}', 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split('\\t')
                if len(parts) == 2:
                    sag_id, contigs_path = parts
                    if os.path.exists(contigs_path):
                        sag_to_contigs[sag_id] = os.path.abspath(contigs_path)
                        print(f"Mapped SAG {sag_id} to contigs: {contigs_path}")
                    else:
                        print(f"Warning: Contigs file not found for {sag_id}: {contigs_path}")
    
    print(f"Successfully loaded {len(sag_to_contigs)} SAG-contigs mappings")
    
    # 构建JSON数据结构
    cluster_data = {}
    
    for _, row in cluster_df.iterrows():
        sag_id = row['SAG_ID']
        cluster_id = str(row['Cluster_ID'])
        cluster_size = row['Cluster_Size']
        
        # 初始化cluster如果不存在
        if cluster_id not in cluster_data:
            cluster_data[cluster_id] = {
                'cluster_id': cluster_id,
                'cluster_size': cluster_size,
                'members': []
            }
        
        # 添加成员信息
        member_info = {
            'sag_id': sag_id
        }
        
        # 添加reads信息（如果在样本表中找到）
        if sag_id in sample_to_reads:
            member_info['read1'] = sample_to_reads[sag_id]['read1']
            member_info['read2'] = sample_to_reads[sag_id]['read2']
        else:
            print(f"Warning: SAG_ID {sag_id} not found in sample sheet")
            member_info['read1'] = None
            member_info['read2'] = None
        
        # 添加contigs文件路径
        if sag_id in sag_to_contigs:
            member_info['individual_contigs'] = sag_to_contigs[sag_id]
            print(f"Found contigs for {sag_id}: {sag_to_contigs[sag_id]}")
        else:
            member_info['individual_contigs'] = None
            print(f"Warning: No contigs found for {sag_id}")
        
        cluster_data[cluster_id]['members'].append(member_info)
    
    # 转换为列表格式
    output_data = {
        'clusters': list(cluster_data.values()),
        'total_clusters': len(cluster_data),
        'total_sags': len(cluster_df)
    }
    
    # 写入JSON文件
    with open('cluster_data.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Successfully created cluster_data.json with {len(cluster_data)} clusters")
    
    # 写入日志
    with open('prepare.log', 'w') as f:
        f.write(f"Cluster preparation completed\\n")
        f.write(f"Total clusters: {len(cluster_data)}\\n")
        f.write(f"Total SAGs: {len(cluster_df)}\\n")
        for cluster_id, cluster_info in cluster_data.items():
            f.write(f"Cluster {cluster_id}: {len(cluster_info['members'])} members\\n")
    """
}