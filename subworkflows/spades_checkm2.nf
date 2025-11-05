#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPADES } from '../modules/spades'
include { CHECKM2 } from '../modules/checkm2'

process PREPARE_CLUSTER_READS {
    tag "cluster_${cluster_id}"
    publishDir "${params.output_structure?.co_assemblies ?: params.outdir}/cluster_reads", mode: 'copy'
    
    input:
    tuple val(cluster_id), val(cluster_info)
    
    output:
    tuple val(cluster_id), val(cluster_info), path("cluster_${cluster_id}_R1.fastq"), path("cluster_${cluster_id}_R2.fastq"), emit: reads
    path "cluster_${cluster_id}_prepare.log", emit: log
    
    script:
    def read1_files = cluster_info.members.collect { "\"${it.read1}\"" }.join(' ')
    def read2_files = cluster_info.members.collect { "\"${it.read2}\"" }.join(' ')
    def member_list = cluster_info.members.collect { it.sag_id }.join(',')
    """
    #!/bin/bash
    
    echo "Preparing reads for cluster ${cluster_id}" > cluster_${cluster_id}_prepare.log
    echo "Cluster size: ${cluster_info.cluster_size}" >> cluster_${cluster_id}_prepare.log
    echo "Number of members: ${cluster_info.members.size()}" >> cluster_${cluster_id}_prepare.log
    echo "Members: ${member_list}" >> cluster_${cluster_id}_prepare.log
    
    # 初始化输出文件
    > cluster_${cluster_id}_R1.fastq
    > cluster_${cluster_id}_R2.fastq
    
    # 合并所有成员的read1文件
    echo "Merging R1 files..." >> cluster_${cluster_id}_prepare.log
    for file in ${read1_files}; do
        if [ -f "\$file" ]; then
            echo "  Adding R1: \$file" >> cluster_${cluster_id}_prepare.log
            cat "\$file" >> cluster_${cluster_id}_R1.fastq
        else
            echo "  WARNING: R1 file not found: \$file" >> cluster_${cluster_id}_prepare.log
        fi
    done
    
    # 合并所有成员的read2文件
    echo "Merging R2 files..." >> cluster_${cluster_id}_prepare.log
    for file in ${read2_files}; do
        if [ -f "\$file" ]; then
            echo "  Adding R2: \$file" >> cluster_${cluster_id}_prepare.log
            cat "\$file" >> cluster_${cluster_id}_R2.fastq
        else
            echo "  WARNING: R2 file not found: \$file" >> cluster_${cluster_id}_prepare.log
        fi
    done
    
    # 统计最终的reads数量
    r1_reads=\$(grep -c "^@" cluster_${cluster_id}_R1.fastq || echo 0)
    r2_reads=\$(grep -c "^@" cluster_${cluster_id}_R2.fastq || echo 0)
    
    echo "Final read counts:" >> cluster_${cluster_id}_prepare.log
    echo "  R1 reads: \$r1_reads" >> cluster_${cluster_id}_prepare.log
    echo "  R2 reads: \$r2_reads" >> cluster_${cluster_id}_prepare.log
    
    if [ \$r1_reads -eq 0 ] || [ \$r2_reads -eq 0 ]; then
        echo "ERROR: No reads found for cluster ${cluster_id}" >> cluster_${cluster_id}_prepare.log
        exit 1
    fi
    """
}

process UPDATE_CLUSTER_JSON {
    tag "update_json"
    container "${params.python.container}"
    publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/prepared_data", mode: 'copy'
    
    input:
    path original_json
    path checkm2_results  // 共组装的CheckM2结果
    path assembly_contigs
    path individual_checkm2_results, stageAs: 'individual_checkm2/*'  // 个体SAG的CheckM2结果
    
    output:
    path "cluster_data_updated.json", emit: updated_json
    path "update.log", emit: log
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import pandas as pd
    import glob
    import os
    from pathlib import Path
    
    # 读取原始JSON
    with open('${original_json}', 'r') as f:
        cluster_data = json.load(f)
    
    # 收集共组装CheckM2结果
    checkm2_files = glob.glob('*_checkm2.tsv')
    checkm2_results = {}
    
    # 收集个体SAG CheckM2结果
    individual_checkm2_files = glob.glob('individual_checkm2/*_checkm2.tsv')
    individual_checkm2_results = {}
    
    # 收集所有组装结果（contigs文件）
    contig_files = glob.glob('*.fasta') + glob.glob('*.fa') + glob.glob('*.contigs.fasta')
    assembly_results = {}
    
    print(f"Found {len(checkm2_files)} co-assembly CheckM2 result files")
    print(f"Found {len(individual_checkm2_files)} individual SAG CheckM2 result files")
    print(f"Found {len(contig_files)} assembly contig files")
    
    # 处理组装结果文件
    for file in contig_files:
        # 从文件名提取cluster_id
        if 'cluster_' in file:
            cluster_id = file.split('cluster_')[1].split('_')[0].split('.')[0]
            assembly_results[cluster_id] = os.path.abspath(file)
            print(f"Found assembly for cluster {cluster_id}: {file}")
    
    # 处理共组装CheckM2结果
    for file in checkm2_files:
        # 从文件名提取cluster_id
        cluster_id = file.replace('cluster_', '').replace('_checkm2.tsv', '')
        
        try:
            df = pd.read_csv(file, sep='\\t')
            if len(df) > 0:
                row = df.iloc[0]  # 取第一行结果
                checkm2_results[cluster_id] = {
                    'completeness': float(row.get('Completeness', 0)),
                    'contamination': float(row.get('Contamination', 0)),
                    'genome_size': int(row.get('Genome_Size', 0)),
                    'gc_content': float(row.get('GC_Content', 0)),
                    'contig_n50': int(row.get('Contig_N50', 0)),
                    'total_coding_sequences': int(row.get('Total_Coding_Sequences', 0))
                }
                print(f"Co-assembly Cluster {cluster_id}: Completeness={row.get('Completeness', 0)}%, Contamination={row.get('Contamination', 0)}%")
            else:
                print(f"Warning: Empty co-assembly CheckM2 result for cluster {cluster_id}")
                checkm2_results[cluster_id] = {
                    'completeness': 0,
                    'contamination': 0,
                    'genome_size': 0,
                    'gc_content': 0,
                    'contig_n50': 0,
                    'total_coding_sequences': 0
                }
        except Exception as e:
            print(f"Error processing co-assembly {file}: {e}")
            checkm2_results[cluster_id] = {
                'completeness': 0,
                'contamination': 0,
                'genome_size': 0,
                'gc_content': 0,
                'contig_n50': 0,
                'total_coding_sequences': 0
            }
    
    # 处理个体SAG CheckM2结果
    for file in individual_checkm2_files:
        # 从文件名提取SAG ID
        sag_id = os.path.basename(file).replace('_checkm2.tsv', '')
        
        try:
            df = pd.read_csv(file, sep='\\t')
            if len(df) > 0:
                row = df.iloc[0]  # 取第一行结果
                individual_checkm2_results[sag_id] = {
                    'completeness': float(row.get('Completeness', 0)),
                    'contamination': float(row.get('Contamination', 0)),
                    'genome_size': int(row.get('Genome_Size', 0)),
                    'gc_content': float(row.get('GC_Content', 0)),
                    'contig_n50': int(row.get('Contig_N50', 0)),
                    'total_coding_sequences': int(row.get('Total_Coding_Sequences', 0))
                }
                print(f"Individual SAG {sag_id}: Completeness={row.get('Completeness', 0)}%, Contamination={row.get('Contamination', 0)}%")
            else:
                print(f"Warning: Empty individual CheckM2 result for SAG {sag_id}")
                individual_checkm2_results[sag_id] = {
                    'completeness': 0,
                    'contamination': 0,
                    'genome_size': 0,
                    'gc_content': 0,
                    'contig_n50': 0,
                    'total_coding_sequences': 0
                }
        except Exception as e:
            print(f"Error processing individual SAG {file}: {e}")
            individual_checkm2_results[sag_id] = {
                'completeness': 0,
                'contamination': 0,
                'genome_size': 0,
                'gc_content': 0,
                'contig_n50': 0,
                'total_coding_sequences': 0
            }
    
    # 更新cluster数据
    for cluster in cluster_data['clusters']:
        cluster_id = cluster['cluster_id']
        
        # 准备round_1结果
        round_1_results = {
            'completeness': 0,
            'contamination': 0,
            'genome_size': 0,
            'gc_content': 0,
            'contig_n50': 0,
            'total_coding_sequences': 0,
            'assembly_status': 'failed',
            'assembly_file': None
        }
        
        # 获取该cluster的CheckM2结果
        cluster_checkm2_data = checkm2_results.get(cluster_id, {
            'completeness': 0,
            'contamination': 0,
            'genome_size': 0,
            'gc_content': 0,
            'contig_n50': 0,
            'total_coding_sequences': 0
        })
        
        # 添加CheckM2结果到cluster级别
        if cluster_id in checkm2_results:
            round_1_results.update({
                'completeness': checkm2_results[cluster_id]['completeness'],
                'contamination': checkm2_results[cluster_id]['contamination'],
                'genome_size': checkm2_results[cluster_id]['genome_size'],
                'gc_content': checkm2_results[cluster_id]['gc_content'],
                'contig_n50': checkm2_results[cluster_id]['contig_n50'],
                'total_coding_sequences': checkm2_results[cluster_id]['total_coding_sequences'],
                'assembly_status': 'completed'
            })
        
        # 添加组装文件路径
        if cluster_id in assembly_results:
            round_1_results['assembly_file'] = assembly_results[cluster_id]
            print(f"Added assembly file for cluster {cluster_id}: {assembly_results[cluster_id]}")
        
        cluster['co_assembly_results'] = {
            'round_1': round_1_results
        }
        
        # 为每个member添加个体SAG的CheckM2质量信息
        if 'members' in cluster:
            for member in cluster['members']:
                sag_id = member.get('sag_id', '')
                
                # 添加个体SAG的质量数据（每个member独有）
                individual_data = individual_checkm2_results.get(sag_id, {
                    'completeness': 0,
                    'contamination': 0
                })
                member['completeness'] = individual_data.get('completeness', 0)
                member['contamination'] = individual_data.get('contamination', 0)
                
                print(f"Added individual CheckM2 quality data to member {sag_id}: completeness={member['completeness']:.1f}%, contamination={member['contamination']:.1f}%")
    
    # 写入更新后的JSON
    with open('cluster_data_updated.json', 'w') as f:
        json.dump(cluster_data, f, indent=2)
    
    # 写入详细日志
    with open('update.log', 'w') as f:
        f.write(f"Updated cluster data with CheckM2 results\\n")
        f.write(f"Total clusters: {len(cluster_data['clusters'])}\\n")
        f.write(f"CheckM2 results processed: {len(checkm2_results)}\\n\\n")
        
        total_members_updated = 0
        for cluster in cluster_data['clusters']:
            cluster_id = cluster['cluster_id']
            members_count = len(cluster.get('members', []))
            total_members_updated += members_count
            
            if 'co_assembly_results' in cluster:
                results = cluster['co_assembly_results']['round_1']
                f.write(f"Cluster {cluster_id}:\\n")
                f.write(f"  - Members: {members_count}\\n")
                f.write(f"  - Completeness: {results['completeness']:.1f}%\\n")
                f.write(f"  - Contamination: {results['contamination']:.1f}%\\n")
                f.write(f"  - Genome size: {results['genome_size']} bp\\n")
                f.write(f"  - Assembly status: {results['assembly_status']}\\n")
                if results['assembly_file']:
                    f.write(f"  - Assembly file: {results['assembly_file']}\\n")
                f.write(f"\\n")
        
        f.write(f"Total members updated with CheckM2 data: {total_members_updated}\\n")
    
    print("JSON update completed successfully")
    print(f"Updated {len(cluster_data['clusters'])} clusters with CheckM2 quality data")
    
    # 统计更新的member数量
    total_members = sum(len(cluster.get('members', [])) for cluster in cluster_data['clusters'])
    individual_results_count = len(individual_checkm2_results)
    
    print(f"Added individual SAG CheckM2 data to {total_members} members")
    print(f"Individual SAG CheckM2 results processed: {individual_results_count}")
    """
}

process PARSE_CLUSTER_JSON {
    container "${params.python.container}"
    input:
    path cluster_json
    
    output:
    path "cluster_*.txt", emit: cluster_files
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    
    # 读取JSON文件
    with open('${cluster_json}', 'r') as f:
        data = json.load(f)
    
    # 为每个cluster创建一个文件
    for cluster in data['clusters']:
        cluster_id = cluster['cluster_id']
        filename = f"cluster_{cluster_id}.txt"
        
        with open(filename, 'w') as f:
            f.write(f"cluster_id\\t{cluster_id}\\n")
            f.write(f"cluster_size\\t{cluster['cluster_size']}\\n")
            f.write("members\\n")
            
            for member in cluster['members']:
                f.write(f"{member['sag_id']}\\t{member['read1']}\\t{member['read2']}\\n")
    """
}

process FILTER_HIGH_CONTAM {
    tag "filter_high_contam"
    container "${params.python.container}"
    publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/prepared_data/high_contam", mode: 'copy'

    
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
    python ${projectDir}/bin/filter_high_contam.py ${cluster_json} output_dir > filter_log.txt 2>&1
    """
}
workflow SPADES_CHECKM2 {
    take:
    cluster_json
    individual_checkm2_results  // 来自main.nf中个体SAG的CheckM2结果

    
    main:
    // 解析JSON文件
    PARSE_CLUSTER_JSON(cluster_json)
    
    // 创建cluster channel
    cluster_ch = PARSE_CLUSTER_JSON.out.cluster_files
        .flatten()
        .map { file ->
            def lines = file.text.split('\n')
            def cluster_id = lines[0].split('\t')[1]
            def cluster_size = lines[1].split('\t')[1] as Integer
            
            def members = []
            for (int i = 3; i < lines.size(); i++) {
                if (lines[i].trim()) {
                    def parts = lines[i].split('\t')
                    if (parts.size() >= 3) {
                        members.add([
                            sag_id: parts[0],
                            read1: parts[1],
                            read2: parts[2]
                        ])
                    }
                }
            }
            
            def cluster_info = [
                cluster_id: cluster_id,
                cluster_size: cluster_size,
                members: members
            ]
            
            tuple(cluster_id, cluster_info)
        }
    
    // 准备每个cluster的reads
    PREPARE_CLUSTER_READS(cluster_ch)
    
    // 对每个cluster进行SPADES组装
    spades_input = PREPARE_CLUSTER_READS.out.reads
        .map { cluster_id, cluster_info, r1, r2 ->
            def meta = [id: "cluster_${cluster_id}"]
            tuple(meta, r1, r2)
        }
    
    SPADES(spades_input)
    
    // 对组装结果进行CheckM2质量评估
    CHECKM2(SPADES.out.contigs)
    
    // 收集所有CheckM2结果和组装结果并更新JSON
    checkm2_results = CHECKM2.out.results
        .map { meta, tsv -> tsv }
        .collect()
    
    assembly_contigs = SPADES.out.contigs
        .map { meta, contigs -> contigs }
        .collect()
    
    UPDATE_CLUSTER_JSON(
        cluster_json,
        checkm2_results,
        assembly_contigs,
        individual_checkm2_results
    )
    FILTER_HIGH_CONTAM (UPDATE_CLUSTER_JSON.out.updated_json)
    emit:
    updated_json = UPDATE_CLUSTER_JSON.out.updated_json
    filtered_json = FILTER_HIGH_CONTAM.out.filtered_clusters  // 新增输出
    assemblies = SPADES.out.contigs
    checkm2_results = CHECKM2.out.results
    logs = UPDATE_CLUSTER_JSON.out.log
}