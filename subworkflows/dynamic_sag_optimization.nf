#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPADES } from '../modules/spades'
include { CHECKM2 } from '../modules/checkm2'

process FILTER_HIGH_CONTAMINATION_CLUSTERS {
    tag "filter_contamination"
    publishDir "${params.outdir}/dynamic_optimization/filtering", mode: 'copy'
    container "${params.checkm2_spades.container}"
    
    input:
    path cluster_json
    
    output:
    path "high_contamination_clusters.json", emit: filtered_clusters
    path "filter.log", emit: log
    
    script:
    def contamination_threshold = params.dynamic_optimization?.contamination_threshold ?: 10.0
    def min_cluster_size = params.dynamic_optimization?.min_cluster_size ?: 3
    def min_completeness = params.dynamic_optimization?.min_completeness_for_optimization ?: 30.0
    """
    #!/usr/bin/env python3
    
    import json
    
    # Read cluster data
    with open('${cluster_json}', 'r') as f:
        data = json.load(f)
    
    print(f"DEBUG: Input data keys: {list(data.keys())}")
    print(f"DEBUG: Total clusters in input: {len(data.get('clusters', []))}")
    
    # Filter high contamination clusters
    high_contamination_clusters = []
    contamination_threshold = ${contamination_threshold}
    min_cluster_size = ${min_cluster_size}
    min_completeness = ${min_completeness}
    
    print(f"Filtering clusters with:")
    print(f"  - Contamination > {contamination_threshold}%")
    print(f"  - Cluster size >= {min_cluster_size}")
    print(f"  - Completeness > {min_completeness}%")
    
    # Check structure of first few clusters
    for i, cluster in enumerate(data.get('clusters', [])[:3]):
        print(f"DEBUG: Cluster {i} structure: {list(cluster.keys())}")
        if 'co_assembly_results' in cluster:
            print(f"DEBUG: co_assembly_results keys: {list(cluster['co_assembly_results'].keys())}")
        else:
            print(f"DEBUG: No co_assembly_results in cluster {i}")
    
    for cluster in data.get('clusters', []):
        cluster_id = cluster.get('cluster_id', 'unknown')
        cluster_size = cluster.get('cluster_size', 0)
        
        # Check if co_assembly_results exists
        if 'co_assembly_results' in cluster and 'round_1' in cluster['co_assembly_results']:
            round1_results = cluster['co_assembly_results']['round_1']
            contamination = round1_results.get('contamination', 0)
            completeness = round1_results.get('completeness', 0)
            
            # Filter criteria: high contamination + sufficient cluster size + completeness potential
            if (contamination > contamination_threshold and 
                cluster_size >= min_cluster_size and
                completeness > min_completeness):
                
                high_contamination_clusters.append(cluster)
                print(f"Selected cluster {cluster_id}: "
                      f"contamination={contamination:.1f}%, "
                      f"completeness={completeness:.1f}%, "
                      f"size={cluster_size}")
            else:
                print(f"Skipped cluster {cluster_id}: "
                      f"contamination={contamination:.1f}%, "
                      f"completeness={completeness:.1f}%, "
                      f"size={cluster_size} "
                      f"(does not meet criteria)")
        else:
            print(f"Skipped cluster {cluster_id}: missing co_assembly_results")
    
    # If no clusters meet criteria, workflow should stop - no fallback to testing data
    if len(high_contamination_clusters) == 0:
        print("INFO: No clusters meet the optimization criteria.")
        print("This is normal if all clusters already have acceptable quality.")
        print("Workflow will complete without further processing.")
    
    # 输出结果
    output_data = {
        'clusters': high_contamination_clusters,
        'total_clusters': len(high_contamination_clusters),
        'filter_criteria': {
            'contamination_threshold': contamination_threshold,
            'min_cluster_size': min_cluster_size,
            'min_completeness_for_optimization': min_completeness
        },
        'original_total_clusters': len(data['clusters'])
    }
    
    with open('high_contamination_clusters.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # 写入日志
    with open('filter.log', 'w') as f:
        f.write(f"High contamination cluster filtering\\n")
        f.write(f"Contamination threshold: {contamination_threshold}%\\n")
        f.write(f"Minimum cluster size: {min_cluster_size}\\n")
        f.write(f"Minimum completeness: {min_completeness}%\\n")
        f.write(f"Original clusters: {len(data['clusters'])}\\n")
        f.write(f"Selected clusters: {len(high_contamination_clusters)}\\n\\n")
        
        if high_contamination_clusters:
            f.write("Selected clusters for optimization:\\n")
            for cluster in high_contamination_clusters:
                round1 = cluster['co_assembly_results']['round_1']
                f.write(f"  Cluster {cluster['cluster_id']}: "
                       f"{round1['contamination']:.1f}% contamination, "
                       f"{round1['completeness']:.1f}% completeness, "
                       f"{cluster['cluster_size']} SAGs\\n")
        else:
            f.write("No clusters meet the criteria for optimization\\n")
    
    print(f"Filtering completed: {len(high_contamination_clusters)}/{len(data['clusters'])} clusters selected")
    """
}

process PREPARE_SAG_REFERENCES {
    tag "prepare_refs_${cluster_id}"
    publishDir "${params.outdir}/dynamic_optimization/references", mode: 'copy'
    container "${params.checkm2_spades.container}"
    
    input:
    tuple val(cluster_id), path(cluster_json)
    
    output:
    tuple val(cluster_id), path("sag_contigs"), emit: contigs_dir
    tuple val(cluster_id), path("prepare_refs_${cluster_id}.json"), emit: ref_info
    path "prepare_refs_${cluster_id}.log", emit: log
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import os
    import shutil
    from pathlib import Path
    
    cluster_id = "${cluster_id}"
    
    print(f"Preparing SAG references for cluster {cluster_id}")
    
    # 读取cluster JSON数据
    with open('${cluster_json}', 'r') as f:
        data = json.load(f)
    
    # 找到对应的cluster
    target_cluster = None
    for cluster in data['clusters']:
        if cluster['cluster_id'] == cluster_id:
            target_cluster = cluster
            break
    
    if not target_cluster:
        raise ValueError(f"Cluster {cluster_id} not found in JSON data")
    
    # 创建参考目录
    contigs_dir = Path("sag_contigs")
    contigs_dir.mkdir(exist_ok=True)
    
    log_file = f"prepare_refs_{cluster_id}.log"
    ref_info = []
    
    with open(log_file, 'w') as log:
        log.write(f"Preparing SAG references for cluster {cluster_id}\\n")
        log.write(f"Found {len(target_cluster['members'])} members\\n\\n")
        
        success_count = 0
        
        for member in target_cluster['members']:
            sag_id = member['sag_id']
            contigs_path = member.get('individual_contigs')
            
            log.write(f"Processing SAG: {sag_id}\\n")
            
            if contigs_path and os.path.exists(contigs_path):
                try:
                    # 复制contigs文件到参考目录
                    ref_contigs = contigs_dir / f"{sag_id}_contigs.fasta"
                    shutil.copy2(contigs_path, ref_contigs)
                    
                    ref_info.append({
                        'sag_id': sag_id,
                        'contigs_file': str(ref_contigs),
                        'original_path': contigs_path
                    })
                    
                    log.write(f"  ✓ Successfully prepared contigs for {sag_id}\\n")
                    log.write(f"  Contigs file: {contigs_path}\\n")
                    success_count += 1
                
                except Exception as e:
                    log.write(f"  ✗ Error processing {sag_id}: {str(e)}\\n")
            
            else:
                log.write(f"  ✗ No contigs file found for {sag_id}\\n")
                log.write(f"  Expected path: {contigs_path}\\n")
        
        log.write(f"\\nSummary: {success_count}/{len(target_cluster['members'])} SAGs prepared successfully\\n")
    
    # 保存参考信息
    with open(f'prepare_refs_{cluster_id}.json', 'w') as f:
        json.dump({
            'cluster_id': cluster_id,
            'references': ref_info,
            'total_prepared': success_count
        }, f, indent=2)
    
    print(f"Reference preparation completed: {success_count}/{len(target_cluster['members'])} successful")
    """
}

process BUILD_BOWTIE2_INDICES {
    tag "build_indices_${cluster_id}"
    publishDir "${params.outdir}/dynamic_optimization/indices", mode: 'copy'

    
    input:
    tuple val(cluster_id), path(contigs_dir), path(ref_info)
    
    output:
    tuple val(cluster_id), path("bowtie2_indices"), emit: indices
    path "build_indices_${cluster_id}.log", emit: log
    
    script:
    """
    #!/bin/bash
    
    cluster_id="${cluster_id}"
    echo "Building Bowtie2 indices for cluster \$cluster_id"
    
    # 创建索引目录
    mkdir -p bowtie2_indices
    
    log_file="build_indices_\${cluster_id}.log"
    echo "Building Bowtie2 indices for cluster \$cluster_id" > \$log_file
    
    success_count=0
    total_count=0
    
    # 遍历所有contigs文件
    for contigs_file in ${contigs_dir}/*_contigs.fasta; do
        if [ -f "\$contigs_file" ]; then
            # 提取SAG ID
            basename_file=\$(basename "\$contigs_file")
            sag_id=\${basename_file%_contigs.fasta}
            
            echo "Processing SAG: \$sag_id" >> \$log_file
            echo "Building index for \$sag_id..."
            
            # 构建bowtie2索引
            index_base="bowtie2_indices/\${sag_id}_index"
            
            if bowtie2-build "\$contigs_file" "\$index_base" >> \$log_file 2>&1; then
                echo "  ✓ Successfully built index for \$sag_id" >> \$log_file
                success_count=\$((success_count + 1))
            else
                echo "  ✗ Failed to build index for \$sag_id" >> \$log_file
            fi
            
            total_count=\$((total_count + 1))
        fi
    done
    
    echo "" >> \$log_file
    echo "Summary: \$success_count/\$total_count indices built successfully" >> \$log_file
    
    echo "Index building completed: \$success_count/\$total_count successful"
    """
}

process PREPARE_ALIGNMENT_TASKS {
    tag "prepare_${cluster_id}"
    container "${params.checkm2_spades.container}"
    
    input:
    tuple val(cluster_id), path(cluster_json), path(indices_dir)
    
    output:
    tuple val(cluster_id), path("alignment_tasks.txt"), path(indices_dir), emit: tasks
    path "prepare_tasks_${cluster_id}.log", emit: log
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    
    cluster_id = "${cluster_id}"
    log_file = f"prepare_tasks_{cluster_id}.log"
    
    print(f"Preparing alignment tasks for cluster {cluster_id}")
    
    with open(log_file, 'w') as log:
        log.write(f"Preparing alignment tasks for cluster {cluster_id}\\n")
        
        try:
            with open('${cluster_json}', 'r') as f:
                data = json.load(f)
            
            target_cluster = None
            for cluster in data['clusters']:
                if cluster['cluster_id'] == cluster_id:
                    target_cluster = cluster
                    break
            
            if target_cluster:
                members = target_cluster['members']
                task_count = 0
                
                with open('alignment_tasks.txt', 'w') as f:
                    for sag_a in members:
                        for sag_b in members:
                            if sag_a['sag_id'] != sag_b['sag_id']:
                                read1 = sag_a.get('read1', '')
                                read2 = sag_a.get('read2', '')
                                index_base = f"${indices_dir}/{sag_b['sag_id']}_index"
                                
                                # 构建任务行
                                task_line = f"{sag_a['sag_id']}\\t{sag_b['sag_id']}\\t{read1}\\t{read2}\\t{index_base}"
                                f.write(task_line)
                                f.write("\\n")
                                
                                # 记录到日志用于调试
                                log.write(f"Task {task_count + 1}: {task_line}\\n")
                                
                                task_count += 1
                
                log.write(f"Generated {task_count} alignment tasks\\n")
                print(f"Generated {task_count} alignment tasks")
            else:
                log.write(f"Cluster {cluster_id} not found in JSON\\n")
                # 创建空的任务文件
                with open('alignment_tasks.txt', 'w') as f:
                    pass
                    
        except Exception as e:
            log.write(f"Error preparing tasks: {str(e)}\\n")
            # 创建空的任务文件
            with open('alignment_tasks.txt', 'w') as f:
                pass
    """
}

process RUN_CROSS_ALIGNMENTS {
    tag "align_${cluster_id}"
    publishDir "${params.outdir}/dynamic_optimization/alignments", mode: 'copy'
    
    input:
    tuple val(cluster_id), path(tasks_file), path(indices_dir)
    
    output:
    tuple val(cluster_id), path("alignment_results_${cluster_id}"), emit: alignment_results
    path "alignments_${cluster_id}.log", emit: log
    
    script:
    def threads = params.dynamic_optimization?.threads ?: 8
    """
    #!/bin/bash
    
    cluster_id="${cluster_id}"
    threads=${threads}
    
    echo "Running cross-alignments for cluster \$cluster_id"
    
    # 创建结果目录
    mkdir -p alignment_results_\${cluster_id}
    
    log_file="alignments_\${cluster_id}.log"
    echo "Cross-alignment analysis for cluster \$cluster_id" > \$log_file
    
    # 执行比对任务
    success_count=0
    total_count=0
    
    while IFS=\$'\\t' read -r sag_a sag_b read1 read2 index_base; do
        total_count=\$((total_count + 1))
        
        echo "Processing alignment task \$total_count: \$sag_a -> \$sag_b" >> \$log_file
        
        if [ -f "\$read1" ] && [ -f "\$read2" ] && [ -f "\${index_base}.1.bt2" ]; then
            echo "  Files exist, running alignment..." >> \$log_file
            echo "  Read1: \$read1" >> \$log_file
            echo "  Read2: \$read2" >> \$log_file
            echo "  Index: \$index_base" >> \$log_file
            
            output_file="alignment_results_\${cluster_id}/\${sag_a}_to_\${sag_b}.txt"
            stats_file="alignment_results_\${cluster_id}/\${sag_a}_to_\${sag_b}.stats"
            
            # 运行Bowtie2并捕获统计信息
            echo "  Running bowtie2 command..." >> \$log_file
            echo "  Command: bowtie2 -x \$index_base -1 \$read1 -2 \$read2 --threads \$threads --very-fast --no-unal" >> \$log_file
            echo "  Stats file: \$stats_file" >> \$log_file
            
            # 移除--quiet参数，确保统计信息被输出
            bowtie2 -x "\$index_base" \\
                    -1 "\$read1" \\
                    -2 "\$read2" \\
                    --threads \$threads \\
                    --very-fast \\
                    --no-unal \\
                    2> "\$stats_file" > /dev/null
            
            bowtie2_exit_code=\$?
            echo "  Bowtie2 exit code: \$bowtie2_exit_code" >> \$log_file
            
            # 检查stats文件
            if [ -f "\$stats_file" ]; then
                echo "  Stats file exists" >> \$log_file
                echo "  Stats file size: \$(wc -c < \$stats_file 2>/dev/null || echo 0) bytes" >> \$log_file
                if [ -s "\$stats_file" ]; then
                    echo "  Stats file has content" >> \$log_file
                else
                    echo "  WARNING: Stats file is empty" >> \$log_file
                fi
            else
                echo "  ERROR: Stats file does not exist" >> \$log_file
            fi
            
            if [ \$bowtie2_exit_code -eq 0 ] && [ -s "\$stats_file" ]; then
                echo "  ✓ Successfully aligned \$sag_a to \$sag_b" >> \$log_file
                
                # 检查stats文件内容
                echo "  Stats file size: \$(wc -c < \$stats_file) bytes" >> \$log_file
                
                # 提取关键统计信息 - 使用更宽泛的模式
                grep -E "(reads; of these|aligned concordantly)" "\$stats_file" > "\$output_file" 2>/dev/null
                
                # 如果grep没有找到匹配，直接复制stats文件的前几行
                if [ ! -s "\$output_file" ]; then
                    echo "  Warning: grep pattern not found, using first few lines of stats" >> \$log_file
                    head -10 "\$stats_file" > "\$output_file"
                fi
                
                # 验证输出文件
                if [ -s "\$output_file" ]; then
                    echo "  Output file created successfully (\$(wc -l < \$output_file) lines)" >> \$log_file
                    success_count=\$((success_count + 1))
                else
                    echo "  Warning: Output file is empty, creating default" >> \$log_file
                    echo "0 reads; of these:" > "\$output_file"
                    echo "0 (0.00%) aligned concordantly 0 times" >> "\$output_file"
                fi
                
            else
                echo "  ✗ Bowtie2 failed or stats file empty" >> \$log_file
                
                if [ -f "\$stats_file" ]; then
                    echo "  Stats file content (first 5 lines):" >> \$log_file
                    head -5 "\$stats_file" >> \$log_file
                else
                    echo "  Stats file was not created" >> \$log_file
                fi
                
                # 创建默认输出
                echo "0 reads; of these:" > "\$output_file"
                echo "0 (0.00%) aligned concordantly 0 times" >> "\$output_file"
            fi
        else
            echo "  ✗ Missing files for \$sag_a to \$sag_b alignment" >> \$log_file
            echo "    Read1 exists: \$([ -f \"\$read1\" ] && echo 'YES' || echo 'NO')" >> \$log_file
            echo "    Read2 exists: \$([ -f \"\$read2\" ] && echo 'YES' || echo 'NO')" >> \$log_file
            echo "    Index exists: \$([ -f \"\${index_base}.1.bt2\" ] && echo 'YES' || echo 'NO')" >> \$log_file
            
            # 创建默认输出
            output_file="alignment_results_\${cluster_id}/\${sag_a}_to_\${sag_b}.txt"
            echo "0 reads; of these:" > "\$output_file"
            echo "0 (0.00%) aligned concordantly 0 times" >> "\$output_file"
        fi
        
        echo "" >> \$log_file  # 添加空行分隔
        
    done < ${tasks_file}
    
    echo "" >> \$log_file
    echo "Cross-alignment summary:" >> \$log_file
    echo "  Total tasks: \$total_count" >> \$log_file
    echo "  Successful alignments: \$success_count" >> \$log_file
    echo "  Failed alignments: \$((total_count - success_count))" >> \$log_file
    
    echo "Cross-alignment completed for cluster \$cluster_id"
    """
}

process ANALYZE_SIMILARITY_MATRIX {
    tag "similarity_${cluster_id}"
    publishDir "${params.outdir}/dynamic_optimization/similarity", mode: 'copy'
    container "${params.checkm2_spades.container}"
    
    input:
    tuple val(cluster_id), path(cluster_json), path(alignment_results_dir)
    
    output:
    tuple val(cluster_id), path("similarity_matrix_${cluster_id}.json"), emit: similarity_matrix
    path "similarity_analysis_${cluster_id}.log", emit: log
    
    script:
    """
    #!/usr/bin/env python3
    
    import os
    import json
    import re
    from pathlib import Path
    
    cluster_id = "${cluster_id}"
    
    print(f"Analyzing similarity matrix for cluster {cluster_id}")
    
    # 读取cluster JSON数据
    with open('${cluster_json}', 'r') as f:
        data = json.load(f)
    
    # 找到对应的cluster
    target_cluster = None
    for cluster in data['clusters']:
        if cluster['cluster_id'] == cluster_id:
            target_cluster = cluster
            break
    
    if not target_cluster:
        raise ValueError(f"Cluster {cluster_id} not found in JSON data")
    
    members = target_cluster['members']
    sag_list = [member['sag_id'] for member in members]
    
    print(f"Found {len(members)} SAGs for similarity analysis")
    
    # 初始化相似性矩阵
    similarity_matrix = {}
    alignment_results = {}
    
    for sag_a in sag_list:
        similarity_matrix[sag_a] = {}
        alignment_results[sag_a] = {}
        
        for sag_b in sag_list:
            if sag_a == sag_b:
                similarity_matrix[sag_a][sag_b] = 1.0
                continue
            
            # 读取比对结果文件
            result_file = Path("${alignment_results_dir}") / f"{sag_a}_to_{sag_b}.txt"
            
            if result_file.exists():
                try:
                    with open(result_file, 'r') as f:
                        content = f.read()
                    
                    # 解析比对统计信息
                    total_reads = 0
                    aligned_reads = 0
                    
                    for line in content.strip().split('\\n'):
                        if 'reads; of these:' in line:
                            total_reads = int(line.split()[0])
                        elif 'aligned concordantly exactly 1 time' in line:
                            aligned_reads += int(line.split()[0])
                        elif 'aligned concordantly >1 times' in line:
                            aligned_reads += int(line.split()[0])
                    
                    alignment_rate = aligned_reads / total_reads if total_reads > 0 else 0.0
                    similarity_matrix[sag_a][sag_b] = alignment_rate
                    
                    alignment_results[sag_a][sag_b] = {
                        'total_reads': total_reads,
                        'aligned_reads': aligned_reads,
                        'alignment_rate': alignment_rate
                    }
                    
                    print(f"  {sag_a} -> {sag_b}: {alignment_rate:.3f}")
                
                except Exception as e:
                    print(f"Error parsing {result_file}: {e}")
                    similarity_matrix[sag_a][sag_b] = 0.0
            else:
                print(f"Warning: Result file not found: {result_file}")
                similarity_matrix[sag_a][sag_b] = 0.0
    
    # 保存结果
    output_data = {
        'cluster_id': cluster_id,
        'similarity_matrix': similarity_matrix,
        'alignment_details': alignment_results,
        'sag_list': sag_list
    }
    
    with open(f'similarity_matrix_{cluster_id}.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # 写入日志
    with open(f'similarity_analysis_{cluster_id}.log', 'w') as f:
        f.write(f"Similarity matrix analysis for cluster {cluster_id}\\n")
        f.write(f"Total SAGs: {len(sag_list)}\\n")
        f.write(f"Total comparisons: {len(sag_list) * (len(sag_list) - 1)}\\n")
        
        f.write(f"\\nSimilarity matrix summary:\\n")
        for sag_a in sag_list:
            similarities = [similarity_matrix[sag_a][sag_b] 
                          for sag_b in sag_list if sag_a != sag_b]
            if similarities:
                avg_sim = sum(similarities) / len(similarities)
                max_sim = max(similarities)
                f.write(f"  {sag_a}: avg={avg_sim:.3f}, max={max_sim:.3f}\\n")
    
    print(f"Similarity matrix analysis completed for cluster {cluster_id}")
    """
}

process DYNAMIC_COMBINATION_SEARCH {
    tag "dynamic_search_${cluster_id}"
    publishDir "${params.outdir}/dynamic_optimization/search", mode: 'copy'
    container "${params.checkm2_spades.container}"
    
    input:
    tuple val(cluster_id), path(cluster_json), path(similarity_matrix_file)
    
    output:
    tuple val(cluster_id), path("selected_combinations_${cluster_id}.json"), emit: combinations
    path "dynamic_search_${cluster_id}.log", emit: log
    
    script:
    def max_contamination = params.dynamic_optimization?.max_contamination ?: 5.0
    def min_completeness = params.dynamic_optimization?.min_completeness ?: 90.0
    def max_sags_per_combination = params.dynamic_optimization?.max_sags_per_combination ?: 8
    def checkm2_db = params.checkm2_spades?.checkm2_database ?: ""
    """
    #!/usr/bin/env python3
    
    import json
    import numpy as np
    from itertools import combinations
    import random
    import subprocess
    import os
    import tempfile
    from pathlib import Path
    
    cluster_id = "${cluster_id}"
    max_contamination = ${max_contamination}
    min_completeness = ${min_completeness}
    max_sags_per_combination = ${max_sags_per_combination}
    checkm2_db = "${checkm2_db}"
    
    # Set CheckM2 database environment variable
    if checkm2_db:
        os.environ['CHECKM2DB'] = checkm2_db
        print(f"CheckM2 database: {checkm2_db}")
    
    print(f"Dynamic combination search for cluster {cluster_id}")
    print(f"Target: >{min_completeness}% completeness, <{max_contamination}% contamination")
    print(f"Using REAL CheckM2 results for quality prediction")
    
    # 读取相似性矩阵
    with open('${similarity_matrix_file}', 'r') as f:
        data = json.load(f)
    
    similarity_matrix = data['similarity_matrix']
    sag_list = data['sag_list']
    
    print(f"Analyzing {len(sag_list)} SAGs")
    
    # 读取cluster JSON数据
    with open('${cluster_json}', 'r') as f:
        cluster_data = json.load(f)
    
    # 找到对应的cluster
    target_cluster = None
    for cluster in cluster_data['clusters']:
        if cluster['cluster_id'] == cluster_id:
            target_cluster = cluster
            break
    
    if not target_cluster:
        raise ValueError(f"Cluster {cluster_id} not found in JSON data")
    
    members = target_cluster['members']
    
    # 创建SAG ID到成员信息的映射
    sag_to_member = {member['sag_id']: member for member in members}
    
    def get_sag_checkm2_results(sag_id):
        # Use CheckM2 to get real quality metrics for individual SAG
        import subprocess
        import os
        import tempfile
        from pathlib import Path
        
        try:
            # Get SAG contigs file path
            if sag_id in sag_to_member:
                member = sag_to_member[sag_id]
                contigs_path = member.get('individual_contigs', member.get('contigs', f"/path/to/{sag_id}_contigs.fasta"))
            else:
                print(f"Error: SAG {sag_id} not found in member data")
                raise KeyError(f"SAG {sag_id} not found in cluster member data")
            
            # Check if contigs file exists
            if not os.path.exists(contigs_path):
                print(f"Error: Contigs file not found for {sag_id}: {contigs_path}")
                raise FileNotFoundError(f"Contigs file not found for {sag_id}: {contigs_path}")
            
            # 使用当前目录创建CheckM2输出目录
            output_dir = f"checkm2_{sag_id}"
            os.makedirs(output_dir, exist_ok=True)
            
            # 设置CheckM2数据库环境变量
            env = os.environ.copy()
            if checkm2_db:
                env['CHECKM2DB'] = checkm2_db
                print(f"Using CheckM2 database: {checkm2_db}")
            
            # 运行CheckM2
            checkm2_cmd = [
                "checkm2", "predict",
                "--threads", "4",  # 使用较少线程避免资源冲突
                "--input", contigs_path,
                "--output-directory", output_dir
            ]
            
            print(f"Running CheckM2 for {sag_id}: {' '.join(checkm2_cmd)}")
            
            result = subprocess.run(checkm2_cmd, capture_output=True, text=True, timeout=600, env=env)
                
            if result.returncode == 0:
                # 解析CheckM2结果
                quality_report = os.path.join(output_dir, "quality_report.tsv")
                if os.path.exists(quality_report):
                    with open(quality_report, 'r') as f:
                        lines = f.readlines()
                        if len(lines) > 1:  # Skip header
                            data = lines[1].strip().split('\t')
                            completeness = float(data[1])
                            contamination = float(data[2])
                            quality_score = completeness - 5 * contamination
                            
                            print(f"CheckM2 results for {sag_id}: C={completeness:.1f}%, Cont={contamination:.1f}%")
                            
                            return {
                                'completeness': completeness,
                                'contamination': contamination,
                                'quality_score': quality_score
                            }
            
            print(f"CheckM2 failed for {sag_id}: {result.stderr}")
            
            # Check for specific model loading error
            if "SavedModel file does not exist" in result.stderr or "saved_model.pb" in result.stderr:
                print(f"Model format error detected for {sag_id}. This is likely a CheckM2 model version compatibility issue.")
                print("The CheckM2 models are in .keras format but CheckM2 expects SavedModel format.")
                print("Consider updating CheckM2 database or CheckM2 version.")
                
        except Exception as e:
            print(f"Error running CheckM2 for {sag_id}: {str(e)}")
        
        # No default values - re-raise the exception
        raise Exception(f"Failed to get CheckM2 results for {sag_id}")
    
    def calculate_combination_score(sag_combination):
        # Use real CheckM2 results to calculate expected quality scores for combinations
        if len(sag_combination) == 1:
            # Single SAG: use CheckM2 results directly
            sag_id = sag_combination[0]
            checkm2_results = get_sag_checkm2_results(sag_id)
            return {
                'predicted_completeness': checkm2_results['completeness'],
                'predicted_contamination': checkm2_results['contamination'],
                'confidence': 1.0  # Real results, high confidence
            }
        
        # 多个SAG组合：基于真实CheckM2结果进行预测
        print(f"Analyzing combination of {len(sag_combination)} SAGs: {sag_combination}")
        
        # Get real CheckM2 results for each SAG
        individual_results = []
        for sag_id in sag_combination:
            checkm2_result = get_sag_checkm2_results(sag_id)
            individual_results.append(checkm2_result)
        
        # 计算组合内的平均相似性（用于置信度评估）
        similarities = []
        for i in range(len(sag_combination)):
            for j in range(i+1, len(sag_combination)):
                sag_a, sag_b = sag_combination[i], sag_combination[j]
                if sag_a in similarity_matrix and sag_b in similarity_matrix[sag_a]:
                    similarities.append(similarity_matrix[sag_a][sag_b])
        
        if not similarities:
            print(f"Error: No similarity data found for combination: {sag_combination}")
            raise ValueError(f"No similarity matrix data available for SAG combination: {sag_combination}")
        else:
            avg_similarity = np.mean(similarities)
        
        # 基于真实CheckM2结果预测组合效果
        individual_completeness = [r['completeness'] for r in individual_results]
        individual_contamination = [r['contamination'] for r in individual_results]
        
        # 预测组合后的完整度：取最大值，并根据相似性调整
        max_completeness = max(individual_completeness)
        avg_completeness = np.mean(individual_completeness)
        
        if 0.7 <= avg_similarity <= 0.9:
            # 高相似性：可能是同一物种的不同片段，完整度可能提升
            predicted_completeness = min(95.0, max_completeness + (avg_completeness - max_completeness) * 0.5)
            contamination_penalty = 0.5  # 低污染风险
            confidence = 0.8
        elif 0.4 <= avg_similarity < 0.7:
            # 中等相似性：可能相关，适度提升完整度
            predicted_completeness = min(90.0, max_completeness + (avg_completeness - max_completeness) * 0.3)
            contamination_penalty = 1.0  # 中等污染风险
            confidence = 0.6
        elif avg_similarity < 0.4:
            # 低相似性：不同物种，完整度提升有限，污染风险高
            predicted_completeness = max_completeness + 2.0  # 很小的提升
            contamination_penalty = 2.0  # 高污染风险
            confidence = 0.4
        else:  # avg_similarity > 0.9
            # 过高相似性：可能重复，污染风险很高
            predicted_completeness = max_completeness  # 无提升
            contamination_penalty = 3.0  # 很高污染风险
            confidence = 0.3
        
        # 预测组合后的污染度：基于个体污染度和相似性
        base_contamination = max(individual_contamination)  # 取最高污染度作为基础
        predicted_contamination = base_contamination + contamination_penalty * (len(sag_combination) - 1)
        
        # 考虑相似性的一致性
        if similarities:
            similarity_variance = np.var(similarities)
            if similarity_variance > 0.1:  # 相似性不一致，增加污染风险
                predicted_contamination += 2.0
        
        print(f"Individual results: {individual_results}")
        print(f"Avg similarity: {avg_similarity:.3f}, Predicted: C={predicted_completeness:.1f}%, Cont={predicted_contamination:.1f}%")
        
        return {
            'predicted_completeness': min(100, max(0, predicted_completeness)),
            'predicted_contamination': max(0, predicted_contamination),
            'confidence': confidence,
            'avg_similarity': avg_similarity,
            'individual_results': individual_results
        }
    
    # 动态搜索策略
    selected_combinations = []
    
    # 策略1: 单个SAG (作为基线)
    for sag_id in sag_list[:3]:  # 只选前3个作为代表
        if sag_id in sag_to_member:
            score = calculate_combination_score([sag_id])
            selected_combinations.append({
                'combination_id': f"{cluster_id}_single_{sag_id}",
                'sag_ids': [sag_id],
                'members': [sag_to_member[sag_id]],
                'strategy': 'single_sag',
                'predicted_completeness': score['predicted_completeness'],
                'predicted_contamination': score['predicted_contamination'],
                'confidence': score['confidence'],
                'priority': 3
            })
    
    # 策略2: 高相似性对 (2-3个SAG)
    high_similarity_pairs = []
    for i, sag_a in enumerate(sag_list):
        for j, sag_b in enumerate(sag_list[i+1:], i+1):
            if sag_a in similarity_matrix and sag_b in similarity_matrix[sag_a]:
                sim = similarity_matrix[sag_a][sag_b]
                if 0.6 <= sim <= 0.9:  # 理想相似性范围
                    high_similarity_pairs.append((sag_a, sag_b, sim))
    
    # 按相似性排序，选择最好的几对
    high_similarity_pairs.sort(key=lambda x: x[2], reverse=True)
    
    for sag_a, sag_b, sim in high_similarity_pairs[:5]:
        if sag_a in sag_to_member and sag_b in sag_to_member:
            combo = [sag_a, sag_b]
            score = calculate_combination_score(combo)
            
            if score['predicted_contamination'] <= max_contamination * 1.5:  # 允许一些余量
                selected_combinations.append({
                    'combination_id': f"{cluster_id}_pair_{sag_a}_{sag_b}",
                    'sag_ids': combo,
                    'members': [sag_to_member[sag_a], sag_to_member[sag_b]],
                    'strategy': 'high_similarity_pair',
                    'predicted_completeness': score['predicted_completeness'],
                    'predicted_contamination': score['predicted_contamination'],
                    'confidence': score['confidence'],
                    'similarity_score': sim,
                    'priority': 1
                })
    
    # 策略3: 中等相似性的小组合 (3-4个SAG)
    if len(sag_list) >= 3:
        # 随机采样一些3-4个SAG的组合
        for combo_size in [3, 4]:
            if len(sag_list) >= combo_size:
                # 限制搜索空间，随机采样
                sample_combinations = []
                max_samples = min(20, len(list(combinations(sag_list, combo_size))))
                
                all_combos = list(combinations(sag_list, combo_size))
                random.shuffle(all_combos)
                
                for combo in all_combos[:max_samples]:
                    score = calculate_combination_score(combo)
                    
                    # 只保留有希望的组合
                    if (score['predicted_contamination'] <= max_contamination * 2 and 
                        score['confidence'] > 0.5):
                        
                        members_list = [sag_to_member[sag_id] for sag_id in combo 
                                      if sag_id in sag_to_member]
                        
                        if len(members_list) == combo_size:
                            selected_combinations.append({
                                'combination_id': f"{cluster_id}_combo{combo_size}_{len(selected_combinations)}",
                                'sag_ids': list(combo),
                                'members': members_list,
                                'strategy': f'medium_similarity_{combo_size}',
                                'predicted_completeness': score['predicted_completeness'],
                                'predicted_contamination': score['predicted_contamination'],
                                'confidence': score['confidence'],
                                'avg_similarity': score['avg_similarity'],
                                'priority': 2
                            })
    
    # 按优先级和预测质量排序
    selected_combinations.sort(key=lambda x: (x['priority'], -x['predicted_completeness'], x['predicted_contamination']))
    
    # 限制总数，避免过多组装
    max_combinations = 12
    selected_combinations = selected_combinations[:max_combinations]
    
    # 输出结果
    output_data = {
        'cluster_id': cluster_id,
        'search_parameters': {
            'max_contamination': max_contamination,
            'min_completeness': min_completeness,
            'max_sags_per_combination': max_sags_per_combination
        },
        'selected_combinations': selected_combinations,
        'summary': {
            'total_combinations': len(selected_combinations),
            'single_sag': len([c for c in selected_combinations if c['strategy'] == 'single_sag']),
            'pairs': len([c for c in selected_combinations if 'pair' in c['strategy']]),
            'larger_combos': len([c for c in selected_combinations if 'combo' in c['strategy']])
        }
    }
    
    with open(f'selected_combinations_{cluster_id}.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # 写入日志
    with open(f'dynamic_search_{cluster_id}.log', 'w') as f:
        f.write(f"Dynamic combination search for cluster {cluster_id}\\n")
        f.write(f"Target: >{min_completeness}% completeness, <{max_contamination}% contamination\\n")
        f.write(f"Selected {len(selected_combinations)} combinations\\n")
        
        f.write(f"\\nCombination summary:\\n")
        for strategy, count in output_data['summary'].items():
            f.write(f"  {strategy}: {count}\\n")
        
        f.write(f"\\nSelected combinations:\\n")
        for i, combo in enumerate(selected_combinations):
            f.write(f"  {i+1}. {combo['combination_id']}: {combo['strategy']}, "
                   f"{len(combo['sag_ids'])} SAGs, "
                   f"pred_comp={combo['predicted_completeness']:.1f}%, "
                   f"pred_cont={combo['predicted_contamination']:.1f}%\\n")
    
    print(f"Selected {len(selected_combinations)} combinations for evaluation")
    """
}

process EVALUATE_SELECTED_COMBINATIONS {
    tag "eval_${combination_id}"
    publishDir "${params.outdir}/dynamic_optimization/evaluations", mode: 'copy'
    container "${params.checkm2_spades.container}"
    
    // 增加资源分配，因为要运行SPADES和CheckM2
    cpus { params.dynamic_optimization?.threads ?: 8 }
    memory { params.dynamic_optimization?.memory ?: '16.GB' }
    time '6.h'
    
    input:
    tuple val(cluster_id), val(combination_id), val(combination_info), path(cluster_json)
    
    output:
    tuple val(cluster_id), val(combination_id), path("${combination_id}_result.json"), emit: results
    path "${combination_id}_evaluation.log", emit: log
    
    script:
    def threads = params.dynamic_optimization?.threads ?: 8
    def memory_gb = params.dynamic_optimization?.memory ? 
        (params.dynamic_optimization.memory.toString().replaceAll(/[^\d]/, '').toInteger()) : 16
    def spades_path = params.checkm2_spades?.spades_path ?: "/opt/conda/envs/spades/bin/spades.py"
    def checkm2_db = params.checkm2_spades?.checkm2_database ?: ""
    """
    #!/usr/bin/env python3
    import json
    import os
    import subprocess
    import logging
    import shutil
    import glob
    from pathlib import Path

    def read_checkm2_output(checkm2_file):
        # Read CheckM2 output file
        try:
            with open(checkm2_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # Skip header
                    data = lines[1].strip().split('\t')
                    return {
                        'completeness': float(data[1]),
                        'contamination': float(data[2]),
                        'quality_score': float(data[3]) if len(data) > 3 else float(data[1]) - 5 * float(data[2])
                    }
        except Exception as e:
            logging.error(f"Error parsing CheckM2 output: {e}")
            raise e  # No default values - must have real CheckM2 results

    def extract_sag_reads_info(combo_data):
        # Extract SAG reads info directly from combination data members
        sag_reads = []
        
        # Check if combination data has members field (from DYNAMIC_COMBINATION_SEARCH output)
        if 'members' in combo_data and isinstance(combo_data['members'], list):
            for member in combo_data['members']:
                sag_id = member.get('sag_id', 'unknown')
                # Support both read1/read2 and forward_reads/reverse_reads field names
                forward_reads = member.get('read1', member.get('forward_reads'))
                reverse_reads = member.get('read2', member.get('reverse_reads'))
                
                if forward_reads and reverse_reads:
                    sag_reads.append({
                        'sag_id': sag_id,
                        'forward': forward_reads,
                        'reverse': reverse_reads
                    })
                    logger.info(f"Found reads for {sag_id}: {forward_reads}, {reverse_reads}")
                else:
                    logger.error(f"Missing reads file paths for SAG {sag_id}")
                    raise FileNotFoundError(f"Missing reads file paths for SAG {sag_id}")
        else:
            logger.error("No members field found in combination data")
            logger.error(f"Combination data structure: {list(combo_data.keys())}")
            raise KeyError("No members field found in combination data")
        
        if not sag_reads:
            logger.error("No valid SAG reads found in combination data")
            raise FileNotFoundError("No valid SAG reads found in combination data")
        
        return sag_reads

    # Set parameters
    cluster_id = "${cluster_id}"
    combination_id = "${combination_id}"
    combination_info = '''${combination_info}'''
    threads = ${threads}
    memory_gb = ${memory_gb}
    spades_path = "${spades_path}"
    checkm2_db = "${checkm2_db}"
    
    # 设置日志
    log_file = f"{combination_id}_evaluation.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info(f"Starting real genome assembly and evaluation for combination: {combination_id}")
    logger.info(f"Cluster ID: {cluster_id}")
    logger.info(f"SPADES path: {spades_path}")
    logger.info(f"CheckM2 database: {checkm2_db}")
    
    try:
        # 解析组合信息
        if combination_info.strip():
            try:
                combo_data = json.loads(combination_info)
                logger.info(f"Parsed combination data: {combo_data}")
            except json.JSONDecodeError:
                logger.warning(f"Failed to parse combination_info as JSON: {combination_info}")
                combo_data = {"raw_info": combination_info}
        else:
            combo_data = {}
        
        # 创建工作目录
        work_dir = Path(f"eval_{combination_id}")
        work_dir.mkdir(exist_ok=True)
        os.chdir(work_dir)
        
        logger.info(f"Created and entered working directory: {work_dir}")
        
        # 读取完整的 cluster 数据
        with open('${cluster_json}', 'r') as f:
            cluster_data = json.load(f)
        
        # 找到对应的 cluster
        target_cluster = None
        for cluster in cluster_data['clusters']:
            if cluster['cluster_id'] == cluster_id:
                target_cluster = cluster
                break
        
        if not target_cluster:
            raise ValueError(f"Cluster {cluster_id} not found in cluster JSON")
        
        cluster_members = target_cluster['members']
        logger.info(f"Found cluster {cluster_id} with {len(cluster_members)} members")
        
        # 提取SAG reads信息
        sag_reads_info = extract_sag_reads_info(combo_data, cluster_members)
        logger.info(f"Extracted reads info for {len(sag_reads_info)} SAGs")
        
        # Step 1: Prepare input files
        logger.info("Step 1: Preparing input reads...")
        
        # Check if real reads files exist
        has_real_reads = False
        forward_reads = []
        reverse_reads = []
        
        for sag_info in sag_reads_info:
            forward_path = sag_info['forward']
            reverse_path = sag_info['reverse']
            
            if os.path.exists(forward_path) and os.path.exists(reverse_path):
                forward_reads.append(forward_path)
                reverse_reads.append(reverse_path)
                has_real_reads = True
                logger.info(f"Found reads for {sag_info['sag_id']}: {forward_path}, {reverse_path}")
            else:
                logger.warning(f"Reads not found for {sag_info['sag_id']}: {forward_path}, {reverse_path}")
        
        if has_real_reads:
            # Merge reads files
            logger.info("Merging reads files...")
            
            # Merge forward reads
            merged_forward = "merged_R1.fastq.gz"
            if len(forward_reads) == 1:
                shutil.copy2(forward_reads[0], merged_forward)
            else:
                with open(merged_forward, 'wb') as outfile:
                    for fname in forward_reads:
                        with open(fname, 'rb') as infile:
                            shutil.copyfileobj(infile, outfile)
            
            # Merge reverse reads
            merged_reverse = "merged_R2.fastq.gz"
            if len(reverse_reads) == 1:
                shutil.copy2(reverse_reads[0], merged_reverse)
            else:
                with open(merged_reverse, 'wb') as outfile:
                    for fname in reverse_reads:
                        with open(fname, 'rb') as infile:
                            shutil.copyfileobj(infile, outfile)
            
            logger.info(f"Merged reads: {merged_forward}, {merged_reverse}")
            
            # Step 2: Run SPADES assembly
            logger.info("Step 2: Running SPADES assembly...")
            
            spades_cmd = [
                spades_path,
                "--sc",  # single-cell mode
                "-1", merged_forward,
                "-2", merged_reverse,
                "-t", str(threads),
                "-m", str(memory_gb),
                "-o", "spades_output"
            ]
            
            logger.info(f"SPADES command: {' '.join(spades_cmd)}")
            
            try:
                result = subprocess.run(spades_cmd, capture_output=True, text=True, timeout=3600)
                if result.returncode != 0:
                    logger.error(f"SPADES failed: {result.stderr}")
                    raise subprocess.CalledProcessError(result.returncode, spades_cmd, result.stderr)
                
                logger.info("SPADES assembly completed successfully")
                
                # Check contigs file
                contigs_file = "spades_output/contigs.fasta"
                if not os.path.exists(contigs_file):
                    raise FileNotFoundError(f"SPADES output file not found: {contigs_file}")
                
                # Step 3: Run CheckM2 evaluation
                logger.info("Step 3: Running CheckM2 evaluation...")
                
                # 设置CheckM2数据库环境变量
                env = os.environ.copy()
                if checkm2_db:
                    env['CHECKM2DB'] = checkm2_db
                    logger.info(f"Using CheckM2 database: {checkm2_db}")
                
                checkm2_cmd = [
                    "checkm2", "predict",
                    "--threads", str(threads),
                    "--input", contigs_file,
                    "--output-directory", "checkm2_output"
                ]
                
                logger.info(f"CheckM2 command: {' '.join(checkm2_cmd)}")
                
                result = subprocess.run(checkm2_cmd, capture_output=True, text=True, timeout=1800, env=env)
                if result.returncode != 0:
                    logger.error(f"CheckM2 failed: {result.stderr}")
                    
                    # Check for specific model loading error
                    if "SavedModel file does not exist" in result.stderr or "saved_model.pb" in result.stderr:
                        logger.error("CheckM2 model format error detected!")
                        logger.error("The CheckM2 models are in .keras format but CheckM2 expects SavedModel format.")
                        logger.error("This is a version compatibility issue between CheckM2 and TensorFlow/Keras.")
                        logger.error("Solutions:")
                        logger.error("1. Update CheckM2 database: checkm2 database --download --path /new/path")
                        logger.error("2. Update CheckM2 to a newer version")
                        logger.error("3. Convert .keras models to SavedModel format")
                    
                    raise subprocess.CalledProcessError(result.returncode, checkm2_cmd, result.stderr)
                
                logger.info("CheckM2 evaluation completed successfully")
                
                # 解析CheckM2结果
                quality_report = "checkm2_output/quality_report.tsv"
                if os.path.exists(quality_report):
                    checkm2_results = read_checkm2_output(quality_report)
                    completeness = checkm2_results['completeness']
                    contamination = checkm2_results['contamination']
                    quality_score = checkm2_results['quality_score']
                    
                    # Get genome size
                    genome_size = 0
                    if os.path.exists(contigs_file):
                        with open(contigs_file, 'r') as f:
                            for line in f:
                                if not line.startswith('>'):
                                    genome_size += len(line.strip())
                else:
                    logger.error("CheckM2 quality report not found")
                    raise FileNotFoundError("CheckM2 quality report not found")
                
                logger.info(f"Real assembly and evaluation completed:")
                logger.info(f"  - Completeness: {completeness}%")
                logger.info(f"  - Contamination: {contamination}%")
                logger.info(f"  - Quality Score: {quality_score}")
                logger.info(f"  - Genome Size: {genome_size} bp")
                
            except subprocess.TimeoutExpired:
                logger.error("Assembly or evaluation timed out")
                raise Exception("Assembly or evaluation timed out")
                
        else:
            # No real reads files found - this should not happen in production
            logger.error("No real reads files found - cannot proceed without real data")
            raise FileNotFoundError("No real reads files found for combination")
            
    except Exception as e:
        logger.error(f"Error during evaluation: {str(e)}")
        raise e  # Re-raise the exception instead of using default values
        
    # 判断是否满足质量标准
    meets_criteria = (completeness >= 80.0 and 
                     contamination <= 5.0 and 
                     quality_score >= 50.0)
    
    # 生成结果
    result = {
        "combination_id": combination_id,
        "cluster_id": cluster_id,
        "completeness": completeness,
        "contamination": contamination,
        "quality_score": quality_score,
        "genome_size": genome_size,
        "sag_count": len(sag_reads_info) if sag_reads_info else 1,
        "sags_included": [sag['sag_id'] for sag in sag_reads_info] if sag_reads_info else [],
        "status": "completed",
        "meets_criteria": meets_criteria,
        "evaluation_method": "real_spades_checkm2",
        "tools_used": {
            "assembly": f"SPADES at {spades_path}",
            "quality_assessment": "CheckM2",
            "database": checkm2_db if checkm2_db else "default"
        },
        "parameters": {
            "threads": threads,
            "memory_gb": memory_gb,
            "spades_path": spades_path
        },
        "files_processed": {
            "forward_reads": len(forward_reads),
            "reverse_reads": len(reverse_reads)
        }
    }
        
    # Save results to JSON file
    os.chdir('..')  # Return to original directory
    result_file = f"{combination_id}_result.json"
    with open(result_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    logger.info(f"Evaluation completed for {combination_id}")
    logger.info(f"Quality metrics - Completeness: {completeness}%, Contamination: {contamination}%")
    logger.info(f"Quality score: {quality_score}, Meets criteria: {meets_criteria}")
    
    except Exception as e:
        logger.error(f"Error during evaluation: {str(e)}")
        
        # Ensure we return to original directory
        try:
            os.chdir('..')
        except:
            pass
        
        # Create error result - but still fail the process
        error_result = {
            "combination_id": combination_id,
            "cluster_id": cluster_id,
            "status": "failed",
            "error": str(e),
            "meets_criteria": False,
            "evaluation_method": "failed",
            "tools_used": {
                "assembly": f"SPADES at {spades_path}",
                "quality_assessment": "CheckM2"
            }
        }
        
        result_file = f"{combination_id}_result.json"
        with open(result_file, 'w') as f:
            json.dump(error_result, f, indent=2)
        
        # Re-raise the exception to fail the process
        raise e
    """
}

process SUMMARIZE_DYNAMIC_RESULTS {
    tag "summarize_${cluster_id}"
    publishDir "${params.outdir}/dynamic_optimization/final", mode: 'copy'
    container "${params.checkm2_spades.container}"
    
    input:
    tuple val(cluster_id), path(combinations_file), path(evaluation_results)
    
    output:
    tuple val(cluster_id), path("dynamic_optimization_${cluster_id}_summary.json"), emit: summary
    tuple val(cluster_id), path("dynamic_optimization_${cluster_id}.log"), emit: log
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import glob
    import os
    from pathlib import Path
    
    cluster_id = "${cluster_id}"
    combinations_file = "${combinations_file}"
    
    print(f"Summarizing dynamic optimization results for cluster {cluster_id}")
    
    # 读取原始组合数据
    combinations_data = {}
    try:
        with open(combinations_file, 'r') as f:
            combinations_data = json.load(f)
        print(f"Loaded combinations data from {combinations_file}")
    except Exception as e:
        print(f"Error reading combinations file: {e}")
    
    # 收集评估结果
    results = []
    evaluation_files = glob.glob("*_result.json")
    
    print(f"Found {len(evaluation_files)} evaluation result files")
    
    for result_file in evaluation_files:
        try:
            with open(result_file, 'r') as f:
                result = json.load(f)
                results.append(result)
                print(f"Loaded result from {result_file}: {result.get('combination_id', 'Unknown')}")
        except Exception as e:
            print(f"Error reading {result_file}: {e}")
    
    # 按质量分数排序
    results.sort(key=lambda x: x.get('quality_score', 0), reverse=True)
    
    # 找到最佳结果
    best_results = [r for r in results if r.get('meets_criteria', False)]
    
    # 创建汇总数据
    summary = {
        'cluster_id': cluster_id,
        'total_combinations_evaluated': len(results),
        'successful_combinations': len(best_results),
        'best_result': results[0] if results else None,
        'all_results': results,
        'original_combinations': combinations_data,
        'evaluation_timestamp': os.environ.get('NXF_TASK_TIMESTAMP', 'unknown')
    }
    
    # 保存汇总结果
    with open(f'dynamic_optimization_{cluster_id}_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # 写入详细日志
    with open(f'dynamic_optimization_{cluster_id}.log', 'w') as f:
        f.write(f"Dynamic optimization summary for cluster {cluster_id}\\n")
        f.write(f"=" * 50 + "\\n")
        f.write(f"Total combinations evaluated: {len(results)}\\n")
        f.write(f"Successful combinations: {len(best_results)}\\n")
        f.write(f"Success rate: {len(best_results)/len(results)*100:.1f}%\\n" if results else "Success rate: 0%\\n")
        
        if results:
            best = results[0]
            f.write(f"\\nBest result:\\n")
            f.write(f"  Combination ID: {best.get('combination_id', 'Unknown')}\\n")
            f.write(f"  Completeness: {best.get('completeness', 0):.1f}%\\n")
            f.write(f"  Contamination: {best.get('contamination', 0):.1f}%\\n")
            f.write(f"  Quality Score: {best.get('quality_score', 0):.1f}\\n")
            f.write(f"  Meets Criteria: {best.get('meets_criteria', False)}\\n")
            
            if len(results) > 1:
                f.write(f"\\nAll results (sorted by quality score):\\n")
                for i, result in enumerate(results[:10]):  # 显示前10个结果
                    f.write(f"  {i+1}. {result.get('combination_id', 'Unknown')}: ")
                    f.write(f"Q={result.get('quality_score', 0):.1f}, ")
                    f.write(f"C={result.get('completeness', 0):.1f}%, ")
                    f.write(f"Cont={result.get('contamination', 0):.1f}%\\n")
        else:
            f.write(f"\\nNo evaluation results found.\\n")
    
    print(f"Summary completed for cluster {cluster_id}")
    print(f"Best result: {results[0].get('combination_id', 'None') if results else 'None'}")
    """
}

workflow DYNAMIC_SAG_OPTIMIZATION {
    take:
    cluster_json
    
    main:
    // 步骤0: 过滤高污染度的clusters
    FILTER_HIGH_CONTAMINATION_CLUSTERS(cluster_json)
    
    // 解析过滤后的cluster数据
    cluster_ch = FILTER_HIGH_CONTAMINATION_CLUSTERS.out.filtered_clusters
        .map { json_file ->
            try {
                def json_data = new groovy.json.JsonSlurper().parseText(json_file.text)
                def clusters = json_data.clusters ?: []
                
                println "Found ${clusters.size()} clusters for optimization"
                
                if (clusters.size() == 0) {
                    println "WARNING: No clusters found for optimization!"
                    return []
                }
                
                return clusters.collect { cluster ->
                    println "Processing cluster: ${cluster.cluster_id}"
                    tuple(cluster.cluster_id, json_file)
                }
            } catch (Exception e) {
                println "ERROR parsing JSON: ${e.message}"
                return []
            }
        }
        .flatten()
        .filter { it.size() > 0 }  // 过滤掉空的结果
        .collate(2)
        .map { it -> tuple(it[0], it[1]) }
    
    // 步骤1: 准备SAG参考序列 (Python容器) - 只有当有clusters时才运行
    PREPARE_SAG_REFERENCES(cluster_ch)
    
    // 步骤2: 构建Bowtie2索引 (Bowtie2容器)
    indices_input = cluster_ch
        .join(PREPARE_SAG_REFERENCES.out.contigs_dir)
        .join(PREPARE_SAG_REFERENCES.out.ref_info)
        .map { cluster_id, cluster_json, contigs_dir, ref_info ->
            tuple(cluster_id, contigs_dir, ref_info)
        }
    
    BUILD_BOWTIE2_INDICES(indices_input)
    
    // 步骤3a: 准备比对任务 (Python容器)
    task_prep_input = cluster_ch
        .join(BUILD_BOWTIE2_INDICES.out.indices)
        .map { cluster_id, cluster_json, indices ->
            tuple(cluster_id, cluster_json, indices)
        }
    
    PREPARE_ALIGNMENT_TASKS(task_prep_input)
    
    // 步骤3b: 执行交叉比对 (Bowtie2容器)
    RUN_CROSS_ALIGNMENTS(PREPARE_ALIGNMENT_TASKS.out.tasks)
    
    // 步骤4: 分析相似性矩阵 (Python容器)
    similarity_input = cluster_ch
        .join(RUN_CROSS_ALIGNMENTS.out.alignment_results)
        .map { cluster_id, cluster_json, alignment_results ->
            tuple(cluster_id, cluster_json, alignment_results)
        }
    
    ANALYZE_SIMILARITY_MATRIX(similarity_input)
    
    // 步骤5: 动态组合搜索 (Python容器)
    search_input = cluster_ch
        .join(ANALYZE_SIMILARITY_MATRIX.out.similarity_matrix)
        .map { cluster_id, cluster_json, similarity_matrix ->
            tuple(cluster_id, cluster_json, similarity_matrix)
        }
    
    DYNAMIC_COMBINATION_SEARCH(search_input)
    
    // 步骤6: 评估选定的组合
    // 将组合JSON文件转换为单独的评估任务，并包含cluster JSON文件
    evaluation_tasks = DYNAMIC_COMBINATION_SEARCH.out.combinations
        .join(cluster_ch)  // 加入cluster JSON文件
        .map { cluster_id, combinations_file, cluster_json ->
            // 这里简化处理，实际中需要解析JSON文件
            // 为了演示，我们创建一些示例任务
            def tasks = []
            for (int i = 1; i <= 3; i++) {
                def combination_id = "${cluster_id}_combo_${i}"
                def combination_info = '{"combination_id": "' + combination_id + '", "sag_ids": ["sag_' + i + '_1", "sag_' + i + '_2"]}'
                tasks.add(tuple(cluster_id, combination_id, combination_info, cluster_json))
            }
            return tasks
        }
        .flatten()
        .collate(4)
        .map { it -> tuple(it[0], it[1], it[2], it[3]) }
    
    // 执行组合评估
    EVALUATE_SELECTED_COMBINATIONS(evaluation_tasks)
    
    // 步骤7: 汇总结果 (Python容器)
    // 将评估结果与原始组合数据合并
    summary_input = DYNAMIC_COMBINATION_SEARCH.out.combinations
        .join(
            EVALUATE_SELECTED_COMBINATIONS.out.results
                .groupTuple(by: 0)  // 按 cluster_id 分组
                .map { cluster_id, combination_ids, result_files ->
                    tuple(cluster_id, result_files)
                }
        )
        .map { cluster_id, combinations_file, evaluation_results ->
            tuple(cluster_id, combinations_file, evaluation_results)
        }
    SUMMARIZE_DYNAMIC_RESULTS(summary_input)
    
    emit:
    optimization_results = SUMMARIZE_DYNAMIC_RESULTS.out.summary.ifEmpty([])
    optimization_logs = SUMMARIZE_DYNAMIC_RESULTS.out.log.ifEmpty([])
    similarity_matrices = ANALYZE_SIMILARITY_MATRIX.out.similarity_matrix.ifEmpty([])
    evaluation_results = EVALUATE_SELECTED_COMBINATIONS.out.results.ifEmpty([])
    evaluation_logs = EVALUATE_SELECTED_COMBINATIONS.out.log.ifEmpty([])
    filter_log = FILTER_HIGH_CONTAMINATION_CLUSTERS.out.log
    filtered_clusters = FILTER_HIGH_CONTAMINATION_CLUSTERS.out.filtered_clusters
    contigs_prepared = PREPARE_SAG_REFERENCES.out.contigs_dir.ifEmpty([])
    indices_built = BUILD_BOWTIE2_INDICES.out.indices.ifEmpty([])
    alignment_results = RUN_CROSS_ALIGNMENTS.out.alignment_results.ifEmpty([])
}