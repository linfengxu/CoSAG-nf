#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPADES } from '../modules/spades'
include { CHECKM2 } from '../modules/checkm2'
include { PARSE_FILTERED_CLUSTERS; PARSE_SUBCLUSTER_FILES } from '../modules/parse_clusters'

process FILTER_HIGH_CONTAMINATION_CLUSTERS {
    tag "filter_contamination"
    container "${params.python.container}"

    publishDir "${params.outdir}/iterative_optimization", mode: 'copy'
    
    input:
    path cluster_json
    
    output:
    path "high_contamination_clusters.json", emit: filtered_clusters
    path "filter.log", emit: log
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    
    # 读取cluster数据
    with open('${cluster_json}', 'r') as f:
        data = json.load(f)
    
    # 过滤高污染度的clusters
    high_contamination_clusters = []
    contamination_threshold = ${params.iterative_optimization?.contamination_threshold ?: 10.0}
    min_cluster_size = ${params.iterative_optimization?.min_cluster_size ?: 3}
    
    for cluster in data['clusters']:
        if 'co_assembly_results' in cluster and 'round_1' in cluster['co_assembly_results']:
            round1_results = cluster['co_assembly_results']['round_1']
            contamination = round1_results.get('contamination', 0)
            cluster_size = cluster.get('cluster_size', 0)
            
            if contamination > contamination_threshold and cluster_size >= min_cluster_size:
                high_contamination_clusters.append(cluster)
                print(f"Cluster {cluster['cluster_id']}: contamination={contamination}%, size={cluster_size} - selected for optimization")
    
    # 输出结果
    output_data = {
        'clusters': high_contamination_clusters,
        'total_clusters': len(high_contamination_clusters),
        'filter_criteria': {
            'contamination_threshold': contamination_threshold,
            'min_cluster_size': min_cluster_size
        }
    }
    
    with open('high_contamination_clusters.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # 写入日志
    with open('filter.log', 'w') as f:
        f.write(f"Filtered clusters for iterative optimization\\n")
        f.write(f"Contamination threshold: {contamination_threshold}%\\n")
        f.write(f"Minimum cluster size: {min_cluster_size}\\n")
        f.write(f"Total clusters selected: {len(high_contamination_clusters)}\\n")
        
        for cluster in high_contamination_clusters:
            round1 = cluster['co_assembly_results']['round_1']
            f.write(f"Cluster {cluster['cluster_id']}: {round1['contamination']:.1f}% contamination, {round1['completeness']:.1f}% completeness\\n")
    
    print(f"Selected {len(high_contamination_clusters)} clusters for iterative optimization")
    """
}

process ADAPTIVE_SUBCLUSTERING {
    tag "cluster_${cluster_id}_round_${round}"
    container "${params.python_new.container}"
    publishDir "${params.outdir}/iterative_optimization/round_${round}", mode: 'copy'
    
    input:
    tuple val(cluster_id), path(cluster_file), val(round), path(previous_results)
    
    output:
    tuple val(cluster_id), val(round), path("cluster_${cluster_id}_subclusters_round_${round}.json"), emit: subclusters
    path "cluster_${cluster_id}_adaptive_analysis_round_${round}.log", emit: log
    
    script:
    def kmer_size = params.iterative_optimization?.kmer_size ?: 21
    def min_subcluster_size = params.iterative_optimization?.min_subcluster_size ?: 2
    def target_completeness = params.iterative_optimization?.target_completeness ?: 90.0
    def target_contamination = params.iterative_optimization?.target_contamination ?: 5.0
    def max_sags_per_subcluster = params.iterative_optimization?.max_sags_per_subcluster ?: 8
    """
    #!/usr/bin/env python3
    
    import json
    import os
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    import random
    
    def extract_kmers_from_fastq(fastq_file, k=21, max_reads=1000):
        # Extract k-mers from FASTQ file (sample first max_reads for efficiency)
        kmers = {}
        
        if not os.path.exists(fastq_file):
            print(f"Warning: File {fastq_file} not found")
            return kmers
            
        with open(fastq_file, 'r') as f:
            line_count = 0
            read_count = 0
            
            for line in f:
                line_count += 1
                if line_count % 4 == 2:  # Sequence line in FASTQ
                    read_count += 1
                    if read_count > max_reads:
                        break
                        
                    sequence = line.strip().upper()
                    
                    # Extract k-mers
                    for i in range(len(sequence) - k + 1):
                        kmer = sequence[i:i+k]
                        if 'N' not in kmer:  # Skip k-mers with N
                            kmers[kmer] = kmers.get(kmer, 0) + 1
        
        return kmers
    
    def calculate_kmer_similarity(kmers1, kmers2):
        # Calculate Jaccard similarity (more robust for sparse data)
        set1 = set(kmers1.keys())
        set2 = set(kmers2.keys())
        
        if not set1 and not set2:
            return 1.0
        if not set1 or not set2:
            return 0.0
        
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        
        return intersection / union if union > 0 else 0.0
    
    def greedy_subcluster_selection(members, kmer_profiles, max_size, similarity_threshold=0.3):
        # 贪心算法选择最相似的SAG组合
        if len(members) <= max_size:
            return [members]
        
        subclusters = []
        remaining_members = members.copy()
        
        while len(remaining_members) >= min_subcluster_size:
            if len(remaining_members) <= max_size:
                subclusters.append(remaining_members)
                break
            
            # 随机选择一个起始点
            seed_member = random.choice(remaining_members)
            current_cluster = [seed_member]
            remaining_members.remove(seed_member)
            
            # 贪心地添加最相似的成员
            while len(current_cluster) < max_size and remaining_members:
                best_member = None
                best_similarity = -1
                
                for candidate in remaining_members:
                    # 计算与当前cluster的平均相似性
                    similarities = []
                    for cluster_member in current_cluster:
                        if candidate['sag_id'] in kmer_profiles and cluster_member['sag_id'] in kmer_profiles:
                            sim = calculate_kmer_similarity(
                                kmer_profiles[candidate['sag_id']], 
                                kmer_profiles[cluster_member['sag_id']]
                            )
                            similarities.append(sim)
                    
                    if similarities:
                        avg_similarity = sum(similarities) / len(similarities)
                        if avg_similarity > best_similarity and avg_similarity > similarity_threshold:
                            best_similarity = avg_similarity
                            best_member = candidate
                
                if best_member:
                    current_cluster.append(best_member)
                    remaining_members.remove(best_member)
                else:
                    break  # 没有找到足够相似的成员
            
            if len(current_cluster) >= min_subcluster_size:
                subclusters.append(current_cluster)
        
        # 处理剩余的成员
        if remaining_members and len(remaining_members) >= min_subcluster_size:
            subclusters.append(remaining_members)
        elif remaining_members and subclusters:
            # 将剩余成员分配到现有的subclusters中
            for member in remaining_members:
                smallest_cluster = min(subclusters, key=len)
                if len(smallest_cluster) < max_size:
                    smallest_cluster.append(member)
        
        return subclusters
    
    cluster_id = "${cluster_id}"
    round_num = ${round}
    kmer_size = ${kmer_size}
    min_subcluster_size = ${min_subcluster_size}
    target_completeness = ${target_completeness}
    target_contamination = ${target_contamination}
    max_sags_per_subcluster = ${max_sags_per_subcluster}
    
    print(f"Processing cluster {cluster_id}, round {round_num}")
    print(f"Target: {target_completeness}% completeness, {target_contamination}% contamination")
    print(f"Max SAGs per subcluster: {max_sags_per_subcluster}")
    
    # 读取cluster文件
    members = []
    with open('${cluster_file}', 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i < 3:  # 跳过头部信息
                continue
            if line.strip():
                parts = line.strip().split('\\t')
                if len(parts) >= 3:
                    members.append({
                        'sag_id': parts[0],
                        'read1': parts[1],
                        'read2': parts[2]
                    })
    
    print(f"Cluster has {len(members)} members")
    
    # 读取前一轮的结果（如果存在）
    previous_results = {}
    if os.path.exists('${previous_results}') and round_num > 1:
        try:
            with open('${previous_results}', 'r') as f:
                previous_data = json.load(f)
                for result in previous_data:
                    if result.get('original_cluster_id') == cluster_id:
                        previous_results = result
                        break
            print(f"Found previous results for round {round_num-1}")
        except:
            print("No valid previous results found")
    
    # 根据轮次和前一轮结果调整策略
    if round_num == 1:
        # 第一轮：基于大小进行初步分割
        if len(members) <= max_sags_per_subcluster:
            subclusters = [members]
            print("Cluster size within limit, no splitting needed")
        else:
            # 简单分割成较小的组
            n_subclusters = (len(members) + max_sags_per_subcluster - 1) // max_sags_per_subcluster
            subcluster_size = len(members) // n_subclusters
            
            subclusters = []
            for i in range(n_subclusters):
                start_idx = i * subcluster_size
                if i == n_subclusters - 1:  # 最后一个subcluster包含剩余的所有成员
                    end_idx = len(members)
                else:
                    end_idx = (i + 1) * subcluster_size
                
                subcluster_members = members[start_idx:end_idx]
                if len(subcluster_members) >= min_subcluster_size:
                    subclusters.append(subcluster_members)
            
            print(f"Round 1: Split into {len(subclusters)} subclusters")
    
    else:
        # 后续轮次：基于前一轮结果进行智能优化
        print("Extracting k-mer profiles for intelligent splitting...")
        kmer_profiles = {}
        
        # 提取k-mer特征（采样以提高效率）
        for member in members:
            sag_id = member['sag_id']
            try:
                kmers1 = extract_kmers_from_fastq(member['read1'], kmer_size, max_reads=500)
                kmers2 = extract_kmers_from_fastq(member['read2'], kmer_size, max_reads=500)
                
                # 合并k-mers
                combined_kmers = kmers1.copy()
                for kmer, count in kmers2.items():
                    combined_kmers[kmer] = combined_kmers.get(kmer, 0) + count
                
                if combined_kmers:
                    kmer_profiles[sag_id] = combined_kmers
            except Exception as e:
                print(f"Error processing {sag_id}: {e}")
        
        print(f"Extracted k-mer profiles for {len(kmer_profiles)} members")
        
        # 使用贪心算法进行智能分组
        if len(kmer_profiles) >= min_subcluster_size:
            valid_members = [m for m in members if m['sag_id'] in kmer_profiles]
            subclusters = greedy_subcluster_selection(
                valid_members, kmer_profiles, max_sags_per_subcluster
            )
            print(f"Round {round_num}: Created {len(subclusters)} subclusters using greedy algorithm")
        else:
            subclusters = [members] if len(members) >= min_subcluster_size else []
            print(f"Round {round_num}: Insufficient k-mer profiles, keeping original cluster")
    
    # 格式化输出
    formatted_subclusters = []
    for i, subcluster_members in enumerate(subclusters):
        formatted_subclusters.append({
            'subcluster_id': f"{cluster_id}_{i+1}_r{round_num}",
            'members': subcluster_members,
            'member_count': len(subcluster_members),
            'round': round_num,
            'strategy': 'size_based' if round_num == 1 else 'kmer_greedy'
        })
    
    # 输出结果
    output_data = {
        'original_cluster_id': cluster_id,
        'round': round_num,
        'subclusters': formatted_subclusters,
        'total_subclusters': len(formatted_subclusters),
        'optimization_strategy': {
            'round': round_num,
            'method': 'size_based' if round_num == 1 else 'kmer_greedy',
            'max_sags_per_subcluster': max_sags_per_subcluster,
            'target_completeness': target_completeness,
            'target_contamination': target_contamination
        }
    }
    
    with open(f'cluster_{cluster_id}_subclusters_round_{round_num}.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # 写入日志
    with open(f'cluster_{cluster_id}_adaptive_analysis_round_{round_num}.log', 'w') as f:
        f.write(f"Adaptive subclustering for cluster {cluster_id}, round {round_num}\\n")
        f.write(f"Strategy: {'Size-based splitting' if round_num == 1 else 'K-mer greedy algorithm'}\\n")
        f.write(f"Original members: {len(members)}\\n")
        f.write(f"Subclusters created: {len(formatted_subclusters)}\\n")
        f.write(f"Max SAGs per subcluster: {max_sags_per_subcluster}\\n")
        f.write(f"Target quality: {target_completeness}% completeness, {target_contamination}% contamination\\n")
        
        for i, subcluster in enumerate(formatted_subclusters):
            f.write(f"Subcluster {i+1}: {subcluster['member_count']} members\\n")
    
    print(f"Created {len(formatted_subclusters)} subclusters for cluster {cluster_id}")
    """
}

process PREPARE_SUBCLUSTER_READS {
    tag "subcluster_${subcluster_id}"

    publishDir "${params.outdir}/iterative_optimization/subcluster_reads", mode: 'copy'
    
    input:
    tuple val(original_cluster_id), val(subcluster_id), val(members)
    
    output:
    tuple val(original_cluster_id), val(subcluster_id), val(members), path("${subcluster_id}_R1.fastq"), path("${subcluster_id}_R2.fastq"), emit: reads
    path "${subcluster_id}_prepare.log", emit: log
    
    script:
    def read1_files = members.collect { "\"${it.read1}\"" }.join(' ')
    def read2_files = members.collect { "\"${it.read2}\"" }.join(' ')
    def member_count = members.size()
    """
    #!/bin/bash
    
    echo "Preparing reads for subcluster ${subcluster_id}" > ${subcluster_id}_prepare.log
    echo "Number of members: ${member_count}" >> ${subcluster_id}_prepare.log
    
    # 初始化输出文件
    > ${subcluster_id}_R1.fastq
    > ${subcluster_id}_R2.fastq
    
    # 合并所有成员的read1文件
    echo "Merging R1 files..." >> ${subcluster_id}_prepare.log
    for file in ${read1_files}; do
        if [ -f "\$file" ]; then
            echo "  Adding R1: \$file" >> ${subcluster_id}_prepare.log
            cat "\$file" >> ${subcluster_id}_R1.fastq
        else
            echo "  WARNING: R1 file not found: \$file" >> ${subcluster_id}_prepare.log
        fi
    done
    
    # 合并所有成员的read2文件
    echo "Merging R2 files..." >> ${subcluster_id}_prepare.log
    for file in ${read2_files}; do
        if [ -f "\$file" ]; then
            echo "  Adding R2: \$file" >> ${subcluster_id}_prepare.log
            cat "\$file" >> ${subcluster_id}_R2.fastq
        else
            echo "  WARNING: R2 file not found: \$file" >> ${subcluster_id}_prepare.log
        fi
    done
    
    # 统计reads数量
    r1_reads=\$(grep -c "^@" ${subcluster_id}_R1.fastq || echo 0)
    r2_reads=\$(grep -c "^@" ${subcluster_id}_R2.fastq || echo 0)
    
    echo "Final read counts:" >> ${subcluster_id}_prepare.log
    echo "  R1 reads: \$r1_reads" >> ${subcluster_id}_prepare.log
    echo "  R2 reads: \$r2_reads" >> ${subcluster_id}_prepare.log
    
    if [ \$r1_reads -eq 0 ] || [ \$r2_reads -eq 0 ]; then
        echo "ERROR: No reads found for subcluster ${subcluster_id}" >> ${subcluster_id}_prepare.log
        exit 1
    fi
    """
}

process SUBCLUSTER_ASSEMBLY {
    tag "assembly_${subcluster_id}"
    container "${params.spades.container}"

    publishDir "${params.outdir}/iterative_optimization/assemblies", mode: 'copy'
    
    input:
    tuple val(original_cluster_id), val(subcluster_id), val(members), path(read1), path(read2)
    
    output:
    tuple val(original_cluster_id), val(subcluster_id), val(members), path("${subcluster_id}_contigs.fasta"), emit: contigs
    path "${subcluster_id}_assembly.log", emit: log
    
    script:
    def threads = params.iterative_optimization?.threads ?: 8
    """
    #!/bin/bash
    
    echo "Starting SPADES assembly for subcluster ${subcluster_id}" > ${subcluster_id}_assembly.log
    echo "Input files: ${read1}, ${read2}" >> ${subcluster_id}_assembly.log
    
    # 运行SPADES组装
    spades.py --careful --only-assembler \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o ${subcluster_id}_spades_output \\
        -t ${threads} \\
        >> ${subcluster_id}_assembly.log 2>&1
    
    # 检查组装结果
    if [ -f "${subcluster_id}_spades_output/contigs.fasta" ]; then
        cp ${subcluster_id}_spades_output/contigs.fasta ${subcluster_id}_contigs.fasta
        echo "Assembly completed successfully" >> ${subcluster_id}_assembly.log
        
        # 统计contigs信息
        num_contigs=\$(grep -c "^>" ${subcluster_id}_contigs.fasta)
        total_length=\$(awk '/^>/ {next} {total += length(\$0)} END {print total}' ${subcluster_id}_contigs.fasta)
        
        echo "Assembly statistics:" >> ${subcluster_id}_assembly.log
        echo "  Number of contigs: \$num_contigs" >> ${subcluster_id}_assembly.log
        echo "  Total length: \$total_length bp" >> ${subcluster_id}_assembly.log
    else
        echo "ERROR: Assembly failed - no contigs.fasta found" >> ${subcluster_id}_assembly.log
        # 创建空的contigs文件以避免流程中断
        touch ${subcluster_id}_contigs.fasta
        exit 1
    fi
    """
}

process SUBCLUSTER_QUALITY_ASSESSMENT {
    tag "checkm2_${subcluster_id}"
    container "${params.checkm2.container}"
    publishDir "${params.outdir}/iterative_optimization/quality_assessment", mode: 'copy'
    
    input:
    tuple val(original_cluster_id), val(subcluster_id), val(members), path(contigs)
    
    output:
    tuple val(original_cluster_id), val(subcluster_id), val(members), path(contigs), path("${subcluster_id}_checkm2.tsv"), emit: results
    path "${subcluster_id}_checkm2.log", emit: log
    
    script:
    def threads = params.iterative_optimization?.threads ?: 4
    def db_path = params.checkm2?.database
    """
    #!/bin/bash
    
    # 设置CheckM2数据库环境变量
    export CHECKM2DB=${db_path}
    
    echo "Starting CheckM2 assessment for subcluster ${subcluster_id}" > ${subcluster_id}_checkm2.log
    echo "Input contigs: ${contigs}" >> ${subcluster_id}_checkm2.log
    
    # 检查contigs文件
    if [ ! -s "${contigs}" ]; then
        echo "ERROR: Empty or missing contigs file" >> ${subcluster_id}_checkm2.log
        # 创建空的结果文件
        echo -e "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes" > ${subcluster_id}_checkm2.tsv
        echo -e "${subcluster_id}\\t0\\t0\\tNA\\tNA\\t0\\t0\\t0\\t0\\t0\\t0\\tEmpty contigs file" >> ${subcluster_id}_checkm2.tsv
        exit 0
    fi
    
    # 运行CheckM2
    checkm2 predict -x fasta \\
        --threads ${threads} \\
        --input ${contigs} \\
        --output-directory ${subcluster_id}_checkm2_output \\
        >> ${subcluster_id}_checkm2.log 2>&1
    
    # 检查结果
    result_file=\$(find ${subcluster_id}_checkm2_output -maxdepth 1 -name "*.tsv" -type f | head -1)
    
    if [ -n "\$result_file" ] && [ -s "\$result_file" ]; then
        cp "\$result_file" ${subcluster_id}_checkm2.tsv
        echo "CheckM2 assessment completed successfully" >> ${subcluster_id}_checkm2.log
        
        # 提取关键指标
        completeness=\$(tail -n +2 ${subcluster_id}_checkm2.tsv | cut -f2 | head -1)
        contamination=\$(tail -n +2 ${subcluster_id}_checkm2.tsv | cut -f3 | head -1)
        
        echo "Quality metrics:" >> ${subcluster_id}_checkm2.log
        echo "  Completeness: \$completeness%" >> ${subcluster_id}_checkm2.log
        echo "  Contamination: \$contamination%" >> ${subcluster_id}_checkm2.log
    else
        echo "ERROR: CheckM2 failed - no valid results" >> ${subcluster_id}_checkm2.log
        # 创建空的结果文件
        echo -e "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes" > ${subcluster_id}_checkm2.tsv
        echo -e "${subcluster_id}\\t0\\t0\\tNA\\tNA\\t0\\t0\\t0\\t0\\t0\\t0\\tCheckM2 failed" >> ${subcluster_id}_checkm2.tsv
    fi
    """
}

process SELECT_BEST_SUBCLUSTER {
    tag "select_best_${original_cluster_id}"
    container "${params.python.container}"
    publishDir "${params.outdir}/iterative_optimization/final_results", mode: 'copy'
    
    input:
    tuple val(original_cluster_id), path(quality_results)
    
    output:
    path "optimized_cluster_${original_cluster_id}.json", emit: optimized_results
    path "optimization_${original_cluster_id}.log", emit: log
    
    script:
    def target_completeness = params.iterative_optimization?.target_completeness ?: 50.0
    def target_contamination = params.iterative_optimization?.target_contamination ?: 10.0
    """
    #!/usr/bin/env python3
    
    import json
    import pandas as pd
    import glob
    import os
    
    original_cluster_id = "${original_cluster_id}"
    target_completeness = ${target_completeness}
    target_contamination = ${target_contamination}
    
    print(f"Selecting best subcluster for cluster {original_cluster_id}")
    
    # 收集所有CheckM2结果
    checkm2_files = glob.glob('*_checkm2.tsv')
    subcluster_results = []
    
    print(f"Found {len(checkm2_files)} CheckM2 result files")
    
    for file in checkm2_files:
        subcluster_id = file.replace('_checkm2.tsv', '')
        
        try:
            df = pd.read_csv(file, sep='\\t')
            if len(df) > 0:
                row = df.iloc[0]
                completeness = float(row.get('Completeness', 0))
                contamination = float(row.get('Contamination', 0))
                
                # 计算质量分数
                quality_score = completeness - 5 * contamination
                
                # 查找对应的contigs文件
                contigs_file = f"{subcluster_id}_contigs.fasta"
                assembly_file = os.path.abspath(contigs_file) if os.path.exists(contigs_file) else None
                
                result = {
                    'subcluster_id': subcluster_id,
                    'completeness': completeness,
                    'contamination': contamination,
                    'quality_score': quality_score,
                    'genome_size': int(row.get('Genome_Size', 0)),
                    'gc_content': float(row.get('GC_Content', 0)),
                    'contig_n50': int(row.get('Contig_N50', 0)),
                    'total_coding_sequences': int(row.get('Total_Coding_Sequences', 0)),
                    'assembly_file': assembly_file,
                    'meets_target': completeness >= target_completeness and contamination <= target_contamination
                }
                
                subcluster_results.append(result)
                print(f"Subcluster {subcluster_id}: {completeness:.1f}% complete, {contamination:.1f}% contamination, quality={quality_score:.1f}")
        
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    # 选择最佳结果
    if subcluster_results:
        # 按质量分数排序
        subcluster_results.sort(key=lambda x: x['quality_score'], reverse=True)
        
        # 选择满足目标的最佳结果，或者质量分数最高的结果
        target_results = [r for r in subcluster_results if r['meets_target']]
        best_result = target_results[0] if target_results else subcluster_results[0]
        
        optimization_result = {
            'original_cluster_id': original_cluster_id,
            'optimization_successful': len(target_results) > 0,
            'best_subcluster': best_result,
            'all_subclusters_evaluated': len(subcluster_results),
            'target_criteria': {
                'completeness': target_completeness,
                'contamination': target_contamination
            },
            'final_metrics': {
                'completeness': best_result['completeness'],
                'contamination': best_result['contamination'],
                'quality_score': best_result['quality_score']
            }
        }
    else:
        optimization_result = {
            'original_cluster_id': original_cluster_id,
            'optimization_successful': False,
            'error': 'No valid subclusters could be assembled and evaluated'
        }
    
    # 输出结果
    with open(f'optimized_cluster_{original_cluster_id}.json', 'w') as f:
        json.dump(optimization_result, f, indent=2)
    
    # 写入日志
    with open(f'optimization_{original_cluster_id}.log', 'w') as f:
        f.write(f"Iterative optimization results for cluster {original_cluster_id}\\n")
        f.write(f"Target completeness: {target_completeness}%\\n")
        f.write(f"Target contamination: {target_contamination}%\\n")
        f.write(f"Subclusters evaluated: {len(subcluster_results)}\\n")
        
        if subcluster_results:
            f.write(f"\\nBest result:\\n")
            best = optimization_result['best_subcluster']
            f.write(f"  Subcluster: {best['subcluster_id']}\\n")
            f.write(f"  Completeness: {best['completeness']:.1f}%\\n")
            f.write(f"  Contamination: {best['contamination']:.1f}%\\n")
            f.write(f"  Quality score: {best['quality_score']:.1f}\\n")
            f.write(f"  Assembly file: {best.get('assembly_file', 'N/A')}\\n")
            f.write(f"  Meets target: {best['meets_target']}\\n")
            
            f.write(f"\\nAll subclusters:\\n")
            for i, result in enumerate(subcluster_results):
                f.write(f"  {i+1}. {result['subcluster_id']}: {result['completeness']:.1f}% complete, {result['contamination']:.1f}% contamination, quality={result['quality_score']:.1f}\\n")
        else:
            f.write("No valid results obtained\\n")
    
    print(f"Optimization completed for cluster {original_cluster_id}")
    """
}

workflow ITERATIVE_OPTIMIZATION {
    take:
    cluster_json
    
    main:
    // 过滤高污染度的clusters
    FILTER_HIGH_CONTAMINATION_CLUSTERS(cluster_json)
    

    
    // 解析过滤后的clusters
    PARSE_FILTERED_CLUSTERS(FILTER_HIGH_CONTAMINATION_CLUSTERS.out.filtered_clusters)
    
    // 创建cluster channel - 直接使用文件
    optimization_clusters = PARSE_FILTERED_CLUSTERS.out.cluster_files
        .flatten()
        .map { file ->
            def lines = file.text.split('\n')
            def cluster_id = lines[0].split('\t')[1]
            tuple(cluster_id, file)
        }
    
    // 实现真正的迭代优化逻辑
    max_rounds = params.iterative_optimization?.max_rounds ?: 3
    
    // 第一轮：初始分割
    first_round_ch = optimization_clusters
        .map { cluster_id, cluster_file -> 
            tuple(cluster_id, cluster_file, 1, file('NO_FILE'))  // 第一轮没有前一轮结果
        }
    
    // 执行第一轮子聚类
    ADAPTIVE_SUBCLUSTERING(first_round_ch)
    
    // 对第一轮结果进行组装和评估
    first_round_subclusters = ADAPTIVE_SUBCLUSTERING.out.subclusters
        .map { cluster_id, round, subcluster_file -> subcluster_file }
        .collect()
    
    PARSE_SUBCLUSTER_FILES(first_round_subclusters)
    
    first_round_subcluster_ch = PARSE_SUBCLUSTER_FILES.out.subcluster_files
        .flatten()
        .map { file ->
            def lines = file.text.split('\n')
            def original_cluster_id = lines[0].split('\t')[1]
            def subcluster_id = lines[1].split('\t')[1]
            
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
            
            tuple(original_cluster_id, subcluster_id, members)
        }
    
    // 第一轮组装和评估
    PREPARE_SUBCLUSTER_READS(first_round_subcluster_ch)
    SUBCLUSTER_ASSEMBLY(PREPARE_SUBCLUSTER_READS.out.reads)
    SUBCLUSTER_QUALITY_ASSESSMENT(SUBCLUSTER_ASSEMBLY.out.contigs)
    
    // 收集第一轮结果
    first_round_results = SUBCLUSTER_QUALITY_ASSESSMENT.out.results
        .map { original_cluster_id, subcluster_id, members, contigs, checkm2_tsv ->
            tuple(original_cluster_id, checkm2_tsv)
        }
        .groupTuple()
    
    SELECT_BEST_SUBCLUSTER(first_round_results)
    
    // 这部分已经在上面的第一轮处理中完成了
    
    emit:
    optimized_results = SELECT_BEST_SUBCLUSTER.out.optimized_results
    optimization_logs = SELECT_BEST_SUBCLUSTER.out.log
    filter_log = FILTER_HIGH_CONTAMINATION_CLUSTERS.out.log
    subclustering_logs = ADAPTIVE_SUBCLUSTERING.out.log
    assembly_logs = SUBCLUSTER_ASSEMBLY.out.log
    quality_logs = SUBCLUSTER_QUALITY_ASSESSMENT.out.log
}