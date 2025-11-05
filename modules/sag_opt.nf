process COSAG_OPTIMIZATION {
    tag "cosag_opt_${json_file.baseName}"
    label 'process_high'
    publishDir "${params.outdir}/sag_optimization2", mode: 'copy'
    
    input:
    path json_file
    
    output:
    path "${json_file.baseName}_output/*", emit: results
    path "${json_file.baseName}_sag_opt.log", emit: log
    path "${json_file.baseName}_output", emit: output_dir, type: 'dir'
    path "${json_file.baseName}_output/*optimized.json", emit: optimized_json, optional: true
    path "${json_file.baseName}_output/all_contig/*.fasta", emit: bset_cosag, optional: true
    
    script:
    """
    # 为每个JSON文件创建独立的输出目录
    mkdir -p ${json_file.baseName}_output
    
    # 运行 SAG 优化脚本
    python ${projectDir}/bin/sag_decontam3.py \\
        --json ${json_file} \\
        --outdir ${json_file.baseName}_output \\
        > ${json_file.baseName}_sag_opt.log 2>&1
    """
}

process COLLECT_ALL_ASSEMBLIES {
    tag "collect_all_assemblies"
    label 'process_low'
    publishDir "${params.outdir}/merged_assemblies", mode: 'copy'
    
    input:
    path bset_fastas
    path json_files  // 改为独立的 path
    path output_dirs  // 改为独立的 path
    
    output:
    path "bset_cosag/*", emit: bset_fastas, optional: true
    path "co_assembly/*", emit: co_assembly_fastas, optional: true
    path "all_accepted_cosgas/*", emit: all_bins  // 新增：所有合格的 bins
    path "filtered_out/*", emit: filtered_out, optional: true
    path "collection_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    import shutil
    from pathlib import Path
    
    # 创建输出目录
    os.makedirs("bset_cosag", exist_ok=True)
    os.makedirs("co_assembly", exist_ok=True)
    os.makedirs("filtered_out", exist_ok=True)
    os.makedirs("all_accepted_cosgas", exist_ok=True)
    
    bset_count = 0
    co_assembly_count = 0
    filtered_count = 0
    
    # 质量过滤标准
    MIN_COMPLETENESS = 50.0
    MAX_CONTAMINATION = 10.0
    
    # 收集 bset_cosag fasta 文件
    bset_files = "${bset_fastas}".split()
    for fasta_file in bset_files:
        if fasta_file and os.path.exists(fasta_file):
            filename = os.path.basename(fasta_file)
            # 复制到 bset_cosag 目录
            shutil.copy(fasta_file, f"bset_cosag/{filename}")
            # 同时复制到 all_bins 目录（用于 dRep）
            shutil.copy(fasta_file, f"all_accepted_cosgas/{filename}")
            bset_count += 1
    
    # 收集 co_assembly fasta 文件（带质量过滤）
    json_files_list = "${json_files}".split()
    filtered_info = []
    passed_info = []
    
    for json_file in json_files_list:
        if not os.path.exists(json_file):
            continue
            
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
        except Exception as e:
            print(f"Error reading {json_file}: {e}")
            continue
        
        # 遍历所有 clusters
        for cluster in data.get("clusters", []):
            cluster_id = cluster.get("cluster_id", "unknown")
            cluster_size = cluster.get("cluster_size", 0)
            co_assembly = cluster.get("co_assembly_results", {})
            
            # 检查 round_1
            round_1 = co_assembly.get("round_1", {})
            assembly_file = round_1.get("assembly_file", "")
            completeness = round_1.get("completeness", 0)
            contamination = round_1.get("contamination", 100)
            genome_size = round_1.get("genome_size", 0)
            contig_n50 = round_1.get("contig_n50", 0)
            
            if not assembly_file or not os.path.exists(assembly_file):
                continue
            
            # 生成文件名
            base_name = os.path.basename(json_file).replace('.json', '')
            new_filename = f"cluster_{cluster_id}_CoSAG.fasta"
            
            # 质量过滤
            if completeness >= MIN_COMPLETENESS and contamination <= MAX_CONTAMINATION:
                # 通过质量标准
                shutil.copy(assembly_file, f"co_assembly/{new_filename}")
                # 同时复制到 all_bins 目录（用于 dRep）
                shutil.copy(assembly_file, f"all_accepted_cosgas/{new_filename}")
                co_assembly_count += 1
                passed_info.append({
                    "file": new_filename,
                    "cluster_id": cluster_id,
                    "cluster_size": cluster_size,
                    "completeness": completeness,
                    "contamination": contamination,
                    "genome_size": genome_size,
                    "contig_n50": contig_n50
                })
            else:
                # 未通过质量标准
                shutil.copy(assembly_file, f"filtered_out/{new_filename}")
                filtered_count += 1
                reasons = []
                
                # 记录未通过的原因
                if completeness < MIN_COMPLETENESS:
                    reasons.append(f"completeness {completeness:.2f} < {MIN_COMPLETENESS}")
                if contamination > MAX_CONTAMINATION:
                    reasons.append(f"contamination {contamination:.2f} > {MAX_CONTAMINATION}")
                
                filtered_info.append({
                    "file": new_filename,
                    "cluster_id": cluster_id,
                    "cluster_size": cluster_size,
                    "completeness": completeness,
                    "contamination": contamination,
                    "genome_size": genome_size,
                    "contig_n50": contig_n50,
                    "reason": reasons
                })
    
    total_bins = bset_count + co_assembly_count
    
    # 生成统计报告
    with open("collection_summary.txt", "w") as f:
        f.write("="*80 + "\\n")
        f.write("Assembly Collection Summary\\n")
        f.write("="*80 + "\\n\\n")
        
        f.write("Quality Filter Criteria:\\n")
        f.write(f"  - Completeness >= {MIN_COMPLETENESS}%\\n")
        f.write(f"  - Contamination <= {MAX_CONTAMINATION}%\\n")
        f.write("\\n")
        
        f.write("Collection Results:\\n")
        f.write(f"  - Total bset_cosag files collected: {bset_count}\\n")
        f.write(f"  - Total co_assembly files collected (passed QC): {co_assembly_count}\\n")
        f.write(f"  - Total bins for dRep: {total_bins}\\n")
        f.write(f"  - Total co_assembly files filtered out: {filtered_count}\\n")
        f.write("\\n")
        
        # Bset_cosag 文件列表
        if bset_count > 0:
            f.write("-"*80 + "\\n")
            f.write(f"Bset_cosag files ({bset_count} files):\\n")
            f.write("-"*80 + "\\n")
            for fasta in sorted(os.listdir("bset_cosag")):
                f.write(f"  - {fasta}\\n")
            f.write("\\n")
        
        # 通过质量标准的 co_assembly 文件
        if co_assembly_count > 0:
            f.write("-"*80 + "\\n")
            f.write(f"Co_assembly files - PASSED QC ({co_assembly_count} files):\\n")
            f.write("-"*80 + "\\n")
            for info in sorted(passed_info, key=lambda x: x['completeness'], reverse=True):
                f.write(f"  - {info['file']}\\n")
                f.write(f"    Cluster ID: {info['cluster_id']} | Cluster size: {info['cluster_size']}\\n")
                f.write(f"    Completeness: {info['completeness']:.2f}% | Contamination: {info['contamination']:.2f}%\\n")
                f.write(f"    Genome size: {info['genome_size']:,} bp | N50: {info['contig_n50']:,} bp\\n")
                f.write("\\n")
        
        # 未通过质量标准的 co_assembly 文件
        if filtered_count > 0:
            f.write("-"*80 + "\\n")
            f.write(f"Co_assembly files - FAILED QC ({filtered_count} files):\\n")
            f.write("-"*80 + "\\n")
            for info in sorted(filtered_info, key=lambda x: x['completeness'], reverse=True):
                f.write(f"  - {info['file']}\\n")
                f.write(f"    Cluster ID: {info['cluster_id']} | Cluster size: {info['cluster_size']}\\n")
                f.write(f"    Completeness: {info['completeness']:.2f}% | Contamination: {info['contamination']:.2f}%\\n")
                f.write(f"    Genome size: {info['genome_size']:,} bp | N50: {info['contig_n50']:,} bp\\n")
                f.write(f"    Reason: {'; '.join(info['reason'])}\\n")
                f.write("\\n")
        
        f.write("="*80 + "\\n")
        f.write("End of Report\\n")
        f.write("="*80 + "\\n")
    
    print(f"\\n{'='*60}")
    print(f"Collection Summary:")
    print(f"{'='*60}")
    print(f"Bset_cosag files: {bset_count}")
    print(f"Co_assembly files (passed QC): {co_assembly_count}")
    print(f"Total bins for dRep: {total_bins}")
    print(f"Co_assembly files (failed QC): {filtered_count}")
    print(f"{'='*60}\\n")
    """
}