#!/usr/bin/env python3
"""
生成Nextflow流程的汇总报告
整合所有分析结果，生成用户友好的汇总表格和报告
"""

import os
import sys
import json
import pandas as pd
import glob
from pathlib import Path
import argparse
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate summary reports for Nextflow pipeline')
    parser.add_argument('--outdir', required=True, help='Pipeline output directory')
    parser.add_argument('--output', default='summary_reports', help='Summary output directory')
    return parser.parse_args()

def collect_individual_assembly_stats(outdir):
    """收集个体组装统计信息"""
    assembly_stats = []
    
    # 查找所有个体组装结果
    assembly_pattern = f"{outdir}/01_individual_assemblies/per_sample/*/spades/contigs/*.fasta"
    assembly_files = glob.glob(assembly_pattern)
    
    for assembly_file in assembly_files:
        sample_id = Path(assembly_file).stem.replace('_contigs', '')
        
        # 基本组装统计
        stats = {
            'Sample_ID': sample_id,
            'Assembly_File': assembly_file,
            'Assembly_Status': 'Completed' if os.path.exists(assembly_file) else 'Failed'
        }
        
        # 读取组装序列统计
        if os.path.exists(assembly_file):
            contigs = []
            current_seq = ""
            
            with open(assembly_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        if current_seq:
                            contigs.append(len(current_seq))
                        current_seq = ""
                    else:
                        current_seq += line.strip()
                if current_seq:
                    contigs.append(len(current_seq))
            
            if contigs:
                stats.update({
                    'Num_Contigs': len(contigs),
                    'Total_Length': sum(contigs),
                    'Max_Contig_Length': max(contigs),
                    'Min_Contig_Length': min(contigs),
                    'N50': calculate_n50(contigs)
                })
            else:
                stats.update({
                    'Num_Contigs': 0,
                    'Total_Length': 0,
                    'Max_Contig_Length': 0,
                    'Min_Contig_Length': 0,
                    'N50': 0
                })
        
        # 读取CheckM2结果
        checkm2_file = f"{outdir}/01_individual_assemblies/per_sample/{sample_id}/checkm2/{sample_id}_checkm2.tsv"
        if os.path.exists(checkm2_file):
            try:
                checkm2_df = pd.read_csv(checkm2_file, sep='\t')
                if len(checkm2_df) > 0:
                    row = checkm2_df.iloc[0]
                    stats.update({
                        'Completeness': row.get('Completeness', 0),
                        'Contamination': row.get('Contamination', 0),
                        'CheckM2_Quality': 'High' if row.get('Completeness', 0) > 70 and row.get('Contamination', 0) < 10 else 'Medium' if row.get('Completeness', 0) > 50 else 'Low'
                    })
            except Exception as e:
                print(f"Warning: Could not read CheckM2 results for {sample_id}: {e}")
                stats.update({
                    'Completeness': 0,
                    'Contamination': 0,
                    'CheckM2_Quality': 'Unknown'
                })
        
        assembly_stats.append(stats)
    
    return pd.DataFrame(assembly_stats)

def calculate_n50(contig_lengths):
    """计算N50值"""
    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length / 2
    
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= target_length:
            return length
    return 0

def collect_clustering_stats(outdir):
    """收集聚类统计信息"""
    clustering_stats = {}
    
    # 读取聚类结果
    clusters_file = f"{outdir}/03_clustering_analysis/hierarchical_clustering/hierarchical_clusters.tsv"
    if os.path.exists(clusters_file):
        clusters_df = pd.read_csv(clusters_file, sep='\t')
        
        clustering_stats = {
            'Total_SAGs': len(clusters_df),
            'Total_Clusters': clusters_df['Cluster_ID'].nunique(),
            'Singleton_Clusters': sum(clusters_df.groupby('Cluster_ID').size() == 1),
            'Multi_SAG_Clusters': sum(clusters_df.groupby('Cluster_ID').size() > 1),
            'Largest_Cluster_Size': clusters_df.groupby('Cluster_ID').size().max(),
            'Mean_Cluster_Size': clusters_df.groupby('Cluster_ID').size().mean()
        }
        
        # 聚类大小分布
        cluster_sizes = clusters_df.groupby('Cluster_ID').size()
        size_distribution = cluster_sizes.value_counts().sort_index()
        clustering_stats['Cluster_Size_Distribution'] = size_distribution.to_dict()
    
    return clustering_stats

def collect_co_assembly_stats(outdir):
    """收集共组装统计信息"""
    co_assembly_stats = []
    
    # 查找共组装结果
    co_assembly_pattern = f"{outdir}/04_co_assemblies/assemblies/cluster_*/contigs/*.fasta"
    co_assembly_files = glob.glob(co_assembly_pattern)
    
    for assembly_file in co_assembly_files:
        cluster_id = Path(assembly_file).parent.parent.name.replace('cluster_', '')
        
        stats = {
            'Cluster_ID': cluster_id,
            'Co_Assembly_File': assembly_file,
            'Co_Assembly_Status': 'Completed' if os.path.exists(assembly_file) else 'Failed'
        }
        
        # 读取组装统计
        if os.path.exists(assembly_file):
            contigs = []
            current_seq = ""
            
            with open(assembly_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        if current_seq:
                            contigs.append(len(current_seq))
                        current_seq = ""
                    else:
                        current_seq += line.strip()
                if current_seq:
                    contigs.append(len(current_seq))
            
            if contigs:
                stats.update({
                    'Num_Contigs': len(contigs),
                    'Total_Length': sum(contigs),
                    'Max_Contig_Length': max(contigs),
                    'N50': calculate_n50(contigs)
                })
        
        # 读取CheckM2结果
        checkm2_file = f"{outdir}/04_co_assemblies/assemblies/cluster_{cluster_id}/checkm2/cluster_{cluster_id}_checkm2.tsv"
        if os.path.exists(checkm2_file):
            try:
                checkm2_df = pd.read_csv(checkm2_file, sep='\t')
                if len(checkm2_df) > 0:
                    row = checkm2_df.iloc[0]
                    stats.update({
                        'Co_Assembly_Completeness': row.get('Completeness', 0),
                        'Co_Assembly_Contamination': row.get('Contamination', 0)
                    })
            except Exception:
                stats.update({
                    'Co_Assembly_Completeness': 0,
                    'Co_Assembly_Contamination': 0
                })
        
        co_assembly_stats.append(stats)
    
    return pd.DataFrame(co_assembly_stats)

def collect_taxonomy_stats(outdir):
    """收集分类统计信息"""
    taxonomy_stats = {}
    
    # 个体SAG分类
    individual_gtdb_file = f"{outdir}/05_taxonomic_classification/individual_sags/gtdbtk.bac120.summary.tsv"
    if os.path.exists(individual_gtdb_file):
        try:
            gtdb_df = pd.read_csv(individual_gtdb_file, sep='\t')
            taxonomy_stats['Individual_SAGs_Classified'] = len(gtdb_df)
            
            # 统计分类层级
            if 'classification' in gtdb_df.columns:
                classifications = gtdb_df['classification'].dropna()
                phyla = [cls.split(';')[1].replace('p__', '') for cls in classifications if ';' in cls]
                taxonomy_stats['Unique_Phyla'] = len(set(phyla))
                taxonomy_stats['Top_Phyla'] = pd.Series(phyla).value_counts().head(5).to_dict()
        except Exception as e:
            print(f"Warning: Could not read individual GTDB results: {e}")
    
    # 共组装分类
    co_assembly_gtdb_file = f"{outdir}/05_taxonomic_classification/co_assemblies/gtdbtk.bac120.summary.tsv"
    if os.path.exists(co_assembly_gtdb_file):
        try:
            gtdb_df = pd.read_csv(co_assembly_gtdb_file, sep='\t')
            taxonomy_stats['Co_Assemblies_Classified'] = len(gtdb_df)
        except Exception as e:
            print(f"Warning: Could not read co-assembly GTDB results: {e}")
    
    return taxonomy_stats

def generate_summary_report(outdir, output_dir):
    """生成综合汇总报告"""
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("Collecting individual assembly statistics...")
    individual_stats = collect_individual_assembly_stats(outdir)
    
    print("Collecting clustering statistics...")
    clustering_stats = collect_clustering_stats(outdir)
    
    print("Collecting co-assembly statistics...")
    co_assembly_stats = collect_co_assembly_stats(outdir)
    
    print("Collecting taxonomy statistics...")
    taxonomy_stats = collect_taxonomy_stats(outdir)
    
    # 保存详细统计表格
    if not individual_stats.empty:
        individual_stats.to_csv(f"{output_dir}/individual_assemblies_summary.tsv", sep='\t', index=False)
        print(f"Saved individual assemblies summary: {len(individual_stats)} samples")
    
    if not co_assembly_stats.empty:
        co_assembly_stats.to_csv(f"{output_dir}/co_assemblies_summary.tsv", sep='\t', index=False)
        print(f"Saved co-assemblies summary: {len(co_assembly_stats)} clusters")
    
    # 生成总体汇总报告
    with open(f"{output_dir}/pipeline_summary_report.txt", 'w') as f:
        f.write("# Nextflow Single-Cell Genome Pipeline Summary Report\n")
        f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Overall Statistics\n")
        f.write(f"Total samples processed: {len(individual_stats) if not individual_stats.empty else 0}\n")
        f.write(f"Successful individual assemblies: {sum(individual_stats['Assembly_Status'] == 'Completed') if not individual_stats.empty else 0}\n")
        f.write(f"Total clusters formed: {clustering_stats.get('Total_Clusters', 0)}\n")
        f.write(f"Co-assemblies generated: {len(co_assembly_stats) if not co_assembly_stats.empty else 0}\n\n")
        
        if clustering_stats:
            f.write("## Clustering Analysis\n")
            f.write(f"Total SAGs clustered: {clustering_stats.get('Total_SAGs', 0)}\n")
            f.write(f"Singleton clusters: {clustering_stats.get('Singleton_Clusters', 0)}\n")
            f.write(f"Multi-SAG clusters: {clustering_stats.get('Multi_SAG_Clusters', 0)}\n")
            f.write(f"Largest cluster size: {clustering_stats.get('Largest_Cluster_Size', 0)}\n")
            f.write(f"Mean cluster size: {clustering_stats.get('Mean_Cluster_Size', 0):.2f}\n\n")
        
        if not individual_stats.empty:
            f.write("## Assembly Quality\n")
            high_quality = sum((individual_stats['Completeness'] > 70) & (individual_stats['Contamination'] < 10))
            medium_quality = sum((individual_stats['Completeness'] > 50) & (individual_stats['Completeness'] <= 70))
            f.write(f"High-quality assemblies (>70% complete, <10% contamination): {high_quality}\n")
            f.write(f"Medium-quality assemblies (50-70% complete): {medium_quality}\n")
            f.write(f"Mean completeness: {individual_stats['Completeness'].mean():.2f}%\n")
            f.write(f"Mean contamination: {individual_stats['Contamination'].mean():.2f}%\n\n")
        
        if taxonomy_stats:
            f.write("## Taxonomic Classification\n")
            f.write(f"Individual SAGs classified: {taxonomy_stats.get('Individual_SAGs_Classified', 0)}\n")
            f.write(f"Co-assemblies classified: {taxonomy_stats.get('Co_Assemblies_Classified', 0)}\n")
            f.write(f"Unique phyla identified: {taxonomy_stats.get('Unique_Phyla', 0)}\n\n")
            
            if 'Top_Phyla' in taxonomy_stats:
                f.write("Top 5 phyla:\n")
                for phylum, count in taxonomy_stats['Top_Phyla'].items():
                    f.write(f"  {phylum}: {count} SAGs\n")
        
        f.write("\n## Output Directory Structure\n")
        f.write("Results are organized in the following directories:\n")
        f.write("- 01_individual_assemblies/: Individual SAG assemblies and quality assessments\n")
        f.write("- 02_similarity_analysis/: MinHash similarity analysis results\n")
        f.write("- 03_clustering_analysis/: Hierarchical clustering results\n")
        f.write("- 04_co_assemblies/: Co-assembly results for each cluster\n")
        f.write("- 05_taxonomic_classification/: GTDB-Tk taxonomic classification\n")
        f.write("- 06_final_results/: Integrated final results and reports\n")
    
    print(f"Summary report generated: {output_dir}/pipeline_summary_report.txt")
    
    # 生成JSON格式的汇总数据
    summary_data = {
        'generation_time': datetime.now().isoformat(),
        'individual_assemblies': len(individual_stats) if not individual_stats.empty else 0,
        'clustering_stats': clustering_stats,
        'co_assemblies': len(co_assembly_stats) if not co_assembly_stats.empty else 0,
        'taxonomy_stats': taxonomy_stats
    }
    
    with open(f"{output_dir}/pipeline_summary.json", 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    print(f"JSON summary generated: {output_dir}/pipeline_summary.json")

def main():
    args = parse_arguments()
    
    if not os.path.exists(args.outdir):
        print(f"Error: Output directory {args.outdir} does not exist")
        sys.exit(1)
    
    output_dir = os.path.join(args.outdir, args.output)
    
    print(f"Generating summary reports for pipeline output: {args.outdir}")
    print(f"Summary reports will be saved to: {output_dir}")
    
    generate_summary_report(args.outdir, output_dir)
    
    print("Summary report generation completed!")

if __name__ == "__main__":
    main()