#!/usr/bin/env python3
"""
脚本功能：将GTDB-Tk分类信息添加到cluster_data_updated.json中
- 匹配sag_id和user_genome（去掉_contigs后缀）
- 添加分类信息到每个成员中

使用方法：
python merge_taxonomy.py <gtdbtk_file> <cluster_file> <output_file>
"""

import json
import pandas as pd
import sys
import argparse
from pathlib import Path

def load_gtdbtk_data(gtdbtk_file):
    """加载GTDB-Tk分类数据"""
    try:
        df = pd.read_csv(gtdbtk_file, sep='\t')
        print(f"成功加载GTDB-Tk数据，共{len(df)}条记录")
        
        # 创建映射字典：去掉_contigs后缀的genome名称 -> 分类信息
        taxonomy_dict = {}
        for _, row in df.iterrows():
            genome_name = row['user_genome']
            # 去掉_contigs后缀
            clean_name = genome_name.replace('_contigs', '') if genome_name.endswith('_contigs') else genome_name
            
            taxonomy_dict[clean_name] = {
                'classification': row['classification'],
                'classification_method': row.get('classification_method', 'N/A'),
                'closest_genome_reference': row.get('closest_genome_reference', 'N/A'),
                'closest_genome_ani': row.get('closest_genome_ani', 'N/A'),
                'msa_percent': row.get('msa_percent', 'N/A'),
                'red_value': row.get('red_value', 'N/A'),
                'warnings': row.get('warnings', 'N/A')
            }
        
        print(f"处理后的分类数据字典包含{len(taxonomy_dict)}条记录")
        return taxonomy_dict
        
    except Exception as e:
        print(f"加载GTDB-Tk数据时出错: {e}")
        return {}

def load_cluster_data(cluster_file):
    """加载cluster数据"""
    try:
        with open(cluster_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        print(f"成功加载cluster数据，共{len(data['clusters'])}个cluster")
        return data
    except Exception as e:
        print(f"加载cluster数据时出错: {e}")
        return None

def merge_taxonomy_info(cluster_data, taxonomy_dict):
    """将分类信息合并到cluster数据中"""
    matched_count = 0
    total_members = 0
    
    for cluster in cluster_data['clusters']:
        for member in cluster['members']:
            total_members += 1
            sag_id = member['sag_id']
            
            # 查找匹配的分类信息
            if sag_id in taxonomy_dict:
                member['taxonomy'] = taxonomy_dict[sag_id]
                matched_count += 1
                print(f"匹配成功: {sag_id} -> {taxonomy_dict[sag_id]['classification']}")
            else:
                # 如果没有找到匹配，添加空的分类信息
                member['taxonomy'] = {
                    'classification': 'Not classified',
                    'classification_method': 'N/A',
                    'closest_genome_reference': 'N/A',
                    'closest_genome_ani': 'N/A',
                    'msa_percent': 'N/A',
                    'red_value': 'N/A',
                    'warnings': 'No GTDB-Tk result'
                }
                print(f"未找到匹配: {sag_id}")
    
    print(f"\n匹配统计:")
    print(f"总成员数: {total_members}")
    print(f"成功匹配: {matched_count}")
    print(f"匹配率: {matched_count/total_members*100:.1f}%")
    
    return cluster_data

def save_updated_data(data, output_file):
    """保存更新后的数据"""
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"成功保存更新后的数据到: {output_file}")
    except Exception as e:
        print(f"保存数据时出错: {e}")

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='将GTDB-Tk分类信息添加到cluster数据中',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python merge_taxonomy.py gtdbtk.bac120.summary.tsv cluster_data_updated.json output.json
        """
    )
    
    parser.add_argument('gtdbtk_file', help='GTDB-Tk分类结果文件 (.tsv)')
    parser.add_argument('cluster_file', help='cluster数据文件 (.json)')
    parser.add_argument('output_file', help='输出文件名 (.json)')
    
    # 如果没有提供参数，显示帮助信息
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    # 检查文件是否存在
    if not Path(args.gtdbtk_file).exists():
        print(f"错误: 找不到文件 {args.gtdbtk_file}")
        sys.exit(1)
    
    if not Path(args.cluster_file).exists():
        print(f"错误: 找不到文件 {args.cluster_file}")
        sys.exit(1)
    
    print("开始处理数据...")
    print(f"GTDB-Tk文件: {args.gtdbtk_file}")
    print(f"Cluster文件: {args.cluster_file}")
    print(f"输出文件: {args.output_file}")
    
    # 加载数据
    taxonomy_dict = load_gtdbtk_data(args.gtdbtk_file)
    if not taxonomy_dict:
        print("无法加载GTDB-Tk数据，退出")
        sys.exit(1)
    
    cluster_data = load_cluster_data(args.cluster_file)
    if not cluster_data:
        print("无法加载cluster数据，退出")
        sys.exit(1)
    
    # 合并数据
    updated_data = merge_taxonomy_info(cluster_data, taxonomy_dict)
    
    # 保存结果
    save_updated_data(updated_data, args.output_file)
    
    print(f"\n处理完成！输出文件: {args.output_file}")

if __name__ == "__main__":
    main()