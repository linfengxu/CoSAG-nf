#!/usr/bin/env python3
import json
import os
import sys
import argparse

def filter_and_export_clusters(input_file, output_dir="filtered_clusters"):
    """
    筛选满足条件的cluster并生成单独的JSON文件
    条件: completeness > 90 且 contamination < 10
    """
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 读取输入JSON文件
    with open(input_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    filtered_count = 0
    
    # 遍历所有clusters
    for cluster in data.get('clusters', []):
        cluster_id = cluster.get('cluster_id')
        co_assembly_results = cluster.get('co_assembly_results', {})
        
        # 检查round_1的结果
        round_1_results = co_assembly_results.get('round_1', {})
        completeness = round_1_results.get('completeness', 0)
        contamination = round_1_results.get('contamination', 100)
        
        # 筛选条件: completeness > 90 且 contamination > 10
        if completeness > 90 and contamination > 10:
            # 创建单独的JSON文件
            output_filename = f"cluster_{cluster_id}.json"
            output_path = os.path.join(output_dir, output_filename)
            
            # 写入单个cluster的数据
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(cluster, f, indent=2, ensure_ascii=False)
            
            filtered_count += 1
            print(f"已导出 cluster_{cluster_id}: completeness={completeness:.2f}, contamination={contamination:.2f}")
    
    print(f"\n总共筛选出 {filtered_count} 个符合条件的clusters")
    print(f"输出目录: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="筛选满足条件的cluster并生成单独的JSON文件")
    parser.add_argument("input_file", help="输入的JSON文件路径")
    parser.add_argument("output_dir", help="输出目录路径")
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input_file):
        print(f"错误: 输入文件 '{args.input_file}' 不存在")
        sys.exit(1)
    
    # 运行筛选
    filter_and_export_clusters(args.input_file, args.output_dir)

if __name__ == "__main__":
    main()