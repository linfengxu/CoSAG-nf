#!/usr/bin/env python3
"""
SAG Cluster Optimizer - 基于tetranucleotide signatures的迭代优化
使用层次聚类识别并移除不相似的SAG，降低co-assembly污染度

使用方法:
python cluster_optimizer.py --json dat/json/144.json --outdir results

目标: Completeness ≥ 90% AND Contamination ≤ 5%
备选: 污染度 < 10% 的最高完整度结果
"""

import json
import numpy as np
import subprocess
import os
from itertools import product
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
import tempfile
import shutil
import logging
from datetime import datetime

class ClusterOptimizer:
    def __init__(self, cluster_json_file, output_dir="optimized_results", 
                 target_contamination=5.0, min_completeness=90.0, max_iterations=10):
        self.cluster_file = cluster_json_file
        self.output_dir = output_dir
        self.target_contamination = target_contamination
        self.min_completeness = min_completeness
        self.max_iterations = max_iterations
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 创建 all_contig 文件夹
        self.all_contig_dir = os.path.join(output_dir, "all_contig")
        os.makedirs(self.all_contig_dir, exist_ok=True)
        
        # 设置日志
        self.setup_logging()
        
        # 加载cluster数据
        with open(cluster_json_file, 'r') as f:
            self.cluster_data = json.load(f)
        
        self.cluster_id = self.cluster_data['cluster_id']
        self.logger.info(f"Optimizing cluster {self.cluster_id} with {self.cluster_data['cluster_size']} SAGs")
        self.logger.info(f"Target: Completeness ≥ {self.min_completeness}% AND Contamination ≤ {self.target_contamination}%")
        self.logger.info(f"Max iterations: {self.max_iterations}")
        
        print(f"Optimizing cluster {self.cluster_id} with {self.cluster_data['cluster_size']} SAGs")
        print(f"Target: Completeness ≥ {self.min_completeness}% AND Contamination ≤ {self.target_contamination}%")
        print(f"Max iterations: {self.max_iterations}")
    
    def setup_logging(self):
        """设置日志记录"""
        log_file = os.path.join(self.output_dir, f"optimization_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        
        # 创建logger
        self.logger = logging.getLogger('ClusterOptimizer')
        self.logger.setLevel(logging.INFO)
        
        # 清除已有的handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)
        
        # 创建文件handler
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.INFO)
        
        # 创建控制台handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # 创建formatter
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        
        # 添加handlers
        self.logger.addHandler(file_handler)
        
        self.logger.info("=== SAG Cluster Optimization Started ===")
        self.logger.info(f"Log file: {log_file}")
        self.logger.info(f"All contigs directory: {self.all_contig_dir}")
        
        print(f"Log file created: {log_file}")
        print(f"All contigs will be saved to: {self.all_contig_dir}")
    
    def calculate_tetranucleotide_signature(self, fasta_file):
        """计算四核苷酸频率特征"""
        nucleotides = ['A', 'T', 'G', 'C']
        tetranucleotides = [''.join(p) for p in product(nucleotides, repeat=4)]
        
        tetra_counts = {tetra: 0 for tetra in tetranucleotides}
        total_tetras = 0
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq = str(record.seq).upper()
                for i in range(len(seq) - 3):
                    tetra = seq[i:i+4]
                    if tetra in tetra_counts and 'N' not in tetra:
                        tetra_counts[tetra] += 1
                        total_tetras += 1
        except Exception as e:
            print(f"Error reading {fasta_file}: {e}")
            return None
        
        if total_tetras == 0:
            return None
            
        # 标准化为频率
        tetra_freq = np.array([tetra_counts[tetra] / total_tetras for tetra in tetranucleotides])
        return tetra_freq
    
    def calculate_signatures_for_all_sags(self):
        """计算所有SAG的tetranucleotide signatures"""
        print("Step 1: Calculating tetranucleotide signatures...")
        
        signatures = {}
        valid_sags = []
        
        for member in self.cluster_data['members']:
            sag_id = member['sag_id']
            contig_file = member['individual_contigs']
            
            print(f"  Processing {sag_id}...")
            signature = self.calculate_tetranucleotide_signature(contig_file)
            
            if signature is not None:
                signatures[sag_id] = signature
                valid_sags.append(member)
            else:
                print(f"  Warning: Skipping {sag_id} due to empty/invalid contigs")
        
        print(f"  Successfully processed {len(signatures)} SAGs")
        return signatures, valid_sags
    
    def perform_hierarchical_clustering(self, signatures):
        """使用层次聚类识别主要亚群"""
        print("Step 2: Performing hierarchical clustering...")
        
        sag_ids = list(signatures.keys())
        signature_matrix = np.array([signatures[sag_id] for sag_id in sag_ids])
        
        # 计算距离矩阵
        distances = pdist(signature_matrix, metric='euclidean')
        linkage_matrix = linkage(distances, method='ward')
        
        # 聚类分析 - 分成2-3个主要群体
        clusters = fcluster(linkage_matrix, t=3, criterion='maxclust')
        
        # 找到最大的群体作为主群
        cluster_counts = {}
        for i, cluster_id in enumerate(clusters):
            if cluster_id not in cluster_counts:
                cluster_counts[cluster_id] = []
            cluster_counts[cluster_id].append(sag_ids[i])
        
        main_cluster = max(cluster_counts.values(), key=len)
        outliers = []
        for cluster_sags in cluster_counts.values():
            if cluster_sags != main_cluster:
                outliers.extend(cluster_sags)
        
        print(f"  Main cluster: {len(main_cluster)} SAGs")
        print(f"  Outliers: {len(outliers)} SAGs - {outliers}")
        
        return main_cluster, outliers
    
    def identify_candidate_sags_to_remove(self, signatures, sag_list, n_candidates=3):
        """识别候选移除的SAG - 返回多个候选而不是单一选择"""
        if len(sag_list) <= 5:  # 保留最少5个SAG
            return []
        
        sag_ids = list(sag_list)
        signature_matrix = np.array([signatures[sag_id] for sag_id in sag_ids])
        
        # 计算每个SAG与其他SAG的平均距离
        avg_distances = {}
        for i, sag_id in enumerate(sag_ids):
            distances = []
            for j, other_sag in enumerate(sag_ids):
                if i != j:
                    dist = np.linalg.norm(signature_matrix[i] - signature_matrix[j])
                    distances.append(dist)
            avg_distances[sag_id] = np.mean(distances)
        
        # 返回距离最大的几个SAG作为候选
        sorted_sags = sorted(avg_distances.items(), key=lambda x: x[1], reverse=True)
        max_candidates = min(n_candidates, len(sag_list) - 5)
        candidates = [sag for sag, _ in sorted_sags[:max_candidates]]
        
        return candidates
    
    def test_sag_combination(self, selected_sags, iteration, test_id=""):
        """测试特定SAG组合的效果"""
        self.logger.info(f"Testing combination {test_id}: {len(selected_sags)} SAGs")
        self.logger.info(f"Selected SAGs: {selected_sags}")
        print(f"    Testing combination {test_id}: {len(selected_sags)} SAGs")
        
        # 进行co-assembly和评估
        r1_file, r2_file, temp_dir = self.create_merged_fastq_files(selected_sags)
        
        try:
            # SPAdes assembly
            output_prefix = f"iter_{iteration}_{test_id}_{len(selected_sags)}sags"
            contigs_file = self.run_spades_assembly(r1_file, r2_file, output_prefix)
            
            if contigs_file:
                # CheckM2 evaluation
                completeness, contamination ,genome_size,contig_n50 = self.run_checkm2_evaluation(contigs_file, output_prefix)
                
                if completeness is not None and contamination is not None:
                    result = {
                        'iteration': iteration,
                        'test_id': test_id,
                        'sag_count': len(selected_sags),
                        'sags': selected_sags.copy(),
                        'completeness': completeness,
                        'contamination': contamination,
                        'genome_size':genome_size,
                        'contig_n50':contig_n50,
                        'contigs_file': contigs_file,
                        'output_prefix': output_prefix
                    }
                    
                    result_msg = f"Results {test_id}: Completeness={completeness:.2f}%, Contamination={contamination:.2f}%"
                    self.logger.info(result_msg)
                    print(f"    {result_msg}")
                    return result
        
        except Exception as e:
            error_msg = f"Error testing combination {test_id}: {e}"
            self.logger.error(error_msg)
            print(f"    {error_msg}")
        
        finally:
            # 清理临时文件
            shutil.rmtree(temp_dir, ignore_errors=True)
        
        return None
    
    def create_merged_fastq_files(self, selected_sags):
        """合并选定SAG的fastq文件"""
        temp_dir = tempfile.mkdtemp()
        
        merged_r1 = os.path.join(temp_dir, f"cluster_{self.cluster_id}_merged_R1.fastq")
        merged_r2 = os.path.join(temp_dir, f"cluster_{self.cluster_id}_merged_R2.fastq")
        
        with open(merged_r1, 'w') as out_r1, open(merged_r2, 'w') as out_r2:
            for member in self.cluster_data['members']:
                if member['sag_id'] in selected_sags:
                    # 合并R1
                    with open(member['read1'], 'r') as r1:
                        shutil.copyfileobj(r1, out_r1)
                    # 合并R2
                    with open(member['read2'], 'r') as r2:
                        shutil.copyfileobj(r2, out_r2)
        
        return merged_r1, merged_r2, temp_dir
    
    def run_spades_assembly(self, r1_file, r2_file, output_prefix):
        """运行SPAdes进行co-assembly"""
        assembly_dir = os.path.join(self.output_dir, f"{output_prefix}_spades")
        
        cmd = [
            "singularity",'exec','-B','/cpfs01',
            '/cpfs01/projects-SSD/cfff-86962b7a8e68_SSD/public/singularity_sif/spades_3.15.5.sif',
            "spades.py",
            "--sc", "--careful",
            "-1", r1_file,
            "-2", r2_file,
            "-o", assembly_dir
        ]
        
        self.logger.info(f"Running SPAdes assembly: {output_prefix}")
        self.logger.info(f"Command: {' '.join(cmd)}")
        print(f"  Running SPAdes: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            contigs_file = os.path.join(assembly_dir, "contigs.fasta")
            
            if os.path.exists(contigs_file):
                contigs_file_abs = os.path.abspath(contigs_file)
                self.logger.info(f"SPAdes assembly completed: {contigs_file}")
                print(f"  ✓ Assembly completed: {contigs_file}")
                
                return contigs_file_abs
            else:
                error_msg = f"contigs.fasta not found in {assembly_dir}"
                self.logger.error(error_msg)
                print(f"  Error: {error_msg}")
                return None
                
        except subprocess.CalledProcessError as e:
            error_msg = f"SPAdes failed: {e}"
            self.logger.error(error_msg)
            self.logger.error(f"stderr: {e.stderr}")
            print(f"  SPAdes failed: {e}")
            print(f"  stderr: {e.stderr}")
            return None
    
    def run_checkm2_evaluation(self, contigs_file, output_prefix):
        """运行CheckM2评估基因组质量"""
        checkm_dir = os.path.join(self.output_dir, f"{output_prefix}_checkm2")
        os.makedirs(checkm_dir, exist_ok=True)
        
        cmd = [
            "singularity",'exec','-B','/cpfs01',
            '/cpfs01/projects-SSD/cfff-86962b7a8e68_SSD/public/singularity_sif/checkm2_1.0.2.sif',
            "checkm2", "predict",
            "--input", contigs_file,
            "--output-directory", checkm_dir,
            "--threads", "20",
            "--database_path",'/cpfs01/projects-SSD/cfff-86962b7a8e68_SSD/public/checkM_db/CheckM2_database/uniref100.KO.1.dmnd',
            "--force"
        ]
        
        self.logger.info(f"Running CheckM2 evaluation: {output_prefix}")
        self.logger.info(f"Command: {' '.join(cmd)}")
        print(f"  Running CheckM2: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 解析CheckM2结果
            quality_file = os.path.join(checkm_dir, "quality_report.tsv")
            if os.path.exists(quality_file):
                with open(quality_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) > 1:
                        data = lines[1].strip().split('\t')
                        completeness = float(data[1])
                        contamination = float(data[2])
                        genome_size = float(data[8])
                        contig_n50 =  float(data[6])
                        self.logger.info(f"CheckM2 evaluation completed: Completeness={completeness:.2f}%, Contamination={contamination:.2f}%")
                        return completeness, contamination,genome_size,contig_n50
            
            error_msg = f"CheckM2 quality report not found or empty: {quality_file}"
            self.logger.error(error_msg)
            return None, None,None, None
            
        except subprocess.CalledProcessError as e:
            error_msg = f"CheckM2 failed: {e}"
            self.logger.error(error_msg)
            print(f"  CheckM2 failed: {e}")
            return None, None,None, None
    
    def optimize_cluster(self):
        """主要的优化流程"""
        print(f"\n=== Optimizing Cluster {self.cluster_id} ===")
        
        # Step 1: 计算tetranucleotide signatures
        signatures, valid_sags = self.calculate_signatures_for_all_sags()
        
        if len(signatures) < 5:
            print("Error: Too few valid SAGs for optimization")
            return
        
        # Step 2: 层次聚类识别主群和离群点
        main_cluster, outliers = self.perform_hierarchical_clustering(signatures)
        
        # Step 3: 从主群开始，移除最不相似的SAG
        current_sags = [sag for sag in signatures.keys()]
        iteration = 1
        best_result = None
        best_under_10_contamination = None  # 污染度<10%的最佳结果
        
        # 用于跟踪完整度趋势的变量
        completeness_history = []
        has_acceptable_result = False          # 是否已经达到可接受结果
        acceptable_max_completeness = 0        # 可接受结果中的最高完整度
        completeness_decline_threshold = 10.0  # 完整度下降超过10%才停止
        min_completeness_threshold = 80.0      # 绝对最低完整度阈值
        
        # 首先移除明显的离群点
        if outliers:
            print(f"\nStep 3: Removing obvious outliers: {outliers}")
            for outlier in outliers:
                if outlier in current_sags:
                    current_sags.remove(outlier)
        
        while len(current_sags) >= 5 and iteration <= self.max_iterations:
            print(f"\n--- Iteration {iteration}: Testing {len(current_sags)} SAGs ---")
            print(f"Current SAGs: {current_sags}")
            
            # Step 4: 测试当前组合
            current_result = self.test_sag_combination(current_sags, iteration, "current")
            
            if current_result:
                # 记录完整度历史
                completeness_history.append(current_result['completeness'])
                
                # 检查是否达到可接受结果 (污染度≤10%)
                is_acceptable = current_result['contamination'] <= 10.0
                if is_acceptable:
                    if not has_acceptable_result:
                        has_acceptable_result = True
                        acceptable_max_completeness = current_result['completeness']
                        print(f"  ✓ First acceptable result achieved (contamination ≤ 10%)")
                    else:
                        # 更新可接受结果中的最高完整度
                        if current_result['completeness'] > acceptable_max_completeness:
                            acceptable_max_completeness = current_result['completeness']
                
                # 检查是否应该停止迭代
                should_stop_due_to_decline = False
                current_completeness = current_result['completeness']
                
                # 只有在已经达到可接受结果的情况下才考虑因完整度下降而停止
                if has_acceptable_result:
                    completeness_drop = acceptable_max_completeness - current_completeness
                    if completeness_drop > completeness_decline_threshold:
                        should_stop_due_to_decline = True
                        print(f"  ⚠ Significant completeness decline from acceptable results:")
                        print(f"    Best acceptable: {acceptable_max_completeness:.1f}% → Current: {current_completeness:.1f}%")
                        print(f"    Drop: {completeness_drop:.1f}% (threshold: {completeness_decline_threshold:.1f}%)")
                        print(f"  → Stopping iteration to preserve genome quality")
                
                # 绝对最低完整度检查（无论是否有可接受结果）
                if current_completeness < min_completeness_threshold:
                    should_stop_due_to_decline = True
                    print(f"  ⚠ Completeness below minimum threshold: {current_completeness:.1f}% < {min_completeness_threshold:.1f}%")
                    print(f"  → Stopping iteration")
                # 更新最佳结果
                is_optimal = (current_result['completeness'] >= self.min_completeness and 
                            current_result['contamination'] <= self.target_contamination)
                
                if is_optimal:
                    if (best_result is None or 
                        not (best_result.get('completeness', 0) >= self.min_completeness and 
                             best_result.get('contamination', 100) <= self.target_contamination) or
                        current_result['completeness'] > best_result.get('completeness', 0)):
                        best_result = current_result
                        print(f"  ✓ New optimal result found!")
                
                # 如果还没有最优结果，保存污染度<10%的最佳结果作为备选
                if current_result['contamination'] <= 10.0:
                    if (best_under_10_contamination is None or 
                        current_result['completeness'] > best_under_10_contamination.get('completeness', 0)):
                        best_under_10_contamination = current_result
                
                # 如果达到最优目标，直接停止迭代
                if is_optimal:
                    print(f"  ✓ Optimal target achieved! Stopping iteration.")
                    break
                
                # 如果完整度显著下降，停止迭代
                if should_stop_due_to_decline:
                    break
            
            # Step 5: 探索不同的SAG移除策略
            if len(current_sags) > 5:
                # 获取候选移除的SAG
                candidates = self.identify_candidate_sags_to_remove(signatures, current_sags, n_candidates=3)
                
                if not candidates:
                    break
                
                print(f"  Testing removal of different SAGs: {candidates}")
                
                # 测试移除不同SAG的效果
                test_results = []
                for i, candidate in enumerate(candidates):
                    test_sags = [sag for sag in current_sags if sag != candidate]
                    test_result = self.test_sag_combination(test_sags, iteration, f"remove_{candidate[:8]}")
                    
                    if test_result:
                        test_results.append((candidate, test_result))
                
                # 检查所有测试结果是否都显示完整度显著下降
                if test_results and has_acceptable_result:
                    max_test_completeness = max(res['completeness'] for _, res in test_results)
                    
                    # 只有在已经有可接受结果的情况下，才考虑预测性停止
                    completeness_drop_prediction = acceptable_max_completeness - max_test_completeness
                    
                    if completeness_drop_prediction > completeness_decline_threshold:
                        print(f"  ⚠ All removal options would cause significant completeness decline:")
                        print(f"    Best acceptable: {acceptable_max_completeness:.1f}% → Best next option: {max_test_completeness:.1f}%")
                        print(f"    Predicted drop: {completeness_drop_prediction:.1f}% (threshold: {completeness_decline_threshold:.1f}%)")
                        
                        # 如果当前结果是可接受的，就停止迭代
                        if current_result and current_result['contamination'] <= 10.0:
                            print(f"  → Stopping iteration (preserving acceptable result)")
                            break
                
                # 选择最佳的移除策略
                if test_results:
                    # 优先选择达到最优目标的结果
                    optimal_results = [(cand, res) for cand, res in test_results 
                                     if res['completeness'] >= self.min_completeness and 
                                        res['contamination'] <= self.target_contamination]
                    
                    if optimal_results:
                        # 选择完整度最高的最优结果
                        best_candidate, best_test_result = max(optimal_results, 
                                                             key=lambda x: x[1]['completeness'])
                        print(f"  ✓ Found optimal combination by removing {best_candidate}")
                        current_sags.remove(best_candidate)
                        
                        # 更新最佳结果
                        best_result = best_test_result
                        print(f"  ✓ Optimal target achieved! Stopping iteration.")
                        break
                    else:
                        # 优先选择污染度低于10%且完整度最高的结果
                        under_10_results = [(cand, res) for cand, res in test_results 
                                          if res['contamination'] <= 10.0]
                        
                        if under_10_results:
                            # 在污染度<10%的结果中选择完整度最高的
                            best_candidate, best_test_result = max(under_10_results, 
                                                                 key=lambda x: x[1]['completeness'])
                            print(f"  → Best option (contamination <10%): remove {best_candidate} "
                                  f"(Completeness={best_test_result['completeness']:.1f}%, "
                                  f"Contamination={best_test_result['contamination']:.1f}%)")
                        else:
                            # 如果没有污染度<10%的结果，选择污染度最低的
                            best_candidate, best_test_result = min(test_results, 
                                                                 key=lambda x: x[1]['contamination'])
                            print(f"  → Best available option: remove {best_candidate} "
                                  f"(Completeness={best_test_result['completeness']:.1f}%, "
                                  f"Contamination={best_test_result['contamination']:.1f}%)")
                        
                        current_sags.remove(best_candidate)
                        
                        # 更新备选结果
                        if best_test_result['contamination'] <= 10.0:
                            if (best_under_10_contamination is None or 
                                best_test_result['completeness'] > best_under_10_contamination.get('completeness', 0)):
                                best_under_10_contamination = best_test_result
                else:
                    # 如果所有测试都失败，使用原来的策略
                    print(f"  All tests failed, removing first candidate: {candidates[0]}")
                    current_sags.remove(candidates[0])
            
            iteration += 1
        
        # 选择最终结果
        final_result = None
        result_type = ""
        
        # 优先选择满足最优条件的结果
        if best_result and (best_result.get('completeness', 0) >= self.min_completeness and 
                           best_result.get('contamination', 100) <= self.target_contamination):
            final_result = best_result
            result_type = "OPTIMAL"
        # 如果没有最优结果，选择污染度<10%的最佳结果
        elif best_under_10_contamination:
            final_result = best_under_10_contamination
            result_type = "BEST_UNDER_10%_CONTAMINATION"
        # 最后选择任何可用的最佳结果
        elif best_result:
            final_result = best_result
            result_type = "BEST_AVAILABLE"
        
        # 输出最终结果
        if final_result:
            print(f"\n=== Optimization Complete ===")
            print(f"Result type: {result_type}")
            print(f"Total iterations: {iteration - 1}")
            print(f"Final result: {final_result['sag_count']} SAGs")
            print(f"Completeness: {final_result['completeness']:.2f}%")
            print(f"Contamination: {final_result['contamination']:.2f}%")
            
            # 评估结果质量
            if (final_result['completeness'] >= self.min_completeness and 
                final_result['contamination'] <= self.target_contamination):
                print("✓ EXCELLENT: Meets both completeness and contamination targets!")
            elif final_result['contamination'] <= 10.0:
                print("✓ GOOD: Low contamination achieved")
            else:
                print("⚠ ACCEPTABLE: Best available result")
            
            print(f"Selected SAGs: {final_result['sags']}")
            
            # 保存优化结果
            result_file = os.path.join(self.output_dir, f"cluster_{self.cluster_id}_optimized.json")
            final_result['result_type'] = result_type
            final_result['optimization_summary'] = {
                'total_iterations': iteration - 1,
                'max_iterations': self.max_iterations,
                'target_completeness': self.min_completeness,
                'target_contamination': self.target_contamination,
                'meets_targets': (final_result['completeness'] >= self.min_completeness and 
                                final_result['contamination'] <= self.target_contamination)
            }
            
            with open(result_file, 'w') as f:
                json.dump(final_result, f, indent=2)
            
            # 保存最佳结果的contigs文件到all_contig文件夹
            self.save_best_contigs(final_result)
            
            # 生成contigs汇总文件
            self.generate_best_contigs_summary(final_result)
            
            # 记录最终结果到日志
            self.logger.info("=== Optimization Complete ===")
            self.logger.info(f"Result type: {result_type}")
            self.logger.info(f"Total iterations: {iteration - 1}")
            self.logger.info(f"Final result: {final_result['sag_count']} SAGs")
            self.logger.info(f"Completeness: {final_result['completeness']:.2f}%")
            self.logger.info(f"Contamination: {final_result['contamination']:.2f}%")
            self.logger.info(f"Selected SAGs: {final_result['sags']}")
            self.logger.info(f"Results saved to: {result_file}")
            
            print(f"Results saved to: {result_file}")
            
            # 显示优化总结
            print(f"\n=== Optimization Summary ===")
            print(f"Completeness trend: {' → '.join([f'{c:.1f}%' for c in completeness_history[-5:]])}")
            if len(completeness_history) >= 2:
                trend = completeness_history[-1] - completeness_history[0]
                print(f"Overall completeness change: {trend:+.1f}%")
                self.logger.info(f"Overall completeness change: {trend:+.1f}%")
            if has_acceptable_result:
                print(f"Best acceptable completeness achieved: {acceptable_max_completeness:.1f}%")
                self.logger.info(f"Best acceptable completeness achieved: {acceptable_max_completeness:.1f}%")
            else:
                print("No acceptable results (contamination ≤ 10%) achieved during optimization")
                self.logger.warning("No acceptable results (contamination ≤ 10%) achieved during optimization")
        else:
            self.logger.error("No valid results obtained during optimization")
            print("No valid results obtained during optimization")
    
    def save_best_contigs(self, final_result):
        """保存最佳结果的contigs文件到all_contig文件夹"""
        if final_result and 'contigs_file' in final_result:
            original_contigs = final_result['contigs_file']
            
            if os.path.exists(original_contigs):
                # 生成最佳结果的文件名
                best_contigs_name = f"cluster_{self.cluster_id}_best_contigs.fasta"
                best_contigs_path = os.path.join(self.all_contig_dir, best_contigs_name)
                
                # 复制最佳结果的contigs文件
                shutil.copy2(original_contigs, best_contigs_path)
                
                # 记录文件信息
                file_size = os.path.getsize(best_contigs_path)
                
                self.logger.info(f"Best contigs file saved: {best_contigs_path}")
                self.logger.info(f"File size: {file_size / (1024 * 1024):.2f} MB")
                
                print(f"✓ Best contigs saved to: {best_contigs_path}")
                print(f"  File size: {file_size / (1024 * 1024):.2f} MB")
                
                return best_contigs_path
            else:
                self.logger.error(f"Best contigs file not found: {original_contigs}")
                print(f"Error: Best contigs file not found: {original_contigs}")
                return None
        else:
            self.logger.warning("No valid final result to save contigs")
            return None
    
    def generate_best_contigs_summary(self, final_result):
        """生成最佳contigs文件汇总"""
        summary_file = os.path.join(self.all_contig_dir, "best_contigs_summary.json")
        
        # 获取最佳contigs文件路径
        best_contigs_name = f"cluster_{self.cluster_id}_best_contigs.fasta"
        best_contigs_path = os.path.join(self.all_contig_dir, best_contigs_name)
        
        summary_data = {
            'cluster_id': self.cluster_id,
            'generation_time': datetime.now().isoformat(),
            'result_type': final_result.get('result_type', 'UNKNOWN'),
            'optimization_summary': {
                'completeness': final_result.get('completeness', 0),
                'contamination': final_result.get('contamination', 0),
                'sag_count': final_result.get('sag_count', 0),
                'selected_sags': final_result.get('sags', []),
                'iteration': final_result.get('iteration', 0),
                'test_id': final_result.get('test_id', ''),
                'output_prefix': final_result.get('output_prefix', '')
            },
            'best_contigs_file': {
                'filename': best_contigs_name,
                'path': best_contigs_path,
                'original_path': final_result.get('contigs_file', ''),
                'exists': os.path.exists(best_contigs_path)
            }
        }
        
        # 添加文件大小信息
        if os.path.exists(best_contigs_path):
            file_size = os.path.getsize(best_contigs_path)
            summary_data['best_contigs_file']['file_size_bytes'] = file_size
            summary_data['best_contigs_file']['file_size_mb'] = round(file_size / (1024 * 1024), 2)
        
        with open(summary_file, 'w') as f:
            json.dump(summary_data, f, indent=2)
        
        # 生成简单的文本汇总
        txt_summary_file = os.path.join(self.all_contig_dir, "README.txt")
        with open(txt_summary_file, 'w') as f:
            f.write(f"Cluster {self.cluster_id} - Best Optimization Result\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Result type: {final_result.get('result_type', 'UNKNOWN')}\n\n")
            
            f.write("Quality Metrics:\n")
            f.write(f"  Completeness: {final_result.get('completeness', 0):.2f}%\n")
            f.write(f"  Contamination: {final_result.get('contamination', 0):.2f}%\n")
            f.write(f"  SAG count: {final_result.get('sag_count', 0)}\n\n")
            
            f.write("Selected SAGs:\n")
            for i, sag in enumerate(final_result.get('sags', []), 1):
                f.write(f"  {i}. {sag}\n")
            
            f.write(f"\nBest contigs file: {best_contigs_name}\n")
            if os.path.exists(best_contigs_path):
                file_size = os.path.getsize(best_contigs_path)
                f.write(f"File size: {file_size / (1024 * 1024):.2f} MB\n")
            
            f.write(f"\nOptimization details:\n")
            f.write(f"  Iteration: {final_result.get('iteration', 0)}\n")
            f.write(f"  Test ID: {final_result.get('test_id', '')}\n")
            f.write(f"  Output prefix: {final_result.get('output_prefix', '')}\n")
        
        self.logger.info(f"Best contigs summary generated: {summary_file}")
        self.logger.info(f"README file generated: {txt_summary_file}")
        
        print(f"✓ Summary files generated:")
        print(f"  - {os.path.basename(summary_file)}")
        print(f"  - {os.path.basename(txt_summary_file)}")
    


def main():
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(description='SAG Cluster Optimizer')
    parser.add_argument('--json', required=True, help='Path to cluster JSON file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.json):
        print(f"Error: File {args.json} not found")
        sys.exit(1)
    
    # 使用指定参数运行优化
    print("=== SAG Cluster Optimizer ===")
    print(f"Input: {args.json}")
    print(f"Output: {args.outdir}")
    print("Target: Completeness ≥ 90% AND Contamination ≤ 5%")
    print("Max iterations: 10")
    print()
    
    optimizer = ClusterOptimizer(
        cluster_json_file=args.json,
        output_dir=args.outdir,
        target_contamination=5.0,
        min_completeness=90.0,
        max_iterations=10
    )
    optimizer.optimize_cluster()

if __name__ == "__main__":
    main()