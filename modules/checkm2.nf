process CHECKM2 {
    tag "${meta.id}"
    label 'process_medium'
    
    container { params.checkm2.container }
    publishDir "${params.output_structure?.individual_assemblies ?: params.outdir}/per_sample/${meta.id}/checkm2",
        mode: 'copy'

    input:
    tuple val(meta), path(contig)
    
    output:
    tuple val(meta), path("*checkm2.tsv"), emit: results, optional: true
    path("*.log"), emit: logs
    path("*_stats.txt"), emit: stats, optional: true
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def threads = params.checkm2?.threads ?: 8
    def db_path = params.checkm2?.database
    def outdir = params.checkm2?.output_dir ?: "checkm2_output"
    def additional_args = params.checkm2?.additional_args ?: ""
    def process_name = task.process
    def min_contig_length = params.checkm2?.min_contig_length ?: 1000
    def min_genome_size = params.checkm2?.min_genome_size ?: 50000
    
    """
    # 设置CheckM2数据库环境变量
    export CHECKM2DB=${db_path}
    
    # 创建输出目录
    mkdir -p ${outdir}
    
    # 检查输入文件
    echo "=== Input file diagnostics ===" > ${prefix}_stats.txt
    echo "Input file: ${contig}" >> ${prefix}_stats.txt
    echo "File size: \$(stat -c%s ${contig}) bytes" >> ${prefix}_stats.txt
    echo "Number of sequences: \$(grep -c '^>' ${contig} || echo 0)" >> ${prefix}_stats.txt
    
    # 统计序列长度
    awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq\$0} END {if (seq) print length(seq)}' ${contig} > seq_lengths.tmp
    
    if [ -s seq_lengths.tmp ]; then
        echo "Sequence length statistics:" >> ${prefix}_stats.txt
        echo "  Min length: \$(sort -n seq_lengths.tmp | head -1)" >> ${prefix}_stats.txt
        echo "  Max length: \$(sort -n seq_lengths.tmp | tail -1)" >> ${prefix}_stats.txt
        echo "  Total length: \$(awk '{sum+=\$1} END {print sum}' seq_lengths.tmp)" >> ${prefix}_stats.txt
        echo "  Sequences >= ${min_contig_length}bp: \$(awk '\$1>=${min_contig_length}' seq_lengths.tmp | wc -l)" >> ${prefix}_stats.txt
        
        # 检查是否有足够长的序列
        long_contigs=\$(awk '\$1>=${min_contig_length}' seq_lengths.tmp | wc -l)
        total_size=\$(awk '{sum+=\$1} END {print sum}' seq_lengths.tmp)
        
        echo "Long contigs (>=${min_contig_length}bp): \$long_contigs" >> ${prefix}_stats.txt
        echo "Total assembly size: \$total_size bp" >> ${prefix}_stats.txt
        
        if [ \$long_contigs -eq 0 ]; then
            echo "WARNING: No contigs >= ${min_contig_length}bp found" >> ${prefix}_stats.txt
            echo "CheckM2 may fail due to insufficient sequence length" >> ${prefix}_stats.txt
        fi
        
        if [ \$total_size -lt ${min_genome_size} ]; then
            echo "WARNING: Total assembly size (\$total_size bp) < ${min_genome_size}bp" >> ${prefix}_stats.txt
            echo "CheckM2 may fail due to insufficient genome size" >> ${prefix}_stats.txt
        fi
    else
        echo "ERROR: No sequences found in input file" >> ${prefix}_stats.txt
        echo "CheckM2 will likely fail" >> ${prefix}_stats.txt
    fi
    
    # 检查数据库
    echo "=== Database check ===" >> ${prefix}_stats.txt
    echo "CHECKM2DB environment variable: \$CHECKM2DB" >> ${prefix}_stats.txt
    if [ -f "\$CHECKM2DB" ]; then
        echo "Database file exists: \$CHECKM2DB" >> ${prefix}_stats.txt
        echo "Database file size: \$(stat -c%s \$CHECKM2DB) bytes" >> ${prefix}_stats.txt
    else
        echo "ERROR: Database file not found: \$CHECKM2DB" >> ${prefix}_stats.txt
    fi
    
    # 运行CheckM2
    echo "=== Running CheckM2 ===" >> ${prefix}_stats.txt
    echo "Command: checkm2 predict -x fasta --threads ${threads} --input ${contig} --output-directory ${outdir} ${additional_args}" >> ${prefix}_stats.txt
    
    set +e  # 允许命令失败而不立即退出
    
    checkm2 predict -x fasta --threads ${threads} \\
        --input ${contig} \\
        --output-directory ${outdir} \\
        ${additional_args} > ${prefix}_checkm2.log 2>&1
    
    exit_code=\$?
    echo "CheckM2 exit code: \$exit_code" >> ${prefix}_stats.txt
    
    set -e  # 重新启用错误时退出
    
    # 检查输出结果
    echo "=== Output check ===" >> ${prefix}_stats.txt
    if [ -d "${outdir}" ]; then
        echo "Output directory contents:" >> ${prefix}_stats.txt
        find ${outdir} -type f -name "*.tsv" >> ${prefix}_stats.txt 2>&1
        
        # 寻找结果文件
        result_file=\$(find ${outdir} -maxdepth 1 -name "*.tsv" -type f | head -1)
        
        if [ -n "\$result_file" ] && [ -s "\$result_file" ]; then
            echo "Found result file: \$result_file" >> ${prefix}_stats.txt
            echo "Result file size: \$(stat -c%s \$result_file) bytes" >> ${prefix}_stats.txt
            cp "\$result_file" ${prefix}_checkm2.tsv
            echo "Results successfully copied" >> ${prefix}_stats.txt
        else
            echo "No valid result file found" >> ${prefix}_stats.txt
            # 创建一个空的结果文件以避免流程中断
            echo -e "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes" > ${prefix}_checkm2.tsv
            echo -e "${prefix}\\t0\\t0\\tNA\\tNA\\t0\\t0\\t0\\t0\\t0\\t0\\tCheckM2 failed - no DIAMOND annotation generated" >> ${prefix}_checkm2.tsv
        fi
    else
        echo "Output directory not created" >> ${prefix}_stats.txt
        # 创建一个空的结果文件
        echo -e "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes" > ${prefix}_checkm2.tsv
        echo -e "${prefix}\\t0\\t0\\tNA\\tNA\\t0\\t0\\t0\\t0\\t0\\t0\\tCheckM2 failed - output directory not created" >> ${prefix}_checkm2.tsv
    fi
    
    # 清理临时文件
    rm -f seq_lengths.tmp
    
    # 生成版本信息
    echo "${process_name}:" > versions.yml
    echo "    checkm2: \$(checkm2 --version 2>&1 | head -n 1 | sed 's/.*checkm2 v//' || echo 'unknown')" >> versions.yml
    """
}
