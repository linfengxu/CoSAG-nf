process GTDBTK_CLASSIFYWF_COSAG {
    tag "gtdb_batch_classify"
    label 'process_high'
    
    container "${params.gtdb.container}"
    containerOptions "--bind ${params.gtdb.database}:/refdata"
    
    // 发布目录配置：保存GTDB-Tk分类结果
    publishDir "${params.output_structure?.taxonomic_classification ?: params.outdir}/co_assemblies",
        mode: 'copy',
        pattern: "gtdb_results/*"
    
    
    input:

    path(bin_fastas)  // All bin fasta files collected together
    
    output:
    path "gtdb_results/*", emit: results
    path "gtdb_results/gtdbtk.bac120.summary.tsv", emit: bac_summary, optional: true
    path "gtdb_results/gtdbtk.ar53.summary.tsv", emit: ar_summary, optional: true
    path "gtdb_results/gtdbtk.bac120.classify.tree", emit: bac_tree, optional: true
    path "gtdb_results/gtdbtk.ar53.classify.tree", emit: ar_tree, optional: true
    path "gtdb_results/gtdbtk.log", emit: gtdb_log
    path "versions.yml", emit: versions                                   // 输出：软件版本信息
    
    script:
    def extension = params.gtdb?.extension ?: "fasta"                        // 基因组文件扩展名
    def args = params.gtdb?.additional_args ?: ""                         // 额外参数
    
    """
    # 创建基因组目录并复制文件 task.cpus
    mkdir -p bin_fastas
    mv ${bin_fastas} bin_fastas/
    mkdir -p gtdb_results
    # 运行GTDB-Tk分类工作流
    gtdbtk classify_wf \\
        --genome_dir bin_fastas \\
        --out_dir gtdb_results \\
        --cpus 30 \\
        --extension ${extension} \\
        --force \\
        --skip_ani_screen \\
        ${args} \\
        > gtdb_results/gtdbtk.log 2>&1
    
    # 生成版本信息
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version 2>&1 | sed 's/gtdbtk: version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

