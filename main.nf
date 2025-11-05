#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    CoSAG-nf: A Scalable Nextflow Pipeline for Co-assembly, Optimization, 
    and Interactive Visualization of High-Throughput Single-Cell Genomes
========================================================================================
    Github   : https://github.com/linfengxu/CoSAG-nf
    Zenodo  : https://your-website.com/CoSAG-nf
    Example_Report : http://www.biostatistics.online/CoSAG/example_report.html
----------------------------------------------------------------------------------------
*/

// Global color definitions
class Colors {
    static final String RESET = "\033[0m"
    static final String BOLD = "\033[1m"
    static final String DIM = "\033[2m"
    static final String RED = "\033[31m"
    static final String GREEN = "\033[32m"
    static final String YELLOW = "\033[33m"
    static final String BLUE = "\033[34m"
    static final String PURPLE = "\033[35m"
    static final String CYAN = "\033[36m"
    static final String WHITE = "\033[37m"
}

def helpMessage() {
    log.info"""
${Colors.CYAN}
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   ${Colors.BOLD}${Colors.WHITE}  ██████╗ ██████╗ ███████╗ █████╗  ██████╗       ███╗   ██╗███████╗  ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE} ██╔════╝██╔═══██╗██╔════╝██╔══██╗██╔════╝       ████╗  ██║██╔════╝  ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE} ██║     ██║   ██║███████╗███████║██║  ███╗█████╗██╔██╗ ██║█████╗    ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE} ██║     ██║   ██║╚════██║██╔══██║██║   ██║╚════╝██║╚██╗██║██╔══╝    ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE} ╚██████╗╚██████╔╝███████║██║  ██║╚██████╔╝      ██║ ╚████║██║       ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}  ╚═════╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝ ╚═════╝       ╚═╝  ╚═══╝╚═╝       ${Colors.RESET}${Colors.CYAN}   ║
║                                                                              ║
║   ${Colors.BOLD}${Colors.YELLOW}A Scalable Nextflow Pipeline for Co-assembly, Optimization,${Colors.RESET}${Colors.CYAN}        ║
║   ${Colors.BOLD}${Colors.YELLOW}and Interactive Visualization of High-Throughput Single-Cell Genomes${Colors.RESET}${Colors.CYAN} ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
${Colors.RESET}

${Colors.BOLD}${Colors.GREEN}📋 DESCRIPTION:${Colors.RESET}
${Colors.WHITE}CoSAG-nf is a comprehensive pipeline designed for processing single-cell amplified 
genomes (SAGs). It performs individual assembly, similarity analysis, hierarchical 
clustering, co-assembly optimization, and taxonomic classification to generate 
high-quality metagenome-assembled genomes (MAGs) from single-cell data.${Colors.RESET}

${Colors.BOLD}${Colors.BLUE}🚀 USAGE:${Colors.RESET}
${Colors.CYAN}nextflow run main.nf -profile singularity[OPTIONS]${Colors.RESET}

${Colors.BOLD}${Colors.RED}📝 REQUIRED PARAMETERS:${Colors.RESET}
${Colors.YELLOW}--samplesheet${Colors.RESET} <file>     Path to sample sheet (TSV format)
${Colors.YELLOW}--outdir${Colors.RESET} <path>          Output directory for results

${Colors.BOLD}${Colors.PURPLE}⚙️  OPTIONAL PARAMETERS:${Colors.RESET}
${Colors.YELLOW}--help${Colors.RESET}                   Show this help message and exit
${Colors.YELLOW}--version${Colors.RESET}               Show pipeline version

${Colors.BOLD}${Colors.GREEN}🔧 PIPELINE CONFIGURATION:${Colors.RESET}
${Colors.YELLOW}--max_memory${Colors.RESET} <string>    Maximum memory (default: ${Colors.DIM}'480.GB'${Colors.RESET})
${Colors.YELLOW}--max_cpus${Colors.RESET} <int>         Maximum CPUs (default: ${Colors.DIM}20${Colors.RESET})
${Colors.YELLOW}--max_time${Colors.RESET} <string>      Maximum time (default: ${Colors.DIM}'24.h'${Colors.RESET})

${Colors.BOLD}${Colors.GREEN}📊 OUTPUT STRUCTURE:${Colors.RESET}
${Colors.WHITE}The pipeline generates a structured output directory:${Colors.RESET}

${Colors.CYAN}results/
├── ${Colors.YELLOW}01_individual_assemblies/${Colors.RESET}    ${Colors.DIM}# Individual SAG assemblies and quality assessment${Colors.RESET}
├── ${Colors.YELLOW}02_similarity_analysis/${Colors.RESET}      ${Colors.DIM}# MinHash similarity analysis results${Colors.RESET}
├── ${Colors.YELLOW}03_clustering_analysis/${Colors.RESET}      ${Colors.DIM}# Hierarchical clustering results${Colors.RESET}
├── ${Colors.YELLOW}04_co_assemblies/${Colors.RESET}           ${Colors.DIM}# Co-assembly results and optimization${Colors.RESET}
├── ${Colors.YELLOW}05_taxonomic_classification/${Colors.RESET} ${Colors.DIM}# GTDB-Tk taxonomic classification${Colors.RESET}
└── ${Colors.YELLOW}06_final_results/${Colors.RESET}           ${Colors.DIM}# Final integrated results and reports${Colors.RESET}

${Colors.BOLD}${Colors.GREEN}💡 EXAMPLE USAGE:${Colors.RESET}

${Colors.DIM}# Basic run${Colors.RESET}
${Colors.CYAN}nextflow run main.nf -profile singularity -c nextflow.config${Colors.RESET}

${Colors.DIM}# Advanced run with custom parameters${Colors.RESET}
${Colors.CYAN}nextflow run main.nf -profile singularity --samplesheet samples.tsv --outdir results --max_cpus 32${Colors.RESET}

${Colors.BOLD}${Colors.YELLOW}📋 SAMPLE SHEET FORMAT:${Colors.RESET}
${Colors.WHITE}Tab-separated file with columns: sampleID, forwardReads, reverseReads${Colors.RESET}

${Colors.BOLD}${Colors.RED}🆘 SUPPORT:${Colors.RESET}
${Colors.CYAN}• GitHub Issues: https://github.com/your-repo/CoSAG-nf/issues${Colors.RESET}
${Colors.CYAN}• Documentation: https://your-docs.com/CoSAG-nf${Colors.RESET}

${Colors.CYAN}╚══════════════════════════════════════════════════════════════════════════════╝${Colors.RESET}
    """.stripIndent()
}

def versionMessage() {
    log.info"""
${Colors.CYAN}
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   ${Colors.BOLD}${Colors.WHITE}██╗   ██╗███████╗██████╗ ███████╗██╗ ██████╗ ███╗   ██╗           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}██║   ██║██╔════╝██╔══██╗██╔════╝██║██╔═══██╗████╗  ██║           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}██║   ██║█████╗  ██████╔╝███████╗██║██║   ██║██╔██╗ ██║           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}╚██╗ ██╔╝██╔══╝  ██╔══██╗╚════██║██║██║   ██║██║╚██╗██║           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE} ╚████╔╝ ███████╗██║  ██║███████║██║╚██████╔╝██║ ╚████║           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}  ╚═══╝  ╚══════╝╚═╝  ╚═╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝           ${Colors.RESET}${Colors.CYAN}   ║
║                                                                              ║
║                        ${Colors.BOLD}${Colors.YELLOW}CoSAG-nf Pipeline Information${Colors.RESET}${Colors.CYAN}                        ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
${Colors.RESET}

${Colors.BOLD}${Colors.GREEN}🔖 PIPELINE VERSION:${Colors.RESET}
${Colors.WHITE}CoSAG-nf${Colors.RESET}                 ${Colors.BOLD}${Colors.YELLOW}v1.0.0${Colors.RESET}

${Colors.BOLD}${Colors.BLUE}⚙️  SYSTEM INFORMATION:${Colors.RESET}
${Colors.WHITE}Nextflow Version${Colors.RESET}        ${Colors.CYAN}${workflow.nextflow.version}${Colors.RESET}
${Colors.WHITE}Container Engine${Colors.RESET}        ${Colors.CYAN}${workflow.containerEngine ?: 'Not specified'}${Colors.RESET}
${Colors.WHITE}Profile${Colors.RESET}                 ${Colors.CYAN}${workflow.profile ?: 'default'}${Colors.RESET}

${Colors.BOLD}${Colors.PURPLE}📁 DIRECTORY INFORMATION:${Colors.RESET}
${Colors.WHITE}Project Directory${Colors.RESET}       ${Colors.DIM}${workflow.projectDir}${Colors.RESET}
${Colors.WHITE}Launch Directory${Colors.RESET}        ${Colors.DIM}${workflow.launchDir}${Colors.RESET}
${Colors.WHITE}Work Directory${Colors.RESET}          ${Colors.DIM}${workflow.workDir}${Colors.RESET}

${Colors.CYAN}╚══════════════════════════════════════════════════════════════════════════════╝${Colors.RESET}
    """.stripIndent()
}
// Initialize default parameter values
params.help = params.help ?: false
params.version = params.version ?: false

// Check for help parameter
if (params.help) {
    helpMessage()
    exit 0
}

// Check for version parameter  
if (params.version) {
    versionMessage()
    exit 0
}
include {SPADES} from "./modules/spades"
include {CHECKM2} from "./modules/checkm2"
include {GTDBTK_CLASSIFYWF} from "./modules/gtdb_classify"
include {GTDBTK_CLASSIFYWF_COSAG} from "./modules/gtdb_classify_cosag"

include {SOURMASH_SKETCH_GENERATION;SOURMASH_COMPARE_MATRIX;SOURMASH_MATRIX_CONVERT;SOURMASH_PROCESS_MATRIX;SOURMASH_QUALITY_FILTER} from "./modules/sourmash"
// Hierarchical clustering modules for SAG binning
include { HIERARCHICAL_CLUSTERING_MATRIX } from './modules/hierarchical_clustering'
include { HIERARCHICAL_CLUSTERING_ANALYSIS } from './modules/hierarchical_clustering'
include { HIERARCHICAL_CLUSTERING_REPORT } from './modules/hierarchical_clustering'
include {PREPARE_CLUSTER_JSON } from './modules/prepare'
include { SPADES_CHECKM2 } from './subworkflows/spades_checkm2'
include { ITERATIVE_OPTIMIZATION } from './subworkflows/iterative_optimization'
include { DYNAMIC_SAG_OPTIMIZATION } from './subworkflows/dynamic_sag_optimization'
include { COSAG_OPTIMIZATION;COLLECT_ALL_ASSEMBLIES } from './modules/sag_opt'
include { UPDATE_OPTIMIZED_JSON;UPDATE_GTDB_JSON;UPDATE_SAG_JSON; FORMAT_REPORTS } from './modules/update_json'

// 工作流开始日志
log.info """
${Colors.CYAN}
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   ${Colors.BOLD}${Colors.WHITE}██████╗ ██╗██████╗ ███████╗██╗     ██╗███╗   ██╗███████╗           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}██╔══██╗██║██╔══██╗██╔════╝██║     ██║████╗  ██║██╔════╝           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}██████╔╝██║██████╔╝█████╗  ██║     ██║██╔██╗ ██║█████╗             ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}██╔═══╝ ██║██╔═══╝ ██╔══╝  ██║     ██║██║╚██╗██║██╔══╝             ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}██║     ██║██║     ███████╗███████╗██║██║ ╚████║███████╗           ${Colors.RESET}${Colors.CYAN}   ║
║   ${Colors.BOLD}${Colors.WHITE}╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝╚═╝╚═╝  ╚═══╝╚══════╝           ${Colors.RESET}${Colors.CYAN}   ║
║                                                                              ║
║                    ${Colors.BOLD}${Colors.YELLOW}CoSAG-nf v1.0.0 - Pipeline Execution${Colors.RESET}${Colors.CYAN}                    ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
${Colors.RESET}

${Colors.BOLD}${Colors.GREEN}🧬 CoSAG-nf: A Scalable Nextflow Pipeline for Co-assembly, 
   Optimization, and Interactive Visualization of 
   High-Throughput Single-Cell Genomes${Colors.RESET}

${Colors.BOLD}${Colors.BLUE}📋 INPUT PARAMETERS:${Colors.RESET}
${Colors.WHITE}Sample sheet${Colors.RESET}     ${Colors.CYAN}${params.samplesheet}${Colors.RESET}
${Colors.WHITE}Output directory${Colors.RESET} ${Colors.CYAN}${params.outdir}${Colors.RESET}
${Colors.WHITE}Max CPUs${Colors.RESET}         ${Colors.YELLOW}${params.max_cpus}${Colors.RESET}
${Colors.WHITE}Max Memory${Colors.RESET}       ${Colors.YELLOW}${params.max_memory}${Colors.RESET}

${Colors.BOLD}${Colors.PURPLE}⚙️  PIPELINE CONFIGURATION:${Colors.RESET}
${Colors.WHITE}Profile${Colors.RESET}          ${Colors.CYAN}${workflow.profile ?: 'default'}${Colors.RESET}
${Colors.WHITE}Nextflow version${Colors.RESET} ${Colors.CYAN}${workflow.nextflow.version}${Colors.RESET}

${Colors.BOLD}${Colors.RED}🚀 Starting pipeline execution...${Colors.RESET}

${Colors.CYAN}╚══════════════════════════════════════════════════════════════════════════════╝${Colors.RESET}
"""


workflow {

    log.info "Setting up input data channel..."
    read_pairs_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def meta = [id: row.sampleID]
            def forward_reads = file(row.forwardReads, checkIfExists: true)
            def reverse_reads = file(row.reverseReads, checkIfExists: true)
            tuple(meta, forward_reads, reverse_reads) 
        }
    
 
   // log.info "Starting SPADES assembly..."
    SPADES(read_pairs_ch) 

    gtdb_input_ch = SPADES.out.contigs
    .map { meta, contig_file -> contig_file }  //
    .collect()  // 


    GTDBTK_CLASSIFYWF(gtdb_input_ch)

    //GTDBTK_CLASSIFYWF.out.sagbac_summary.view()

    sag_contigs_mapping = SPADES.out.contigs
        .map { meta, contigs -> 
            // SAG_ID\tcontigs_path
            tuple(meta.id, contigs)
        }
        .collectFile(name: 'sag_contigs_mapping.tsv', newLine: true) { sag_id, contigs_path ->
            "${sag_id}\t${contigs_path}"
        }
    
    CHECKM2(SPADES.out.contigs)
    
    // MinHash signature
    //log.info "Starting MinHash signature generation..."
    SOURMASH_SKETCH_GENERATION(
        SPADES.out.contigs,
        params.sourmash?.ksize ?: 51,
        params.sourmash?.scaled ?: 1000
    )
    // Collect all signatures and compute similarity matrix
    signatures_ch = SOURMASH_SKETCH_GENERATION.out.signature
        .map { meta, signature -> signature }  // Extract signature files
        .collect()  // Collect all signatures
    
    SOURMASH_COMPARE_MATRIX(
        signatures_ch,
        params.sourmash?.ksize ?: 51,
        params.sourmash?.scaled ?: 1000
    )
    
    // Convert .npy matrix to labeled CSV format
    SOURMASH_MATRIX_CONVERT(
        SOURMASH_COMPARE_MATRIX.out.similarity_matrix,
        SOURMASH_COMPARE_MATRIX.out.labels,
        signatures_ch,
        params.sourmash?.ksize ?: 51,
        params.sourmash?.scaled ?: 1000
    )
    
    // Process the raw similarity matrix with Python
    SOURMASH_PROCESS_MATRIX(
        SOURMASH_MATRIX_CONVERT.out.raw_matrix,
        SOURMASH_MATRIX_CONVERT.out.log,
        params.sourmash?.ksize ?: 51,
        params.sourmash?.scaled ?: 1000
    )
// Apply quality filtering based on connectivity (optional)
    def enable_quality_filter = params.sourmash?.enable_quality_filter ?: true
    
    if (enable_quality_filter) {
        SOURMASH_QUALITY_FILTER(
            SOURMASH_PROCESS_MATRIX.out.similarity_matrix,
            SOURMASH_PROCESS_MATRIX.out.comparison_details,
            params.sourmash?.min_similarity_threshold ?: 0.05,
            params.sourmash?.min_connections ?: 1
        )
        
        // Use filtered matrix for clustering
        similarity_matrix_for_clustering = SOURMASH_QUALITY_FILTER.out.filtered_matrix
    } else {
        // Use original matrix for clustering
        similarity_matrix_for_clustering = SOURMASH_PROCESS_MATRIX.out.similarity_matrix
    }

    HIERARCHICAL_CLUSTERING_MATRIX(
        similarity_matrix_for_clustering,
        params.hierarchical_clustering?.distance_metric ?: 'euclidean'
    )
    
    HIERARCHICAL_CLUSTERING_ANALYSIS(
        HIERARCHICAL_CLUSTERING_MATRIX.out.distance_matrix,
        params.hierarchical_clustering?.linkage_method ?: 'complete',
        params.hierarchical_clustering?.criterion ?: 'inconsistent',
        params.hierarchical_clustering?.threshold ?: 0.95
    )
    
    HIERARCHICAL_CLUSTERING_REPORT(
        HIERARCHICAL_CLUSTERING_ANALYSIS.out.clusters,
        similarity_matrix_for_clustering,
        HIERARCHICAL_CLUSTERING_ANALYSIS.out.report,
        HIERARCHICAL_CLUSTERING_ANALYSIS.out.validation
    )

     // 
    PREPARE_CLUSTER_JSON(
        HIERARCHICAL_CLUSTERING_ANALYSIS.out.clusters,
        params.samplesheet,
        sag_contigs_mapping 
    )

    individual_checkm2_results = CHECKM2.out.results
        .map { meta, tsv -> tsv }
        .collect()

    SPADES_CHECKM2(
        PREPARE_CLUSTER_JSON.out.cluster_json,
        individual_checkm2_results
    )

  
    sag_input = SPADES_CHECKM2.out.filtered_json
        .flatten()  // 
    
    COSAG_OPTIMIZATION(
        sag_input
    )

    COLLECT_ALL_ASSEMBLIES(
    COSAG_OPTIMIZATION.out.bset_cosag.collect().ifEmpty([]),
    SPADES_CHECKM2.out.updated_json.collect(),
    COSAG_OPTIMIZATION.out.output_dir.collect()
)
    GTDBTK_CLASSIFYWF_COSAG(COLLECT_ALL_ASSEMBLIES.out.all_bins)
    
    UPDATE_OPTIMIZED_JSON(SPADES_CHECKM2.out.updated_json,COSAG_OPTIMIZATION.out.optimized_json.collect())
    UPDATE_GTDB_JSON(UPDATE_OPTIMIZED_JSON.out.update_optimized_json,GTDBTK_CLASSIFYWF_COSAG.out.bac_summary)
    UPDATE_SAG_JSON(GTDBTK_CLASSIFYWF.out.sagbac_summary,UPDATE_GTDB_JSON.out.update_gtdb_json)
    FORMAT_REPORTS(UPDATE_SAG_JSON.out.update_sag_json)

}


// 工作流完成时的日志和汇总报告生成
workflow.onComplete {
    log.info """
${Colors.CYAN}
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   ${Colors.BOLD}${Colors.WHITE}██████╗ ██████╗ ███╗   ███╗██████╗ ██╗     ███████╗████████╗███████╗${Colors.RESET}${Colors.CYAN} ║
║   ${Colors.BOLD}${Colors.WHITE}██╔════╝██╔═══██╗████╗ ████║██╔══██╗██║     ██╔════╝╚══██╔══╝██╔════╝${Colors.RESET}${Colors.CYAN} ║
║   ${Colors.BOLD}${Colors.WHITE}██║     ██║   ██║██╔████╔██║██████╔╝██║     █████╗     ██║   █████╗  ${Colors.RESET}${Colors.CYAN} ║
║   ${Colors.BOLD}${Colors.WHITE}██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ██║     ██╔══╝     ██║   ██╔══╝  ${Colors.RESET}${Colors.CYAN} ║
║   ${Colors.BOLD}${Colors.WHITE}╚██████╗╚██████╔╝██║ ╚═╝ ██║██║     ███████╗███████╗   ██║   ███████╗${Colors.RESET}${Colors.CYAN} ║
║   ${Colors.BOLD}${Colors.WHITE} ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝   ╚═╝   ╚══════╝${Colors.RESET}${Colors.CYAN} ║
║                                                                              ║
║                        ${Colors.BOLD}${Colors.YELLOW}CoSAG-nf Pipeline Completed!${Colors.RESET}${Colors.CYAN}                         ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
${Colors.RESET}

${Colors.BOLD}${Colors.BLUE}📊 EXECUTION SUMMARY:${Colors.RESET}
${Colors.WHITE}Pipeline${Colors.RESET}        ${Colors.BOLD}${Colors.CYAN}CoSAG-nf v1.0.0${Colors.RESET}
${Colors.WHITE}Status${Colors.RESET}          ${workflow.success ? "${Colors.BOLD}${Colors.GREEN}✅ SUCCESS${Colors.RESET}" : "${Colors.BOLD}${Colors.RED}❌ FAILED${Colors.RESET}"}
${Colors.WHITE}Duration${Colors.RESET}        ${Colors.YELLOW}${workflow.duration}${Colors.RESET}
${Colors.WHITE}Exit status${Colors.RESET}     ${workflow.success ? "${Colors.GREEN}${workflow.exitStatus}${Colors.RESET}" : "${Colors.RED}${workflow.exitStatus}${Colors.RESET}"}

${Colors.BOLD}${Colors.PURPLE}📁 DIRECTORIES:${Colors.RESET}
${Colors.WHITE}Work directory${Colors.RESET}  ${Colors.DIM}${workflow.workDir}${Colors.RESET}
${Colors.WHITE}Output directory${Colors.RESET} ${Colors.CYAN}${params.outdir}${Colors.RESET}

${workflow.success ? 
    "${Colors.BOLD}${Colors.GREEN}🎉 RESULTS READY!${Colors.RESET}\n${Colors.WHITE}Check the structured output directory for your analysis results.${Colors.RESET}\n\n${Colors.BOLD}${Colors.CYAN}📋 Key Results:${Colors.RESET}\n${Colors.GREEN}•${Colors.RESET} ${Colors.WHITE}Individual assemblies: 01_individual_assemblies/${Colors.RESET}\n${Colors.GREEN}•${Colors.RESET} ${Colors.WHITE}Final results: 06_final_results/${Colors.RESET}\n${Colors.GREEN}•${Colors.RESET} ${Colors.WHITE}HTML reports: 06_final_results/html_reports/${Colors.RESET}" : 
    "${Colors.BOLD}${Colors.RED}💥 PIPELINE FAILED!${Colors.RESET}\n${Colors.WHITE}Please check the error logs and try again.${Colors.RESET}\n${Colors.RED}Error: ${workflow.errorReport ?: 'No specific error reported'}${Colors.RESET}"}

${Colors.CYAN}╚══════════════════════════════════════════════════════════════════════════════╝${Colors.RESET}
    """.stripIndent()
    
    // 生成汇总报告
    if (workflow.success && params.generate_summary) {
        log.info "Generating pipeline summary reports..."
        
        def summary_script = "${projectDir}/bin/generate_summary_reports.py"
        def summary_cmd = "python3 ${summary_script} --outdir ${params.outdir} --output 06_final_results/summary_tables"
        
        try {
            def proc = summary_cmd.execute()
            proc.waitFor()
            
            if (proc.exitValue() == 0) {
                log.info "Summary reports generated successfully!"
            } else {
                log.warn "Failed to generate summary reports: ${proc.err.text}"
            }
        } catch (Exception e) {
            log.warn "Error generating summary reports: ${e.message}"
        }
    }
    
    // 输出结果目录结构信息
    log.info """
${Colors.CYAN}╔══════════════════════════════════════════════════════════════════════════════╗${Colors.RESET}
${Colors.CYAN}║${Colors.RESET}                        ${Colors.BOLD}${Colors.YELLOW}Results Directory Structure${Colors.RESET}                         ${Colors.CYAN}║${Colors.RESET}
${Colors.CYAN}╚══════════════════════════════════════════════════════════════════════════════╝${Colors.RESET}
${Colors.WHITE}Main output: ${Colors.CYAN}${params.outdir}${Colors.RESET}

${Colors.BOLD}${Colors.GREEN}Key directories:${Colors.RESET}
${Colors.GREEN}📁${Colors.RESET} ${Colors.YELLOW}01_individual_assemblies/${Colors.RESET}  ${Colors.DIM}- Individual SAG assemblies & quality${Colors.RESET}
${Colors.GREEN}📁${Colors.RESET} ${Colors.YELLOW}02_similarity_analysis/${Colors.RESET}    ${Colors.DIM}- MinHash similarity matrices${Colors.RESET}
${Colors.GREEN}📁${Colors.RESET} ${Colors.YELLOW}03_clustering_analysis/${Colors.RESET}    ${Colors.DIM}- Hierarchical clustering results${Colors.RESET}
${Colors.GREEN}📁${Colors.RESET} ${Colors.YELLOW}04_co_assemblies/${Colors.RESET}         ${Colors.DIM}- Co-assembly results per cluster${Colors.RESET}
${Colors.GREEN}📁${Colors.RESET} ${Colors.YELLOW}05_taxonomic_classification/${Colors.RESET} ${Colors.DIM}- GTDB-Tk classification${Colors.RESET}
${Colors.GREEN}📁${Colors.RESET} ${Colors.YELLOW}06_final_results/${Colors.RESET}         ${Colors.DIM}- Integrated results & reports${Colors.RESET}

${Colors.BOLD}${Colors.CYAN}Key files to check:${Colors.RESET}
${Colors.GREEN}📄${Colors.RESET} ${Colors.WHITE}06_final_results/json_reports/final_result.json${Colors.RESET}
${Colors.GREEN}📄${Colors.RESET} ${Colors.WHITE}06_final_results/html_reports/embedded_report.html${Colors.RESET}
${Colors.GREEN}📄${Colors.RESET} ${Colors.WHITE}06_final_results/summary_tables/pipeline_summary_report.txt${Colors.RESET}
${Colors.CYAN}╚══════════════════════════════════════════════════════════════════════════════╝${Colors.RESET}
    """.stripIndent()
}