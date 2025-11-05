process UPDATE_OPTIMIZED_JSON {
    tag "update_optimized_json"
    label 'process_high'
    publishDir "${params.output_structure?.final_results ?: params.outdir}/json_reports", mode: 'copy'
    container "${params.python.container}"

    input:
    path updated_json                    // cluster_data_updated.json
    path optimized_jsons                 // cluster_*_optimized.json 
    
    output:
    path "output/update_optimized.json", emit: update_optimized_json
    
    script:
    """
    mkdir -p output
    

    python3 ${projectDir}/bin/update_cluster_optimized2.py \\
        ${updated_json} \\
        ${optimized_jsons} \\
        -o update_optimized.json
    
    mv update_optimized.json output/ 2>/dev/null || true
    """
}

process UPDATE_GTDB_JSON {
    tag "update_GTDB_json"
    label 'process_high'
    publishDir "${params.output_structure?.final_results ?: params.outdir}/json_reports", mode: 'copy'
    container "${params.python.container}"

    input:
    path updated_json                    // cluster_data_updated.json
    path gtdb_jsons                 // cluster_*_optimized.json 
    
    output:
    path "output/update_GTDB.json", emit: update_gtdb_json
    
    script:
    """
    mkdir -p output
    

    python3 ${projectDir}/bin/update_gtdbtk_classification2.py \\
        ${updated_json} \\
        ${gtdb_jsons} \\
        -o update_GTDB.json
    
    mv update_GTDB.json output/ 2>/dev/null || true
    """
}

process UPDATE_SAG_JSON {
    tag "update_SAG_json"
    label 'process_high'
    publishDir "${params.output_structure?.final_results ?: params.outdir}/json_reports", mode: 'copy'
    container "${params.python.container}"

    input:
    path updated_json                    // cluster_data_updated.json
    path gtdb_jsons                 // cluster_*_optimized.json 
    
    output:
    path "output/final_result.json", emit: update_sag_json
    
    script:
    """
    mkdir -p output
    

    python3 ${projectDir}/bin/merge_taxonomy.py \\
        ${updated_json} \\
        ${gtdb_jsons} \\
        final_result.json
    
    mv final_result.json output/ 2>/dev/null || true
    """
}

process FORMAT_REPORTS {
    tag "format_reports"
    label 'process_high'
    publishDir "${params.output_structure?.final_results ?: params.outdir}/html_reports", mode: 'copy'
    container "${params.python.container}"

    input:
    path updated_json                    // cluster_data_updated.json
    
    output:
    path "output/embedded_report.html", emit: update_sag_json
    
    script:
    """
    mkdir -p output
    

    python3 ${projectDir}/bin/foramt_report.py \\
        ${updated_json} \\
    
    mv embedded_report.html output/ 2>/dev/null || true
    """
}