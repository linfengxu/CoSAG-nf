process COLLECT_ASSEMBLIES {
    tag "collect_assemblies"
    label 'process_low'

   // container 'python:3.11-slim'
    container 'quay.io/xulf2022/python3.8_bio:v1'


    input:
    path(optimized_jsons)     // all cluster_*_optimized.json collected
    path(coassembly_json)     // cluster_data_updated.json

    output:
    path("collection/bset_cosag/*.fasta"),         emit: bset_fastas,      optional: true
    path("collection/co_assembly/*.fasta"),         emit: coasm_fastas,     optional: true
    path("collection/all_accepted_cosags/*.fasta"), emit: all_bins
    path("collection/filtered_out/*.fasta"),        emit: filtered_out,     optional: true
    path("collection/collection_summary.txt"),      emit: summary
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def ca_arg    = coassembly_json.name != 'NO_FILE'
                    ? "--coassembly_json ${coassembly_json}" : ''
    """
    collect_assemblies.py \\
        --optimized_jsons ${optimized_jsons} \\
        ${ca_arg} \\
        --out_dir           collection \\
        --min_completeness  ${params.min_completeness} \\
        --max_contamination ${params.max_contamination} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
