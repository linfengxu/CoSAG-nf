process ASSESS_HQ_MAG_FINAL {
    tag 'hq_mag_final'
    label 'process_low'

    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    tuple path(cluster_json), path(barrnap_gffs)

    output:
    path('cluster_data_gtdbtk.json'), emit: updated_json
    path('versions.yml'),            emit: versions

    when:
    params.hq_mag_final_assessment && (task.ext.when == null || task.ext.when)

    script:
    def gt = params.hq_mag_completeness_gt
    def lt = params.hq_mag_contamination_lt
    def barr = barrnap_gffs instanceof List ? barrnap_gffs : [barrnap_gffs]
    """
    mkdir -p _barr
    ${barr.collect { "cp '${it}' _barr/ 2>/dev/null || true" }.join('\n')}

    assess_hq_mag_final.py \\
        --json ${cluster_json} \\
        --barrnap_dir _barr \\
        --completeness_gt ${gt} \\
        --contamination_lt ${lt} \\
        --require_trna false \\
        --out_json cluster_data_gtdbtk.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
