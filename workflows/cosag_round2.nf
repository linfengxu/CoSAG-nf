#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    workflows/cosag_round2.nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Standalone round-2 workflow:
      Round1 updated_clusters.json
        -> build round2 contigs from clusters[*].coassembly_contigs
        -> MINHASH_CLUSTER (round=2)
        -> COASSEMBLY_ROUND (round=2)
        -> QUALITY_TAXONOMY
        -> final cluster_data_gtdbtk.json
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { MINHASH_CLUSTER } from '../subworkflows/local/minhash_cluster/main'
include { COASSEMBLY_ROUND } from '../subworkflows/local/coassembly_round/main'
include { QUALITY_TAXONOMY } from '../subworkflows/local/quality_taxonomy/main'
include { FORMAT_REPORT    } from '../modules/local/format_report/main'

workflow COSAG_ROUND2 {

    take:
    round1_updated_json_in
    sag_bac_summary_in
    sag_ar_summary_in

    main:
    ch_round1_updated_json = round1_updated_json_in
    ch_sag_bac_summary = sag_bac_summary_in
    ch_sag_ar_summary = sag_ar_summary_in
    def round2_max_contam = (params.max_contamination ?: 10) as Double

    // Build round2 MinHash input from round1 JSON coassembly_contigs.
    // Exclude high-contamination clusters from entering round2.
    ch_round2_contigs = ch_round1_updated_json
        .splitJson(path: 'clusters')
        .filter { c ->
            if (!c?.coassembly_contigs) {
                return false
            }
            def contam = c?.checkm2_coassembly?.contamination
            return contam == null || (contam as Double) <= round2_max_contam
        }
        .map { c ->
            def meta = [id: c.cluster_id, cluster_id: c.cluster_id, round: 2]
            tuple(meta, file(c.coassembly_contigs))
        }

    MINHASH_CLUSTER(
        ch_round2_contigs,
        2,
        [],
        [],
        ch_round1_updated_json,
        file("NO_FILE"),
    )

    COASSEMBLY_ROUND(
        MINHASH_CLUSTER.out.cluster_json,
        2,
    )

    QUALITY_TAXONOMY(
        COASSEMBLY_ROUND.out.updated_json,
        COASSEMBLY_ROUND.out.opt_jsons,
        ch_sag_bac_summary,
        ch_sag_ar_summary,
    )

    FORMAT_REPORT(
        QUALITY_TAXONOMY.out.final_json
    )

    emit:
    round2_cluster_json = MINHASH_CLUSTER.out.cluster_json
    round2_updated_json = COASSEMBLY_ROUND.out.updated_json
    final_json          = QUALITY_TAXONOMY.out.final_json
    report_html         = FORMAT_REPORT.out.report
    bac_summary         = QUALITY_TAXONOMY.out.bac_summary
    ar_summary          = QUALITY_TAXONOMY.out.ar_summary
}

