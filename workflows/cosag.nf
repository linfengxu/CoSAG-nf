/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    linfengxu/cosag-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Single-cell amplified genome co-assembly and quality filtering pipeline
----------------------------------------------------------------------------------------
*/

include { SAG_ASSEMBLY     } from '../subworkflows/local/sag_assembly/main'
include { MINHASH_CLUSTER  } from '../subworkflows/local/minhash_cluster/main'
include { COASSEMBLY       } from '../subworkflows/local/coassembly/main'
include { QUALITY_TAXONOMY } from '../subworkflows/local/quality_taxonomy/main'
include { FORMAT_REPORT    } from '../modules/local/format_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COSAG {

    take:
    samplesheet   // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ] ]

    main:
    ch_versions = Channel.empty()

    /*
    ============================================================
        MODULE 1 — Individual SAG Assembly & QC
    ============================================================
        Input  : samplesheet [ meta, [r1, r2] ]
        Output : contigs     [ meta, fasta ]
                 checkm2     path(quality_report.tsv)
                 bac_summary path(gtdbtk.bac120.summary.tsv)
                 ar_summary  path(gtdbtk.ar53.summary.tsv)
    */
    ch_samplesheet_path = Channel.fromPath(params.input, checkIfExists: true)

    SAG_ASSEMBLY ( samplesheet )
    ch_versions = ch_versions.mix( SAG_ASSEMBLY.out.versions )

    // Materialize summary channels once, then duplicate with mix() to avoid
    // fragile channel reuse patterns that can trigger runtime recursion.
    ch_sag_bac_base = SAG_ASSEMBLY.out.bac_summary.map { meta, summary -> summary }
    ch_sag_ar_base  = SAG_ASSEMBLY.out.ar_summary.map  { meta, summary -> summary }
    ch_sag_bac_qt   = ch_sag_bac_base.mix(Channel.empty())
    ch_sag_bac_emit = ch_sag_bac_base.mix(Channel.empty())
    ch_sag_ar_qt    = ch_sag_ar_base.mix(Channel.empty())
    ch_sag_ar_emit  = ch_sag_ar_base.mix(Channel.empty())

    /*
    ============================================================
        MODULE 2 — MinHash Clustering
    ============================================================
        Input  : contigs     [ meta, fasta ]
                 samplesheet path
        Output : cluster_json path(cluster_data.json)
    */
    MINHASH_CLUSTER(
    SAG_ASSEMBLY.out.contigs,
    1,
    SAG_ASSEMBLY.out.pass_tsv,
    file(params.input),
    [],
    file("NO_FILE")
)
    ch_versions = ch_versions.mix( MINHASH_CLUSTER.out.versions )

    /*
    ============================================================
        MODULE 3 — Co-Assembly, TNF Optimization & JSON Update
    ============================================================
        Input  : cluster_json       path(cluster_data.json)
                 individual_checkm2 path(quality_report.tsv)
        Output : cosag_contigs  [ meta, fasta ]
                 updated_json   path(cluster_data_updated.json)
                 opt_jsons      [ path(*_tnf_optimized.json) ] collected
    */
    COASSEMBLY (
        MINHASH_CLUSTER.out.cluster_json
    )
    ch_versions = ch_versions.mix( COASSEMBLY.out.versions )

    /*
    ============================================================
        MODULE 4 — CoSAG Taxonomic Classification & Final JSON
    ============================================================
        Input  : updated_json  path(cluster_data_updated.json)
                 opt_jsons     [ path(*_tnf_optimized.json) ] collected
        Output : final_json    path(cluster_data_gtdbtk.json)
    */
    QUALITY_TAXONOMY (
        COASSEMBLY.out.updated_json,
        COASSEMBLY.out.opt_jsons,
        ch_sag_bac_qt,
        ch_sag_ar_qt,
    )
    ch_versions = ch_versions.mix( QUALITY_TAXONOMY.out.versions )

    FORMAT_REPORT (
        QUALITY_TAXONOMY.out.final_json
    )
    ch_versions = ch_versions.mix( FORMAT_REPORT.out.versions )

    emit:
    final_json   = QUALITY_TAXONOMY.out.final_json
    report_html  = FORMAT_REPORT.out.report
    updated_json = COASSEMBLY.out.updated_json
    sag_bac_summary = ch_sag_bac_emit
    sag_ar_summary  = ch_sag_ar_emit
    bac_summary  = QUALITY_TAXONOMY.out.bac_summary
    ar_summary   = QUALITY_TAXONOMY.out.ar_summary
    versions     = ch_versions
}