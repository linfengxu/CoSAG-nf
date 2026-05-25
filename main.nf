#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    linfengxu/CoSAG-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/linfengxu/CoSAG-nf
    Docs   : https://github.com/linfengxu/CoSAG-nf/tree/main/docs
----------------------------------------------------------------------------------------
*/

include { COSAG                   } from './workflows/cosag'
include { COSAG_ROUND2            } from './workflows/cosag_round2'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_cosag_pipeline/main'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_cosag_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COSAG_NF {

    take:
    samplesheet

    main:
    COSAG ( samplesheet )

    emit:
    final_json  = COSAG.out.final_json
    updated_json = COSAG.out.updated_json
    sag_bac_summary = COSAG.out.sag_bac_summary
    sag_ar_summary  = COSAG.out.sag_ar_summary
    bac_summary  = COSAG.out.bac_summary
    ar_summary   = COSAG.out.ar_summary
    versions    = COSAG.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    def run_round = (params.round ?: 1) as Integer
    def base_outdir = params.outdir

    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    COSAG_NF (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    if (run_round == 2) {
        if (!base_outdir.endsWith('/round2')) {
            params.outdir = "${base_outdir}/round2"
        }
        COSAG_ROUND2(
            COSAG_NF.out.updated_json,
            COSAG_NF.out.sag_bac_summary,
            COSAG_NF.out.sag_ar_summary,
        )
    }

    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        [],
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
