include { GTDBTK_CLASSIFYWF as GTDBTK_COSAG } from '../../../modules/local/gtdbtk_classifywf/main'
include { UPDATE_CLUSTER_OPTIMIZED           } from '../../../modules/local/update_cluster_optimized/main'
include { EXTRACT_COSAG_CONTIGS              } from '../../../modules/local/extract_cosag_contigs/main'
include { BARRNAP_COSAG                      } from '../../../modules/local/barrnap_cosag/main'
include { UPDATE_GTDBTK_RESULTS as UPDATE_GTDBTK_RESULTS_SAG   } from '../../../modules/local/update_gtdbtk_results/main'
include { UPDATE_GTDBTK_RESULTS as UPDATE_GTDBTK_RESULTS_COSAG } from '../../../modules/local/update_gtdbtk_results/main'
include { ASSESS_HQ_MAG_FINAL                          } from '../../../modules/local/assess_hq_mag_final/main'

workflow QUALITY_TAXONOMY {

    take:
    cluster_json    // path — cluster_data_updated.json
    opt_jsons       // path — cluster_*_tnf_optimized.json collected
    sag_bac_summary // path — gtdbtk.bac120.summary.tsv from all_sags
    sag_ar_summary  // path — gtdbtk.ar53.summary.tsv from all_sags

    main:
    ch_versions = Channel.empty()

    // ── Step 1: 更新 JSON — 写入 TNF 优化结果 ────────────────────────────
    UPDATE_CLUSTER_OPTIMIZED (
        cluster_json,
        opt_jsons.ifEmpty( file("NO_OPT_JSON") ),
    )
    ch_versions = ch_versions.mix( UPDATE_CLUSTER_OPTIMIZED.out.versions )

    // ── Step 2: 从 cluster_data_tnf.json 提取并重命名通过阈值的 contigs ──
    EXTRACT_COSAG_CONTIGS (
        UPDATE_CLUSTER_OPTIMIZED.out.updated_json
    )
    ch_versions = ch_versions.mix( EXTRACT_COSAG_CONTIGS.out.versions )

    // ── Step 2b: 对抽取的 CoSAG FASTA 跑 barrnap（rRNA）──────────────────────
    ch_cosag_rrna_in = EXTRACT_COSAG_CONTIGS.out.fastas
        .flatMap { files ->
            def lst = files instanceof List ? files : [files]
            lst.collect { f ->
                tuple([id: f.baseName.toString()], f)
            }
        }

    BARRNAP_COSAG ( ch_cosag_rrna_in )
    ch_versions = ch_versions.mix( BARRNAP_COSAG.out.versions )

    // ── Step 3: GTDB-Tk 分类（CoSAG）────────────────────────────────────
    EXTRACT_COSAG_CONTIGS.out.fastas
        | collect
        | filter { fastas -> fastas && fastas.size() > 0 }
        | map { fastas -> [ [id: 'all_cosags'], fastas ] }
        | GTDBTK_COSAG
    ch_versions = ch_versions.mix( GTDBTK_COSAG.out.versions )

    // ── Step 4a: 更新 JSON — 先写入 GTDB-Tk SAG(member) 分类结果 ─────────
    ch_sag_bac_summary = sag_bac_summary
        .map { v -> v instanceof List && v.size() > 1 ? v[1] : v }
        .ifEmpty( file("NO_BAC_SUMMARY") )
    ch_sag_ar_summary = sag_ar_summary
        .map { v -> v instanceof List && v.size() > 1 ? v[1] : v }
        .ifEmpty( file("NO_AR_SUMMARY") )

    UPDATE_GTDBTK_RESULTS_SAG (
        UPDATE_CLUSTER_OPTIMIZED.out.updated_json,
        ch_sag_bac_summary,
        ch_sag_ar_summary,
        "sag",
    )
    ch_versions = ch_versions.mix( UPDATE_GTDBTK_RESULTS_SAG.out.versions )

    // ── Step 4b: 更新 JSON — 再叠加 GTDB-Tk CoSAG 分类结果 ────────────────
    ch_cosag_bac_summary = GTDBTK_COSAG.out.bac_summary
        .map { meta, summary -> summary }
        .ifEmpty( file("NO_BAC_SUMMARY") )
    ch_cosag_ar_summary = GTDBTK_COSAG.out.ar_summary
        .map { meta, summary -> summary }
        .ifEmpty( file("NO_AR_SUMMARY") )

    UPDATE_GTDBTK_RESULTS_COSAG (
        UPDATE_GTDBTK_RESULTS_SAG.out.updated_json,
        ch_cosag_bac_summary,
        ch_cosag_ar_summary,
        "coassembly",
    )
    ch_versions = ch_versions.mix( UPDATE_GTDBTK_RESULTS_COSAG.out.versions )

    // ── Step 5（仅最终评估）：写入 hq_mag_assessment，不改上游抽取/GTDB 逻辑 ──
    // Dummy channels must live in workflow scope (not inside if) — NF avoids shadowing Channel in branch bodies.
    ch_dummy_barr_for_assess = Channel.of([[ file("${projectDir}/assets/hq_mag_dummy/cluster_DUMMY_barrnap.gff", checkIfExists: true) ]])
    ch_barr_for_assess = params.run_cosag_rrna_annotation ? BARRNAP_COSAG.out.gff.map { _meta, g -> g }.collect() : ch_dummy_barr_for_assess

    ch_final_json = UPDATE_GTDBTK_RESULTS_COSAG.out.updated_json
    if (params.hq_mag_final_assessment) {
        ch_assess_in = UPDATE_GTDBTK_RESULTS_COSAG.out.updated_json.merge(ch_barr_for_assess)
            .map { row ->
                def r = row instanceof List ? row : [row]
                tuple(r[0], r[1])
            }
        ASSESS_HQ_MAG_FINAL(ch_assess_in)
        ch_versions = ch_versions.mix(ASSESS_HQ_MAG_FINAL.out.versions)
        ch_final_json = ASSESS_HQ_MAG_FINAL.out.updated_json
    }

    emit:
    final_json   = ch_final_json
    bac_summary  = GTDBTK_COSAG.out.bac_summary
    ar_summary   = GTDBTK_COSAG.out.ar_summary
    versions     = ch_versions
}
