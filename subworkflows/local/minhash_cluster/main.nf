/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/minhash_cluster/main.nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Round 1: contigs 来自 SAG_ASSEMBLY，最后调 PREPARE_CLUSTER_JSON
    Round 2+: contigs 来自上一轮 coassembly final_contigs，最后调 PREPARE_ROUND_JSON
              checkm2_pass / samplesheet 传 []（不使用）
----------------------------------------------------------------------------------------
*/

include { SOURMASH_SKETCH                  } from '../../../modules/local/sourmash/sketch/main'
include { SOURMASH_COMPARE                 } from '../../../modules/local/sourmash/compare/main'
include { SOURMASH_PROCESS_MATRIX          } from '../../../modules/local/sourmash/process_matrix/main'
include { SOURMASH_MATRIX_CONVERT          } from '../../../modules/local/sourmash/matrix_convert/main'
include { SOURMASH_QUALITY_FILTER          } from '../../../modules/local/sourmash/quality_filter/main'
include { HIERARCHICAL_CLUSTER_MATRIX      } from '../../../modules/local/hierarchical_cluster/matrix/main'
include { HIERARCHICAL_CLUSTERING_ANALYSIS } from '../../../modules/local/hierarchical_cluster/analysis/main'
include { HIERARCHICAL_CLUSTER_REPORT      } from '../../../modules/local/hierarchical_cluster/report/main'
include { PREPARE_CLUSTER_JSON             } from '../../../modules/local/prepare_cluster_json/main'
include { PREPARE_ROUND_JSON               } from '../../../modules/local/prepare_round_json/main'


workflow MINHASH_CLUSTER {

    take:
    contigs       // channel: [ val(meta), path(fasta) ]
                  //   Round 1: SAG 个体 contigs（来自 SAG_ASSEMBLY）
                  //   Round 2+: coassembly contigs（来自 COASSEMBLY_ROUND.out.final_contigs）
    round         // val: 当前轮次编号（1, 2, 3...）
    // ── Round 1 专用（Round 2+ 传 [] 占位）──
    checkm2_pass  // path: sags_pass.tsv
    samplesheet   // path: paired_end_qc_samples.tsv
    // ── Round 2+ 专用（Round 1 传 [] 占位）──
    prev_json     // path: 上一轮 updated_clusters.json
    // ── Round 3+ 专用（Round 1/2 传 NO_FILE 占位）──
    round1_json   // path: Round 1 updated JSON，供 PREPARE_ROUND_JSON 溯源 SAG reads

    main:
    ch_versions = Channel.empty()

    // ── contig 路径映射（两种轮次格式一致：meta.id \t fasta 路径）──────────
    // Round 1: meta.id = SAG_ID      → sag_contigs_mapping.tsv
    // Round 2+: meta.id = cluster_id → cluster_contigs_mapping.tsv
    contigs
        | map { meta, fasta -> "${meta.id}\t${fasta.toRealPath()}" }
        | collectFile(
            name:    "round${round}_contigs_mapping.tsv",
            newLine: true,
            sort:    true,
        )
        | set { ch_contig_mapping }

    // ── Step 1: MinHash sketch ────────────────────────────────────────────
    SOURMASH_SKETCH ( contigs )
    ch_versions = ch_versions.mix( SOURMASH_SKETCH.out.versions )

    SOURMASH_SKETCH.out.signature
        | map { meta, sig -> sig }
        | collect
        | set { ch_signatures }

    // ── Step 2: 全局相似度矩阵 ───────────────────────────────────────────
    SOURMASH_COMPARE(
        ch_signatures,
        params.sourmash_ksize  ?: 51,
        params.sourmash_scaled ?: 1000,
    )
    ch_versions = ch_versions.mix( SOURMASH_COMPARE.out.versions )

    // ── Step 3: 矩阵转换 ─────────────────────────────────────────────────
    SOURMASH_MATRIX_CONVERT(
        SOURMASH_COMPARE.out.similarity_matrix,
        SOURMASH_COMPARE.out.labels,
        ch_signatures,
        params.sourmash_ksize  ?: 51,
        params.sourmash_scaled ?: 1000,
    )

    // ── Step 4: 矩阵处理 ─────────────────────────────────────────────────
    SOURMASH_PROCESS_MATRIX(
        SOURMASH_MATRIX_CONVERT.out.raw_matrix,
        SOURMASH_MATRIX_CONVERT.out.log,
        params.sourmash_ksize  ?: 51,
        params.sourmash_scaled ?: 1000,
    )
    ch_versions = ch_versions.mix( SOURMASH_PROCESS_MATRIX.out.versions )

    // ── Step 5: 可选质量过滤 ─────────────────────────────────────────────
    if ( params.sourmash_enable_quality_filter ) {
        SOURMASH_QUALITY_FILTER ( SOURMASH_PROCESS_MATRIX.out.similarity_matrix )
        ch_versions      = ch_versions.mix( SOURMASH_QUALITY_FILTER.out.versions )
        ch_sim_matrix    = SOURMASH_QUALITY_FILTER.out.filtered_matrix
    } else {
        ch_sim_matrix    = SOURMASH_PROCESS_MATRIX.out.similarity_matrix
    }

    // ── Step 6: 相似度矩阵 → 距离矩阵 ───────────────────────────────────
    HIERARCHICAL_CLUSTER_MATRIX(
        ch_sim_matrix,
        params.distance_metric ?: 'ani',
        params.sourmash_ksize
    )
    ch_versions = ch_versions.mix( HIERARCHICAL_CLUSTER_MATRIX.out.versions )

    // ── Step 7: 层次聚类 ─────────────────────────────────────────────────
    HIERARCHICAL_CLUSTERING_ANALYSIS(
        HIERARCHICAL_CLUSTER_MATRIX.out.distance_matrix,
        params.cluster_linkage_method ?: 'complete',
        params.cluster_criterion      ?: 'inconsistent',
        params.cluster_threshold      ?: 0.95,
    )
    ch_versions = ch_versions.mix( HIERARCHICAL_CLUSTERING_ANALYSIS.out.versions )

    // ── Step 8: 聚类报告 ─────────────────────────────────────────────────
    HIERARCHICAL_CLUSTER_REPORT(
        HIERARCHICAL_CLUSTERING_ANALYSIS.out.clusters,
        ch_sim_matrix,
        HIERARCHICAL_CLUSTERING_ANALYSIS.out.report,
        HIERARCHICAL_CLUSTERING_ANALYSIS.out.validation,
    )
    ch_versions = ch_versions.mix( HIERARCHICAL_CLUSTER_REPORT.out.versions )

    // ── Step 9: 生成 JSON（按轮次分支）──────────────────────────────────
    //
    //  Round 1 → PREPARE_CLUSTER_JSON
    //    输入: clusters_tsv / samplesheet / contig_mapping / checkm2_pass
    //    同时合并 reads，输出 round1_clusters.json
    //
    //  Round 2+ → PREPARE_ROUND_JSON
    //    输入: prev_json / clusters_tsv / round编号
    //    从 prev_json 溯源 leaf_sags → 合并原始 reads，输出 roundN_clusters.json

    if ( round == 1 ) {
        PREPARE_CLUSTER_JSON(
            HIERARCHICAL_CLUSTERING_ANALYSIS.out.clusters,
            samplesheet,
            ch_contig_mapping,
            checkm2_pass,
        )
        ch_versions  = ch_versions.mix( PREPARE_CLUSTER_JSON.out.versions )
        ch_out_json  = PREPARE_CLUSTER_JSON.out.json
    } else {
        PREPARE_ROUND_JSON(
            prev_json,
            HIERARCHICAL_CLUSTERING_ANALYSIS.out.clusters,
            round1_json,
            round,
        )
        ch_versions  = ch_versions.mix( PREPARE_ROUND_JSON.out.versions )
        ch_out_json  = PREPARE_ROUND_JSON.out.json
    }

    emit:
    cluster_json = ch_out_json    // → COASSEMBLY_ROUND
    versions     = ch_versions
}
