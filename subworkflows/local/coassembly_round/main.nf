/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/coassembly_round/main.nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    单轮共组装（Round N）:
        cluster_json (roundN_clusters.json，含 merged_reads 完整路径)
            → SPADES_COASSEMBLY
            → CHECKM2
            → UPDATE_CLUSTER_JSON  (填 coassembly_contigs + checkm2_coassembly)
            → branch: needs_opt / direct
            → COSAG_OPTIMIZATION   (仅满足条件的 cluster，直接产出优化 contig)
            → emit: final_contigs, updated_json, opt_jsons

    触发 TNF 优化的条件（均满足时才优化）:
        params.cosag_optimize               = true
        checkm2_coassembly.contamination    > params.cosag_opt_max_contamination (默认 10.0)
        checkm2_coassembly.completeness     > params.cosag_opt_min_completeness  (默认 90.0)
----------------------------------------------------------------------------------------
*/

include { SPADES as SPADES_COASSEMBLY } from '../../../modules/local/spades/main'
include { CHECKM2                     } from '../../../modules/local/checkm2/main'
include { UPDATE_CLUSTER_JSON         } from '../../../modules/local/update_cluster_json/main'
include { COSAG_OPTIMIZATION          } from '../../../modules/local/cosag_optimization/main'


workflow COASSEMBLY_ROUND {

    take:
    cluster_json    // path: roundN_clusters.json
    round           // val:  轮次编号

    main:
    ch_versions = Channel.empty()

    def do_optimize = params.cosag_optimize               ?: false
    def min_comp    = (params.cosag_opt_min_completeness  ?: 90.0) as Double
    def max_contam  = (params.cosag_opt_max_contamination ?: 10.0) as Double

    /*
    ── Step 1: 解析 JSON → [meta, [r1, r2]] ──
    */
    cluster_json
        | splitJson( path: 'clusters' )
        | map { cluster ->
            def meta = [
                id:         cluster.cluster_id,
                cluster_id: cluster.cluster_id,
                round:      round,
            ]
            tuple( meta, [
                file(cluster.merged_reads.read1),
                file(cluster.merged_reads.read2),
            ])
        }
        | set { ch_spades_input }

    /*
    ── Step 2: SPAdes 共组装 ──
    */
    SPADES_COASSEMBLY ( ch_spades_input )
    ch_versions = ch_versions.mix( SPADES_COASSEMBLY.out.versions )

    /*
    ── Step 3: CheckM2 批量评估 ──
    */
    SPADES_COASSEMBLY.out.contigs
        | map { meta, fasta -> fasta }
        | collect
        | CHECKM2
    ch_versions = ch_versions.mix( CHECKM2.out.versions )

    /*
    ── Step 4: contig 路径列表  cluster_id <TAB> /abs/path ──
    */
    SPADES_COASSEMBLY.out.contigs
        | map { meta, fasta ->
            "${meta.cluster_id}\t${fasta.toRealPath()}"
        }
        | collectFile(
            name:    "round${round}_contig_list.txt",
            newLine: true,
            sort:    true,
        )
        | set { ch_contig_list }

    /*
    ── Step 5: 更新 JSON ──
    输出:
      updated_json      → 完整 JSON（供下一轮 PREPARE_ROUND_JSON 使用）
      filtered_clusters → per-cluster JSON（供 branch 分流）
    */
    UPDATE_CLUSTER_JSON (
        cluster_json,
        CHECKM2.out.results,
        ch_contig_list,
    )
    ch_versions = ch_versions.mix( UPDATE_CLUSTER_JSON.out.versions )

    /*
    ── Step 6: 分流 ──
    needs_opt: do_optimize && contam > max_contam && comp > min_comp
    direct:    其余所有 cluster
    */
    UPDATE_CLUSTER_JSON.out.filtered_clusters
        | flatten
        | map { json ->
            def c      = new groovy.json.JsonSlurper().parse(json)
            def comp   = c.checkm2_coassembly?.completeness  ?: 90.0
            def contam = c.checkm2_coassembly?.contamination ?: 10.0
            def meta   = [
                id:         c.cluster_id,
                cluster_id: c.cluster_id,
                round:      round,
            ]
            def needs_opt = do_optimize && (contam > max_contam) && (comp > min_comp)
            tuple( meta, json, needs_opt )
        }
        | branch { meta, json, needs_opt ->
            optimize: needs_opt == true
            direct:   true
        }
        | set { ch_split }

    /*
    ── Step 7: TNF 优化 ──
    */
    ch_split.optimize
        | map { meta, json, _ -> tuple( meta, json ) }
        | COSAG_OPTIMIZATION
    ch_versions = ch_versions.mix( COSAG_OPTIMIZATION.out.versions )

    /*
    ── Step 8: 直接使用优化流程产出的最佳 contig（不再二次 SPAdes）──
    */
    COSAG_OPTIMIZATION.out.best_contigs
        | set { ch_optimized_contigs }

    /*
    ── 汇总最终 contigs（优化后 + 直接通过）──
    直接通过的 cluster 从 per-cluster JSON 取 coassembly_contigs 路径
    */
    ch_split.direct
        | map { meta, json, _ ->
            def c = new groovy.json.JsonSlurper().parse(json)
            tuple( meta, file(c.coassembly_contigs) )
        }
        | set { ch_direct_contigs }

    ch_final_contigs = ch_optimized_contigs
        | mix( ch_direct_contigs )

    COSAG_OPTIMIZATION.out.optimized_json
        | map { meta, json -> json }
        | collect
        | set { ch_opt_jsons }

    emit:
    final_contigs   = ch_final_contigs                      // [ val(meta), path ] → 下一轮 or QUALITY_TAXONOMY
    updated_json    = UPDATE_CLUSTER_JSON.out.updated_json  // path → 下一轮 PREPARE_ROUND_JSON
    opt_jsons       = ch_opt_jsons
    checkm2_results = CHECKM2.out.results
    versions        = ch_versions
}
