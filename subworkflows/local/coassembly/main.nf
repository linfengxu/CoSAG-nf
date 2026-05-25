/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/coassembly/main.nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COASSEMBLY_ROUND as COASSEMBLY_ROUND_1 } from '../coassembly_round/main'
include { COASSEMBLY_ROUND as COASSEMBLY_ROUND_2 } from '../coassembly_round/main'
include { COASSEMBLY_ROUND as COASSEMBLY_ROUND_3 } from '../coassembly_round/main'
include { MINHASH_CLUSTER  as MINHASH_CLUSTER_R2 } from '../minhash_cluster/main'
include { MINHASH_CLUSTER  as MINHASH_CLUSTER_R3 } from '../minhash_cluster/main'


workflow COASSEMBLY {

    take:
    round1_json     // path: round1_clusters.json

    main:
    ch_versions = Channel.empty()
    def n_rounds = (params.cosag_rounds ?: 1) as Integer

    /*
    ════════════════════════════════════════
        Round 1
    ════════════════════════════════════════
    */
    COASSEMBLY_ROUND_1 ( round1_json, 1 )
    ch_versions = ch_versions.mix( COASSEMBLY_ROUND_1.out.versions )

    ch_last_contigs = COASSEMBLY_ROUND_1.out.final_contigs
    ch_last_json    = COASSEMBLY_ROUND_1.out.updated_json

    // Round 1 JSON 保留，供 Round 3+ 溯源 SAG reads
    ch_round1_json  = round1_json

    /*
    ════════════════════════════════════════
        Round 2
    ════════════════════════════════════════
    */
    if ( n_rounds >= 2 ) {

        ch_last_contigs
            | map { meta, fasta ->
                tuple( [ id: meta.cluster_id ], fasta )
            }
            | set { ch_r2_contigs }

        MINHASH_CLUSTER_R2(
            ch_r2_contigs,
            2,
            [],
            [],
            ch_last_json,
            file("NO_FILE"),  // round1_json，Round 2 不需要
        )
        ch_versions = ch_versions.mix( MINHASH_CLUSTER_R2.out.versions )

        COASSEMBLY_ROUND_2 ( MINHASH_CLUSTER_R2.out.cluster_json, 2 )
        ch_versions = ch_versions.mix( COASSEMBLY_ROUND_2.out.versions )

        ch_last_contigs = COASSEMBLY_ROUND_2.out.final_contigs
        ch_last_json    = COASSEMBLY_ROUND_2.out.updated_json
    }

    /*
    ════════════════════════════════════════
        Round 3
        prev_json = Round 2 updated（member_type=cluster）
        需要传入 round1_json 溯源原始 SAG reads
    ════════════════════════════════════════
    */
    if ( n_rounds >= 3 ) {

        ch_last_contigs
            | map { meta, fasta ->
                tuple( [ id: meta.cluster_id ], fasta )
            }
            | set { ch_r3_contigs }

        MINHASH_CLUSTER_R3(
            ch_r3_contigs,
            3,
            [],
            [],
            ch_last_json,
            ch_round1_json,     // ← 传入 Round 1 JSON 供 PREPARE_ROUND_JSON 溯源
        )
        ch_versions = ch_versions.mix( MINHASH_CLUSTER_R3.out.versions )

        COASSEMBLY_ROUND_3 ( MINHASH_CLUSTER_R3.out.cluster_json, 3 )
        ch_versions = ch_versions.mix( COASSEMBLY_ROUND_3.out.versions )

        ch_last_contigs = COASSEMBLY_ROUND_3.out.final_contigs
        ch_last_json    = COASSEMBLY_ROUND_3.out.updated_json
    }

    // opt_jsons：所有轮次的优化结果合并
    ch_all_opt_jsons = COASSEMBLY_ROUND_1.out.opt_jsons
        | mix( n_rounds >= 2 ? COASSEMBLY_ROUND_2.out.opt_jsons : Channel.empty() )
        | mix( n_rounds >= 3 ? COASSEMBLY_ROUND_3.out.opt_jsons : Channel.empty() )

    emit:
    final_contigs   = ch_last_contigs
    updated_json    = ch_last_json
    opt_jsons       = ch_all_opt_jsons
    versions        = ch_versions
}
