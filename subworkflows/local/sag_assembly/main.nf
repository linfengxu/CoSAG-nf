/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/sag_assembly/main.nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULE 1: Individual SAG assembly, quality assessment, and taxonomic classification

    Flow:
        reads [ val(meta), path(reads) ]
            → SPADES              (个体 SAG 组装)
            → CHECKM2             (批量质量评估)
            → FILTER_SAGS         (污染度过滤，默认只卡 contamination > 10%)
            → GTDBTK_CLASSIFYWF   (GTDB-Tk 分类，只对通过过滤的 SAG)
            → emit: contigs, checkm2, bac_summary, ar_summary
----------------------------------------------------------------------------------------
*/

include { SPADES            } from '../../../modules/local/spades/main'
include { CHECKM2           } from '../../../modules/local/checkm2/main'
include { FILTER_SAGS       } from '../../../modules/local/filter_sags/main'
include { GTDBTK_CLASSIFYWF } from '../../../modules/local/gtdbtk_classifywf/main'


workflow SAG_ASSEMBLY {

    take:
    reads   // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ] ]

    main:
    ch_versions = Channel.empty()

    /*
    ────────────────────────────────────────────────────────────
        Step 1: 个体 SAG 组装
    ────────────────────────────────────────────────────────────
    */
    SPADES ( reads )
    ch_versions = ch_versions.mix( SPADES.out.versions )

    /*
    ────────────────────────────────────────────────────────────
        Step 2: 批量 CheckM2 质量评估
    ────────────────────────────────────────────────────────────
    所有 SAG contigs 收集后批量提交，比逐个提交效率高
    */
    SPADES.out.contigs
        | map { meta, contig -> contig }
        | collect
        | CHECKM2
    ch_versions = ch_versions.mix( CHECKM2.out.versions )

    /*
    ────────────────────────────────────────────────────────────
        Step 3: 过滤高污染 SAG
    ────────────────────────────────────────────────────────────
    默认只卡 contamination > params.max_contamination (10%)
    completeness 不在此阶段过滤（MDA 扩增偏差导致完整度普遍偏低）
    */
    FILTER_SAGS ( CHECKM2.out.results )
    ch_versions = ch_versions.mix( FILTER_SAGS.out.versions )

    /*
    ────────────────────────────────────────────────────────────
        Step 4: 筛选通过过滤的 contigs
    ────────────────────────────────────────────────────────────
    sags_pass.tsv 的 Name 列格式为 <sag_id>_contigs
    需要去掉 _contigs 后缀才能和 meta.id 匹配
    */
    FILTER_SAGS.out.pass_tsv
        | splitCsv( header: true, sep: '\t' )
        | map { row -> row.Name.replaceAll('_contigs$', '') }
        | set { ch_pass_names }

    SPADES.out.contigs
        | map { meta, contig -> [ meta.id, meta, contig ] }
        | join( ch_pass_names | map { name -> [ name, name ] } )
        | map { id, meta, contig, name -> [ meta, contig ] }
        | set { ch_pass_contigs }

    /*
    ────────────────────────────────────────────────────────────
        Step 5: GTDB-Tk 分类（只对通过过滤的 SAG）
    ────────────────────────────────────────────────────────────
    个体 SAG GTDB 分类用于：
      1. 验证下游 MinHash 聚类合理性（同 cluster 应同物种）
      2. 背景多样性描述（论文方法部分）
    */
    ch_pass_contigs
        | map { meta, contig -> contig }
        | collect
        | map { contigs -> [ [ id: 'all_sags' ], contigs ] }
        | GTDBTK_CLASSIFYWF
    ch_versions = ch_versions.mix( GTDBTK_CLASSIFYWF.out.versions )

    emit:
    contigs     = ch_pass_contigs                   // [ val(meta), path(fasta) ] → MINHASH_CLUSTER
    all_contigs = SPADES.out.contigs                // [ val(meta), path(fasta) ] 全部，备用
    pass_tsv    = FILTER_SAGS.out.pass_tsv          // path(sags_pass.tsv)
    checkm2     = CHECKM2.out.results               // path(quality_report.tsv) → COASSEMBLY
    bac_summary = GTDBTK_CLASSIFYWF.out.bac_summary // path(gtdbtk.bac120.summary.tsv)
    ar_summary  = GTDBTK_CLASSIFYWF.out.ar_summary  // path(gtdbtk.ar53.summary.tsv)
    versions    = ch_versions
}