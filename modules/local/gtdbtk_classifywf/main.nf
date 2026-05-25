process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_1'
    //container '/mnt/data/xulf/software/sif/gtdbtk_with_ps.sif'
    containerOptions "--bind ${params.gtdb_database}:/refdata"

    input:
    tuple val(meta), path(bin_fastas)

    output:
    tuple val(meta), path("gtdb_results/*"),                                    emit: results
    tuple val(meta), path("gtdb_results/gtdbtk.bac120.summary.tsv"), optional: true, emit: bac_summary
    tuple val(meta), path("gtdb_results/gtdbtk.ar53.summary.tsv"),   optional: true, emit: ar_summary
    tuple val(meta), path("gtdb_results/gtdbtk.bac120.classify.tree"), optional: true, emit: bac_tree
    tuple val(meta), path("gtdb_results/gtdbtk.ar53.classify.tree"),  optional: true, emit: ar_tree
    tuple val(meta), path("gtdb_results/gtdbtk.log"),                           emit: gtdb_log
    path "versions.yml",                                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def skip_ani   = params.gtdb_skip_ani_screen ? '--skip_ani_screen' : ''
    def min_perc   = params.gtdb_min_perc_aa
    """
    mkdir -p bin_fastas gtdb_results

    # collect all fasta 
    for f in ${bin_fastas}; do
        cp \$f bin_fastas/
    done

    gtdbtk classify_wf \\
        --genome_dir bin_fastas \\
        --out_dir gtdb_results \\
        --cpus ${task.cpus} \\
        --extension fasta \\
        --min_perc_aa ${min_perc} \\
        --force \\
        ${skip_ani} \\
        ${args} \\
        > gtdb_results/gtdbtk.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version 2>&1 | sed 's/gtdbtk: version //; s/ .*//')
    END_VERSIONS
    """
}