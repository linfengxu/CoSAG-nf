process CHECKM2 {
    tag "checkm2_batch"
    label 'process_high'

    container 'quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0'
    //container '/mnt/data/xulf/software/sif/checkm2_1.0.2.sif'

    input:
    path(contigs)   // Collect SAG  contigs.fasta

    output:
    path("checkm2_results/quality_report.tsv"), emit: results
    path("checkm2_results.log"),                emit: log
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def add_args = params.checkm2_additional_args ?: ''
    """
    # collect
    mkdir -p input_bins
    for f in ${contigs}; do
        cp \$f input_bins/
    done

    export CHECKM2DB=${params.checkm2_db}

    checkm2 predict \\
        --input input_bins \\
        --output-directory checkm2_results \\
        --extension fasta \\
        --threads ${task.cpus} \\
        ${add_args} \\
        ${args} \\
        > checkm2_results.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version 2>&1 | sed 's/.*checkm2 v//')
    END_VERSIONS
    """
}