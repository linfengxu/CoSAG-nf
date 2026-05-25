process SPADES {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    //container '/mnt/data/xulf/software/sif/spades_3.15.5.sif'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_contigs.fasta"),        emit: contigs
    tuple val(meta), path("${meta.id}_assembly_graph.gfa"),   emit: graphs,   optional: true
    tuple val(meta), path("${meta.id}_spades.log"),           emit: logs
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = meta.id
    def memory      = task.memory.toGiga()
    def kmer_sizes  = params.spades_kmers
    def sc          = params.spades_sc      ? '--sc'                  : ''
    def careful     = params.spades_careful  ? '--careful'             : ''
    def add_args    = params.spades_additional_args ?: ''
    """
    spades.py \\
        -t ${task.cpus} \\
        -m ${memory} \\
        -k ${kmer_sizes} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --sc \\
        --careful \\
        --disable-rr \\
        --disable-gzip-output \\
        -o ./assembly \\
        ${add_args} \\
        ${args}

    cp assembly/contigs.fasta          ${prefix}_contigs.fasta
    cp assembly/assembly_graph*.gfa    ${prefix}_assembly_graph.gfa 2>/dev/null || true
    cp assembly/spades.log             ${prefix}_spades.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}
