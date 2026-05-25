process BARRNAP_COSAG {
    tag "${meta.id}"
    label 'process_low'

    container "quay.io/biocontainers/barrnap:0.9--0"
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_barrnap.gff"), emit: gff
    path("${meta.id}.16S.fa"), emit: rna16
    path("${meta.id}_barrnap.log"), emit: log
    path("${meta.id}_versions_barrnap.yml"), emit: versions

    when:
    params.run_cosag_rrna_annotation && (task.ext.when == null || task.ext.when)

    script:
    def args = task.ext.args ?: ''
    """
    barrnap \\
        --threads ${task.cpus} \\
        --outseq ${meta.id}_barrnap.rrna.fa \\
        ${args} \\
        ${fasta} \\
        > ${meta.id}_barrnap.gff \\
        2> ${meta.id}_barrnap.log

    # Extract 16S from barrnap --outseq FASTA (headers like >16S_rRNA::...)
    : > ${meta.id}.16S.fa
    if [ -s ${meta.id}_barrnap.rrna.fa ]; then
        awk '/^>/{keep=0} /^>16S_rRNA/{keep=1} keep' ${meta.id}_barrnap.rrna.fa >> ${meta.id}.16S.fa
    fi

    cat <<-END_VERSIONS > ${meta.id}_versions_barrnap.yml
    "${task.process}":
        barrnap: \$(barrnap --version 2>&1 | head -1 | sed 's/barrnap //')
    END_VERSIONS
    """
}
