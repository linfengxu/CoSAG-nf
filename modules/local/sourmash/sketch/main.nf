process SOURMASH_SKETCH {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0'
    //container '/mnt/data/xulf/software/sif/sourmash_4.9.2--hdfd78af_1.sif'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_k${params.sourmash_ksize}.sig"), emit: signature, optional: true
    path("${meta.id}_sketch.log"),                                      emit: log
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = meta.id
    def filter_enabled = params.sourmash_filter_small_sketches ?: false
    """
    sourmash sketch dna \\
        -p k=${params.sourmash_ksize},scaled=${params.sourmash_scaled} \\
        --name ${prefix} \\
        -o ${prefix}_k${params.sourmash_ksize}.sig \\
        ${fasta} \\
        ${args} \\
        > ${prefix}_sketch.log 2>&1

    if [ -f "${prefix}_k${params.sourmash_ksize}.sig" ]; then
        if [ "${filter_enabled}" == "true" ]; then
            n_hashes=\$(sourmash sig describe ${prefix}_k${params.sourmash_ksize}.sig 2>&1 \\
                | grep -i "num hashes\\|n hashes\\|hashes" \\
                | head -1 \\
                | grep -oE '[0-9]+' \\
                | tail -1 || echo 0)
            echo "n_hashes=\${n_hashes}" >> ${prefix}_sketch.log
            if [ "\${n_hashes}" -lt "${params.sourmash_min_sketch_size}" ]; then
                echo "WARNING: ${prefix} has only \${n_hashes} hashes (min: ${params.sourmash_min_sketch_size})" \\
                    >> ${prefix}_sketch.log
                echo "WARNING: ${prefix} may have unreliable similarity estimates" \\
                    >> ${prefix}_sketch.log
                rm -f ${prefix}_k${params.sourmash_ksize}.sig
            fi
        else
            echo "INFO: sketch size filtering disabled, keeping all signatures" \\
                >> ${prefix}_sketch.log
        fi
    else
        echo "ERROR: Failed to create signature for ${prefix}" >> ${prefix}_sketch.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS
    """
}