process SOURMASH_COMPARE_NF {
    tag "sourmash_compare"
    label 'process_high'

    container 'quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0'
    //container '/mnt/data/xulf/software/sif/sourmash_4.9.2--hdfd78af_1.sif'
    input:
    path(signatures)   // all .sig files collected

    output:
    path("similarity_matrix.npy"),            emit: matrix_npy
    path("similarity_matrix.csv"),            emit: matrix_csv
    path("similarity_matrix.npy.labels.txt"), emit: labels
    path("sourmash_compare.log"),             emit: log
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    n_sigs=\$(ls *.sig 2>/dev/null | wc -l)
    echo "Found \${n_sigs} signature files" | tee sourmash_compare.log

    if [ "\${n_sigs}" -lt 2 ]; then
        echo "ERROR: Need at least 2 signatures for comparison" | tee -a sourmash_compare.log
        exit 1
    fi

    sourmash compare \\
        *.sig \\
        --ksize ${params.sourmash_ksize} \\
        -o similarity_matrix.npy \\
        --csv similarity_matrix.csv \\
        ${args} \\
        2>&1 | tee -a sourmash_compare.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS
    """
}


process SOURMASH_COMPARE {
    tag "sourmash_compare"
    label 'process_high'
    container '/mnt/data/xulf/software/sif/sourmash_4.9.2--hdfd78af_1.sif'


    input:
    path(signatures)
    val(ksize)
    val(scaled)

    output:
    path("similarity_matrix.npy"), emit: similarity_matrix
    path("similarity_matrix.npy.labels.txt"), emit: labels
    path("sourmash_compare.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    # Get all signature files
    SIG_FILES=(\$(ls *.sig))
    echo "Found \${#SIG_FILES[@]} signature files" | tee sourmash_compare.log
    
    # Check if we have enough signatures
    if [ \${#SIG_FILES[@]} -lt 2 ]; then
        echo "Error: Need at least 2 valid signatures for comparison" | tee -a sourmash_compare.log
        # Create empty files for pipeline continuation
        touch similarity_matrix.npy
        touch similarity_matrix.npy.labels.txt
        exit 0
    fi
    
    # Use sourmash compare to generate binary matrix
    echo "Running sourmash compare on \${#SIG_FILES[@]} signatures..." | tee -a sourmash_compare.log
    echo "Command: sourmash compare *.sig -o similarity_matrix.npy --ksize ${ksize}" | tee -a sourmash_compare.log
    
    sourmash compare *.sig -o similarity_matrix.npy --ksize ${ksize} 2>&1 | tee -a sourmash_compare.log
    
    # Check if sourmash compare succeeded
    if [ ! -f "similarity_matrix.npy" ]; then
        echo "Error: sourmash compare failed to create .npy file" | tee -a sourmash_compare.log
        touch similarity_matrix.npy
        touch similarity_matrix.npy.labels.txt
        exit 1
    fi
    
    # Check if labels file was created
    if [ ! -f "similarity_matrix.npy.labels.txt" ]; then
        echo "Warning: Labels file not created by sourmash" | tee -a sourmash_compare.log
        touch similarity_matrix.npy.labels.txt
    fi
    
    echo "sourmash compare completed successfully" | tee -a sourmash_compare.log
    echo "Matrix file size: \$(ls -lh similarity_matrix.npy | awk '{print \$5}')" | tee -a sourmash_compare.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS
    """
}