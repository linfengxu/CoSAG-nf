process SPADES {
    tag "${meta.id}"
    label 'process_high'
    
    container {params.spades.container}
    
   // publishDir "${params.outdir}/per_sample/${meta.id}/spades",
   //     mode: 'copy',
   //     saveAs: { filename ->
    //        if (filename.endsWith('.fasta')) "contigs/$filename"
   //         else if (filename.endsWith('.log')) "logs/$filename"
     //       else if (filename.endsWith('.gfa')) "graphs/$filename"
   //         else null
   //     }
    publishDir "${params.output_structure?.individual_assemblies ?: params.outdir}/per_sample/${meta.id}/spades",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.fasta')) "contigs/$filename"
            else if (filename.endsWith('.log')) "logs/$filename"
            else if (filename.endsWith('.gfa')) "graphs/$filename"
            else null
        }
    input:
    tuple val(meta), path(forward), path(reverse)
    
    output:
    tuple val(meta), path("*contigs.fasta"), emit: contigs
    tuple val(meta), path("*assembly_graph*.gfa"), optional: true, emit: graphs
    path("*.log"), emit: logs
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def memory = task.memory.toGiga()
    def kmer_sizes = params.spades?.kmers ?: "21,33,55,77,99"
    def threads = task.cpus
    def single_cell_mode = params.spades?.sc ? "--sc" : "--sc"
    def careful_mode = params.spades?.careful ? "--careful" : "--careful"
    def disable_rr = params.spades?.disable_rr ? "--disable-rr" : "--disable-rr"
    def disable_gzip = params.spades?.disable_gzip ? "--disable-gzip-output" : "--disable-gzip-output"
    def additional_args = params.spades?.additional_args ?: ""
    
    """
    spades.py \\
        -t ${threads} \\
        ${single_cell_mode} \\
        ${careful_mode} \\
        ${disable_rr} \\
        ${disable_gzip} \\
        -m ${memory} \\
        -k ${kmer_sizes} \\
        -1 ${forward} \\
        -2 ${reverse} \\
        -o ./assembly \\
        ${additional_args}
    
    # Copy only contigs.fasta file to current directory with sample prefix
    cp assembly/contigs.fasta ${prefix}_contigs.fasta
    
    # Copy assembly graph if available
    cp assembly/assembly_graph*.gfa ${prefix}_assembly_graph.gfa 2>/dev/null || true
    
    # Copy log file
    cp assembly/spades.log ${prefix}_spades.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}
