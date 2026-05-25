process PREPARE_CLUSTER_READS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), val(read1s), val(read2s)
    // read1s / read2s: comma-separated paths from parse_cluster_json.py

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: reads
    path "${meta.id}_prepare.log",                                                     emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    echo "Preparing reads for ${prefix}" > ${prefix}_prepare.log

    stream_file() {
        local f="\$1"
        if [[ "\$f" == *.gz ]]; then zcat "\$f"; else cat "\$f"; fi
    }

    # Merge R1
    {
        for f in ${read1s.replace(',', ' ')}; do
            if [ -f "\$f" ]; then
                echo "  R1: \$f" >> ${prefix}_prepare.log
                stream_file "\$f"
            else
                echo "  WARNING: R1 not found: \$f" >> ${prefix}_prepare.log
            fi
        done
    } | gzip > ${prefix}_R1.fastq.gz

    # Merge R2
    {
        for f in ${read2s.replace(',', ' ')}; do
            if [ -f "\$f" ]; then
                echo "  R2: \$f" >> ${prefix}_prepare.log
                stream_file "\$f"
            else
                echo "  WARNING: R2 not found: \$f" >> ${prefix}_prepare.log
            fi
        done
    } | gzip > ${prefix}_R2.fastq.gz

    r1_reads=\$(zcat ${prefix}_R1.fastq.gz | awk 'NR%4==1' | wc -l || echo 0)
    r2_reads=\$(zcat ${prefix}_R2.fastq.gz | awk 'NR%4==1' | wc -l || echo 0)

    echo "R1 reads: \$r1_reads" >> ${prefix}_prepare.log
    echo "R2 reads: \$r2_reads" >> ${prefix}_prepare.log

    if [ "\$r1_reads" -eq 0 ] || [ "\$r2_reads" -eq 0 ]; then
        echo "ERROR: No reads merged for ${prefix}" >> ${prefix}_prepare.log
        exit 1
    fi
    """
}
