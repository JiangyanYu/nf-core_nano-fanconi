process SAMTOOLS_IMPORT {
    label 'process_medium'

    conda "bioconda::samtools=1.16.1"
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        tuple val(meta), path(fastq_files)

    output:
        tuple val(meta), path("*.unmapped.bam"), emit: bam
        path("versions.yml"), emit: versions

    script:
        def fq = fastq_files instanceof List ? fastq_files.join(' ') : fastq_files
        """
        samtools import \\
            ${fq} \\
            > ${meta.sample}.unmapped.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: $(samtools --version | head -n1 | sed 's/^samtools //')
        END_VERSIONS
        """
}
