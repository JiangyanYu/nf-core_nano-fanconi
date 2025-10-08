process FAIDX_REFERENCE {
    label 'process_medium'

    conda "bioconda::samtools=1.16"
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        path fasta
        
    output:
        path("genome.fa"), emit: fasta
        path("genome.fa.fai"), emit: fasta_index
        path("versions.yml"), emit: versions

    script:
        """
        # Link the input FASTA as genome.fa
        ln -sf ${fasta} genome.fa
    
        # Create FAI index
        samtools faidx genome.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
            wget: \$(samtools --version | grep '^samtools' | head -n1)
        END_VERSIONS
        """
}
