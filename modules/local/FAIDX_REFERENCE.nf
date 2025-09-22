process PREPARE_REFERENCES {
    label 'process_medium'

    conda "bioconda::samtools=1.16"
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        val(fasta_path)
        
    output:
        path("genome.fa"), emit: fasta
        path("genome.fa.fai"), emit: fai
        path("versions.yml"), emit: versions

    script:
        """

        # Create FAI index
        samtools faidx genome.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: $(samtools --version | head -n1 | sed 's/^samtools //')
            wget: $(wget --version | head -n1)
        END_VERSIONS
        """
}
