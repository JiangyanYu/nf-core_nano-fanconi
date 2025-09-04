process PREPARE_REFERENCES {
    label 'process_medium'

    conda "bioconda::samtools=1.16"
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        val(profile)
        val(fasta_url)
        val(fasta_index_url)

    output:
        path("genome.fa"), emit: fasta
        path("genome.fa.fai"), emit: fai
        path("versions.yml"), emit: versions

    script:
        """
        # Download FASTA
        wget -O genome.fa "${fasta_url}"
        
        # Create FAI index
        samtools faidx genome.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: $(samtools --version | head -n1 | sed 's/^samtools //')
            wget: $(wget --version | head -n1)
        END_VERSIONS
        """
}
