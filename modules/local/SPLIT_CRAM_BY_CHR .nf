process SPLIT_VCF_BY_CHR {
    label 'process_medium'

    conda "bioconda::samtools=1.16"
    container "quay.io/biocontainers/samtools:1.16--hfe4b78e_1"

    input:
        tuple val(meta), path(cram_file)            
        val(chr)
        path(fasta)


    output:
        tuple val(meta), path("chr*.cram"), emit: split_crams
        path("versions.yml"), emit: versions

    script:
     """
        samtools view -O cram -o ${meta.sample}.${chr}.cram ${cram_file} ${chr}
        
        samtools index ${meta.sample}.${chr}.cram ${meta.sample}.${chr}.cram.crai
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: $(samtools --version | head -n1 | sed 's/^samtools //')
        END_VERSIONS
        """
}
