process SPLIT_CRAM_BY_CHROM {
    tag "$meta.id:$chrom"
    maxForks 24  // Limits the number of concurrent executions of this process to 24
    label 'process_medium'

    conda "bioconda::samtools=1.16"
    container "quay.io/biocontainers/samtools:1.16--hfe4b78e_1"

    input:
        tuple val(meta), path(cram), path(crai), val(chrom)
        path(fasta)

    output:
        tuple val(meta), path("${meta.id}.${chrom}.cram"), emit: cram
        tuple val(meta), path("${meta.id}.${chrom}.cram.crai"), emit: crai
        path("versions.yml"), emit: versions

    script:
     """
        samtools view \
            -@ ${task.cpus} \\
            -O cram \\
            -o ${meta.id}.${chrom}.cram \\
            ${cram} \\
            ${chrom}
            
        samtools index -@ ${task.cpus} ${meta.id}.${chrom}.cram 
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
        END_VERSIONS
        """
}
