process SPLIT_VCF_BY_CHROM {
    tag "$meta.id: $caller:$chrom"
    maxForks 24  // Limits the number of concurrent executions of this process to 24
    label 'process_medium'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
        tuple val(meta), path(vcf), val(caller), val(chrom)

    output:
        tuple val(meta), path("${meta.id}.${caller}.${chrom}.vcf.gz"), emit: vcf
        tuple val(meta), path("${meta.id}.${caller}.${chrom}.vcf.gz.tbi"), emit: tbi
        path("versions.yml"), emit: versions

    script:
        """
        tabix ${vcf}

        bcftools view \\
            --threads ${task.cpus} \\
            -r ${chrom} \\
            -O z \\
            -o ${meta.id}.${caller}.${chrom}.vcf.gz \\
            ${vcf}

        tabix ${meta.id}.${caller}.${chrom}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
        END_VERSIONS
        """
}
