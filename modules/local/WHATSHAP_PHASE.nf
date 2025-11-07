process WHATSHAP_PHASE {
    tag "$meta.id:$chrom"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com/repository/docker/jiangyanyu/docker-whatshap/' :
        'jiangyanyu/docker-whatshap:v240302' }"

    input:
        tuple val(meta), path(split_cram), path(split_crai), path(split_vcf), path(split_tbi), val(caller), val(chrom)
        path(fasta)

    output:
        tuple val(meta), path("${meta.id}.${caller}.${chrom}.phased.vcf.gz"), path("${meta.id}.${caller}.${chrom}.phased.vcf.gz.tbi"), val(caller), val(chrom), emit: vcf_tbi_caller_chrom
        path  ("versions.yml")                                       , emit: versions

    script:

    """
    whatshap phase -o ${meta.id}.${caller}.${chrom}.phased.vcf.gz \\
        --reference=${fasta} \\
        ${meta.id}.${caller}.${chrom}.vcf.gz \\
        ${meta.id}.${chrom}.cram


    tabix ${meta.id}.${caller}.${chrom}.phased.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    END_VERSIONS

    """
}
