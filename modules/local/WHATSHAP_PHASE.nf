process WHATSHAP_PHASE {
    tag "$meta.id:$chrom"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com/repository/docker/jiangyanyu/docker-whatshap/' :
        'jiangyanyu/docker-whatshap:v240302' }"

    input:
        // tuple val(meta), path(cram_file), path(cram_crai_file), val(chrom)
        // tuple val(meta), path(split_vcfs), path(split_vcfs_tbi), val(caller), val(chrom)
        tuple val(meta), path(cram_file), path(cram_crai_file), path(split_vcfs), path(split_vcfs_tbi), val(caller), val(chrom)
        path(fasta)
        path(fasta_index)


    output:
        tuple val(meta), path("${meta.id}.${caller}.${chrom}.phased.vcf.gz")           , emit: vcf
        path  ("versions.yml")                                       , emit: versions

    script:

    """
    whatshap phase -o ${meta.id}.${caller}.${chrom}.phased.vcf.gz \\
        --reference=${fasta} \\
        ${meta.id}.${caller}.${chrom}.vcf.gz \\
        ${meta.id}.${chrom}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    END_VERSIONS

    """
}
