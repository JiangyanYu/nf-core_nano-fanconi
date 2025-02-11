process WHATSHAP_PHASE {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com/repository/docker/jiangyanyu/docker-whatshap/' :
        'jiangyanyu/docker-whatshap:v240302' }"

    input:
    tuple val(meta), path(bam_file), path(bam_bai_file)
    path(reference_fasta)
    path(index)


    output:
        tuple val(meta), path("${meta.sample}_phased.vcf")           , emit: phased_vcf
        path  ("versions.yml")                                       , emit: versions

    script:

    """

    #whatshap phase -o ${meta.sample}_phased.vcf \\
    #    --reference=${reference_fasta} \\
    #    ${meta.sample}.vcf.gz ${meta.sample}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    END_VERSIONS

    """
}
