process WHATSHAP {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/docker-whatshap:v240302' :
        'jiangyanyu/docker-whatshap:v240302' }"

    input:
        tuple val(meta), path(bam_bai_vcf_files), path(sniffles_vcf), path(sniffles_vcf_tbi)
        path(reference_fasta)
        path(index)

    output:
        tuple val(meta), path("${meta.sample}*_phased.vcf")          , emit: phased_vcf
        path  ("versions.yml")                                       , emit: versions

    script:
    // def vcf_file = phased_vcf.name != 'NO_FILE.vcf' ? "$phased_vcf" : "${meta.sample}.phased.vcf.gz"
    def vcf_file = sniffles_vcf.name != 'test.vcf' ? "$sniffles_vcf" : "${meta.sample}.vcf.gz"
    """
    whatshap phase \\
    -o ${meta.sample}_phased.vcf \\
    --reference ${reference_fasta} \\
    ${meta.sample}.vcf.gz ${meta.sample}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}
