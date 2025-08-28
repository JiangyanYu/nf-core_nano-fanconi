process WHATSHAP_HAPLOTAG {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/docker-whatshap:v240302' :
        'jiangyanyu/docker-whatshap:v240302' }"

    input:
        tuple val(meta), path(cram_file), path(crai_file)
        tuple val(meta), path(phased_merged_vcf), path(phased_merged_tbi)
        path(fasta)
        path(fasta_index)

    output:
        tuple val(meta), path("${meta.sample}*.haplotagged.cram")     , emit: cram
        tuple val(meta), path("${meta.sample}*.haplotagged.cram.crai") , emit: crai
        path  ("versions.yml")                                       , emit: versions

    script:
    // def vcf_file = phased_merged_vcf.name != 'NO_FILE.vcf' ? "$phased_merged_vcf" : "${meta.sample}.phased.vcf.gz"
    // def vcf_file = phased_merged_vcf.name != 'test.vcf' ? "$phased_merged_vcf" : "${meta.sample}.vcf.gz"
    """

    whatshap haplotag --tag-supplementary --ignore-read-groups --output-threads=${task.cpus} \\
    -o ${meta.sample}.haplotagged.cram --reference ${fasta} ${meta.sample}.vcf.gz ${meta.sample}.cram && \\
    samtools index -@ ${task.cpus} ${meta.sample}.haplotagged.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}
