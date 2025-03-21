process WHATSHAP_HAPLOTAG {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/docker-whatshap:v240302' :
        'jiangyanyu/docker-whatshap:v240302' }"

    input:
        tuple val(meta), path(bam_file), path(bam_bai_file)
        tuple val(meta), path(phased_vcf), path(phased_tbi)
        path(reference_fasta)
        path(index)

    output:
        tuple val(meta), path("${meta.sample}*.haplotagged.bam")     , emit: bam
        tuple val(meta), path("${meta.sample}*.haplotagged.bam.bai") , emit: bai
        path  ("versions.yml")                                       , emit: versions

    script:
    // def vcf_file = phased_vcf.name != 'NO_FILE.vcf' ? "$phased_vcf" : "${meta.sample}.phased.vcf.gz"
    // def vcf_file = phased_vcf.name != 'test.vcf' ? "$phased_vcf" : "${meta.sample}.vcf.gz"
    """

    whatshap haplotag --tag-supplementary --ignore-read-groups --output-threads=${task.cpus} \\
    -o ${meta.sample}.haplotagged.bam --reference ${reference_fasta} ${meta.sample}.vcf.gz ${meta.sample}.sorted.bam && \\
    samtools index ${meta.sample}.haplotagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}
