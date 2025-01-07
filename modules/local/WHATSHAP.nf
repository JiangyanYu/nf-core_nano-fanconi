process WHATSHAP {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-whatshap:240302' :
        'ghcr.io/dhslab/docker-whatshap:240302' }"

    input:
        val(meta)
	path(sniffle_vcf_file)
	path(sorted_bam)
	path(reference_fasta)
	path(index)

    output:
        tuple val(meta), path("${meta.sample}*.haplotagged.bam")     , emit: bam
        tuple val(meta), path("${meta.sample}*.haplotagged.bam.bai") , emit: bai
        path  ("versions.yml")                                       , emit: versions

    script:
    def vcf_file = sniffle_vcf_file.name != 'NO_FILE.vcf' ? "$sniffle_vcf_file" : "${meta.sample}.vcf.gz"
    """
    whatshap haplotag --tag-supplementary --ignore-read-groups --output-threads=${task.cpus} \\
    -o ${meta.sample}.haplotagged.bam --reference ${reference_fasta} $vcf_file ${sorted_bam} && \\
    samtools index ${meta.sample}.haplotagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version)
    END_VERSIONS
    """
}
