process WHATSHAP_HAPLOTAG {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/docker-whatshap:v240302' :
        'jiangyanyu/docker-whatshap:v240302' }"

    input:
        tuple val(meta), path(split_crams), path(split_crais)
        tuple val(meta), path(phased_split_vcfs), path(phased_split_tbis)
        val(chr)
        path(fasta)

    output:
        tuple val(meta), path("${meta.id}*.${chr}.haplotagged.cram")     , emit: cram
        tuple val(meta), path("${meta.id}*.${chr}.haplotagged.cram.crai") , emit: crai
        path  ("versions.yml")                                       , emit: versions

    script:
    // def vcf_file = phased_merged_vcf.name != 'NO_FILE.vcf' ? "$phased_merged_vcf" : "${meta.id}.phased.vcf.gz"
    // def vcf_file = phased_merged_vcf.name != 'test.vcf' ? "$phased_merged_vcf" : "${meta.id}.vcf.gz"
    """

    # Filter by MG>=95
    samtools view --reference ${fasta} -h -e '[mg] && [mg]>=95' ${meta.id}.cram | \\

    whatshap haplotag --tag-supplementary --ignore-read-groups --output-threads=${task.cpus} \\
    -o ${meta.id}.haplotagged.cram --reference ${fasta} ${meta.id}.vcf.gz /dev/stdin

    samtools view --reference ${fasta} -h -e '[mg] && [mg]<95' -O cram -o ${meta.id}.not_haplotagged.cram ${meta.id}.cram
    
    samtools merge -@ ${task.cpus} -O cram -o ${meta.id}.haplotagged_merged.cram ${meta.id}.haplotagged.cram ${meta.id}.not_haplotagged.cram

    rm ${meta.id}.not_haplotagged.cram ${meta.id}.haplotagged.cram

    mv ${meta.id}.haplotagged_merged.cram ${meta.id}.haplotagged.cram

    samtools index -@ ${task.cpus} ${meta.id}.haplotagged.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}
