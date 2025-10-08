process MERGE_BASECALL {
    label 'process_medium'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
        tuple val(meta), path(input_bams)

    output:
        tuple val(meta), path("${meta.id ?: meta.sample}.unaligned.bam"), emit: merged_bam
        path  ("versions.yml"), emit: versions

    script:
    def prefix = meta.id ?: meta.sample
    def is_single_file = input_bams instanceof Path
    def merge_cmd = is_single_file ? 
        "ln -s ${input_bams} ${prefix}.unaligned.bam" :
        "samtools merge -f -@ ${task.cpus} ${prefix}.unaligned.bam ${input_bams.join(' ')}"

    """
        ${merge_cmd}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
        END_VERSIONS
    """
}
