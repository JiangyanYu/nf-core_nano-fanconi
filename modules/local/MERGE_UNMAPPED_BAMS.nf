process MERGE_UNMAPPED_BAMS {
    label 'process_medium'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
        tuple val(meta), path(unmapped_bams)

    output:
        tuple val(meta), path("${meta.id ?: meta.sample}.unaligned.bam"), emit: merged_unmapped_bam
        path  ("versions.yml"), emit: versions

    script:
    def prefix = meta.id ?: meta.sample
    def is_single_file = unmapped_bams.size() == 1
    def merge_cmd = is_single_file ? 
        "ln -s ${unmapped_bams[0]} ${prefix}.unaligned.bam" :
        "samtools merge -f -@ ${task.cpus} ${prefix}.unaligned.bam ${unmapped_bams.join(' ')}"

// Print the meta information for debugging
 println "Merging unmapped BAMs for sample: ${meta.sample}, ID: ${meta.id ?: 'N/A'}"

    """
        ${merge_cmd}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
        END_VERSIONS
    """
}
