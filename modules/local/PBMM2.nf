process PBMM2 {
    maxForks 8  // Limits the number of concurrent executions of this process to 8
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.2' :
        'jiangyanyu/pacbio_wgs:v1.2' }"

    input:

        tuple val(meta), path (bam_paths) 
        path (fasta)

    output:
        tuple val(meta), path ("*.cram")       , emit: cram
        path "versions.yml"                   , emit: versions

    script:
        def args = task.ext.args ?: ''
        """      
        samtools concat -@ ${task.cpus} ${bam_paths} | \\
        pbmm2 align \\
                --num-threads ${task.cpus} \\
                --preset CCS \\
                ${args} \\
                ${fasta} \\
                /dev/stdin | \\
        samtools sort -@ ${task.cpus} -O cram -o ${meta.sample}.cram
        samtools index -@ ${task.cpus} ${meta.sample}.cram

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pbmm2: \$(pbmm2 --version 2>&1)
        END_VERSIONS
        """
}
