process PBMM2 {
    maxForks 2  // Limits the number of concurrent executions of this process to 2
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.2' :
        'jiangyanyu/pacbio_wgs:v1.2' }"

    input:

        tuple val(meta), path (reads_paths) 
        path (index)

    output:
        tuple val(meta), path ("*.bam")       , emit: bam
        path "versions.yml"                   , emit: versions

    script:
        def args = task.ext.args ?: ''
        """        
        pbmm2 align \\
                --num-threads ${task.cpus} \\
                --preset CCS \\
                --sort \\
                ${index} \\
                ${reads_paths} \\
                ${meta.sample}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pbmm2: \$(pbmm2 --version 2>&1)
        END_VERSIONS
        """
}
