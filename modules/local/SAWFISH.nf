process SAWFISH {
    maxForks 2  // Limits the number of concurrent executions of this process to 2
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.2' :
        'jiangyanyu/pacbio_wgs:v1.2' }"

    input:

        tuple val(meta), path(bam_file), path(bam_bai_file) 
        path (index)

    output:
        path ("joint-call/*alignment*")                  , emit: bam
        path ("joint-call/genotyped.sv.vcf.gz")          , emit: vcf
        path ("joint-call/genotyped.sv.vcf.gz.tbi")      , emit: tbi
        path "versions.yml"                                             , emit: versions

    script:
        def args = task.ext.args ?: ''
        """   
        source /opt/condaetc/etc/profile.d/conda.sh  
        conda activate sawfish

        sawfish discover \\
                --threads ${task.cpus} \\
                --ref ${index} \\
                --bam ${meta.sample}.sorted.bam \\
                --output-dir discover

        sawfish joint-call
                --threads ${task.cpus} \\
                --sample discover \\
                --output-dir joint-call

        cat <<-END_VERSIONS > versions.yml

        "${task.process}":
            sawfish: \$(sawfish --version 2>&1)
        END_VERSIONS
        """
}
