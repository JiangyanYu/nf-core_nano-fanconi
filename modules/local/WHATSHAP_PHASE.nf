// Function to determine the label
def determineLabel() {
    return params.use_gpu ? 'process_gpu_long' : 'process_high'
}
def processLabel = determineLabel()

process WHATSHAP_PHASE {
    maxForks 2  // Limits the number of concurrent executions of this process to 2
    label processLabel

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/docker-dorado:v241016' :
        'jiangyanyu/docker-dorado:v241016' }"

    input:
        //tuple val(meta), path(bam_bai_vcf_files), path(phased_vcf), path(phased_vcf_tbi)
        path(reference_fasta)
        path(index)

    output:
        tuple val(meta), path("${meta.sample}*_phased.vcf")          , emit: phased_vcf
        path  ("versions.yml")                                       , emit: versions

    script:
        """
        echo "Simple test of the pipeline."

        """
}
