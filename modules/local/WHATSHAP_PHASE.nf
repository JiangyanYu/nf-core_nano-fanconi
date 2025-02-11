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
    
        tuple val(meta), path (bam_file), path (bam_bai_file), path (sniffles_vcf), path (sniffles_vcf_tbi)
        path(reference_fasta)
        path(index)

    output:
        tuple val(meta), path("${meta.sample}*_phased.vcf")          , emit: phased_vcf
        path  ("versions.yml")                                       , emit: versions

    script:
    // def vcf_file = phased_vcf.name != 'NO_FILE.vcf' ? "$phased_vcf" : "${meta.sample}.phased.vcf.gz"
    // def vcf_file = sniffles_vcf.name != 'test.vcf' ? "$sniffles_vcf" : "${meta.sample}.vcf.gz"
    """
    echo "WORKING DIRECTORY: \$(pwd)"
    echo "META: ${meta}"
    echo "BAM FILE: $bam_file"
    echo "BAI FILE: $bam_bai_file"
    echo "VCF FILE: $sniffles_vcf"
    echo "TBI FILE: $sniffles_vcf_tbi"
    echo "REFERENCE FASTA: $reference_fasta"
    echo "INDEX FILE: $index"

    ls -lh  # List all files in the work directory

    #whatshap phase \\
    #-o ${meta.sample}_phased.vcf \\
    #--reference ${reference_fasta} \\
    #${sniffles_vcf} ${bam_file} ## ${meta.sample}.vcf.gz ${meta.sample}.sorted.bam

    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    whatshap: \$(whatshap --version |sed 's/^.*Version: //')
    #END_VERSIONS
    """
}
