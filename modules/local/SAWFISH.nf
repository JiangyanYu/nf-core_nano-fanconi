process SAWFISH {
    tag "$meta.id"
    maxForks 8  // Limits the number of concurrent executions of this process to 8
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.2' :
        'jiangyanyu/pacbio_wgs:v1.2' }"

    input:
        tuple val(meta), path(cram)
        tuple val(meta), path(crai)
        path (fasta)
        path (fasta_index)

    output:
        tuple val(meta), path ("${meta.id}.joint-call/genotyped.sv.vcf.gz"), emit: vcf
        tuple val(meta), path ("${meta.id}.joint-call/genotyped.sv.vcf.gz.tbi"), emit: tbi
        path "versions.yml", emit: versions

    script:
        """   
        source /opt/conda/etc/profile.d/conda.sh  
        conda activate sawfish

        sawfish discover \\
                --threads ${task.cpus} \\
                --ref ${fasta} \\
                --bam ${cram} \\
                --output-dir ${meta.id}.discover

        sawfish joint-call \\
                --threads ${task.cpus} \\
                --sample ${meta.id}.discover \\
                --output-dir ${meta.id}.joint-call

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sawfish: \$(sawfish --version 2>&1)
        END_VERSIONS
        """
}
        // tuple val(meta), path ("joint-call/*alignment*"), emit: bam