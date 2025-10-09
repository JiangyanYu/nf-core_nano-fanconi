process PBMM2_FROM_BAM {
    maxForks 8  // Limits the number of concurrent executions of this process to 8
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.2' :
        'jiangyanyu/pacbio_wgs:v1.2' }"

    input:
        tuple val(meta), path (unmapped_bams, stageAs: "input_bam_??.bam") 
        path (fasta)
        path (fasta_index)

    output:
        tuple val(meta), path ("*.cram"), emit: cram
        path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: ''
        """
        echo "${unmapped_bams}" | \\
        sed 's/ /\\n/g' | \\
        cat > ${meta.id}.fofn \\

        pbmm2 align \\
                --num-threads ${task.cpus} \\
                --preset CCS \\
                ${args} \\
                ${fasta} \\
                ${meta.id}.fofn | \\
        samtools sort -@ ${task.cpus} --reference ${fasta} /dev/stdin | \\
        samtools addreplacerg -@ ${task.cpus} -r "ID:${meta.id}\\tSM:${meta.id}" -O cram -o ${meta.id}.cram /dev/stdin \\
        samtools index -@ ${task.cpus} ${meta.id}.cram

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | sed 's/^samtools //')
            pbmm2: \$(pbmm2 --version 2>&1 | head -n 1)
        END_VERSIONS
        """
}

// 
        