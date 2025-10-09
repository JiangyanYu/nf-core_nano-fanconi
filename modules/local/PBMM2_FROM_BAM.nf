process PBMM2_FROM_BAM {
    maxForks 8  // Limits the number of concurrent executions of this process to 8
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.2' :
        'jiangyanyu/pacbio_wgs:v1.2' }"

    input:
        tuple val(meta), path (unmapped_bams)
        path (fasta)
        path (fasta_index)

    output:
        tuple val(meta), path ("*.cram"), emit: cram
        path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: ''
        """
        echo "${unmapped_bams}" | \\
        sed 's/ /\\n/g' > ${meta.id}.fofn

        pbmm2 align \\
                ${fasta} \\
                ${meta.id}.fofn \\
                ${meta.id}.bam \\
                --num-threads ${task.cpus} \\
                --preset CCS

        samtools sort -@ ${task.cpus} ${meta.id}.bam | \\
        samtools addreplacerg -@ ${task.cpus} --reference ${fasta} -r "ID:${meta.id}\\tSM:${meta.id}" -O cram -o ${meta.id}.cram /dev/stdin

        samtools index -@ ${task.cpus} ${meta.id}.cram
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | sed 's/^samtools //'),
            pbmm2: \$(pbmm2 --version | head -n 1 | sed 's/^pbmm2 //g')
        END_VERSIONS
        """
}

//         | \\
        // #cat > ${meta.id}.bam

        // #samtools sort -@ ${task.cpus} | \\
        // #samtools addreplacerg -@ ${task.cpus} --reference ${fasta} -r "ID:${meta.id}\\tSM:${meta.id}" -O cram -o ${meta.id}.cram /dev/stdin \\

        // #samtools index -@ ${task.cpus} ${meta.id}.cram
        