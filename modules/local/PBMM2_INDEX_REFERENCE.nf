process PBMM2_INDEX_REFERENCE {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.2' :
        'jiangyanyu/pacbio_wgs:v1.2' }"

    input:
        path (fasta)
        path (fasta_index)

    output:
        path ("*.mmi"), emit: mmi
        path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: ''
        """
        pbmm2 index \\
                ${fasta} \\
                ${fasta}.mmi \\
                --num-threads ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pbmm2: \$(pbmm2 --version | head -n 1 | sed 's/^pbmm2 //g')
        END_VERSIONS
        """
}
