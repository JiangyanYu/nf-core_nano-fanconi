process SNIFFLES {
    tag "$meta.sample"
    label 'process_high'

    conda "bioconda::sniffles=1.0.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:1.0.12--h8b12597_1' :
        'quay.io/biocontainers/sniffles:1.0.12--h8b12597_1' }"

    input:
    tuple val(meta), path(input)


    output:
    tuple val(meta), path("${meta.sample}*_sniffles.vcf"), emit: sv_calls
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sniffles \
        -m  $input \
        -v ${meta.sample}_sniffles.vcf \
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --help 2>&1 | grep Version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}

