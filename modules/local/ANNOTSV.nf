process ANNOTSV {
    tag "$meta.sample"
    label 'process_high'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/annotsv:3.4.4--py312hdfd78af_0' :
        'quay.io/biocontainers/annotsv:3.4.4--py312hdfd78af_0' }"
        
    containerOptions {
        "-v ${params.annotsvinput}:${params.annotsvinput} -v ${params.annotsvAnnotationsDir}:${params.annotsvAnnotationsDir}"
    }

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path("${prefix}*.tsv"),  emit: tsv
    path "versions.yml"                    ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.sample}"
    
    // Apply annotation mode flag to command
    def mode = params.annotsvMode
    outputFile = "${prefix}_AnnotSV.tsv"
    
    //Pass any additional flags to the AnnotSV 
    //def extraArgs = params.extraAnnotsvFlags ?: ''
    """
    AnnotSV \\
        -SVinputFile ${params.annotsvinput} \\
        -annotationsDir ${params.annotsvAnnotationsDir} \\
        -bedtools bedtools \\
        -bcftools bcftools \\
        -annotationMode ${params.annotsvMode} \\
        -genomeBuild ${params.annotsvGenomeBuild} \\
        -includeCI 1 \\
        -overwrite 1 \\
        -outputFile ${outputFile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(AnnotSV --version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}
