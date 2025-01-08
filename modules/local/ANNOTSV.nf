process ANNOTSV {
    tag "$meta.sample"
    label 'process_high'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/annotsv:3.4.4--py312hdfd78af_0' :
        'quay.io/biocontainers/annotsv:3.4.4--py312hdfd78af_0' }"

	input:
	tuple val(meta), path(vcf_file)
	path(annotsvDir)
	val(annotsvMode)

	output:
    tuple val(meta), path("${prefix}*.tsv"),  emit: tsv
    path "versions.yml"                    ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.sample}"
	
	// Change output file name based on annotation mode
	def outputFile = null
	    if (mode == 'full') {
               outputFile = "${prefix}_full_AnnotSV.tsv"
            } else if (mode == 'split') {
               outputFile = "${prefix}_split_AnnotSV.tsv"
            } else if (mode == 'both') {
               outputFile = "${prefix}_both_AnnotSV.tsv"
            } else {
               throw new RuntimeException("Invalid option for --annotSV: ${mode}")}
	
	//Pass any additional flags to the AnnotSV 
	def extraArgs = params.extraAnnotsvFlags ?: ''
	"""
	AnnotSV \\
		-SVinputFile ${vcf_file} \\
		-annotationsDir ${annotsvDir} \\
		-bedtools bedtools -bcftools bcftools \\
		-annotationMode ${annotsvMode} \\
		-genomeBuild ${params.genomeBuild} \\
		-includeCI 1 \\
		-overwrite 1 \\
		-outputFile ${outputFile} ${extraArgs}
		
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    
	"""
}
