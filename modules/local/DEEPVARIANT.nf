// Function to determine the label
def determineLabel() {
    return params.use_gpu ? 'process_gpu_long' : 'process_high'
}
def processLabel = determineLabel()

process DEEPVARIANT {
    tag "$meta.sample"
    label processLabel

    container "google/deepvariant:1.8.0-gpu"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(phased_bam_file), path(phased_bai_file)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("${prefix}.unfiltered.vcf.gz")  ,  emit: vcf
    tuple val(meta), path("${prefix}.unfiltered.g.vcf.gz"),  emit: gvcf
    path "versions.yml"                        ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.sample}"
    //def regions = intervals ? "--regions ${intervals}" : ""

    """

    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${fasta} \\
        --reads=${prefix}.sorted.bam \\
        --output_vcf=${prefix}.unfiltered.vcf.gz \\
        --output_gvcf=${prefix}.unfiltered.g.vcf.gz \\
        ${args} \\
        --num_shards=${task.cpus} 

    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
