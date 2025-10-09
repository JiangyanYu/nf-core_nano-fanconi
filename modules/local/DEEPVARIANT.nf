// Function to determine the label
def determineLabel() {
    return params.use_gpu ? 'process_gpu_long' : 'process_high'
}
def processLabel = determineLabel()

process DEEPVARIANT {
    tag "$meta.id"
    label processLabel

    container "google/deepvariant:1.9.0-gpu"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
        tuple val(meta), path(cram)
        tuple val(meta), path(crai)
        path(fasta)
        path(fai)

    output:
        tuple val(meta), path("${meta.id}.deepvariant.unfiltered.vcf.gz"), emit: vcf
        tuple val(meta), path("${meta.id}.deepvariant.unfiltered.vcf.gz.tbi"), emit: tbi
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        /opt/deepvariant/bin/run_deepvariant \\
            --model_type=ONT_R104 \\
            --ref=${fasta} \\
            --reads=${cram} \\
            --output_vcf=${meta.id}.deepvariant.unfiltered.vcf.gz \\
            --num_shards=${task.cpus} 

        tabix ${meta.id}.deepvariant.unfiltered.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
        END_VERSIONS
        """
}
