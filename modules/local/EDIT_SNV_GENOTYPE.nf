process EDIT_SNV_GENOTYPE {
    // tag "$meta.sample"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.3' :
        'jiangyanyu/pacbio_wgs:v1.3' }"

    input:
        // tuple val(meta), path(snv_vcf_file), path(snv_tbi_file),path(sv_vcf_file), path(sv_tbi_file)
        tuple val(meta), path(sv_vcf_file), path(sv_tbi_file)

    output:
        tuple val(meta), path("${meta.sample}_gt.converted.vcf")       , emit: vcf
        path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core_nano-fanconi/bin/


    """
    mkdir test

    python SNV_modify_GT.py \\
        --snv_vcf ${meta.sample}_filtered.vcf.gz  \\
        --sv_vcf genotyped.sv.vcf.gz \\
        --output_vcf ${meta.sample}_gt.converted.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
