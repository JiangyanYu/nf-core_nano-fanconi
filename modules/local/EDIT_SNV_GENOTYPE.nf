process EDIT_SNV_GENOTYPE {
    tag "$meta.id:$chrom"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.3' :
        'jiangyanyu/pacbio_wgs:v1.3' }"

    // conda "bioconda::bcftools=1.16"
    // container "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"


    input:
        tuple val(meta), path(snv_vcf_file), path(snv_tbi_file), path(sv_vcf_file), path(sv_tbi_file), val(chrom)
        
    output:
        tuple val(meta), path("${meta.id}.deepvariant.${chrom}.edited_gt.vcf"), val("deepvariant"), val(chrom), emit: vcf
        path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core_nano-fanconi/bin/


    """
    SNV_modify_GT.py \\
        --snv_vcf ${meta.id}.deepvariant.${chrom}.vcf \\
        --sv_vcf ${meta.id}.sawfish.${chrom}.vcf \\
        --output_vcf ${meta.id}.deepvariant.${chrom}.edited_gt.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
