process EDIT_SNV_GENOTYPE {
    tag "$meta.id:$chrom"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.3' :
        'jiangyanyu/pacbio_wgs:v1.3' }"

    // conda "bioconda::bcftools=1.16"
    // container "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"


    input:
        tuple val(meta), path(snv_vcf), path(snv_tbi), val(snv_caller), path(sv_vcf), path(sv_tbi), val(sv_caller), val(chrom)
        path regions_csv
        
    output:
        tuple val(meta), path("${meta.id}.${snv_caller}.${chrom}.edited_gt.vcf"), val(snv_caller), val(chrom), emit: vcf
        path "versions.yml"                                                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core_nano-fanconi/bin/


    """
    SNV_modify_GT.py \\
        --snv_vcf ${snv_vcf} \\
        --sv_vcf ${sv_vcf} \\
        --regions_csv ${regions_csv} \\
        --output_vcf ${meta.id}.deepvariant.${chrom}.edited_gt.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
