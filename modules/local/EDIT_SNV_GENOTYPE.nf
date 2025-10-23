process EDIT_SNV_GENOTYPE {
    tag "$meta.id:$chrom"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'jiangyanyu/pacbio_wgs:v1.3' :
        'jiangyanyu/pacbio_wgs:v1.3' }"

    input:
        tuple val(meta), path(snv_vcf_file)
        tuple val(meta), path(snv_tbi_file)
        tuple val(meta), path(sv_vcf_file),
        tuple val(meta), path(sv_tbi_file)

    output:
        tuple val(meta), path("${meta.id}.${caller}.${chrom}.edited_gt.vcf.gz"), emit: vcf
        tuple val(meta), path("${meta.id}.${caller}.${chrom}.edited_gt.vcf.gz.tbi"), emit: tbi
        path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core_nano-fanconi/bin/


    """
    SNV_modify_GT.py \\
        --snv_vcf ${snv_vcf_file} \\
        --sv_vcf ${sv_vcf_file} \\
        --output_vcf ${meta.id}.${caller}.${chrom}.edited_gt.vcf


    bcftools view \\
        --threads ${task.cpus} \\
        -O z \\
        -o ${meta.id}.${caller}.${chrom}.edited_gt.vcf.gz \\
        ${meta.id}.${caller}.${chrom}.edited_gt.vcf


    tabix ${meta.id}.${caller}.${chrom}.edited_gt.vcf.gz
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
