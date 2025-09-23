process SPLIT_VCF_BY_CHR {
    label 'process_medium'

    conda "bioconda::bcftools=1.16"
    container "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"

    input:
        tuple val(meta), path(deepvariant_vcf)
        val(chr)

    output:
        tuple val(meta), path("chr*.vcf.gz"), emit: split_vcfs
        path("versions.yml"), emit: versions

    script:
        """
        bcftools view -r ${chr} ${deepvariant_vcf} -O z -o ${meta.sample}.${chr}.vcf.gz
        
        bcftools index ${meta.sample}.${chr}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: $(bcftools --version | head -n1 | sed 's/^bcftools //')
        END_VERSIONS
        """
}
