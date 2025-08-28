process SPLIT_VCF_BY_CHR {
    label 'process_medium'

    conda "bioconda::bcftools=1.16"
    container "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"

    input:
        tuple val(meta), path(deepvariant_vcf)

    output:
        tuple val(meta), path("chr*.vcf.gz"), emit: split_vcfs
        path("versions.yml"), emit: versions

    script:
        def chroms = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]
        def split_cmds = chroms.collect { chr ->
            "bcftools view -r ${chr} ${deepvariant_vcf} -O z -o ${meta.sample}.${chr}.vcf.gz"
        }.join("\n")
        """
        ${split_cmds}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: $(bcftools --version | head -n1 | sed 's/^bcftools //')
        END_VERSIONS
        """
}
