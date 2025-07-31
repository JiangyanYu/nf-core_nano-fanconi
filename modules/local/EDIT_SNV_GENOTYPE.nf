process EDIT_SNV_GENOTYPE {
    tag "$meta.sample"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'mdivr/snpdep:v1' :
        'mdivr/snpdep:v1' }"

    input:
    path SNV
    path SV

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core_nano-fanconi/bin/
    """
    // first use SV.vcf to define the large deletion breakpoint genomic positions
    // then python editsnvgt.py 
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
