/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNanofanconi.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

 ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
 ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
 ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
 ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


include { FAST5_TO_POD5                                 } from '../modules/local/FAST5_TO_POD5'
include { DORADO_BASECALLER                             } from '../modules/local/DORADO_BASECALLER'
include { MERGE_BASECALL as MERGE_BASECALL_ID           } from '../modules/local/MERGE_BASECALL'
include { MERGE_BASECALL as MERGE_BASECALL_SAMPLE       } from '../modules/local/MERGE_BASECALL'
include { DORADO_BASECALL_SUMMARY                       } from '../modules/local/DORADO_BASECALL_SUMMARY'
include { PYCOQC                                        } from '../modules/local/PYCOQC'
include { DORADO_ALIGNER                                } from '../modules/local/DORADO_ALIGNER'
include { SAMTOOLS_SORT                                 } from '../modules/local/SAMTOOLS_SORT'
include { SAMTOOLS_INDEX                                } from '../modules/local/SAMTOOLS_INDEX'
include { SAMTOOLS_STATS                                } from '../modules/local/SAMTOOLS_STATS.nf'
include { SNIFFLES                                      } from '../modules/local/SNIFFLES.nf'
include { BCFTOOLS_SORT as SNIFFLES_SORT_VCF            } from '../modules/nf-core/bcftools/sort/main.nf'
include { TABIX_BGZIP as SNIFFLES_BGZIP_VCF             } from '../modules/nf-core/tabix/bgzip/main.nf'
include { TABIX_TABIX as SNIFFLES_TABIX_VCF             } from '../modules/nf-core/tabix/tabix/main.nf'
include { WHATSHAP                                      } from '../modules/local/WHATSHAP.nf'
include { MOSDEPTH                                      } from '../modules/local/MOSDEPTH.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                                       } from '../modules/local/MULTIQC'
// include { PEPPER                                        } from '../modules/local/PEPPER'
// include { MODKIT                                        } from '../modules/local/MODKIT'
// include { MODKIT_TO_BW                                  } from '../modules/local/MODKIT_TO_BW'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NANOFANCONI {

    ch_versions = Channel.empty()
    

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOFANCONI: input check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_phased_vcf = INPUT_CHECK.out.reads.map{ meta, files -> [[sample: meta.sample],meta.fast5_path, meta.vcf, meta.vcf_tbi] }.dump(tag: "ch_phased_vcf")

if (params.reads_format == 'bam' ) {
    INPUT_CHECK
    .out
    .reads
    .flatMap { meta, files -> 
        def bam_path = meta.fast5_path
        def bam_files = []
        if (file(bam_path).isDirectory()) {
            bam_files = file("${bam_path}/*.bam")
        } else if (bam_path.endsWith('.bam')) {
            bam_files = [file(bam_path)]
        }
        bam_files.collect { [[sample: meta.sample], it] }  // Create a list of [meta, file] pairs
    }
    .groupTuple(by: 0) // group bams by meta (i.e sample) which is zero-indexed
    // .dump(tag: 'basecall_sample', pretty: true)
    .set { ch_basecall_sample_merged_bams } // set channel name
}

    MERGE_BASECALL_SAMPLE (
        ch_basecall_sample_merged_bams
    )
    ch_versions = ch_versions.mix(MERGE_BASECALL_SAMPLE.out.versions)


    SAMTOOLS_SORT (
         MERGE_BASECALL_SAMPLE.out.merged_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    
    

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOFANCONI: Sniffles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (params.run_sniffles) {

        /*
         * Call structural variants with sniffles
         */
        SNIFFLES( SAMTOOLS_SORT.out.bam )
        ch_versions = ch_versions.mix(SNIFFLES.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        SNIFFLES_SORT_VCF( SNIFFLES.out.sv_calls )
        ch_sv_calls_vcf = SNIFFLES_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(SNIFFLES_SORT_VCF.out.versions)

        /*
         * Index sniffles vcf.gz
         */
        SNIFFLES_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = SNIFFLES_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(SNIFFLES_TABIX_VCF.out.versions)

    }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOFANCONI: whatshap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    if (params.run_whatshap) {
        //
        // MODULE: Index PEPPER bam
        //
        sample_meta = INPUT_CHECK.out.reads.map{ meta, files -> [[sample: meta.sample]] }.dump(tag: "sample_meta")
        WHATSHAP (
            sample_meta, SNIFFLES_SORT_VCF.out.vcf, SAMTOOLS_SORT.out.bam, file(params.fasta), file(params.fasta_index)
        )
        ch_versions = ch_versions.mix(WHATSHAP.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOFANCONI: whatshap depth calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

        //
        // MODULE: MOSDEPTH for depth calculation
        //
        ch_mosdepth_input = WHATSHAP.out.bam.mix(WHATSHAP.out.bai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
        MOSDEPTH (
            ch_mosdepth_input
        )
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
        
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOFANCONI: currently remove methylation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
        //
        // MODULE: MODKIT to extract methylation data
        //
        
        if (params.extract_methylation) {
            ch_modkit_input = WHATSHAP.out.bam.mix(WHATSHAP.out.bai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
            MODKIT (
                ch_modkit_input,
                file(params.fasta)
            )
            ch_versions = ch_versions.mix(MODKIT.out.versions)

            ch_modkit_to_bw_input = MODKIT.out.hap1_bed.join(MODKIT.out.hap2_bed).join(MODKIT.out.combined_bed)
            MODKIT_TO_BW (
                ch_modkit_to_bw_input,
                file(params.fasta_index)
            )
            ch_versions = ch_versions.mix(MODKIT_TO_BW.out.versions)

            SAMTOOLS_STATS (
                WHATSHAP.out.bam
            )
            ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
        }
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOFANCONI: CUSTOM_DUMPSOFTWAREVERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowNanofanconi.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowNanofanconi.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    if (params.reads_format == 'fast5' || params.reads_format == 'pod5') {
        ch_multiqc_files = ch_multiqc_files.mix(PYCOQC.out.json.collect{it[1]}.ifEmpty([]))
    }

    if (params.run_whatshap) {
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_bed.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_csi.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_bed.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_csi.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/