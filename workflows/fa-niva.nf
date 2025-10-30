nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
// WorkflowFaniva.initialise(params, log)

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

include { FAIDX_REFERENCE                                      } from '../modules/local/FAIDX_REFERENCE.nf'
include { PBMM2_INDEX_REFERENCE                                } from '../modules/local/PBMM2_INDEX_REFERENCE.nf'
include { SAMTOOLS_BGZIP                                       } from '../modules/nf-core/samtools/bgzip.nf'
include { SAMTOOLS_FAIDX                                       } from '../modules/nf-core/samtools/faidx.nf'
include { FAST5_TO_POD5                                        } from '../modules/local/FAST5_TO_POD5.nf'
include { DORADO_BASECALLER                                    } from '../modules/local/DORADO_BASECALLER.nf'
include { DORADO_BASECALL_SUMMARY                              } from '../modules/local/DORADO_BASECALL_SUMMARY.nf'
include { PYCOQC                                               } from '../modules/local/PYCOQC.nf'
include { PBMM2_FROM_BAM                                       } from '../modules/local/PBMM2_FROM_BAM.nf'
// include { SAMTOOLS_STATS                                } from '../modules/local/SAMTOOLS_STATS.nf'
include { EXTRACT_LOW_MG_FROM_CRAM                             } from '../modules/local/EXTRACT_LOW_MG_FROM_CRAM.nf'
include { SPLIT_CRAM_BY_CHROM                                  } from '../modules/local/SPLIT_CRAM_BY_CHROM.nf'
include { DEEPVARIANT                                          } from '../modules/local/DEEPVARIANT.nf'
include { SAWFISH                                              } from '../modules/local/SAWFISH.nf'
include { SPLIT_VCF_BY_CHROM as SPLIT_VCF_BY_CHROM_DEEPVARIANT } from '../modules/local/SPLIT_VCF_BY_CHROM.nf'
include { SPLIT_VCF_BY_CHROM as SPLIT_VCF_BY_CHROM_SAWFISH     } from '../modules/local/SPLIT_VCF_BY_CHROM.nf'
include { EDIT_SNV_GENOTYPE                                    } from '../modules/local/EDIT_SNV_GENOTYPE.nf'
// include { BCFTOOLS_SORT as SNIFFLES_SORT_VCF            } from '../modules/nf-core/bcftools/sort/main.nf'
// include { TABIX_BGZIP as SNIFFLES_BGZIP_VCF             } from '../modules/nf-core/tabix/bgzip/main.nf'
// include { TABIX_TABIX as SNIFFLES_TABIX_VCF             } from '../modules/nf-core/tabix/tabix/main.nf'
// include { ANNOTSV_SAWFISH                               } from '../modules/local/ANNOTSV_SAWFISH.nf'
// include { ANNOTSV_DEEPVARIANT                           } from '../modules/local/ANNOTSV_DEEPVARIANT.nf'
// include { BCFTOOLS_FILTER as DEEPVARIANT_FILTER_VCF     } from '../modules/nf-core/bcftools/filter/main.nf'
// include { TABIX_BGZIP as EDIT_SNV_GENOTYPE_BGZIP_VCF    } from '../modules/nf-core/tabix/bgzip/main.nf'
// include { TABIX_TABIX as EDIT_SNV_GENOTYPE_TABIX_VCF    } from '../modules/nf-core/tabix/tabix/main.nf'
// include { WHATSHAP_PHASE                                } from '../modules/local/WHATSHAP_PHASE.nf'
// include { WHATSHAP_HAPLOTAG                             } from '../modules/local/WHATSHAP_HAPLOTAG.nf'
// include { BCFTOOLS_SORT as PHASE_SORT_VCF               } from '../modules/nf-core/bcftools/sort/main.nf'
// include { TABIX_BGZIP as PHASE_BGZIP_VCF                } from '../modules/nf-core/tabix/bgzip/main.nf'
// include { TABIX_TABIX as PHASE_TABIX_VCF                } from '../modules/nf-core/tabix/tabix/main.nf'
// include { MOSDEPTH                                      } from '../modules/local/MOSDEPTH.nf'
// include { TABIX_TABIX as DEEPVARIANT_TABIX_VCF          } from '../modules/nf-core/tabix/tabix/main.nf'
// include { TABIX_TABIX as DEEPVARIANT_TABIX_GVCF         } from '../modules/nf-core/tabix/tabix/main.nf'
// include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
// include { MULTIQC                                       } from '../modules/local/MULTIQC.nf'
// include { PEPPER                                        } from '../modules/local/PEPPER'
// include { MODKIT                                        } from '../modules/local/MODKIT'
// include { MODKIT_TO_BW                                  } from '../modules/local/MODKIT_TO_BW'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Define chromosome names
// This should be extracted from the reference fai index file
def chroms = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]
        

// Info required for completion email and summary
def multiqc_report = []

workflow FANIVA {

    ch_versions = Channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FANIVA: FAIDX_REFERENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FAIDX_REFERENCE(

        file(params.fasta).toRealPath()
    )
    ch_fasta = FAIDX_REFERENCE.out.fasta
    ch_fasta_index = FAIDX_REFERENCE.out.fasta_index
    ch_versions = ch_versions.mix(FAIDX_REFERENCE.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: pbmm2 index reference
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    PBMM2_INDEX_REFERENCE(

        ch_fasta,
        ch_fasta_index
    )
    ch_fasta_mmi = PBMM2_INDEX_REFERENCE.out.mmi
    ch_versions = ch_versions.mix(PBMM2_INDEX_REFERENCE.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: input check
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (

        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_phased_vcf = INPUT_CHECK.out.reads.map{ meta, files -> [[sample: meta.sample],meta.vcf] }.dump(tag: "ch_phased_vcf")


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FANIVA: fast5-pod5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // fast5 input
    if (params.reads_format == 'fast5') {
        INPUT_CHECK
        .out
        .reads
        .map { meta, files -> 
            def fast5_path = meta.input_path
    
            // Check if fast5_path is null or empty
            if (!fast5_path) {
                throw new IllegalArgumentException("fast5_path is null or empty")
            }
        
            def fast5_files = []
            
            // TO DO: provide raw.github link to download files automatically
    
            if (file(fast5_path).isDirectory()) {
                fast5_files = file("${fast5_path}/*.fast5")
            } else if (fast5_path.endsWith('.fast5')) {
                fast5_files = [file(fast5_path)]
            }
            
            [meta, fast5_files]
        }
        .flatMap { meta, files ->
            def chunks = files.toList().collate(params.dorado_files_chunksize)  // chunk files into groups of 2
            def chunkList = []
            for (int i = 0; i < chunks.size(); i++) {
                def newMeta = meta.clone()  // clone the meta to avoid modifying the original
                newMeta.chunkNumber = i + 1  // add chunk number, starting from 1
                chunkList << [newMeta, chunks[i]]
            }
            return chunkList
        }
        // .dump(tag: 'input', pretty: true)
        .set { ch_fast5 }

    FAST5_TO_POD5 (
        ch_fast5
    )

    FAST5_TO_POD5
    .out
    .pod5
    .set { ch_pod5 } 

    ch_versions = ch_versions.mix(FAST5_TO_POD5.out.versions)

    } else if (params.reads_format == 'pod5') {
    INPUT_CHECK
    .out
    .reads
    .map { meta, files -> 
        def pod5_path = meta.input_path

        // Check if fast5_path is null or empty
            if (!pod5_path) {
                throw new IllegalArgumentException("pod5_path is null or empty")
            }

        def pod5_files = []

        if (file(pod5_path).isDirectory()) {
            pod5_files = file("${pod5_path}/*.pod5")
        } else if (pod5_path.endsWith('.pod5')) {
            pod5_files = [file(pod5_path)]
        }
        [meta, pod5_files]
    }
    .flatMap { meta, files ->
        def chunks = files.toList().collate(params.dorado_files_chunksize)  // chunk files into groups of 2
        def chunkList = []
        for (int i = 0; i < chunks.size(); i++) {
            def newMeta = meta.clone()  // clone the meta to avoid modifying the original
            newMeta.chunkNumber = i + 1  // add chunk number, starting from 1
            chunkList << [newMeta, chunks[i]]
        }
        return chunkList
    }
    // .dump(tag: 'input_pod5', pretty: true)
    .set { ch_pod5 }
    }


    if (params.reads_format == 'pod5' || params.reads_format == 'fast5') {
        DORADO_BASECALLER (
            ch_pod5
        )
        ch_versions = ch_versions.mix(DORADO_BASECALLER.out.versions)
        DORADO_BASECALLER
        .out
        .bam
        .map { meta, bam -> [[id: meta.id, sample: meta.sample, flowcell: meta.flowcell, batch: meta.batch, kit: meta.kit] , bam]} // make sample name the only meta (remove flow cell and other info)
        .groupTuple(by: 0) // group bams by meta (i.e sample) which zero indexed
        // .dump(pretty: true)
        //.set { ch_basecall_single_bams }
        .set { ch_unmapped_bams }

        // Dorado basecall summary
        DORADO_BASECALL_SUMMARY (
            
            //ch_basecall_single_bams
            ch_unmapped_bams
        )

        DORADO_BASECALL_SUMMARY
        .out
        .summary
        // .dump(pretty: true)
        .set { ch_basecall_summary }


        // MODULE: PycoQC (QC from Basecall results)
        PYCOQC (
            ch_basecall_summary
        )
        ch_versions = ch_versions.mix(PYCOQC.out.versions)

    }

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: pbmm2 alignment from BAM
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // PBMM2_FROM_BAM (

    //     ch_basecall_single_bams,
    //     ch_fasta,
    //     ch_fasta_index,
    //     ch_fasta_mmi
    // )
    // ch_pbmm2_cram = PBMM2_FROM_BAM.out.cram
    // ch_pbmm2_crai = PBMM2_FROM_BAM.out.crai
    // ch_versions = ch_versions.mix(PBMM2_FROM_BAM.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: Manage if reads_format is bam
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    if (params.reads_format == 'bam' ) {
        INPUT_CHECK
        .out
        .reads
        .flatMap { meta, files -> 
            def bam_path = meta.input_path
            def bam_files = []
            if (file(bam_path).isDirectory()) {
                bam_files = file("${bam_path}/*.bam")
            } else if (bam_path.endsWith('.bam')) {
                bam_files = [file(bam_path)]
            }
            bam_files.collect { [[id: meta.id], it] }  // Create a list of [meta, file] pairs
        }
        .groupTuple(by: 0) // group bams by meta (i.e sample) which is zero-indexed
        .map { meta, files -> 
            // Remove duplicates from the file list
            def unique_files = files.unique()
            [meta, unique_files]
        }
        .set { ch_unmapped_bams } // set channel name

    }


    // /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: pbmm2 alignment from BAM
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    PBMM2_FROM_BAM (

        ch_unmapped_bams,
        ch_fasta,
        ch_fasta_index,
        ch_fasta_mmi
    )
    // ch_pbmm2_cram = PBMM2_FROM_BAM.out.cram
    // ch_pbmm2_crai = PBMM2_FROM_BAM.out.crai
    ch_pbmm2_cram_crai = PBMM2_FROM_BAM.out.cram_crai
    ch_versions = ch_versions.mix(PBMM2_FROM_BAM.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: EXTRACT_LOW_MG_FROM_CRAM
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // Extract low mapping quality reads from CRAM
    EXTRACT_LOW_MG_FROM_CRAM(

        ch_pbmm2_cram_crai,
        ch_fasta,
    )
    ch_pbmm2_low_mg_cram_crai = EXTRACT_LOW_MG_FROM_CRAM.out.cram_crai
    ch_versions = ch_versions.mix(EXTRACT_LOW_MG_FROM_CRAM.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: SPLIT_CRAM_BY_CHROM
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // Split CRAM by chromosome
    ch_pbmm2_cram_crai
        .flatMap { meta, cram, crai ->
            chroms.collect { chr ->
                [meta, cram, crai, chr]
            }
        }
        .set { ch_pbmm2_cram_crai_chrom }


    // Split CRAM by chromosome
    SPLIT_CRAM_BY_CHROM(

        ch_pbmm2_cram_crai_chrom,
        ch_fasta
    )
    ch_split_by_chrom_cram_crai = SPLIT_CRAM_BY_CHROM.out.cram_crai
    ch_versions = ch_versions.mix(SPLIT_CRAM_BY_CHROM.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: DeepVariant
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
               
    DEEPVARIANT( 

        ch_pbmm2_cram_crai,
        ch_fasta,
        ch_fasta_index
    )  
    ch_deepvariant_vcf_tbi  = DEEPVARIANT.out.vcf_tbi
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: Sawfish
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    SAWFISH(

        ch_pbmm2_cram_crai,
        ch_fasta,
        ch_fasta_index
    )
    ch_sawfish_vcf_tbi  = SAWFISH.out.vcf_tbi
    ch_versions = ch_versions.mix(SAWFISH.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: AnnotSV-Sawfish
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

//         if (params.run_annotsv) {
    
//             ANNOTSV_SAWFISH (
//                 SAWFISH.out.vcf
//             )

//             ch_versions = ch_versions.mix(ANNOTSV_SAWFISH.out.versions)
        
//         }
        
//     }
    

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: SPLIT_VCF_BY_CHROM for deepvariant
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // Split VCF by chromosome
    ch_deepvariant_vcf_tbi
        .flatMap { meta, vcf, tbi, caller ->
            chroms.collect { chr ->
                [meta, vcf, tbi, caller, chr]
            }
        }
        .set { ch_deepvariant_vcf_tbi_chrom }


    // Split VCF by chromosome
    SPLIT_VCF_BY_CHROM_DEEPVARIANT(

        ch_deepvariant_vcf_tbi_chrom
    )
    ch_deepvariant_split_by_chrom_vcf_tbi = SPLIT_VCF_BY_CHROM_DEEPVARIANT.out.vcf_tbi
    ch_versions = ch_versions.mix(SPLIT_VCF_BY_CHROM_DEEPVARIANT.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: SPLIT_VCF_BY_CHROM for sawfish
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // Split VCF by chromosome
    ch_sawfish_vcf_tbi
        .flatMap { meta, vcf, tbi, caller ->
            chroms.collect { chr ->
                [meta, vcf, tbi, caller, chr]
            }
        }
        .set { ch_sawfish_vcf_tbi_chrom }


    // Split VCF by chromosome
    SPLIT_VCF_BY_CHROM_SAWFISH(

        ch_sawfish_vcf_tbi_chrom
    )
    ch_sawfish_split_by_chrom_vcf_tbi = SPLIT_VCF_BY_CHROM_SAWFISH.out.vcf_tbi
    ch_versions = ch_versions.mix(SPLIT_VCF_BY_CHROM_SAWFISH.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: EDIT_SNV_GENOTYPE - Match DeepVariant and Sawfish by chromosome (FIXED)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // Debug the input channels first
    ch_deepvariant_split_by_chrom_vcf_tbi
        .view { "DeepVariant input: ${it}" }
        .set { ch_deepvariant_debug }
    
    ch_sawfish_split_by_chrom_vcf_tbi
        .view { "Sawfish input: ${it}" }
        .set { ch_sawfish_debug }

    // Match DeepVariant and Sawfish VCF files by sample and chromosome
    ch_deepvariant_debug
        .map { meta, vcf, tbi, caller, chrom -> 
            // Create join key: [sample_id, chromosome]
            [[meta.id, chrom], [meta, vcf, tbi, caller, chrom]] 
        }
        .join(
            ch_sawfish_debug.map { meta, vcf, tbi, caller, chrom -> 
                // Create matching join key: [sample_id, chromosome]
                [[meta.id, chrom], [vcf, tbi, caller]] 
            }, 
            by: 0
        )
        .map { key, deepvariant_data, sawfish_data ->
            // Reconstruct tuple: extract data from nested lists
            def meta = deepvariant_data[0]
            def snv_vcf = deepvariant_data[1]
            def snv_tbi = deepvariant_data[2] 
            def snv_caller = deepvariant_data[3]
            def chrom = deepvariant_data[4]
            
            def sv_vcf = sawfish_data[0]
            def sv_tbi = sawfish_data[1]
            def sv_caller = sawfish_data[2]
            
            [meta, snv_vcf, snv_tbi, snv_caller, sv_vcf, sv_tbi, sv_caller, chrom]
        }
        .view { meta, snv_vcf, snv_tbi, snv_caller, sv_vcf, sv_tbi, sv_caller, chrom ->
            "DEBUG MATCHED: Processing ${meta.id} ${chrom} - DV:${snv_vcf.name} SF:${sv_vcf.name}"
        }
        .set { ch_matched_vcf_by_chrom }

    // Count how many matched tuples we have
    ch_matched_vcf_by_chrom
        .count()
        .view { "Total matched chromosomes for EDIT_SNV_GENOTYPE: $it" }

    // Create the channel again for actual processing (since count() consumes it)
    ch_deepvariant_debug
        .map { meta, vcf, tbi, caller, chrom -> 
            [[meta.id, chrom], [meta, vcf, tbi, caller, chrom]] 
        }
        .join(
            ch_sawfish_debug.map { meta, vcf, tbi, caller, chrom -> 
                [[meta.id, chrom], [vcf, tbi, caller]] 
            }, 
            by: 0
        )
        .map { key, deepvariant_data, sawfish_data ->
            def meta = deepvariant_data[0]
            def snv_vcf = deepvariant_data[1]
            def snv_tbi = deepvariant_data[2] 
            def snv_caller = deepvariant_data[3]
            def chrom = deepvariant_data[4]
            
            def sv_vcf = sawfish_data[0]
            def sv_tbi = sawfish_data[1]
            def sv_caller = sawfish_data[2]
            
            [meta, snv_vcf, snv_tbi, snv_caller, sv_vcf, sv_tbi, sv_caller, chrom]
        }
        .set { ch_matched_vcf_by_chrom_final }

    // Load SNV_modify_regions csv file
    ch_SNV_modify_regions = Channel.of(file(params.SNV_modify_regions))

    EDIT_SNV_GENOTYPE(
        ch_matched_vcf_by_chrom_final,
        ch_SNV_modify_regions
    )
    ch_deepvariant_vcf_chrom_edited_gt = EDIT_SNV_GENOTYPE.out.vcf
    ch_versions = ch_versions.mix(EDIT_SNV_GENOTYPE.out.versions)

}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: WHATSHAP_PHASE
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
//     // Run WHATSHAP_PHASE on each split VCF
//     WHATSHAP_PHASE(
//         ch_split_vcfs.map { meta, vcf_path -> [meta, vcf_path, ch_pbmm2_cram, file(params.fasta), file(params.fasta_index)] }
//     )

//     // Collect phased VCFs from WHATSHAP_PHASE
//     ch_phased_chr_vcfs = WHATSHAP_PHASE.out.phased_vcf.collect()

//     // Merge chromosome VCFs into a single VCF
//     MERGE_VCF(
//         tuple(val('meta'), ch_phased_chr_vcfs)
//     )
//     MERGE_VCF
//     .out
//     .merged_vcf
//     .set { ch_merged_phased_vcf }
//     // Merge chromosome VCFs into a single VCF
//     MERGE_VCF(
//         [ 'meta', ch_phased_chr_vcfs ]
//     )
//     MERGE_VCF
//     .out
//     .merged_vcf
//     .set { ch_merged_phased_vcf }
//     ch_multiqc_files = ch_multiqc_files.mix(ch_final_haplotagged_vcf.collect{it[1]}.ifEmpty([]))

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: WHATSHAP_HAPLOTAG
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
//     // Run WHATSHAP_HAPLOTAG on the final merged, haplotagged VCF
//     WHATSHAP_HAPLOTAG(
//         ch_pbmm2_cram,
//         ch_final_haplotagged_vcf.map{ meta, vcf -> vcf },
//         ch_fasta,
//         ch_fasta_index
//     )
//     ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: SAMTOOLS_STATS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
//     SAMTOOLS_STATS (
//         WHATSHAP_HAPLOTAG.out.cram
//     )
//     ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: AnnotSV-Deepvariant
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

//         if (params.run_annotsv) {
    
//             ANNOTSV_DEEPVARIANT (
//                DEEPVARIANT.out.vcf
//             )

//         ch_versions = ch_versions.mix(ANNOTSV_DEEPVARIANT.out.versions)
        
//         }
        
//     }




// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: whatshap
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */


//     // if (params.run_whatshap) {
//     //     //
//     //     // MODULE: whatshap for phasing
//     //     //


//     //     if (params.joint_SNV_SV_phasing) {
//     //         //
//     //         // Edit SNV genotype when large deletion occurs
//     //         //

//     //         add_meta1 = INPUT_CHECK.out.reads
//     //             .map{ meta, files -> [[sample: meta.sample], meta.vcf] }
//     //             .dump(tag: "add_meta1")

//     //         ch_snv_vcf = DEEPVARIANT_FILTER_VCF.out.filteredvcf
//     //             .mix(DEEPVARIANT_TABIX_VCF.out.tbi)
//     //             .groupTuple(size:2)
//     //             .map{ meta, files -> [ meta, files.flatten() ]}
//     //          deepvariant_vcf = ch_snv_vcf.join(add_meta1).dump(tag: "joined")


//     //         add_meta2 = INPUT_CHECK.out.reads
//     //             .map{ meta, files -> [[sample: meta.sample], meta.vcf_tbi] }
//     //             .dump(tag: "add_meta2")

//     //         ch_sv_vcf = SAWFISH.out.vcf
//     //             .mix(SAWFISH.out.tbi)
//     //             .groupTuple(size:2)
//     //             .map{ meta, files -> [ meta, files.flatten() ]}
//     //         sawfish_vcf = ch_sv_vcf.join(add_meta2).dump(tag: "joined")


//     //         EDIT_SNV_GENOTYPE (
//     //             deepvariant_vcf,
//     //             sawfish_vcf
//     //         )

//     //         ch_versions = ch_versions.mix(EDIT_SNV_GENOTYPE.out.versions)

//     //         /*
//     //         * Zip and Index edited vcf file
//     //         */
//     //         EDIT_SNV_GENOTYPE_BGZIP_VCF ( EDIT_SNV_GENOTYPE.out.vcf )
//     //         ch_versions = ch_versions.mix( EDIT_SNV_GENOTYPE_BGZIP_VCF.out.versions)

//     //         EDIT_SNV_GENOTYPE_TABIX_VCF ( EDIT_SNV_GENOTYPE_BGZIP_VCF.out.output )
//     //         ch_versions = ch_versions.mix( EDIT_SNV_GENOTYPE_TABIX_VCF.out.versions)



//     //         //
//     //         // Input converted snv.vcf for phasing
//     //         //


//     //         add_meta3 = INPUT_CHECK.out.reads
//     //             .map{ meta, files -> [[sample: meta.sample],meta.vcf_tbi] }
//     //             .dump(tag: "add_meta3")


//     //         ch_phase_vcf = EDIT_SNV_GENOTYPE_BGZIP_VCF.out.output
//     //             .mix(EDIT_SNV_GENOTYPE_TABIX_VCF.out.tbi)
//     //             .groupTuple(size:2)
//     //             .map{ meta, files -> [ meta, files.flatten() ]}
//     //         phase_vcf = ch_phase_vcf.join(add_meta3).dump(tag: "joined")

                       
//     //     } else {
            


//     //         add_meta4 = INPUT_CHECK.out.reads
//     //             .map{ meta, files -> [[sample: meta.sample],meta.vcf_tbi] }
//     //             .dump(tag: "add_meta4")


//     //         ch_phase_vcf = DEEPVARIANT_FILTER_VCF.out.filteredvcf
//     //             .mix(DEEPVARIANT_TABIX_VCF.out.tbi)
//     //             .groupTuple(size:2)
//     //             .map{ meta, files -> [ meta, files.flatten() ]}
//     //         phase_vcf = ch_phase_vcf.join(add_meta4).dump(tag: "joined")

//     //     }

//     //     ch_phase_cram = SAMTOOLS_SORT.out.cram
//     //         .mix(SAMTOOLS_SORT.out.crai)
//     //         .groupTuple(size:2)
//     //         .map{ meta, files -> [ meta, files.flatten() ]}

//     //     phase_cram = ch_phase_cram.join(ch_phased_vcf).dump(tag: "joined")

        

//     //      WHATSHAP_PHASE (
//     //          phase_cram,
//     //          phase_vcf,
//     //          file(params.fasta),
//     //          file(params.fasta_index)
//     //      )
//     //     ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)

//     //     /*
//     //      * Sort phased variants with bcftools
//     //      */
//     //     PHASE_SORT_VCF( WHATSHAP_PHASE.out.phased_vcf )
//     //     ch_sv_phase_vcf = PHASE_SORT_VCF.out.vcf
//     //     ch_versions = ch_versions.mix(PHASE_SORT_VCF.out.versions)

//     //     /*
//     //      * Index phased vcf.gz
//     //      */
//     //     PHASE_TABIX_VCF( ch_sv_phase_vcf )
//     //     ch_sv_calls_tbi  = PHASE_TABIX_VCF.out.tbi
//     //     ch_versions = ch_versions.mix( PHASE_TABIX_VCF.out.versions)

//     //     //
//     //     // MODULE: whatshap for haplotag
//     //     //
//     //     ch_haplotag_cram = SAMTOOLS_SORT.out.cram
//     //         .mix(SAMTOOLS_SORT.out.crai)
//     //         .groupTuple(size:2)
//     //         .map{ meta, files -> [ meta, files.flatten() ]}

//     //     haplotag_cram = ch_haplotag_cram.join(ch_phased_vcf).dump(tag: "joined")


//     //     add_meta = INPUT_CHECK.out.reads
//     //     .map{ meta, files -> [[sample: meta.sample],meta.vcf_tbi] }
//     //     .dump(tag: "add_meta")


//     //     ch_haplotag_vcf = PHASE_SORT_VCF.out.vcf
//     //         .mix(PHASE_TABIX_VCF.out.tbi)
//     //         .groupTuple(size:2)
//     //         .map{ meta, files -> [ meta, files.flatten() ]}
            
//     //     haplotag_vcf = ch_haplotag_vcf.join(add_meta).dump(tag: "joined")
     
//     //      WHATSHAP_HAPLOTAG (
//     //          haplotag_cram,
//     //          haplotag_vcf,
//     //          file(params.fasta),
//     //          file(params.fasta_index)
//     //      )
         
//         //ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: whatshap depth calculation
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

//         //
//         // MODULE: MOSDEPTH for depth calculation
//         //
//     ch_mosdepth_input = WHATSHAP_HAPLOTAG.out.cram.mix(WHATSHAP_HAPLOTAG.out.crai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
//         MOSDEPTH (
//             ch_mosdepth_input
//         )
//         ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
//     // ...existing code...
    
    
// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: currently remove methylation
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// /*
//         //
//         // MODULE: MODKIT to extract methylation data
//         //
        
//         if (params.extract_methylation) {
//             ch_modkit_input = WHATSHAP.out.cram.mix(WHATSHAP.out.crai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
//             MODKIT (
//                 ch_modkit_input,
//                 file(params.fasta)
//             )
//             ch_versions = ch_versions.mix(MODKIT.out.versions)

//             ch_modkit_to_bw_input = MODKIT.out.hap1_bed.join(MODKIT.out.hap2_bed).join(MODKIT.out.combined_bed)
//             MODKIT_TO_BW (
//                 ch_modkit_to_bw_input,
//                 file(params.fasta_index)
//             )
//             ch_versions = ch_versions.mix(MODKIT_TO_BW.out.versions)

//             SAMTOOLS_STATS (
//                 WHATSHAP.out.cram
//             )
//             ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
//         }
// */

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     FANIVA: CUSTOM_DUMPSOFTWAREVERSIONS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

//     CUSTOM_DUMPSOFTWAREVERSIONS (
//         ch_versions.unique().collectFile(name: 'collated_versions.yml')
//     )

//     //
//     // MODULE: MultiQC
//     //
//     workflow_summary    = WorkflowFaniva.paramsSummaryMultiqc(workflow, summary_params)
//     ch_workflow_summary = Channel.value(workflow_summary)

//     methods_description    = WorkflowFaniva.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
//     ch_methods_description = Channel.value(methods_description)

//     ch_multiqc_files = Channel.empty()
//     ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//     ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
//     ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

//     if (params.reads_format == 'fast5' || params.reads_format == 'pod5') {
//         ch_multiqc_files = ch_multiqc_files.mix(PYCOQC.out.json.collect{it[1]}.ifEmpty([]))
//     }

//     if (params.run_whatshap) {
//         ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]}.ifEmpty([]))
//         ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]}.ifEmpty([]))
//         ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.collect{it[1]}.ifEmpty([]))
//         ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_bed.collect{it[1]}.ifEmpty([]))
//         ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_csi.collect{it[1]}.ifEmpty([]))
//         ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_bed.collect{it[1]}.ifEmpty([]))
//         ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_csi.collect{it[1]}.ifEmpty([]))
//     }

//     MULTIQC (
//         ch_multiqc_files.collect(),
//         ch_multiqc_config.toList(),
//         ch_multiqc_custom_config.toList(),
//         ch_multiqc_logo.toList()
//     )
//     multiqc_report = MULTIQC.out.report.toList()
// }

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     COMPLETION EMAIL AND SUMMARY
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */



