/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMutectplatypus.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta,
    params.fasta_fai,
    params.dict,
    params.germline_resource,
    params.germline_resource_idx,
    params.haplotype_map,
    params.drivers,
    params.ngscheckmate_bed,
    params.mappability,
    params.mappability_tbi
    ]

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope

fasta                 = params.fasta                 ? Channel.fromPath(params.fasta).collect()                 : Channel.empty()
fasta_fai             = params.fasta_fai             ? Channel.fromPath(params.fasta_fai).collect()             : Channel.empty()
dict                  = params.dict                  ? Channel.fromPath(params.dict).collect()                  : Channel.empty()
germline_resource     = params.germline_resource     ? Channel.fromPath(params.germline_resource).collect()     : Channel.empty()
germline_resource_idx = params.germline_resource_idx ? Channel.fromPath(params.germline_resource_idx).collect() : Channel.empty()
haplotype_map         = params.haplotype_map         ? Channel.fromPath(params.haplotype_map).collect()         : Channel.empty()
drivers               = params.drivers               ? Channel.fromPath(params.drivers).collect()               : Channel.empty()
//mappability_bw        = params.mappability_bw        ? Channel.fromPath(params.mappability_bw).collect()        : Channel.empty()
mappability           = params.mappability           ? Channel.fromPath(params.mappability).collect()           : Channel.empty()
mappability_tbi       = params.mappability_tbi       ? Channel.fromPath(params.mappability_tbi).collect()       : Channel.empty()
ngscheckmate_bed      = params.ngscheckmate_bed      ? Channel.fromPath(params.ngscheckmate_bed).collect()                  : Channel.empty()
intervals_ch          = params.intervals             ? Channel.fromPath(params.intervals).collect()             : []

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
vep_cache_version    = params.vep_cache_version     ?: Channel.empty()
vep_genome           = params.vep_genome            ?: Channel.empty()
vep_species          = params.vep_species           ?: Channel.empty()
vep_cache            = params.vep_cache             ? Channel.fromPath(params.vep_cache).collect()               : []

// Initialise file channels based on params not defined in the params.genomes[params.genome]
pon                   = params.pon                   ? Channel.fromPath(params.pon).collect()                    : []
pon_idx               = params.pon_idx               ? Channel.fromPath(params.pon_idx).collect()                : []

// Initialise value channels not defined in params.genomes[params.genome]

seqz_het             = params.seqz_het ? Channel.value(params.seqz_het) : Channel.empty()
bin                  = params.bin                    ?: Channel.empty()
gender               = params.gender                 ?: Channel.empty()
ploidy               = params.ploidy                 ?: Channel.empty()
pan                  = params.pan                    ?: Channel.empty()
tef                  = params.tef                    ?: Channel.empty()

// initialise sequenza gamma values

switch (params.sequenza_gamma) {
                                case 5 : seq_gam = 280
                                break;
                                case 4 : seq_gam = 200
                                break;
                                case 3 : seq_gam = 150
                                break;
                                case 2 : seq_gam = 100
                                break;
                                case 1 : seq_gam = 50
                                break;
    }
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

// Sample QC on BAM files
include { BAM_SAMPLEQC } from '../subworkflows/local/bam_sampleqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Installed from modules/local
//

include { CREATE_INTERVALS_BED            } from '../modules/local/create_intervals_bed/main'
include { BUILD_INTERVALS                 } from '../modules/local/build_intervals/main'
//include { PLATYPUS_CALLVARIANTS           } from '../modules/local/platypus/main'
//include { PLATYPUS_FILTER                 } from '../modules/local/filter_platypus/main'
include { VCF_SPLIT                       } from '../subworkflows/local/vcf_split.nf'
include { ZIP_VCF as ZIP_MUTECT_ANN_VCF   } from '../modules/local/zip_vcf/main'
//include { ZIP_VCF as ZIP_PLATYPUS_ANN_VCF } from '../modules/local/zip_vcf/main'
include { ZIP_VCF as ZIP_MUTECT_MONO_VCF  } from '../modules/local/zip_vcf/main'
include { VCF2MAF                         } from '../modules/local/vcf2maf/main'
include { CONCAT_VCF as CONCAT_MUTECT     } from '../modules/local/concat_vcf/main'
//include { CONCAT_VCF as CONCAT_PLATYPUS   } from '../modules/local/concat_vcf/main'
include { SEQUENZAUTILS_MERGESEQZ         } from '../modules/local/sequenzautils/mergeseqz/main'
include { SEQUENZAUTILS_BINNING           } from '../modules/local/sequenzautils/seqzbin/main'
include { SEQUENZAUTILS_RSEQZ             } from '../modules/local/sequenzautils/seqz_R/main.nf'
//include { MAPPABILITY                     } from '../modules/local/mappability/main'
include { BCFTOOLS_MAPPABILITY            } from '../modules/nf-core/bcftools/annotate/main'
include { EVOVERSE_CNAQC                  } from '../modules/local/evoverse/main'

//
// MODULE: Installed directly from nf-core/modules
//

include { PICARD_CROSSCHECKFINGERPRINTS   } from '../modules/nf-core/picard/crosscheckfingerprints/main'
include { GATK4_MUTECT2                   } from '../modules/nf-core/gatk4/mutect2/main'
include { GATK4_LEARNREADORIENTATIONMODEL } from '../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS          } from '../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_GETPILEUPSUMMARIES as GET_PS_TUMOUR } from '../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GET_PS_NORM } from '../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATHER_PS_TUMOUR } from '../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATHER_PS_NORM } from '../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    } from '../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         } from '../modules/nf-core/gatk4/filtermutectcalls/main'
//include { ENSEMBLVEP as PLAT_VEP          } from '../modules/nf-core/ensemblvep/main'
include { ENSEMBLVEP                      } from '../modules/nf-core/ensemblvep/main'
include { SEQUENZAUTILS_GCWIGGLE          } from '../modules/nf-core/sequenzautils/gcwiggle/main'
include { SEQUENZAUTILS_BAM2SEQZ          } from '../modules/nf-core/sequenzautils/bam2seqz/main'
//include { SEQUENZAUTILS_HETSNPS           } from '../modules/nf-core/modules/sequenzautils/hetsnps/main'
//include { ADD_MAPPABILITY                 } from '../modules/local/add_mappability/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Function to split input channel into tumour and normal for mutect and collect samples from a single patient

def make_mutect_input(input) {
    return input
        .map { meta, files -> [ meta.patient, meta.id, meta.status, files[0],files[1]] }
        .groupTuple()
        .map { patient, id, status, bam, bai -> [ patient, id[status.findIndexValues { it ==~ /tumour/ }], id[status.findIndexValues { it ==~ /normal/ }], bam, bai ]}
}

// Info required for completion email and summary
def multiqc_report = []

workflow MUTECT_PLATYPUS {

    ch_versions = Channel.empty()
    result_intervals = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    //INPUT_CHECK.out.bams.view{"input check out: $it"}
    INPUT_CHECK.out.bams
        .map { meta, files -> [ meta.patient, meta.id, meta.status, files[0], files[1] ] }
        .groupTuple()
        .map { patient, id, status, bams, bais -> [ patient, id[status.findIndexValues { it ==~ /tumour/ }], id[status.findIndexValues { it ==~ /normal/ }], bams, bais ]}
        .set{ mutect_input }

//    mutect_input = make_mutect_input(INPUT_CHECK.out.bams)
//    mutect_input.view{"mutect input: $it"}
    
    //mutect_input.view{"mutect input: $it"}
    //crosscheck_input = make_crosscheck_imput(INPUT_CHECK.out.bams)
    //pair_input = crosscheck_input.control.combine(crosscheck_input.tumour, by:0)


    //
    // create intervals to split jobs.
    //

    if (!params.intervals) {
        result_intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS(fasta_fai))
        } else {
        result_intervals = CREATE_INTERVALS_BED(intervals_ch)
        }

    result_intervals = result_intervals.flatten()
        .map{ intervalFile ->
            def duration = 0.0
            for (line in intervalFile.readLines()) {
                final fields = line.split('\t')
                if (fields.size() >= 5) duration += fields[4].toFloat()
                else {
                    start = fields[1].toInteger()
                    end = fields[2].toInteger()
                    duration += (end - start) / params.nucleotides_per_second
                }
            }
            [duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
        .flatten().collate(2)
        .map{duration, intervalFile -> intervalFile}
    

    mutect_input.combine(result_intervals).map{ patient, which_tumour, which_norm, bams, bais, intervals ->
        [patient, intervals.baseName + "_" + patient, which_tumour, which_norm, bams, bais, intervals]
        }
        .set{bam_intervals}


    //PICARD_CROSSCHECKFINGERPRINTS ( mutect_input, haplotype_map )
    INPUT_CHECK.out.bams
                    .map{ meta, files -> [ meta, files[0], files[1] ] }
                    .set { ngscheckmate_input }

    BAM_SAMPLEQC(ngscheckmate_input, ngscheckmate_bed, fasta)

    GATK4_MUTECT2(
        bam_intervals,
        fasta,
        fasta_fai,
        dict,
        pon,
        pon_idx,
        germline_resource,
        germline_resource_idx,
    )

    orientation_input = GATK4_MUTECT2.out.f1r2.groupTuple()

    GATK4_LEARNREADORIENTATIONMODEL ( orientation_input )

    concat_input = GATK4_MUTECT2.out.vcf.groupTuple()

    if (!params.intervals) {
        CONCAT_MUTECT ( concat_input, fasta_fai, [] )
        } else {
        CONCAT_MUTECT ( concat_input, fasta_fai, intervals_ch )
        }

    merge_stats_input = GATK4_MUTECT2.out.stats
                                        .groupTuple()

    GATK4_MERGEMUTECTSTATS ( merge_stats_input )

    // split input to create tables for each tumour
    INPUT_CHECK.out.bams
        .branch {   tumour: it[0]["status"] == "tumour"
                    normal: it[0]["status"] == "normal"
            }
        .set{pileup}

    pileup.tumour.combine(result_intervals)
                .map{ meta, files, intervals ->
                [ meta, meta.id + "_" + intervals.baseName, files, intervals] }
                .set{ pileup_tumour_intervals }

    GET_PS_TUMOUR (
        pileup_tumour_intervals,
        fasta,
        fasta_fai,
        dict,
        germline_resource,
        germline_resource_idx
    )
    gather_pileup_tumour_input = GET_PS_TUMOUR.out.table
                                    .groupTuple()
                                    .map{ meta, interval_ids, table -> [meta,table]}

    GATHER_PS_TUMOUR ( gather_pileup_tumour_input, dict )

    pileup.normal.combine(result_intervals)
                .map{ meta, files, intervals ->
                [ meta, meta.id + "_" + intervals.baseName, files, intervals] }
                .set{ pileup_normal_intervals }

    GET_PS_NORM (
        pileup_normal_intervals,
        fasta,
        fasta_fai,
        dict,
        germline_resource,
        germline_resource_idx
    )

    // remove intervals for join with GATHER_PS_NORM    
    GET_PS_NORM.out.table
                    .groupTuple()
                    .map{ meta, interval_ids, table -> [meta,table]}
                    .set{ gather_pileup_normal_input }

    GATHER_PS_NORM ( gather_pileup_normal_input, dict )

    GATHER_PS_TUMOUR.out.table
                        .map{ meta, table ->
                        [ meta.patient, [ sample:meta.sample, status:meta.status, id:meta.id], table ] }
                        .set{ tumour_gatherpileup_input }

    GATHER_PS_NORM.out.table
                        .map{ meta, table ->
                        [ meta.patient, [ sample:meta.sample, status:meta.status, id:meta.id ], table ] }
                        .set{ normal_gatherpileup_input  }
    
    tumour_gatherpileup_input.combine(normal_gatherpileup_input, by:0)
                        .map{ patient, meta1, table, meta2, table2 ->
                                    [ meta1 + [patient:patient], table, table2] }
                        .set{ calculate_contamination_input }

    GATK4_CALCULATECONTAMINATION( calculate_contamination_input )

    GATK4_CALCULATECONTAMINATION.out.contamination
                        .map{ meta, contamination_table ->
                                    [ meta.patient , [sample:meta.sample,status:meta.status, id:meta.id],
                                    contamination_table ]}
                        .groupTuple()
                        //.view{"contamination input: $it"}
                        .set{ contamination_input }

    GATK4_CALCULATECONTAMINATION.out.segmentation
                        .map { meta, segmentation_table ->
                                    [meta.patient, [sample:meta.sample,status:meta.status, id:meta.id],
                                    segmentation_table]}
                        .groupTuple()
                        //.view{"segmentation input: $it"}
                        .set{ segmentation_input }

    contamination_input.join(segmentation_input)
                      //  .view{"joint $it"}
                        .set{ for_filter }


    CONCAT_MUTECT.out.vcf.join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior)         
                            .join(for_filter)
                            .set { filter_input }


    filter_input.join(GATK4_MERGEMUTECTSTATS.out.stats)
                .map{patient, vcf, tbi, prior, meta1, contamination_table, meta2, segementation_table, stats ->
                    [patient, vcf, tbi, prior, contamination_table, segementation_table, stats] }
                .set{ filter_mutect }

    GATK4_FILTERMUTECTCALLS ( filter_mutect, fasta, fasta_fai, dict )
    
    ENSEMBLVEP (
        GATK4_FILTERMUTECTCALLS.out.vcf,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        []
    )
    //ENSEMBLVEP.out.vcf.view{"ensembl out $it"}

    BCFTOOLS_MAPPABILITY(ENSEMBLVEP.out.vcf, mappability, mappability_tbi)

    BCFTOOLS_MAPPABILITY.out.vcf
              //  .view{ "bcftools mappability input: $it"}
                .map { patient, file -> [ patient , "spacer", "spacer2", file ] }
                .set { PATIENT_VCF }

    ZIP_MUTECT_ANN_VCF ( PATIENT_VCF )

    pileup.normal
                .map{ meta, files -> [meta.patient,meta, files]}
                .set{normal_for_split}

    pileup.tumour
                .map{ meta, files -> [meta.patient,meta, files]}
                .set{tumour_for_split}

    BCFTOOLS_MAPPABILITY.out.vcf
                .set{mappability_for_split}

    normal_for_split.combine(tumour_for_split, by:0)
                .map{ patient, meta_control, files_control, meta_tumour, files_tumour ->
                [ patient, meta_control, meta_tumour] }
                .combine(mappability_for_split, by:0)
                .set{samples_for_split}

    
    VCF_SPLIT(samples_for_split)

    VCF_SPLIT.out.vcf
            .map { meta_control, meta_tumour, vcf ->
                    [[patient:meta_control.patient, control:meta_control.id, tumour:meta_tumour.id],vcf]}
            .set{vcf2maf_input}

    VCF2MAF ( vcf2maf_input, fasta )
    
    vcf2maf_input
            .map { meta, vcf ->
                  [meta.patient, meta.control, meta.tumour, vcf ]} 
            .set { zip_mono_mutect_input }

    ZIP_MUTECT_MONO_VCF ( zip_mono_mutect_input )
    //gatk_filter_out = GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi)

    //platypus_input = gatk_filter_out.combine(bam_intervals, by:0)

    //PLATYPUS_CALLVARIANTS(  platypus_input,
    //                        fasta,
    //                        fasta_fai )

    //concat_platypus_input = PLATYPUS_CALLVARIANTS.out.vcf.groupTuple()

    /*
    if (!params.intervals) {
        CONCAT_PLATYPUS ( concat_platypus_input, fasta_fai, [] )
        } else {
        CONCAT_PLATYPUS ( concat_platypus_input, fasta_fai, intervals_ch )
        }

    filter_platypus_input = CONCAT_PLATYPUS.out.vcf.join(bam_intervals)
                                .map{patient,vcf,tbi,region,which_tumour,which_norm,bam,bai,bed ->
                                [patient, vcf, which_norm]}

    PLATYPUS_FILTER ( filter_platypus_input, tef )

    PLAT_VEP (
        PLATYPUS_FILTER.out.vcf,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        []
    )

    PLAT_VEP.out.vcf
                .map { patient, file -> [ patient , "spacer", "spacer2", file ] }
                .set { PLAT_PATIENT_VCF }

    ZIP_PLATYPUS_ANN_VCF ( PLAT_PATIENT_VCF )
    */
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    SEQUENZAUTILS_GCWIGGLE(fasta)

    // Gather built wig.gz file or get them from the params
    wiggle = file(params.wiggle).exists() ? Channel.fromPath(params.wiggle).collect() : SEQUENZAUTILS_GCWIGGLE.out.wiggle

// Function to split input channel into tumour and normal for sequenza
    INPUT_CHECK.out.bams
                    //.map{ meta, files ->
                    //[meta.patient,meta.sample,meta.status,meta.id, files] }
                    .branch { meta, files ->
                        tumour: meta.status == "tumour"
                        control: meta.status == "normal"
                    }
                    .set{seq_split}
    
    seq_split.tumour
                    .map { meta , files -> 
                         [meta.patient, [ id:meta.id , status:meta.status], files]}
                    .set{tumour_bam}
    seq_split.control
                    .map { meta , files -> 
                         [meta.patient, [ id:meta.id, status:meta.status], files]}
                    .set { control_bam }
    
    tumour_bam
            .combine(control_bam, by:0)
            .set{seq_input_pair}

    seqz_chr = fasta_fai
        .splitCsv(sep: "\t")
        .map{ chr -> chr[0][0] }
        .filter( ~/^chr\d+|^chr[X,Y]|^\d+|[X,Y]/ )
        .collect()


    seq_input_pair
            .map { patient, meta1, files1, meta2, files2 ->
 //                   tumour_id  meta1.id 
                    [meta1 + [patient:patient], files1, files2]
                    }
            .set{ seq_input_matched }

    SEQUENZAUTILS_BAM2SEQZ(seq_input_matched,
                            fasta,
                            fasta_fai,
                            seqz_het,
                            wiggle,
                            seqz_chr)

    SEQUENZAUTILS_BAM2SEQZ.out.seqz
                            .groupTuple()
                            .set{merge_seqz_input}

    SEQUENZAUTILS_MERGESEQZ (merge_seqz_input)

//    SEQUENZAUTILS_HETSNPS(SEQUENZAUTILS_MERGESEQZ.out.concat_seqz)
//    SEQUENZAUTILS_MERGESEQZ.out.concat_seqz.view()
    SEQUENZAUTILS_BINNING(SEQUENZAUTILS_MERGESEQZ.out.concat_seqz, bin)

    if ( params.sequenza_tissue_type == "PDO" ) {
        purity = Channel.of(["PDO_1",0.95],["PDO_2",0.98])
        rseqz_input = SEQUENZAUTILS_BINNING.out.seqz_bin.combine(purity) 
        SEQUENZAUTILS_RSEQZ(    rseqz_input,
                                gender,
                                ploidy,
                                seq_gam )
    } else if ( params.sequenza_tissue_type == "TISSUE" ) {
        purity = Channel.of(["TISSUE_1",0.1],["TISSUE_2",0.3],["TISSUE_3",0.5],["TISSUE_4", 0.7],["TISSUE_5",0.9])
        rseqz_input = SEQUENZAUTILS_BINNING.out.seqz_bin.combine(purity)
        SEQUENZAUTILS_RSEQZ(    rseqz_input,
                                gender,
                                ploidy,
                                seq_gam )
    }

    
   // SEQUENZAUTILS_RSEQZ.out.rseqz.view{"RSEQZ: $it"} 
    ZIP_MUTECT_MONO_VCF.out.vcf
                .map{ patient, control, tumour, vcf_gz, vcf_tbi ->
                      [patient, tumour, vcf_gz, vcf_tbi]}
                .set{ vcf_for_combine }

    SEQUENZAUTILS_RSEQZ.out.rseqz
        .map { meta, out_dir ->
                [ meta.patient, meta.id, out_dir ] }
        .set{ rseqz_for_combine }


        vcf_for_combine
            .combine(rseqz_for_combine,by:1)
            .map{ id, patient, vcf_gz, vcf_tbi, patient2, segment_dir ->
                    [ [patient:patient, id:id ], [vcf_gz,vcf_tbi], segment_dir]}
            .view{"evoInput $it"}
            .set{ evo_input }
    
    EVOVERSE_CNAQC(evo_input, ploidy, drivers )

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowMytest.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ENSEMBLVEP.out.report.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_FILTERMUTECTCALLS.out.stats.collect().ifEmpty([]))

    //multiqc_mutect_input = mutect_input.collect()
    multiqc_mutect_input = mutect_input
        .map{ patient, which_tumour, which_norm, bam, bai -> patient }

    MULTIQC (
        ch_multiqc_files.collect(), multiqc_mutect_input
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
