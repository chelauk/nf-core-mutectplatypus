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
    params.drivers
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
drivers               = params.drivers               ? Channel.fromPath(params.drivers).collect()               : Channel.empty()
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

seqz_het             = params.seqz_het               ?: Channel.empty()
bin                  = params.bin                    ?: Channel.empty()
gender               = params.gender                 ?: Channel.empty()
ploidy               = params.ploidy                 ?: Channel.empty()

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
include { PLATYPUS_CALLVARIANTS           } from '../modules/local/platypus/main'
include { PLATYPUS_FILTER                 } from '../modules/local/filter_platypus/main'
include { ZIP_VCF as ZIP_MUTECT_ANN_VCF   } from '../modules/local/zip_vcf/main'
include { ZIP_VCF as ZIP_PLATYPUS_ANN_VCF } from '../modules/local/zip_vcf/main'

//
// MODULE: Installed directly from nf-core/modules
//

include { GATK4_MUTECT2                   } from '../modules/nf-core/modules/gatk4/mutect2/main'
include { GATK4_LEARNREADORIENTATIONMODEL } from '../modules/nf-core/modules/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS          } from '../modules/nf-core/modules/gatk4/mergemutectstats/main'
include { CONCAT_VCF                      } from '../modules/nf-core/modules/concat_vcf/main'
include { CONCAT_VCF as CONCAT_PLATYPUS   } from '../modules/nf-core/modules/concat_vcf/main'
include { GATK4_GETPILEUPSUMMARIES        } from '../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GET_PS_NORM } from '../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     } from '../modules/nf-core/modules/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATHER_PS_NORM } from '../modules/nf-core/modules/gatk4/gatherpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    } from '../modules/nf-core/modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         } from '../modules/nf-core/modules/gatk4/filtermutectcalls/main'
include { ENSEMBLVEP as PLAT_VEP          } from '../modules/nf-core/modules/ensemblvep/main'
include { ENSEMBLVEP                      } from '../modules/nf-core/modules/ensemblvep/main'
include { SEQUENZAUTILS_GCWIGGLE          } from '../modules/nf-core/modules/sequenzautils/gcwiggle/main'
include { SEQUENZAUTILS_BAM2SEQZ          } from '../modules/nf-core/modules/sequenzautils/bam2seqz/main'
include { SEQUENZAUTILS_MERGESEQZ         } from '../modules/nf-core/modules/sequenzautils/mergeseqz/main'
//include { SEQUENZAUTILS_HETSNPS           } from '../modules/nf-core/modules/sequenzautils/hetsnps/main'
include { SEQUENZAUTILS_BINNING           } from '../modules/nf-core/modules/sequenzautils/seqzbin/main'
include { SEQUENZAUTILS_RSEQZ             } from '../modules/nf-core/modules/sequenzautils/seqz_R/main.nf'
include { EVOVERSE_CNAQC                  } from '../modules/local/evoverse/main'
include { MULTIQC                         } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

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
    mutect_input = make_mutect_input(INPUT_CHECK.out.bams)

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


    mutect_input.combine(result_intervals).map{ patient, which_tumour, which_norm, bam, bai, intervals ->
        [patient, intervals.baseName + "_" + patient, which_tumour, which_norm, bam, bai, intervals]
        }
        .set{bam_intervals}

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
        CONCAT_VCF ( concat_input, fasta_fai, [] )
        } else {
        CONCAT_VCF ( concat_input, fasta_fai, intervals_ch )
        }

    merge_stats_input = GATK4_MUTECT2.out.stats.groupTuple()

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

    GATK4_GETPILEUPSUMMARIES (
        pileup_tumour_intervals,
        fasta,
        fasta_fai,
        dict,
        germline_resource,
        germline_resource_idx
    )
    gather_pileup_tumour_input = GATK4_GETPILEUPSUMMARIES.out.table
                                    .groupTuple()
                                    .map{ meta, interval_ids, table -> [meta,table]}

    GATK4_GATHERPILEUPSUMMARIES ( gather_pileup_tumour_input, dict )

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
    gather_pileup_normal_input = GET_PS_NORM.out.table
                                    .groupTuple()
                                    .map{ meta, interval_ids, table -> [meta,table]}

    GATHER_PS_NORM ( gather_pileup_normal_input, dict )


    tumour_gatherpileup_input = GATK4_GATHERPILEUPSUMMARIES.out.table
                                    .map{ meta, table ->
                                    [ meta.patient, meta.sample, meta.status, meta.id, table ] }

    normal_gatherpileup_input = GATHER_PS_NORM.out.table
                                    .map{ meta, table ->
                                    [ meta.patient, meta.sample, meta.status, meta.id, table ] }

    contamination_input = tumour_gatherpileup_input.combine(normal_gatherpileup_input, by:0)
                                    .map{ patient, sample1, status1, id1, table, sample2, status2, id2, table2 ->
                                    [patient, sample1, table, table2] }

    GATK4_CALCULATECONTAMINATION( contamination_input )

    contamination_input = GATK4_CALCULATECONTAMINATION.out.contamination.groupTuple()

    segmentation_input = GATK4_CALCULATECONTAMINATION.out.segmentation.groupTuple()

    for_filter = contamination_input.join(segmentation_input)
                                    .map{ patient, samples1, contamination, samples2, segmentation ->
                                    [patient,contamination,segmentation] }

    filter_input = CONCAT_VCF.out.vcf.join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior)
    filter_input = filter_input.join(for_filter)
    filter_input = filter_input.join(GATK4_MERGEMUTECTSTATS.out.stats)

    GATK4_FILTERMUTECTCALLS ( filter_input, fasta, fasta_fai, dict )

    ENSEMBLVEP (
        GATK4_FILTERMUTECTCALLS.out.vcf,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        []
    )

    ZIP_MUTECT_ANN_VCF ( ENSEMBLVEP.out.vcf )

    gatk_filter_out = GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi)

    platypus_input = gatk_filter_out.combine(bam_intervals, by:0)

    PLATYPUS_CALLVARIANTS(  platypus_input,
                            fasta,
                            fasta_fai )

    concat_platypus_input = PLATYPUS_CALLVARIANTS.out.vcf.groupTuple()


    if (!params.intervals) {
        CONCAT_PLATYPUS ( concat_platypus_input, fasta_fai, [] )
        } else {
        CONCAT_PLATYPUS ( concat_platypus_input, fasta_fai, intervals_ch )
        }

    filter_platypus_input = CONCAT_PLATYPUS.out.vcf.join(bam_intervals)
                                .map{patient,vcf,tbi,region,which_tumour,which_norm,bam,bai,bed ->
                                [patient, vcf, which_norm]}

    PLATYPUS_FILTER ( filter_platypus_input )

    PLAT_VEP (
        PLATYPUS_FILTER.out.vcf,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        []
    )

    ZIP_PLATYPUS_ANN_VCF ( PLAT_VEP.out.vcf )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    SEQUENZAUTILS_GCWIGGLE(fasta)

    // Gather built wig.gz file or get them from the params
    wiggle = file(params.wiggle).exists() ? Channel.fromPath(params.wiggle).collect() : SEQUENZAUTILS_GCWIGGLE.out.wiggle

// Function to split input channel into tumour and normal for sequenza
    INPUT_CHECK.out.bams
                    .map{ meta, files ->
                    [meta.patient,meta.sample,meta.status,meta.id, files] }
                    .branch {
                        tumour: it[2] == "tumour"
                        normal: it[2] == "normal"
                    }
                    .set{seq_split}

    seq_split.tumour.combine(seq_split.normal, by:0)
                    .set{seq_input_pair}

    seqz_chr = fasta_fai
        .splitCsv(sep: "\t")
        .map{ chr -> chr[0][0] }
        .filter( ~/^chr\d+|^chr[X,Y]|^\d+|[X,Y]/ )

    seq_input_pair.combine(seqz_chr)
                .map{ patient, sample1, status1, id1, files1, sample2, status2, id2, files2, chr ->
                [ patient, id1, chr, files1, files2] }
                .set{ seq_input_chr }

    SEQUENZAUTILS_BAM2SEQZ(seq_input_chr,
                            fasta,
                            seqz_het,
                            wiggle)

    SEQUENZAUTILS_BAM2SEQZ.out.seqz
                            .groupTuple(by:[0,1])
                            .set{merge_seqz_input}

    SEQUENZAUTILS_MERGESEQZ (merge_seqz_input)

//    SEQUENZAUTILS_HETSNPS(SEQUENZAUTILS_MERGESEQZ.out.concat_seqz)

    SEQUENZAUTILS_BINNING(SEQUENZAUTILS_MERGESEQZ.out.concat_seqz, bin)

    SEQUENZAUTILS_RSEQZ(SEQUENZAUTILS_BINNING.out.seqz_bin, gender, ploidy)

    //ZIP_MUTECT_ANN_VCF.out.vcf.view()
    evo_input = SEQUENZAUTILS_RSEQZ.out.rseqz.combine(ZIP_MUTECT_ANN_VCF.out.vcf, by:0 )
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

    MULTIQC (
        ch_multiqc_files.collect()
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
