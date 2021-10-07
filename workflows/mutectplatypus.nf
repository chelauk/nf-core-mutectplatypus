/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMutectplatypus.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta,
    params.fasta_fai,
    params.dict,
    params.germline_resource,
    params.germline_resource_idx
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

//
// MODULE: platypus module
//
def platypusvariant_options   = modules['platypusvariant']
include { PLATYPUSVARIANT } from '../modules/local/platypusvariant' addParams( options: platypusvariant_options )

//
// MODULE: bcftools module
//

def bcftools_options          = modules['bcftools']
include { BCFTOOLS_CONCAT } from '../modules/nf-core/modules/bcftools/concat/main' addParams( options: bcftools_options)

//
// MODULE: filter platypus module
//

def filter_platypus_options    = modules['filter_platypus']
include { FILTERPLATYPUS } from '../modules/local/filterplatypus.nf' addParams( options: filter_platypus_options)

//
// MODULE: bgzip module
//

def bgzip_options              = modules['bgzip_vcfs']
include { BGZIP } from '../modules/local/bgzip' addParams( options: bgzip_options)

//
// MODULE: create_bed_intervals module
//
def bed_intervals_options      = modules['bed_intervals']

// Stage dummy file to be used as an optional input where required
include { CREATE_INTERVALS_BED }   from '../modules/local/create_intervals_bed/main' addParams(options: bed_intervals_options)
include { BUILD_INTERVALS }        from '../modules/local/build_intervals/main'

include { GATK4_MUTECT2 }          from '../modules/local/gatk/mutect2'

include { GATK4_GETPILEUPSUMMARIES } from '../modules/nf-core/modules/gatk4/getpileupsummaries/main'

include { GATK4_GATHERPILEUPSUMMARIES } from '../modules/nf-core/modules/gatk4/gatherpileupsummaries/main'

include { GATK4_CALCULATECONTAMINATION } from '../modules/nf-core/modules/gatk4/calculatecontamination/main'
ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope

fasta                 = params.fasta                 ? Channel.fromPath(params.fasta).collect()                 : ch_dummy_file
fasta_fai             = params.fasta_fai             ? Channel.fromPath(params.fasta_fai).collect()             : ch_dummy_file
dict                  = params.dict                  ? Channel.fromPath(params.dict).collect()                  : ch_dummy_file
germline_resource     = params.germline_resource     ? Channel.fromPath(params.germline_resource).collect()     : ch_dummy_file
germline_resource_idx = params.germline_resource_idx ? Channel.fromPath(params.germline_resource_idx).collect() : ch_dummy_file

// Initialise input sample
csv_file = file(params.input)
input_samples  = extract_csv(csv_file)

mutect_input = make_mutect_input(input_samples)

def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.patient && row.sample)) log.warn "Missing or unknown field in csv file header"
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing
        def meta = [:]

        //TODO since it is mandatory: error/warning if not present?
        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()
        meta.id                       = meta.patient + "_" +  meta.sample

        // If no gender specified, gender is not considered
        // gender is only mandatory for somatic CNV
        if (row.gender) meta.gender = row.gender.toString()
        else meta.gender = "NA"

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toString()
        else meta.status = "NA"

        if (row.bam && row.bai) {
            def bam        = file(row.bam, checkIfExists: true)
            def bai        = file(row.bai, checkIfExists: true)
            return [meta, [bam, bai]]
        // recalibration
        }
    }
}

def make_mutect_input(input) {
    return input
        .map { meta, files -> [ meta.patient, meta.id, meta.status, files[0],files[1]] }
		.groupTuple()
		.map { patient, id, status, bam, bai -> [ patient, id[status.findIndexValues { it ==~ /tumour/ }], id[status.findIndexValues { it ==~ /control/ }], bam, bai ]}
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MUTECTPLATYPUS {

    ch_software_versions = Channel.empty()

	result_intervals = Channel.empty()
    if (!params.intervals)
        result_intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS(fasta_fai))
    else
        result_intervals = CREATE_INTERVALS_BED(file(params.intervals))

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
        [patient, patient + "_" + intervals.baseName, which_tumour, which_norm, bam, bai, intervals]
    }.set{bam_intervals}

	GATK4_MUTECT2(
	    bam_intervals,
        fasta,
        fasta_fai,
        dict,
        germline_resource,
        germline_resource_idx,
    )


    // split input to create tables for each tumour
    input_samples.branch{ tumour: it[0]["status"] == "tumour"
	                      normal: it[0]["status"] == "control"
						 }
                  .set{pileup}

    pileup.tumour
        .combine(result_intervals)
        .map { meta, files, intervals -> [ meta.patient, meta.id, meta.id + "_" + intervals.baseName, meta.status, files[0],files[1], intervals] }
        .set{pileuptumour_intervals}

    GATK4_GETPILEUPSUMMARIES(
        pileuptumour_intervals,
        germline_resource,
        germline_resource_idx
        )

    gather_pileups = GATK4_GETPILEUPSUMMARIES.out.table
	  .groupTuple(by: 1)
	  .map { patient, id, id_intervals, status, pileup_tables -> [ patient.unique()[0], id, id_intervals, status.unique()[0], pileup_tables ] }


  GATK4_GATHERPILEUPSUMMARIES (
        gather_pileups,
 		dict
    )
  GATK4_CALCULATECONTAMINATION (
        GATK4_GATHERPILEUPSUMMARIES.out.gathered_table
  )
	//
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMutectplatypus.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
