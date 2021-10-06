#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/mutectplatypus
========================================================================================
    Github : https://github.com/nf-core/mutectplatypus
    Website: https://nf-co.re/mutectplatypus
    Slack  : https://nfcore.slack.com/channels/mutectplatypus
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.dbsnp                 = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi             = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.dict                  = WorkflowMain.getGenomeAttribute(params, 'dict')
params.fasta                 = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai             = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.germline_resource     = WorkflowMain.getGenomeAttribute(params, 'germline_resource')
params.germline_resource_idx = WorkflowMain.getGenomeAttribute(params, 'germline_resource_idx')
params.intervals             = WorkflowMain.getGenomeAttribute(params, 'intervals')
params.known_indels          = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi      = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { MUTECTPLATYPUS } from './workflows/mutectplatypus'

//
// WORKFLOW: Run main nf-core/mutectplatypus analysis pipeline
//
workflow NFCORE_MUTECTPLATYPUS {
    MUTECTPLATYPUS ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_MUTECTPLATYPUS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
