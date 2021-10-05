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

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

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
