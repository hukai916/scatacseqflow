#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/scatacseqflow
========================================================================================
 nf-core/scatacseqflow Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/scatacseqflow
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

log.info Utils.logo(workflow, params.monochrome_logs)

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    // TODO nf-core: Update typical command used to run pipeline
    def command = "nextflow run main.nf --input samplesheet.csv --preprocess 10xgenomics --genome GRCh37 -profile docker"
    log.info NfcoreSchema.paramsHelp(workflow, params, json_schema, command)
    log.info Workflow.citation(workflow)
    log.info Utils.dashedLine(params.monochrome_logs)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params, json_schema)
log.info NfcoreSchema.paramsSummaryLog(workflow, params, json_schema)
log.info Workflow.citation(workflow)
log.info Utils.dashedLine(params.monochrome_logs)

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////

// Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////
include { PREPROCESS } from './workflows/pipeline' addParams( summary_params: summary_params )
include { DOWNSTREAM } from './workflows/pipeline' addParams( summary_params: summary_params )

workflow  SCATACSEQFLOW {
  if (params.preprocess) {
    log.info "Running preprocess ..."
    PREPROCESS ()

    if (params.preprocess == "default") {
      // DOWNSTREAM (PREPROCESS.out.bam, PREPROCESS.out.fragments)
    } else if (params.preprocess == "10xgenomics") {

    } else if (params.preprocess == "biorad") {

    }
    // DOWNSTREAM ()
  } else {
    DOWNSTREAM ()
  }
}

workflow {
  SCATACSEQFLOW ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
