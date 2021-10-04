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

def modules = params.modules.clone()

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////

// Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////
include { PREPROCESS } from './workflows/pipeline' addParams( summary_params: summary_params )
include { DOWNSTREAM } from './workflows/pipeline' addParams( summary_params: summary_params )
include { SPLIT_BED  } from './modules/local/split_bed' addParams( options: modules['split_bed'] )
include { SPLIT_BAM  } from './modules/local/split_bam' addParams( options: modules['split_bam'] )

// Parse samplesheet:
if (params.input_preprocess) {
  Channel
  .from(file(params.input_preprocess, checkIfExists: true))
  .splitCsv(header: true, sep: ",", strip: true)
  .map {
    row ->
      [ row.sample_name, row.path_fastq_1, row.path_fastq_2, row.path_barcode ]
  }
  .unique()
  .set { ch_samplesheet }
}
// Parse ArchR samplesheet:
if (params.input_archr) {
  Channel
  .from(file(params.input_archr, checkIfExists: true))
  .splitCsv(header: true, sep: ",", strip: true)
  .map {
    row ->
      [ row.sample_name, row.file_path ]
  }
  .unique()
  .set { ch_samplesheet_archr }

  // log.info "ch_samplesheet_archr: " + ch_samplesheet_archr.view()
}

workflow  SCATACSEQFLOW {
  if (params.preprocess) {
    log.info "Running preprocess ..."
    PREPROCESS (ch_samplesheet)
    // PREPROCESS.out.view() // if put view() here, error?

    log.info "Running downstream analysis with ArchR ..."
    if (params.preprocess == "default") {
      // if PREPROCESS emits multiple output, must use .out[index].view()
      // if PREPROCESS emits only one output, use .out.view() is fine.
      // PREPROCESS.out[0].view() result files for multiQC
      // PREPROCESS.out[1].view() remove_duplicated bam files for splitting
      // PREPROCESS.out[2].view() fragement file for ArchR
      log.info "TEST preprocess output"
      PREPROCESS.out[3].view()
      log.info "TEST HERE"

      DOWNSTREAM (PREPROCESS.out[3])
      SPLIT_BED(DOWNSTREAM.out[1])
      // SPLIT_BAM(ch_samplesheet_archr, DOWNSTREAM.out[2].collect(), PREPROCESS.out[1].collect(), "[^:]*")
    } else if (params.preprocess == "10xgenomics") {
      // DOWNSTREAM (PREPROCESS.out[1], PREPROCESS.out[2])
      // SPLIT_BED(DOWNSTREAM.out[1])
      // SPLIT_BAM(PREPROCESS.out[bam_filter].collect(), DOWNSTREAM.out[1].collect(), "NA")
    } else if (params.preprocess == "biorad") {

    }
  } else {
    DOWNSTREAM (ch_samplesheet_archr)
    // DOWNSTREAM.out[2].collect().view()
    // DOWNSTREAM.out[2].collect().flatten().filter( ~/^.*\.tsv$/ ).view()
    SPLIT_BED(DOWNSTREAM.out[1])
    ch_test = Channel.fromPath( '/Users/kaihu/Projects/workflow/test_data/10x_genomics_5k/remove_duplicate/*.bam' )
    if (!(params.barcode_regex)) {
      SPLIT_BAM(ch_samplesheet_archr, DOWNSTREAM.out[2].collect(), PREPROCESS.out[1].collect(), "NA")
    } else {
      SPLIT_BAM(ch_samplesheet_archr, DOWNSTREAM.out[2].collect(), PREPROCESS.out[1].collect(), params.barcode_regex, params.barcode_regex)
    }
  }
}

workflow {
  SCATACSEQFLOW ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
