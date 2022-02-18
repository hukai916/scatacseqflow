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

// Validate supplied genome:
include { get_genome_ucsc } from './modules/local/genome_ucsc'
include { get_genome_ensembl } from './modules/local/genome_ensembl'
if (params.ref_fasta_ucsc) {
  genome_ucsc_list = get_genome_ucsc()

  if (!(genome_ucsc_list.contains(params.ref_fasta_ucsc))) {
    exit 1, "Supplied UCSC genome not supported: " + params.ref_fasta_ucsc + " !"
  }
}
if (params.ref_fasta_ensembl) {
  genome_ensembl_list = get_genome_ensembl()

  if (!(genome_ensembl_list.contains(params.ref_fasta_ensembl))) {
    exit 1, "Supplied Ensembl genome not supported: " + params.ref_fasta_ensembl + " !"
  }
}
// List supported genomes if --support_genome flag:
if (params.support_genome) {
  genome_ucsc_list = get_genome_ucsc()
  genome_ensembl_list = get_genome_ensembl()

  log.info "Supported UCSC genomes: \n" + genome_ucsc_list.join(", ")
  log.info "\nSupported Ensembl genomes: \n" + genome_ensembl_list.join(", ")
  exit 1
}

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
  .set { ch_samplesheet_preprocess }
} else if (params.input_archr) { // Parse ArchR samplesheet:
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
} else {
  exit 1, "Must specify eitehr --input_archr or --input_preprocess!"
}

// Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////
include { PREPROCESS_DEFAULT } from './workflows/preprocess_default' addParams( summary_params: summary_params )
include { PREPROCESS_10XGENOMICS } from './workflows/preprocess_10xgenomics' addParams( summary_params: summary_params )
include { DOWNSTREAM_ARCHR } from './workflows/downstream_archr' addParams( summary_params: summary_params )
include { SPLIT_BED  } from './modules/local/split_bed' addParams( options: modules['split_bed'] )
include { SPLIT_BAM  } from './modules/local/split_bam' addParams( options: modules['split_bam'] )
include { MULTIQC    } from './modules/local/multiqc' addParams( options: modules['multiqc'] )
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)

workflow SCATACPIPE {
  take:
    input_archr
    input_preprocess
    ch_samplesheet

  main:
    if (input_archr) {
      log.info "Running DOWNSTREAM_ARCHR ..."

      DOWNSTREAM_ARCHR (ch_samplesheet, "preprocess_null")
      SPLIT_BED (DOWNSTREAM_ARCHR.out[1])
      MULTIQC (DOWNSTREAM_ARCHR.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
    } else if (input_preprocess) {
      log.info "Running PREPROCESS ..."

      if (input_preprocess == "default") {
        PREPROCESS_DEFAULT (ch_samplesheet)
        DOWNSTREAM_ARCHR (PREPROCESS_DEFAULT.out[2], "preprocess_default")
        SPLIT_BED (DOWNSTREAM_ARCHR.out[1]) // take a tuple (sample_name, fragment_path, tsv_path) as input
        SPLIT_BAM (PREPROCESS_DEFAULT.out[3], DOWNSTREAM_ARCHR.out[2].collect(), PREPROCESS_DEFAULT.out[4].collect(), "[^:]*") // input: sample_name, all_bams, all_fragments, barcode_regex
        MULTIQC(PREPROCESS_DEFAULT.out[0].mix(DOWNSTREAM_ARCHR.out[0].ifEmpty([])).mix(Channel.from(ch_multiqc_config)).collect())
      } else if (input_preprocess == "10xgenomics") {
        PREPROCESS_10XGENOMICS (ch_samplesheet)
        DOWNSTREAM_ARCHR (PREPROCESS_10XGENOMICS.out[2], "preprocess_10xgenomics")
        SPLIT_BED (DOWNSTREAM_ARCHR.out[1])
        SPLIT_BAM (PREPROCESS_10XGENOMICS.out[3], DOWNSTREAM_ARCHR.out[2].collect(), PREPROCESS_10XGENOMICS.out[4].collect(), "NA")
        MULTIQC (DOWNSTREAM_ARCHR.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
      } else if (input_preprocess == "biorad") {
        exit 1, "biorad to be added"
      } else {
        exit 1, "must supply valid preprocess option"
      }
    } else {
      exit 1, "Pls supply either --input_archr or --input_preprocess"
    }
}

workflow {
  if (params.input_archr) {
    SCATACPIPE (params.input_archr, params.input_preprocess, ch_samplesheet_archr)
  } else if (params.input_preprocess) {
    SCATACPIPE (params.input_archr, params.input_preprocess, ch_samplesheet_preprocess)
  }
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

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
