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

// Print suppot genome:
// if (params.support_genome) {
//
//   include { DOWNLOAD_TEST } from './modules/local/download_test'
// }

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
  .set { ch_samplesheet }
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
  // exit 1, "Must specify eitehr --input_archr or --input_preprocess!"
}

// Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////
// include { PREPROCESS } from './workflows/pipeline' addParams( summary_params: summary_params )
include { PREPROCESS_DEFAULT } from './workflows/pipeline' addParams( summary_params: summary_params )
include { PREPROCESS_10XGENOMICS } from './workflows/pipeline' addParams( summary_params: summary_params )
include { DOWNSTREAM } from './workflows/pipeline' addParams( summary_params: summary_params )
include { SPLIT_BED  } from './modules/local/split_bed' addParams( options: modules['split_bed'] )
include { SPLIT_BAM  } from './modules/local/split_bam' addParams( options: modules['split_bam'] )
include { MULTIQC    } from './modules/local/multiqc' addParams( options: modules['multiqc'] )
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)

// For test, to be deleted:
include { DOWNLOAD_FROM_UCSC } from './modules/local/download_from_ucsc'
include { DOWNLOAD_FROM_UCSC_GTF } from './modules/local/download_from_ucsc_gtf'
include { DOWNLOAD_FROM_ENSEMBL } from './modules/local/download_from_ensembl'
include { DOWNLOAD_FROM_ENSEMBL_GTF } from './modules/local/download_from_ensembl_gtf'
include { BUILD_BSGENOME } from './modules/local/build_bsgenome'
include { BUILD_TXDB } from './modules/local/build_txdb'
// include { GENEID_TO_SYMBOL } from './modules/local/geneid_to_symbol'
include { PREP_GENOME } from './modules/local/prep_genome'
include { PREP_GTF } from './modules/local/prep_gtf'

include { BUILD_GENE_ANNOTATION } from './modules/local/build_gene_annotation' addParams( options: modules['build_gene_annotation'] )
include { BUILD_GENOME_ANNOTATION } from './modules/local/build_genome_annotation' addParams( options: modules['build_genome_annotation'] )
//
// if (params.support_genome) {
//   log.info "Genomes that can be downloaded:\n"
//
//   DOWNLOAD_TEST()
//   exit 1
// }

workflow  SCATACSEQFLOW {
  if (params.preprocess) {
    log.info "Running preprocess ..."
    // PREPROCESS (ch_samplesheet)

    if (params.preprocess == "default") {
      // if PREPROCESS emits multiple output, must use .out[index].view()
      // if PREPROCESS emits only one output, use .out.view() is fine.
      PREPROCESS_DEFAULT (ch_samplesheet)
      DOWNSTREAM (PREPROCESS_DEFAULT.out[2], "preprocess_default")
      SPLIT_BED (DOWNSTREAM.out[1]) // take a tuple (sample_name, fragment_path, tsv_path) as input
      SPLIT_BAM (PREPROCESS_DEFAULT.out[3], DOWNSTREAM.out[2].collect(), PREPROCESS_DEFAULT.out[4].collect(), "[^:]*") // input: sample_name, all_bams, all_fragments, barcode_regex
      log.info "HERE: downstream_res: " + DOWNSTREAM.out[0].view()
      // Add MultiQC module here:
      MULTIQC(PREPROCESS_DEFAULT.out[0].mix(DOWNSTREAM.out[0].ifEmpty([])).mix(Channel.from(ch_multiqc_config)).collect())
    } else if (params.preprocess == "10xgenomics") {
      PREPROCESS_10XGENOMICS (ch_samplesheet)
      DOWNSTREAM (PREPROCESS_10XGENOMICS.out[2], "preprocess_10xgenomics")
      SPLIT_BED (DOWNSTREAM.out[1])
      SPLIT_BAM (PREPROCESS_10XGENOMICS.out[3], DOWNSTREAM.out[2].collect(), PREPROCESS_10XGENOMICS.out[4].collect(), "NA")
      MULTIQC (DOWNSTREAM.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
    } else if (params.preprocess == "biorad") {
      exit 1, "biorad to be added"
    } else {
      exit 1, "must supply valid preproess option"
    }
  } else {
    DOWNSTREAM (ch_samplesheet_archr, "preprocess_null")
    SPLIT_BED (DOWNSTREAM.out[1])
    MULTIQC (DOWNSTREAM.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
    // ch_test = Channel.fromPath( '/Users/kaihu/Projects/workflow/test_data/10x_genomics_5k/remove_duplicate/*.bam' )
  }
}

workflow {
  if (0) { // for quick testing
    log.info "eeee"

    // exit 1, "just here"

    // DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
    // DOWNLOAD_FROM_UCSC_GTF (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))

    // DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
    // DOWNLOAD_FROM_ENSEMBL_GTF (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
    // ADD_CHR_GENOME()
    // ADD_CHR_GTF()
    // PPRE_GENOME_ARCHR(): add "chr", retrive "primary"
    //
    PREP_GENOME(Channel.fromPath(params.test_fasta), "custom_genome")
    PREP_GTF(PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, Channel.fromPath(params.test_gtf))

    BUILD_BSGENOME(PREP_GENOME.out.genome_fasta)
    // // ARCHR_CREATE_GENOME_ANNOTATION()
    BUILD_TXDB (BUILD_BSGENOME.out.bsgenome, PREP_GTF.out.gtf)
    BUILD_GENE_ANNOTATION(BUILD_TXDB.out.txdb, PREP_GTF.out.gtf, params.species_latin_name)

    if (params.archr_blacklist_bed) {
      BUILD_GENOME_ANNOTATION(BUILD_BSGENOME.out.bsgenome, BUILD_GENE_ANNOTATION.out.gene_annotation, params.archr_blacklist_bed)
    } else {
      BUILD_GENOME_ANNOTATION(BUILD_BSGENOME.out.bsgenome, BUILD_GENE_ANNOTATION.out.gene_annotation, "$projectDir/assets/file_token.txt")
    }

    // if (!params.archr_blacklist_bed || (params.archr_blacklist_bed == '')) {
    //   BUILD_GENOME_ANNOTATION(BUILD_GENOME_ANNOTATION.out.bsgenome, "$projectDir/assets/file_token.txt", params.species_latin_name)
    // } else {
    //   BUILD_GENOME_ANNOTATION(BUILD_BSGENOME.out.bsgenome, params.archr_blacklist_bed, params.species_latin_name)
    // }

    // if (params.species_latin_name) {
    //   GENEID_TO_SYMBOL (Channel.fromPath(params.test_gtf), params.species_latin_name)
    // } else {
    //   log.info "Param: Pls also supply --species_latin_name."
    //   // exit 1, "exit code here"// exit info wont be displayed anyway
    // }

  } else {
    SCATACSEQFLOW ()
  }
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
