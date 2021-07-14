////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Validate input parameters
// Workflow.validateWorkflowParams(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
// def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
//
// def multiqc_options   = modules['multiqc']
// multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''
//
// Modules: local
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['csv':'']] )
include { GET_10XGENOMICS_FASTQ } from '../modules/local/get_10xgenomics_fastq'   addParams( options: modules['get_10xgenomics_fastq'] )
include { GET_BIORAD_FASTQ      } from '../modules/local/get_biorad_fastq'        addParams( options: modules['get_biorad_fastq'] )

include { CELLRANGER_ATAC_COUNT } from '../modules/local/cellranger_atac_count'   addParams( options: modules['cellranger_atac_count'] )
include { CORRECT_BARCODE       } from '../modules/local/correct_barcode'         addParams( options: modules['correct_barcode'] )
include { MATCH_READS           } from '../modules/local/match_reads'             addParams( options: modules['match_reads'] )
include { FASTQC                } from '../modules/local/fastqc'                  addParams( options: modules['fastqc'] )

include { BIORAD_FASTQC         } from '../modules/local/biorad_fastqc'           addParams( options: modules['biorad_fastqc'] )
include { BIORAD_TRIM_READS     } from '../modules/local/biorad_atac_seq_trim_reads'       addParams( options: modules['biorad_atac_seq_trim_reads'] )

// // Modules: nf-core/modules
// include { FASTQC                } from '../modules/nf-core/software/fastqc/main'  addParams( options: modules['fastqc']            )
// include { MULTIQC               } from '../modules/nf-core/software/multiqc/main' addParams( options: multiqc_options              )
//
// // Subworkflows: local
// include { INPUT_CHECK           } from '../subworkflows/local/input_check'        addParams( options: [:]                          )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
// parse samplesheet.csv
if (params.input) {
  Channel
  .from(file(params.input, checkIfExists: true))
  .splitCsv(header: true, sep: ",", strip: true)
  .map {
    row ->
      [ row.sample_name, row.path_fastq_1, row.path_fastq_2, row.path_barcode ]
  }
  .unique()
  .set { ch_samplesheet }
}

// ch_samplesheet.view()
// Info required for completion email and summary
// def multiqc_report = []


workflow PREPROCESS {
  // Module: prepare 10xgenomics folder structure
  if (params.preprocess == "default") {
    log.info "INFO: --preprocess: default"
    GET_10XGENOMICS_FASTQ (ch_samplesheet)
    // module: fastQC
    FASTQC (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)

    // module: barcode correction: correct barcode fastq given whitelist and barcode fastq file
    if (!(params.barcode_whitelist)) {
      log.info "NOTICE: --barcode_whitelist: not supplied, skip barcode correction!"
    } else {
      CORRECT_BARCODE (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, params.barcode_whitelist)
      // module: match read1 and read2
      MATCH_READS (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.corrected_barcode, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
    }

    // module: debarcode: add barcode sequence to the beginning of the fastq sequence identifier with sinto

    // module: trimming

    // module: mapping with bwa or minimap2: mark duplicate

    // module: generate fragement file with sinto

  } else if (params.preprocess == "10xgenomics") {
    log.info "INFO: --preprocess: 10xgenomics"
    if (params.ref_cellranger == "") {
      exit 1, 'Parameter --ref_cellranger: pls supply full path to reference!'
    }

    // Module: prepare fastq folder
    GET_10XGENOMICS_FASTQ (ch_samplesheet)
    CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.fastq_folder, params.ref_cellranger)
    // sample_name = PARSEUMI.out.umi.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()
    // sample_fastq_folder = GET_10XGENOMICS_FASTQ.out.fastq.to

  } else if (params.preprocess == "biorad") {
    log.info "INFO: --preprocess: biorad"
    log.info "INFO: must use biorad compatible sequencing results!"
    GET_BIORAD_FASTQ (ch_samplesheet)
    BIORAD_FASTQC (GET_BIORAD_FASTQ.out.sample_name, GET_BIORAD_FASTQ.out.fastq_folder)
    BIORAD_TRIM_READS (BIORAD_FASTQC.out.sample_name, BIORAD_FASTQC.out.fastqc_results)

  } else {
    log.info "ERROR: for parameter --preprocess, choose from default, 10xgenomics, biorad."
  }
}

workflow SCATACSEQFLOW {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    // INPUT_CHECK (
    //     ch_input
    // )

    /*
     * MODULE: Run FastQC
     */
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))


    // /*
    //  * MODULE: Pipeline reporting
    //  */
    // // Get unique list of files containing version information
    // ch_software_versions
    //     .map { it -> if (it) [ it.baseName, it ] }
    //     .groupTuple()
    //     .map { it[1][0] }
    //     .flatten()
    //     .collect()
    //     .set { ch_software_versions }
    // GET_SOFTWARE_VERSIONS (
    //     ch_software_versions
    // )

    /*
     * MODULE: MultiQC
     */
    // workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)
    //
    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report       = MULTIQC.out.report.toList()
    // ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

// workflow.onComplete {
//     Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
//     Completion.summary(workflow, params, log)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
