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
if (params.preprocess) {
  if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
} else {
  // will perform downstream analysis
  // if (params.sample_name) {} else { exit 1, "Input sample_name must be specified!" }
  // if (params.fragment) { ch_fragment = file(params.fragment) } else { exit 1, "Input fragment must be provided!" }
  if (params.input_archr) {} else { exit 1, "--input_archr samplesheet must be specified!"}

}

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
include { BIORAD_ATAC_SEQ_DEBARCODE } from '../modules/local/biorad_atac_seq_debarcode'           addParams( options: modules['biorad_atac_seq_debarcode'] )
include { BIORAD_ATAC_SEQ_TRIM_READS } from '../modules/local/biorad_atac_seq_trim_reads'       addParams( options: modules['biorad_atac_seq_trim_reads'] )
include { BIORAD_ATAC_SEQ_BWA   } from '../modules/local/biorad_atac_seq_bwa'     addParams( options: modules['biorad_atac_seq_bwa'] )
include { BIORAD_ATAC_SEQ_ALIGNMENT_QC   } from '../modules/local/biorad_atac_seq_alignment_qc'     addParams( options: modules['biorad_atac_seq_alignment_qc'] )
include { BIORAD_ATAC_SEQ_FILTER_BEADS   } from '../modules/local/biorad_atac_seq_filter_beads'     addParams( options: modules['biorad_atac_seq_filter_beads'] )
include { ADD_BARCODE_TO_READS       } from '../modules/local/add_barcode_to_reads'    addParams( options: modules['add_barcode_to_reads'] )
include { CUTADAPT         } from '../modules/local/cutadapt'    addParams( options: modules['cutadapt'] )

include { DOWNLOAD_FROM_UCSC        } from '../modules/local/download_from_ucsc'    addParams( options: modules['download_from_ucsc'] )
include { DOWNLOAD_FROM_ENSEMBL     } from '../modules/local/download_from_ensembl'    addParams( options: modules['download_from_ensembl'] )
include { GET_PRIMARY_GENOME        } from '../modules/local/get_primary_genome'    addParams( options: modules['get_primary_genome'] )
include { BWA_INDEX        } from '../modules/local/bwa_index'    addParams( options: modules['bwa_index'] )
include { BWA_MAP          } from '../modules/local/bwa_map'    addParams( options: modules['bwa_map'] )

include { MINIMAP2_INDEX   } from '../modules/local/minimap2_index'    addParams( options: modules['minimap2_index'] )
include { MINIMAP2_MAP     } from '../modules/local/minimap2_map'    addParams( options: modules['minimap2_map'] )

include { QUALIMAP         } from '../modules/local/qualimap'    addParams( options: modules['qualimap'] )
include { GET_FRAGMENTS         } from '../modules/local/get_fragments'    addParams( options: modules['get_fragments'] )

include { DOWNLOAD_FROM_UCSC_GTF } from '../modules/local/download_from_ucsc_gtf'    addParams( options: modules['download_from_ucsc_gtf'] )
include { FIX_UCSC_GTF } from '../modules/local/fix_ucsc_gtf'    addParams( options: modules['fix_ucsc_gtf'] )
include { DOWNLOAD_FROM_ENSEMBL_GTF } from '../modules/local/download_from_ensembl_gtf'    addParams( options: modules['download_from_ensembl_gtf'] )
include { CELLRANGER_INDEX } from '../modules/local/cellranger_index'             addParams( options: modules['cellranger_index'] )

// For ArchR functions:
include { ARCHR_CREATE_ARROWFILES } from '../modules/local/archr_create_arrowfiles' addParams( options: modules['archr_create_arrowfiles'] )
include { ARCHR_ADD_DOUBLETSCORES } from '../modules/local/archr_add_doubletscores' addParams( options: modules['archr_add_doubletscores'] )
include { ARCHR_ARCHRPROJECT } from '../modules/local/archr_archrproject' addParams( options: modules['archr_archrproject'] )
include { ARCHR_ARCHRPROJECT_QC } from '../modules/local/archr_archrproject_qc' addParams( options: modules['archr_archrproject_qc'] )


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
  // Check to see if parameter modules are specified or not
  if (params.preprocess == "default") {
    log.info "INFO: --preprocess: default"
  } else if (params.preprocess == "10xgenomics") {
    log.info "INFO: --preprocess: 10xgenomics"
  } else if (params.preprocess == "biorad") {
    log.info "INFO: --preprocess: biorad"
    log.info "INFO: must use biorad compatible sequencing data!"
  } else {
    log.info "ERROR: for parameter --preprocess, choose from default, 10xgenomics, biorad."
  }

  if (params.preprocess == "default") {
    if (!(params.barcode_whitelist)) {
      log.info "NOTICE: --barcode_whitelist: not supplied, skip barcode correction!"
    }
  }

  // Module: prepare 10xgenomics folder structure
  if (params.preprocess == "default") {
    // log.info "INFO(2): --preprocess: default"
    GET_10XGENOMICS_FASTQ (ch_samplesheet)
    // module: fastQC
    FASTQC (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)

    // module: barcode correction (optional) and add barcode: correct barcode fastq given whitelist and barcode fastq file
    if (!(params.barcode_whitelist)) {
      // log.info "NOTICE(2): --barcode_whitelist: not supplied, skip barcode correction!"

      ADD_BARCODE_TO_READS (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
    } else {
      CORRECT_BARCODE (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, params.barcode_whitelist, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
      // module: match read1 and read2
      // MATCH_READS (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.corrected_barcode, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
      // Note that the above might be problematic, since MATCH_READS would take inputs from two channels, the instance of samples may not match.

      MATCH_READS (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.corrected_barcode, CORRECT_BARCODE.out.read1_fastq, CORRECT_BARCODE.out.read2_fastq)

      ADD_BARCODE_TO_READS (MATCH_READS.out.sample_name, MATCH_READS.out.barcode_fastq, MATCH_READS.out.read1_fastq, MATCH_READS.out.read2_fastq)
    }

    // module: trimming off adapter
    CUTADAPT (ADD_BARCODE_TO_READS.out.sample_name, ADD_BARCODE_TO_READS.out.read1_fastq, ADD_BARCODE_TO_READS.out.read2_fastq, params.read1_adapter, params.read2_adapter)

    // module: mapping with bwa or minimap2: mark duplicate
    // bwa or minimap2
    if (params.mapper == 'bwa') {
      log.info "INFO: --mapper: bwa"

      if (!params.ref_bwa_index) {
        log.info "INFO: --ref_bwa_index not provided, checking --ref_fasta and --ref_fasta_ucsc/--ref_fasta_ensembl ..."

        if (params.ref_fasta) {
          log.info "INFO: --ref_fasta provided, use it for building bwa index."

          // module : bwa_index
          BWA_INDEX (params.ref_fasta)
          // mapping with the built index
        } else if (params.ref_fasta_ucsc) {
          // exit 1, 'WARNING: --ref_fasta_ucsc is not supported yet, pls use --ref_fasta_ensembl!'
          log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build bwa index, and map with bwa ..."

          // module : download_from_ucsc
          DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc)
          // module : extract primary sequence
          GET_PRIMARY_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta)
          // module : bwa_index
          BWA_INDEX (GET_PRIMARY_GENOME.out.genome_fasta)
        } else if (params.ref_fasta_ensembl) {
          log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."

          // module : download_from_ucsc
          DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, params.ensembl_release)
          // module : bwa_index
          BWA_INDEX (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta)
        } else {
          exit 1, 'Parameter --ref_fasta_ucsc/--ref_fasta_ensembl: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
        }
        // module : bwa_map
        BWA_MAP (CUTADAPT.out.sample_name, CUTADAPT.out.trimed_read1_fastq, CUTADAPT.out.trimed_read2_fastq, BWA_INDEX.out.bwa_index_folder)
      } else {
        // use user provided bwa index for mapping
        // module : bwa_map
        BWA_MAP (CUTADAPT.out.sample_name, CUTADAPT.out.trimed_read1_fastq, CUTADAPT.out.trimed_read2_fastq, params.ref_bwa_index)
      }
    } else if (params.mapper == "minimap2") {
      log.info "INFO: --mapper: minimap2"

      if (!params.ref_minimap2_index) {
        log.info "INFO: --ref_minimap2_index not provided, check --ref_fasta and --ref_fasta_uscs/--ref_fasta_ensembl ..."

        if (params.ref_fasta) {
          log.info "INFO: --ref_fasta provided, use it to build minimap2 index."

          // module : bwa_index
          MINIMAP2_INDEX (params.ref_fasta)
          // mapping with the built index
        } else if (params.ref_fasta_ucsc) {
          log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build minimap2 index, and map with minimap2 ..."

          // module : download_from_ucsc
          DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc)
          // module : get_primary_genome
          GET_PRIMARY_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta)
          // module : bwa_index
          MINIMAP2_INDEX (GET_PRIMARY_GENOME.out.genome_fasta)
        } else if (params.ref_fasta_ensembl) {
          log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."

          // module : download_from_ensembl
          DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, params.ensembl_release)
          // module : bwa_index
          MINIMAP2_INDEX (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta)
        } else {
          exit 1, 'Parameter --ref_fasta_ucsc/--ref_fasta_ensembl: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
        }
        // module : minimap2_map
        MINIMAP2_MAP (CUTADAPT.out.sample_name, CUTADAPT.out.trimed_read1_fastq, CUTADAPT.out.trimed_read2_fastq, MINIMAP2_INDEX.out.minimap2_index)
      } else {
        // use user provided bwa index for mapping
        // module : minimap2_map
        MINIMAP2_MAP (CUTADAPT.out.sample_name, CUTADAPT.out.trimed_read1_fastq, CUTADAPT.out.trimed_read2_fastq, params.ref_minimap2_index)
      }
    } else {
      exit 1, 'Parameter --mapper: pls supply a mapper to use, eiter bwa or minimap2!'
    }

    // module: bamqc with qualimap
    if (params.mapper == 'bwa') {
      QUALIMAP (BWA_MAP.out.sample_name, BWA_MAP.out.bam)
    } else if (params.mapper == "minimap2") {
      QUALIMAP (MINIMAP2_MAP.out.sample_name, MINIMAP2_MAP.out.bam)
    }

    // module: generate fragment file with sinto
    if (params.mapper == 'bwa') {
      GET_FRAGMENTS (BWA_MAP.out.sample_name, BWA_MAP.out.bam)
    } else if (params.mapper == "minimap2") {
      GET_FRAGMENTS (MINIMAP2_MAP.out.sample_name, MINIMAP2_MAP.out.bam)
    }

    // module: generate fragement file with sinto
  } else if (params.preprocess == "10xgenomics") {
    // log.info "INFO: --preprocess: 10xgenomics"
    // Module: prepare fastq folder
    GET_10XGENOMICS_FASTQ (ch_samplesheet)

    if (params.ref_cellranger == "") {
      log.info "Parameter --ref_cellranger not supplied, checking --ref_cellranger_ucsc/--ref_cellranger_ensembl!"

      if (params.ref_cellranger_ucsc) {
        log.info "Parameter --ref_cellranger_ucsc provided, will download genome, gtf, and build index with cellranger-atac."

        // Module: download ucsc genome
        DOWNLOAD_FROM_UCSC (params.ref_cellranger_ucsc)
        // Module: download ucsc gtf
        DOWNLOAD_FROM_UCSC_GTF (params.ref_cellranger_ucsc)
        // Module: fix gtf
        FIX_UCSC_GTF (DOWNLOAD_FROM_UCSC_GTF.out.gtf)
        // Module: extract primary genome
        GET_PRIMARY_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta)
        // Module: prepare cellranger index
        CELLRANGER_INDEX (GET_PRIMARY_GENOME.out.genome_fasta, FIX_UCSC_GTF.out.gtf, DOWNLOAD_FROM_UCSC.out.genome_name)
      } else if (params.ref_cellranger_ensembl) {
        // Module: download ensembl genome
        DOWNLOAD_FROM_ENSEMBL (params.ref_cellranger_ensembl, params.ensembl_release)
        // Module: download ensembl gtf
        DOWNLOAD_FROM_ENSEMBL_GTF (params.ref_cellranger_ensembl, params.ensembl_release)
        // Module: prepare cellranger index
        CELLRANGER_INDEX (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL_GTF.out.gtf, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
      }
      // Module: run cellranger-atac count
      CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.fastq_folder, CELLRANGER_INDEX.out.index)
    } else {
      log.info "Parameter --ref_cellranger supplied, will use it as index folder."
      CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.fastq_folder, params.ref_cellranger)
      // sample_name = PARSEUMI.out.umi.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()
      // sample_fastq_folder = GET_10XGENOMICS_FASTQ.out.fastq.to
    }

} else if (params.preprocess == "biorad") {
    // log.info "INFO: --preprocess: biorad"
    // log.info "INFO: must use biorad compatible sequencing data!"
    if (params.ref_bwa_index == "") {
      exit 1, 'Parameter --ref_bwa_index: pls supply full path to bwa index folder!'
    }
    // if (params.ref_bwa_fasta == "") {
    if (params.ref_fasta == "") {
      exit 1, 'Parameter --ref_fasta: pls supply full path to reference fasta file!'
    }
    if (params.biorad_genome == "") {
      exit 1, 'Parameter --biorad_genome: pls choose from "hg19", "hg38","mm10", or "hg19-mm10"!'
    }

    GET_BIORAD_FASTQ (ch_samplesheet)
    BIORAD_FASTQC (GET_BIORAD_FASTQ.out.sample_name, GET_BIORAD_FASTQ.out.fastq_folder)
    BIORAD_ATAC_SEQ_DEBARCODE (GET_BIORAD_FASTQ.out.sample_name, GET_BIORAD_FASTQ.out.fastq_folder)
    // Note that BIORAD_ATAC_SEQ_TRIM_READS must be performed after debarcode.
    BIORAD_ATAC_SEQ_TRIM_READS (BIORAD_ATAC_SEQ_DEBARCODE.out.sample_name, BIORAD_ATAC_SEQ_DEBARCODE.out.debarcoded_reads)
    BIORAD_ATAC_SEQ_BWA (BIORAD_ATAC_SEQ_TRIM_READS.out.sample_name, BIORAD_ATAC_SEQ_TRIM_READS.out.trimmed_reads, params.ref_bwa_index)
    BIORAD_ATAC_SEQ_ALIGNMENT_QC (BIORAD_ATAC_SEQ_BWA.out.sample_name, BIORAD_ATAC_SEQ_BWA.out.alignments, params.ref_fasta)
    BIORAD_ATAC_SEQ_FILTER_BEADS (BIORAD_ATAC_SEQ_BWA.out.sample_name, BIORAD_ATAC_SEQ_BWA.out.alignments, params.biorad_genome)
  } else {
    // log.info "ERROR: for parameter --preprocess, choose from default, 10xgenomics, biorad."
  }
}

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
}

workflow DOWNSTREAM {

    ch_software_versions = Channel.empty()
    log.info "INFO: --downstream: ArchR"
    // Module: create ArrowFile
    // ARCHR_CREATE_ARROWFILES(params.sample_name, params.fragment, params.archr_genome, params.archr_thread)
    ARCHR_CREATE_ARROWFILES(ch_samplesheet_archr, params.archr_genome, params.archr_thread)

    // Module: add DoubletScores
    ARCHR_ADD_DOUBLETSCORES(ARCHR_CREATE_ARROWFILES.out.sample_name, ARCHR_CREATE_ARROWFILES.out.arrowfile)

    // Module: create ArchRProject
    // ARCHR_ARCHRPROJECT(ARCHR_ADD_DOUBLETSCORES.out.sample_name, params.archr_genome, params.archr_thread, ARCHR_ADD_DOUBLETSCORES.out.arrowfile) // Note, ARCH_ADD_DOUBLETSCORES will modify arrowfile in place, therefore, in ARCHR_ARCHRPROJECT, must use the arrowfile generated from ARCHR_ADD_DOUBLETSCORES, otherwise, ARCHR_CREATE_ARROWFILES generate arrowfile will be updated each time ARCHR_ADD_DOUBLETSCORES, so that the -resume won't work for ARCHR_ARCHRPROJECT as long as ARCHR_ADD_DOUBLETSCORES runs.
    // ARCHR_ARCHRPROJECT(ARCHR_ADD_DOUBLETSCORES.out.sample_name, params.archr_genome, params.archr_thread, ARCHR_ADD_DOUBLETSCORES.out.arrowfile)
    ch_samplename_list = ARCHR_ADD_DOUBLETSCORES.out.sample_name.toSortedList( { a, b -> a.getName() <=> b.getName() })
    ch_arrowfile_list = ARCHR_ADD_DOUBLETSCORES.out.arrowfile.toSortedList( { a, b -> a.getName() <=> b.getName() })
    //
    ch_samplename_list.view()
    ch_arrowfile_list.view()

    ch_arrowfile_list.join(",").view()

    // ARCHR_ADD_DOUBLETSCORES.out.arrowfile.collect().toSortedList

    // ARCHR_ARCHRPROJECT()

    // NOTE: should use collect to collect all samples, and merge them into one ArchRProject.

    // Module: ArchRProject QC
    // ARCHR_ARCHRPROJECT_QC(ARCHR_ARCHRPROJECT.out.sample_name, ARCHR_ARCHRPROJECT.out.archr_project)






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
