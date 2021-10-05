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
log.info "TEST: " + params.preprocess + params.input_preprocess

if (params.preprocess) {
  if (params.input_preprocess) {
    ch_input = file(params.input_preprocess)
  } else {
      exit 1, 'Input samplesheet not specified!'
  }
} else {
    if (params.input_archr) {
    } else {
        exit 1, "--input_archr samplesheet must be specified!"
    }
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
include { get_bsgenome } from '../modules/local/functions'

// Modules: local
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['csv':'']] )
include { GET_10XGENOMICS_FASTQ } from '../modules/local/get_10xgenomics_fastq'   addParams( options: modules['get_10xgenomics_fastq'] )
include { GET_BIORAD_FASTQ      } from '../modules/local/get_biorad_fastq'        addParams( options: modules['get_biorad_fastq'] )

include { CELLRANGER_ATAC_COUNT } from '../modules/local/cellranger_atac_count'   addParams( options: modules['cellranger_atac_count'] )
include { CORRECT_BARCODE       } from '../modules/local/correct_barcode'         addParams( options: modules['correct_barcode'] )
include { CORRECT_BARCODE_PHENIQS } from '../modules/local/correct_barcode_pheniqs' addParams( options: modules['correct_barcode_pheniqs'] )
include { MATCH_READS           } from '../modules/local/match_reads'             addParams( options: modules['match_reads'] )
include { MATCH_READS_TRIMMED   } from '../modules/local/match_reads_trimmed'     addParams( options: modules['match_reads_trimmed'] )
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

include { BAM_FILTER       } from '../modules/local/bam_filter'    addParams( options: modules['bam_filter'] )
include { REMOVE_DUPLICATE } from '../modules/local/remove_duplicate'    addParams( options: modules['remove_duplicate'] )
include { QUALIMAP         } from '../modules/local/qualimap'    addParams( options: modules['qualimap'] )
include { GET_FRAGMENTS    } from '../modules/local/get_fragments'    addParams( options: modules['get_fragments'] )

include { DOWNLOAD_FROM_UCSC_GTF } from '../modules/local/download_from_ucsc_gtf'    addParams( options: modules['download_from_ucsc_gtf'] )
include { FIX_UCSC_GTF } from '../modules/local/fix_ucsc_gtf'    addParams( options: modules['fix_ucsc_gtf'] )
include { DOWNLOAD_FROM_ENSEMBL_GTF } from '../modules/local/download_from_ensembl_gtf'    addParams( options: modules['download_from_ensembl_gtf'] )
include { CELLRANGER_INDEX } from '../modules/local/cellranger_index'             addParams( options: modules['cellranger_index'] )

// For ArchR functions:
include { ARCHR_GET_ANNOTATION } from '../modules/local/archr_get_annotation' addParams( options: modules['archr_get_annotation'] )
include { ARCHR_GET_ANNOTATION_CUSTOM } from '../modules/local/archr_get_annotation_custom' addParams( options: modules['archr_get_annotation_custom'] )
include { ARCHR_CREATE_ARROWFILES } from '../modules/local/archr_create_arrowfiles' addParams( options: modules['archr_create_arrowfiles'] )
include { ARCHR_CREATE_ARROWFILES_ANNOTATION } from '../modules/local/archr_create_arrowfiles_annotation' addParams( options: modules['archr_create_arrowfiles_annotation'] )
include { ARCHR_ADD_DOUBLETSCORES } from '../modules/local/archr_add_doubletscores' addParams( options: modules['archr_add_doubletscores'] )
include { ARCHR_ARCHRPROJECT } from '../modules/local/archr_archrproject' addParams( options: modules['archr_archrproject'] )
include { ARCHR_ARCHRPROJECT_ANNOTATION } from '../modules/local/archr_archrproject_annotation' addParams( options: modules['archr_archrproject_annotation'] )
include { ARCHR_ARCHRPROJECT_QC } from '../modules/local/archr_archrproject_qc' addParams( options: modules['archr_archrproject_qc'] )
include { ARCHR_DIMENSION_REDUCTION } from '../modules/local/archr_dimension_reduction' addParams( options: modules['archr_dimension_reduction'] )
include { ARCHR_BATCH_CORRECTION } from '../modules/local/archr_batch_correction' addParams( options: modules['archr_batch_correction'] )
include { ARCHR_CLUSTERING } from '../modules/local/archr_clustering' addParams( options: modules['archr_clustering'] )
include { ARCHR_EMBEDDING } from '../modules/local/archr_embedding' addParams( options: modules['archr_embedding'] )
include { ARCHR_MARKER_GENE } from '../modules/local/archr_marker_gene' addParams( options: modules['archr_marker_gene'] )
include { ARCHR_SCRNASEQ_UNCONSTRAINED } from '../modules/local/archr_scrnaseq_unconstrained' addParams( options: modules['archr_scrnaseq_unconstrained'] )
include { ARCHR_SCRNASEQ_CONSTRAINED } from '../modules/local/archr_scrnaseq_constrained' addParams( options: modules['archr_scrnaseq_constrained'] )
include { ARCHR_PSEUDO_BULK } from '../modules/local/archr_pseudo_bulk' addParams( options: modules['archr_pseudo_bulk'] )
include { ARCHR_PSEUDO_BULK_CLUSTERS } from '../modules/local/archr_pseudo_bulk_clusters' addParams( options: modules['archr_pseudo_bulk_clusters'] )
include { ARCHR_PSEUDO_BULK_CLUSTERS2 } from '../modules/local/archr_pseudo_bulk_clusters2' addParams( options: modules['archr_pseudo_bulk_clusters2'] )
include { ARCHR_CALL_PEAKS } from '../modules/local/archr_call_peaks' addParams( options: modules['archr_call_peaks'] )
include { ARCHR_CALL_PEAKS_CLUSTERS } from '../modules/local/archr_call_peaks_clusters' addParams( options: modules['archr_call_peaks_clusters'] )
include { ARCHR_CALL_PEAKS_CLUSTERS2 } from '../modules/local/archr_call_peaks_clusters2' addParams( options: modules['archr_call_peaks_clusters2'] )
include { ARCHR_GET_MARKER_PEAKS_CLUSTERS } from '../modules/local/archr_get_marker_peaks_clusters' addParams( options: modules['archr_get_marker_peaks_clusters'] )
include { ARCHR_GET_MARKER_PEAKS_CLUSTERS2 } from '../modules/local/archr_get_marker_peaks_clusters2' addParams( options: modules['archr_get_marker_peaks_clusters2'] )
include { ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS } from '../modules/local/archr_marker_peaks_in_tracks_clusters' addParams( options: modules['archr_marker_peaks_in_tracks_clusters'] )
include { ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2 } from '../modules/local/archr_marker_peaks_in_tracks_clusters2' addParams( options: modules['archr_marker_peaks_in_tracks_clusters2'] )
include { ARCHR_PAIRWISE_TEST_CLUSTERS } from '../modules/local/archr_pairwise_test_clusters' addParams( options: modules['archr_pairwise_test_clusters'] )
include { ARCHR_PAIRWISE_TEST_CLUSTERS2 } from '../modules/local/archr_pairwise_test_clusters2' addParams( options: modules['archr_pairwise_test_clusters2'] )
include { ARCHR_MOTIF_ENRICHMENT_CLUSTERS } from '../modules/local/archr_motif_enrichment_clusters' addParams( options: modules['archr_motif_enrichment_clusters'] )
include { ARCHR_MOTIF_ENRICHMENT_CLUSTERS2 } from '../modules/local/archr_motif_enrichment_clusters2' addParams( options: modules['archr_motif_enrichment_clusters2'] )
include { ARCHR_MOTIF_DEVIATIONS_CLUSTERS } from '../modules/local/archr_motif_deviations_clusters' addParams( options: modules['archr_motif_deviations_clusters'] )
include { ARCHR_MOTIF_DEVIATIONS_CLUSTERS2 } from '../modules/local/archr_motif_deviations_clusters2' addParams( options: modules['archr_motif_deviations_clusters2'] )
include { ARCHR_FOOTPRINTING_CLUSTERS } from '../modules/local/archr_footprinting_clusters' addParams( options: modules['archr_footprinting_clusters'] )
include { ARCHR_FOOTPRINTING_CLUSTERS2 } from '../modules/local/archr_footprinting_clusters2' addParams( options: modules['archr_footprinting_clusters2'] )
include { ARCHR_COACCESSIBILITY_CLUSTERS } from '../modules/local/archr_coaccessibility_clusters' addParams( options: modules['archr_coaccessibility_clusters'] )
include { ARCHR_COACCESSIBILITY_CLUSTERS2 } from '../modules/local/archr_coaccessibility_clusters2' addParams( options: modules['archr_coaccessibility_clusters2'] )
include { ARCHR_PEAK2GENELINKAGE_CLUSTERS2 } from '../modules/local/archr_peak2genelinkage_clusters2' addParams( options: modules['archr_peak2genelinkage_clusters2'] )
include { ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS } from '../modules/local/archr_get_positive_tf_regulator_clusters' addParams( options: modules['archr_get_positive_tf_regulator_clusters'] )
include { ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2 } from '../modules/local/archr_get_positive_tf_regulator_clusters2' addParams( options: modules['archr_get_positive_tf_regulator_clusters2'] )
include { ARCHR_TRAJECTORY_CLUSTERS2 } from '../modules/local/archr_trajectory_clusters2' addParams( options: modules['archr_trajectory_clusters2'] )
include { ARCHR_GET_CLUSTERING_TSV } from '../modules/local/archr_get_clustering_tsv' addParams( options: modules['archr_get_clustering_tsv'] )

// // Modules: nf-core/modules
// include { FASTQC                } from '../modules/nf-core/software/fastqc/main'  addParams( options: modules['fastqc']            )
// include { MULTIQC               } from '../modules/nf-core/software/multiqc/main' addParams( options: multiqc_options              )
//
// // Subworkflows: local
// include { INPUT_CHECK           } from '../subworkflows/local/input_check'        addParams( options: [:]                          )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow PREPROCESS {
  take:
    ch_samplesheet

  main:
    // Log: check and see which parameter modules are specified
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

    // Perform preprocess accordingly
    if (params.preprocess == "default") {
      if (!(params.barcode_whitelist)) {
        log.info "NOTICE: --barcode_whitelist: not supplied, skip barcode correction!"
      }
      // log.info "INFO(2): --preprocess: default"
      GET_10XGENOMICS_FASTQ (ch_samplesheet)
      // module: fastQC
      FASTQC (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)

      // module: barcode correction (optional) and add barcode: correct barcode fastq given whitelist and barcode fastq file
      if (!(params.barcode_whitelist)) {
        // log.info "NOTICE(2): --barcode_whitelist: not supplied, skip barcode correction!"
        ADD_BARCODE_TO_READS (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
      } else {
        // Allow users to choose from barcode_correction.R or Pheniqs:
        if (params.barcode_correction == "naive") {
          CORRECT_BARCODE (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, params.barcode_whitelist, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
        // MATCH_READS (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.corrected_barcode, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
        // Note that the above might be problematic, since MATCH_READS would take inputs from two channels, the instance of samples may not match.
          MATCH_READS (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.corrected_barcode, CORRECT_BARCODE.out.read1_fastq, CORRECT_BARCODE.out.read2_fastq)

          ADD_BARCODE_TO_READS (MATCH_READS.out.sample_name, MATCH_READS.out.barcode_fastq, MATCH_READS.out.read1_fastq, MATCH_READS.out.read2_fastq)
        } else {
          // use pheniqs:
          CORRECT_BARCODE_PHENIQS (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, params.barcode_whitelist, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)

          MATCH_READS (CORRECT_BARCODE_PHENIQS.out.sample_name, CORRECT_BARCODE_PHENIQS.out.corrected_barcode, CORRECT_BARCODE_PHENIQS.out.read1_fastq, CORRECT_BARCODE_PHENIQS.out.read2_fastq)
        }

}

      // module: trimming off adapter
      if ((params.barcode_whitelist) && (params.barcode_correction == "pheniqs")) {
        CUTADAPT (MATCH_READS.out.sample_name, MATCH_READS.out.read1_fastq, MATCH_READS.out.read2_fastq, params.read1_adapter, params.read2_adapter)
      } else {
        CUTADAPT (ADD_BARCODE_TO_READS.out.sample_name, ADD_BARCODE_TO_READS.out.read1_fastq, ADD_BARCODE_TO_READS.out.read2_fastq, params.read1_adapter, params.read2_adapter)
      }

      // module: MATCH_READS_TRIMMED: in case user choose to trim based on quality and read pair gets unbalanced.
      MATCH_READS_TRIMMED (CUTADAPT.out.sample_name, CUTADAPT.out.trimed_read1_fastq, CUTADAPT.out.trimed_read2_fastq)

      // TODO: fragment generation should take input from deduplcated bam files so that fragments will be garanteed unique. Make the duplicate removal as a default cause it won't hurt. Otherwise visualizaiton is problematic.
      // For fragement files, sinto will take care of the duplication because each line is unique? (need confirm).
      // After mapping, if we need to split the bam files based on clusters, we need to remove duplicates, which is determined by the barcode information: PICARD or SAMTOOLS.
      // Duplication should be determined: barcode (cell), start and end coordination, should also shift the soft-trimming starting position.

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
          BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, BWA_INDEX.out.bwa_index_folder)
        } else {
          // use user provided bwa index for mapping, module : bwa_map
          BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, params.ref_bwa_index)
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
          MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, MINIMAP2_INDEX.out.minimap2_index)
        } else {
            // use user provided bwa index for mapping
            // module : minimap2_map
            MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, params.ref_minimap2_index)
        }
      } else {
          exit 1, 'Parameter --mapper: pls supply a mapper to use, eiter bwa or minimap2!'
      }

      // module: filter out poorly mapped reads
      if (params.mapper == 'bwa') {
        BAM_FILTER (BWA_MAP.out.sample_name, BWA_MAP.out.bam, params.filter)
      } else if (params.mapper == "minimap2") {
          BAM_FILTER (MINIMAP2_MAP.out.sample_name, MINIMAP2_MAP.out.bam, params.filter)
      }

      // module: remove duplicates based on cell barcode, start, end
      REMOVE_DUPLICATE(BAM_FILTER.out.sample_name, BAM_FILTER.out.bam)

      // module: bamqc with qualimap for raw bam files
      // TODO: Run Qualimap on the final filtered deduplicated bam file.
      QUALIMAP (REMOVE_DUPLICATE.out.sample_name, REMOVE_DUPLICATE.out.bam)

      // module: generate fragment file with sinto
      // use raw bam file since ArchR may take advantage of the duplication info.
      GET_FRAGMENTS (BAM_FILTER.out.sample_name, BAM_FILTER.out.bam)
    } else if (params.preprocess == "10xgenomics") {
        if (!params.ref_cellranger) {
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
            // Module: prepare fastq folder
            GET_10XGENOMICS_FASTQ (ch_samplesheet)
            // Module: run cellranger-atac count
            CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.fastq_folder, CELLRANGER_INDEX.out.index_folder)
          } else if (params.ref_cellranger_ensembl) {
            // Module: download ensembl genome
            DOWNLOAD_FROM_ENSEMBL (params.ref_cellranger_ensembl, params.ensembl_release)
            // Module: download ensembl gtf
            DOWNLOAD_FROM_ENSEMBL_GTF (params.ref_cellranger_ensembl, params.ensembl_release)
            // Module: prepare cellranger index
            CELLRANGER_INDEX (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL_GTF.out.gtf, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
            // Module: prepare fastq folder
            GET_10XGENOMICS_FASTQ (ch_samplesheet)
            // Module: run cellranger-atac count
            CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.fastq_folder, CELLRANGER_INDEX.out.index_folder)
            } else {
                exit 1, "--ref_cellranger_ucsc/--ref_cellranger_ensembl must be specified!"
              }
        } else {
          log.info "Parameter --ref_cellranger supplied, will use it as index folder."
          GET_10XGENOMICS_FASTQ (ch_samplesheet)
          CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.fastq_folder, params.ref_cellranger)
          // sample_name = PARSEUMI.out.umi.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()
          // sample_fastq_folder = GET_10XGENOMICS_FASTQ.out.fastq.to
        }
    } else if (params.preprocess == "biorad") {
      // log.info "INFO: --preprocess: biorad"
      // log.info "INFO: must use biorad compatible sequencing data!"
      if (!params.ref_bwa_index) {
        exit 1, 'Parameter --ref_bwa_index: pls supply full path to bwa index folder!'
      }
      if (!params.ref_fasta) {
        exit 1, 'Parameter --ref_fasta: pls supply full path to reference fasta file!'
      }
      if (!params.biorad_genome) {
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
    }

    // Collect all output results for MultiQC report:
    res_files = Channel.empty()
    // res_files = res_files.mix(Channel.from(ch_multiqc_config))
    // res_files = res_files.mix(Channel.from(ch_multiqc_custom_config).collect().ifEmpty([]))

    // Use try-catch since if certain module is not run, module.out becomes undefined.
    // FASTQC module:
    try {
      res_files = res_files.mix(FASTQC.out.zip.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // CORRECT_BARCODE module:
    try {
      res_files = res_files.mix(CORRECT_BARCODE.out.corrected_barcode_summary.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // CORRECT_BARCODE_PHENIQS module:
    try {
      res_files = res_files.mix(CORRECT_BARCODE_PHENIQS.out.corrected_barcode_summary.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // CUTADAPT module:
    try {
      res_files = res_files.mix(CUTADAPT.out.log.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // QUALIMAP module:
    try {
      res_files = res_files.mix(QUALIMAP.out.bamqc.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // CELLRANGER_INDEX module:
    try {
      res_files.mix(CELLRANGER_INDEX.out.bwa_index_folder.collect().ifEmpty([]))
    } catch (Exception ex) {
    }
    // CELLRANGER_ATAC_COUNT module:
    try {
      res_files = res_files.mix(CELLRANGER_ATAC_COUNT.out.cellranger_atac_count.collect().ifEmpty([]))
    } catch (Exception ex) {
    }

    // if (params.preprocess == "default") {
    //     res_files = res_files.mix(FASTQC.out.zip.collect().ifEmpty([]))
    //   if (params.barcode_correction == "pheniqs") {
    //     res_files = res_files.mix(CORRECT_BARCODE_PHENIQS.out.corrected_barcode_summary.collect().ifEmpty([]))
    //   } else if (params.barcode_correction == "naive") {
    //     res_files = res_files.mix(CORRECT_BARCODE.out.corrected_barcode_summary.collect().ifEmpty([]))
    //   }
    //   res_files = res_files.mix(CUTADAPT.out.log.collect().ifEmpty([]))
    //   res_files = res_files.mix(QUALIMAP.out.bamqc.collect().ifEmpty([]))
    // } else if (params.preprocess == "10xgenomics") {
    //   res_files = res_files.mix(CELLRANGER_ATAC_COUNT.out.cellranger_atac_count.collect().ifEmpty([]))
    //   if (!params.ref_cellranger) {
    //     res_files.mix(CELLRANGER_INDEX.out.bwa_index_folder.collect().ifEmpty([]))
    //   }
    // }

  emit:
    res_files // out[0]: res folders for MultiQC report
    GET_FRAGMENTS.out.fragments // out[1]: for split bed
    GET_FRAGMENTS.out.ch_fragment // out[2]: fragment ch for ArchR
    REMOVE_DUPLICATE.out.sample_name // out[3]: for split bam
    REMOVE_DUPLICATE.out.bam // out[4]: for split bam
}

workflow DOWNSTREAM {
  take:
    ch_samplesheet_archr

  main:
    ch_software_versions = Channel.empty()
    log.info "INFO: --downstream: ArchR"
    // Module: check if ArchR genome matches with preprocess genome, and create custome Genome if needed.
    (bsgenome, genome_status) = get_bsgenome(params.archr_genome, params.archr_custom_genome, params.archr_txdb, params.archr_org, params.archr_bsgenome, params.ref_fasta_ucsc, params.ref_fasta_ensembl, params.ref_cellranger_ucsc, params.ref_cellranger_ensembl)

    // Module: createArrowFile and addDoubletScores
    if (["custom"].contains(genome_status)) {
      log.info "INFO: ArchR will build gene/genomeAnnotation files with custom TxDb, Org, and BSgenome files supplied by user."

      ARCHR_GET_ANNOTATION_CUSTOM(params.archr_txdb, params.archr_org, params.archr_bsgenome)
      ARCHR_CREATE_ARROWFILES_ANNOTATION(ch_samplesheet_archr, ARCHR_GET_ANNOTATION_CUSTOM.out.geneAnnotation, ARCHR_GET_ANNOTATION_CUSTOM.out.genomeAnnotation, ARCHR_GET_ANNOTATION_CUSTOM.out.user_rlib, params.archr_thread)
      // Module: add DoubletScores
      ARCHR_ADD_DOUBLETSCORES(ARCHR_CREATE_ARROWFILES_ANNOTATION.out.sample_name, ARCHR_CREATE_ARROWFILES_ANNOTATION.out.arrowfile)
      ch_samplename_list = ARCHR_ADD_DOUBLETSCORES.out.sample_name.toSortedList()
      ch_arrowfile_list = ARCHR_ADD_DOUBLETSCORES.out.arrowfile.toSortedList( { a, b -> a.getName() <=> b.getName() })
    } else if (["ready", "ready_ucsc", "ready_ensembl"].contains(genome_status)) {
      // "ready" means ArchR natively supported genome
      if (genome_status == "ready_ensembl") {
        // log.info "INFO: ArchR will use natively supported ArchR genome (though ensembl genome supplied): " + bsgenome
        exit 1, "INFO: ensembl genome supplied, need to build ArchR genome first, pls supply corresponding TxDb, Org, and BSgenome objects via --txdb, --org, and --bsgenome parameters."
      } else {
        log.info "INFO: ArchR will use natively supported ArchR genome: " + bsgenome
      }
      ARCHR_CREATE_ARROWFILES(ch_samplesheet_archr, bsgenome, params.archr_thread)
      // Module: add DoubletScores
      ARCHR_ADD_DOUBLETSCORES(ARCHR_CREATE_ARROWFILES.out.sample_name, ARCHR_CREATE_ARROWFILES.out.arrowfile)
      ch_samplename_list = ARCHR_ADD_DOUBLETSCORES.out.sample_name.toSortedList()
      ch_arrowfile_list = ARCHR_ADD_DOUBLETSCORES.out.arrowfile.toSortedList( { a, b -> a.getName() <=> b.getName() })
    } else if (["need_build", "need_build_ucsc", "need_build_ensembl"].contains(genome_status)) {
        // "need_build" means ArchR needs to build gene/genomeAnnotation files first
        if (genome_status == "need_build_ensembl") {
          // log.info "INFO: ArchR will build gene/genomeAnnotation files with (though ensembl genome supplied): " + bsgenome
          exit 1, "INFO: ensembl genome supplied, need to build ArchR genome first, pls supply corresponding TxDb, Org, and BSgenome objects via --txdb, --org, and --bsgenome parameters."
        } else {
          log.info "INFO: ArchR will build gene/genomeAnnotation files with: " + bsgenome
        }
        ARCHR_GET_ANNOTATION(params.archr_genome)
        ARCHR_CREATE_ARROWFILES_ANNOTATION(ch_samplesheet_archr, ARCHR_GET_ANNOTATION.out.geneAnnotation, ARCHR_GET_ANNOTATION.out.genomeAnnotation, ARCHR_GET_ANNOTATION.out.user_rlib, params.archr_thread)
        // Module: add DoubletScores
        ARCHR_ADD_DOUBLETSCORES(ARCHR_CREATE_ARROWFILES_ANNOTATION.out.sample_name, ARCHR_CREATE_ARROWFILES_ANNOTATION.out.arrowfile)
        ch_samplename_list = ARCHR_ADD_DOUBLETSCORES.out.sample_name.toSortedList()
        ch_arrowfile_list = ARCHR_ADD_DOUBLETSCORES.out.arrowfile.toSortedList( { a, b -> a.getName() <=> b.getName() })
    } else {
        exit 1, "Must supply ArchR genomes!"
    }

    // Module: create ArchRProject and ArchRProjectQC
    if (["custom"].contains(genome_status)) {
      ARCHR_ARCHRPROJECT_ANNOTATION(ch_arrowfile_list, ARCHR_GET_ANNOTATION_CUSTOM.out.geneAnnotation, ARCHR_GET_ANNOTATION_CUSTOM.out.genomeAnnotation, ARCHR_GET_ANNOTATION_CUSTOM.out.user_rlib)
      ARCHR_ARCHRPROJECT_QC(ARCHR_ARCHRPROJECT_ANNOTATION.out.archr_project, params.archr_filter_ratio)
    } else if (["ready", "ready_ucsc"].contains(genome_status)) {
      ARCHR_ARCHRPROJECT(ch_arrowfile_list, bsgenome, params.archr_thread)
      ARCHR_ARCHRPROJECT_QC(ARCHR_ARCHRPROJECT.out.archr_project, params.archr_filter_ratio)
    } else if (["need_build", "need_build_ucsc"].contains(genome_status)) {
      ARCHR_ARCHRPROJECT_ANNOTATION(ch_arrowfile_list, ARCHR_GET_ANNOTATION.out.geneAnnotation, ARCHR_GET_ANNOTATION.out.genomeAnnotation, ARCHR_GET_ANNOTATION.out.user_rlib)
      ARCHR_ARCHRPROJECT_QC(ARCHR_ARCHRPROJECT_ANNOTATION.out.archr_project, params.archr_filter_ratio)
    }

    // Module: dimension reduction
    ARCHR_DIMENSION_REDUCTION(ARCHR_ARCHRPROJECT_QC.out.archr_project)

    // Module: batch correction with harmony
    ARCHR_BATCH_CORRECTION(ARCHR_DIMENSION_REDUCTION.out.archr_project)

    // Module: clustering with Seurat's FindClusters() function
    ARCHR_CLUSTERING(ARCHR_BATCH_CORRECTION.out.archr_project)

    // Module: single-cell embeddings
    ARCHR_EMBEDDING(ARCHR_CLUSTERING.out.archr_project)

    // Module: find marker gene
    ARCHR_MARKER_GENE(ARCHR_EMBEDDING.out.archr_project)

    // Module: integrate with matching scRNAseq data
    if (!(params.archr_scrnaseq)) {
      params.groupby_cluster = "Clusters"
      log.info "NOTICE: --archr_scrnaseq: not supplied, skip integrative analysis with scRNA-seq!"
      // ARCHR_PSEUDO_BULK(ARCHR_MARKER_GENE.out.archr_project, params.groupby_cluster)
      ARCHR_PSEUDO_BULK_CLUSTERS(ARCHR_MARKER_GENE.out.archr_project)
      // For each Arrorproject, you can have only one set of peak set unless you copy arrow files and create another arrowproject. That is why we implemented ARCHR_PSEUDO_BULK_CLUSTERS and ARCHR_PSEUDO_BULK_CLUSTERS2
    } else {
        params.groupby_cluster = "Clusters2"
        log.info "NOTICE: --archr_scrnaseq: supplied, will perform integrative analysis with scRNA-seq!"
        ARCHR_PSEUDO_BULK_CLUSTERS(ARCHR_MARKER_GENE.out.archr_project)
        ARCHR_SCRNASEQ_UNCONSTRAINED(ARCHR_MARKER_GENE.out.archr_project, params.archr_scrnaseq)
        // log.info "INFO: use the following cluster names to define --archr_scrnaseq_grouplist."
        ARCHR_SCRNASEQ_UNCONSTRAINED.out.cell_type_scrna
          .splitText()
          .subscribe onNext: { String str -> println "Cluster name from scRNAseq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --archr_scrnaseq_grouplist ***\n"}

        if ((!params.archr_scrnaseq_grouplist)) {
          log.info "NOTICE: --archr_scrnaseq_grouplist: not supplied, skip constrained integration!"
          // ARCHR_PSEUDO_BULK(ARCHR_SCRNASEQ_UNCONSTRAINED.out.archr_project, params.groupby_cluster)
          ARCHR_PSEUDO_BULK_CLUSTERS2(ARCHR_SCRNASEQ_UNCONSTRAINED.out.archr_project)
        } else {
            log.info "NOTICE: --archr_scrnaseq_grouplist: supplied, will perform constrained integration!"
            ARCHR_SCRNASEQ_CONSTRAINED(ARCHR_SCRNASEQ_UNCONSTRAINED.out.archr_project, params.archr_scrnaseq, params.archr_scrnaseq_grouplist)
            // ARCHR_PSEUDO_BULK(ARCHR_SCRNASEQ_CONSTRAINED.out.archr_project, params.groupby_cluster)
            ARCHR_PSEUDO_BULK_CLUSTERS2(ARCHR_SCRNASEQ_CONSTRAINED.out.archr_project)
        }
    }

    // Module: call peaks
    if (params.groupby_cluster == "Clusters") {
      ARCHR_CALL_PEAKS_CLUSTERS(ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_project)
    } else if (params.groupby_cluster == "Clusters2") {
      ARCHR_CALL_PEAKS_CLUSTERS(ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_project)
      ARCHR_CALL_PEAKS_CLUSTERS2(ARCHR_PSEUDO_BULK_CLUSTERS2.out.archr_project)
    }

    // Module: identify marker peaks and perform MA/Volcano plots
    if (params.groupby_cluster == "Clusters") {
      ARCHR_GET_MARKER_PEAKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project)
      ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.group_names
        .splitText()
        .subscribe onNext: { String str -> println "Group name from scATAC-seq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --pairwise_test_clusters_1/2 and --marker_peak_clusters***\n"}
    } else if (params.groupby_cluster == "Clusters2") {
      ARCHR_GET_MARKER_PEAKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project)
      ARCHR_GET_MARKER_PEAKS_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project)
      ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.group_names
        .splitText()
        .subscribe onNext: { String str -> println "Group name from scATAC-seq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --pairwise_test_clusters_1/2 ***\n"}
      ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.group_names
        .splitText()
        .subscribe onNext: { String str -> println "Group name from scRNA-seq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --pairwise_test_clusters2_1/2 and --marker_peak_clusters2***\n"}
    }

    // Module: plot peaks in browser tracks
    if (params.groupby_cluster == "Clusters") {
      if (!(params.marker_peak_geneSymbol && params.marker_peak_clusters)) {
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: not supplied, skip marker peak plotting on browser tracks!"
      } else {
        // Perform plotting
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: supplied, will perform marker peak plotting on browser tracks!"
        ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.marker_peak_geneSymbol, params.marker_peak_clusters)
      }
    } else if (params.groupby_cluster == "Clusters2") {
      if (!(params.marker_peak_geneSymbol && params.marker_peak_clusters)) {
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: not supplied, skip marker peak plotting on browser tracks!"
      } else {
          // Perform plotting
          log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: supplied, will perform marker peak plotting on browser tracks!"
          ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.marker_peak_geneSymbol, params.marker_peak_clusters)
      }

      if (!(params.marker_peak_geneSymbol && params.marker_peak_clusters2)) {
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters2: not supplied, skip marker peak plotting on browser tracks!"
      } else {
        // Perform plotting
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters2: supplied, will perform marker peak plotting on browser tracks!"
        ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project, ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.archr_marker_peaks, params.marker_peak_geneSymbol, params.marker_peak_clusters2)
      }
    }

    // Module: perform pairwise test
    if (params.groupby_cluster == "Clusters") {
      if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
        log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip pairwise plotting!"
      } else {
        // Perform plotting
        log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform pairwise plotting!"
        ARCHR_PAIRWISE_TEST_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2)
      }
    } else if (params.groupby_cluster == "Clusters2") {
        if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
          log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip pairwise plotting!"
      } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform pairwise plotting!"
          ARCHR_PAIRWISE_TEST_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2)
      }

      if (!(params.pairwise_test_clusters2_1 && params.pairwise_test_clusters2_2)) {
        log.info "NOTICE: --pairwise_test_clusters2_1/2: not supplied, skip pairwise plotting!"
      } else {
        // Perform plotting
        log.info "NOTICE: --pairwise_test_clusters2_1/2: supplied, perform pairwise plotting!"
        ARCHR_PAIRWISE_TEST_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project, params.pairwise_test_clusters2_1, params.pairwise_test_clusters2_2)
      }
    }

    // Module: motif enrichment: note that ARCHR_MOTIF_ENRICHMENT_CLUSTERS and ARCHR_MOTIF_ENRICHMENT_CLUSTERS2 are exactly the same except for the outdir name.
    if (params.groupby_cluster == "Clusters") {
      if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
        log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip motif enrichment!"
      } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform motif enrichment!"
          ARCHR_MOTIF_ENRICHMENT_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_PAIRWISE_TEST_CLUSTERS.out.archr_marker_test, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2, params.custom_peaks)
      }
    } else if (params.groupby_cluster == "Clusters2") {
        if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
          log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip motif enrichment!"
        } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform motif enrichment!"
          ARCHR_MOTIF_ENRICHMENT_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_PAIRWISE_TEST_CLUSTERS.out.archr_marker_test, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2, params.custom_peaks)
        }

      if (!(params.pairwise_test_clusters2_1 && params.pairwise_test_clusters2_2)) {
        log.info "NOTICE: --pairwise_test_clusters2_1/2: not supplied, skip motif enrichment!"
      } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters2_1/2: supplied, perform motif enrichment!"
          ARCHR_MOTIF_ENRICHMENT_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project, ARCHR_PAIRWISE_TEST_CLUSTERS2.out.archr_marker_test, ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.archr_marker_peaks, params.pairwise_test_clusters2_1, params.pairwise_test_clusters2_2, params.custom_peaks)
      }
    }

    if ((params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2) || (params.pairwise_test_clusters2_1 && params.pairwise_test_clusters2_2)) {
      // Module: motif deviation,require motif enrichment result
      if (params.groupby_cluster == "Clusters") {
        ARCHR_MOTIF_DEVIATIONS_CLUSTERS(ARCHR_MOTIF_ENRICHMENT_CLUSTERS.out.archr_project, params.custom_peaks)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_MOTIF_DEVIATIONS_CLUSTERS(ARCHR_MOTIF_ENRICHMENT_CLUSTERS.out.archr_project, params.custom_peaks)
          ARCHR_MOTIF_DEVIATIONS_CLUSTERS2(ARCHR_MOTIF_ENRICHMENT_CLUSTERS2.out.archr_project, params.custom_peaks)
      }

      // Module: footprinting
      if (params.groupby_cluster == "Clusters") {
        ARCHR_FOOTPRINTING_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project, ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_dir)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_FOOTPRINTING_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project, ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_dir)
          ARCHR_FOOTPRINTING_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project, ARCHR_PSEUDO_BULK_CLUSTERS2.out.archr_dir)
      }

      // Module: integrative analysis
      // Below are for integrative analysis: co-accessibility; peak2genelinkage; positive TF regulators.
      // Module: co-accessibility (for both clusters and clusters2)
      if (params.groupby_cluster == "Clusters") {
        ARCHR_COACCESSIBILITY_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_COACCESSIBILITY_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
          ARCHR_COACCESSIBILITY_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project)
      }

      // Module: peak2genelinkage: for clusters2 only
      if (params.groupby_cluster == "Clusters2") {
        ARCHR_PEAK2GENELINKAGE_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project)
      }

      // Module: identify "positive" TF-regulators
      if (params.groupby_cluster == "Clusters") {
        ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
          ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project)
      }

      // Module: trajectory analysis: for Clusters2 only
      // TODO: Module: trajectory analysis: for Clusters using Gene Score Matrix
      if (params.groupby_cluster == "Clusters2") {
        if (!params.trajectory_groups) {
          log.info "Parameter --trajectory_groups not supplied, checking trajectory analysis!"
        } else {
            log.info "Parameter --trajectory_groups supplied, will perform trajectory analysis!"
            ARCHR_TRAJECTORY_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project, params.trajectory_groups)
        }
      } else {
        log.info "Parameter --scrnaseq not supplied, skip trajectory analysis!"
      }
    }


    // Module: prepare clustering tsv file for spliting using sinto fragment
    if (params.groupby_cluster == "Clusters") {
      ARCHR_GET_CLUSTERING_TSV(ARCHR_CLUSTERING.out.archr_project, ch_samplesheet_archr, "Clusters")
    } else if (params.groupby_cluster == "Clusters2") {
      ARCHR_GET_CLUSTERING_TSV(ARCHR_PSEUDO_BULK_CLUSTERS2.out.archr_project, ch_samplesheet_archr, "Clusters2")
    }

    // Module: split fragment bed files:
    // ch_tsv = ARCHR_GET_CLUSTERING_TSV.out.tsv // each instance may output more than 1 clustering tsv file.
    // SPLIT_BED(ARCHR_GET_CLUSTERING_TSV.out.tsv, ARCHR_GET_CLUSTERING_TSV.out.fragment)

    // Module: split bam files if available:
    // SPLIT_BAM()

    // Collect all output results for MultiQC report:
    res_files = Channel.empty()
    // res_files = res_files.mix(Channel.from(ch_multiqc_config))
    // res_files = res_files.mix(Channel.from(ch_multiqc_custom_config).collect().ifEmpty([]))

    // ARCHR_CREATE_ARROWFILES module:
    try {
      res_files = res_files.mix(ARCHR_CREATE_ARROWFILES.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_ADD_DOUBLETSCORES
    try {
      res_files = res_files.mix(ARCHR_ADD_DOUBLETSCORES.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_ARCHRPROJECT_QC:
    try {
      res_files = res_files.mix(ARCHR_ARCHRPROJECT_QC.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_CLUSTERING:
    try {
      res_files = res_files.mix(ARCHR_CLUSTERING.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_EMBEDDING:
    try {
      res_files = res_files.mix(ARCHR_EMBEDDING.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MARKER_GENE:
    try {
      res_files = res_files.mix(ARCHR_MARKER_GENE.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_SCRNASEQ_UNCONSTRAINED:
    try {
      res_files = res_files.mix(ARCHR_SCRNASEQ_UNCONSTRAINED.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_SCRNASEQ_CONSTRAINED:
    try {
      res_files = res_files.mix(ARCHR_SCRNASEQ_CONSTRAINED.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_MARKER_PEAKS_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_MARKER_PEAKS_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_PAIRWISE_TEST_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_PAIRWISE_TEST_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_PAIRWISE_TEST_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_PAIRWISE_TEST_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MOTIF_DEVIATIONS_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MOTIF_DEVIATIONS_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_COACCESSIBILITY_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_COACCESSIBILITY_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_COACCESSIBILITY_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_COACCESSIBILITY_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_TRAJECTORY_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_TRAJECTORY_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}

    // Note some module may not run and therefore may not have out and therefore erro
    // res_folders = res_folders.mix(ARCHR_PEAK2GENELINKAGE_CLUSTERS2.out.res_dir.collect().ifEmpty([]))
    // res_folders = res_folders.mix(ARCHR_TRAJECTORY_CLUSTERS2.out.res_dir.collect().ifEmpty([]))

  emit:
    res_files.collect()
    ARCHR_GET_CLUSTERING_TSV.out.res // Here if using collect(), only the first element will be used for split_bed module. For split bed
    ARCHR_GET_CLUSTERING_TSV.out.tsv // for split bam
    // [ARCHR_GET_CLUSTERING_TSV.out.sample_name, ARCHR_GET_CLUSTERING_TSV.out.tsv]
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
