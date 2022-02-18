/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
// NfcoreSchema pops compilation error.

// Validate input parameters
// WorkflowScatacpipe.initialise(params, log)

// TODO: add required input list here:
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

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
// include { get_bsgenome } from '../modules/local/functions'

// include the supported genomes
include { get_genome_ucsc } from '../modules/local/genome_ucsc'
include { get_genome_ensembl } from '../modules/local/genome_ensembl'
genome_ensembl_list = get_genome_ensembl()
genome_ucsc_list = get_genome_ucsc()
if (params.ref_fasta_ensembl) {
  if (!genome_ensembl_list.contains(params.ref_fasta_ensembl)) {
    exit 1, "Pls use --support_genome to show a list of supported genomes!"
  }
}
if (params.ref_fasta_ucsc) {
  if (!genome_ucsc_list.contains(params.ref_fasta_ucsc)) {
    exit 1, "Pls use --support_genome to show a list of supported genomes!"
  }
}
if (params.archr_genome) {
  if (!genome_ensembl_list.contains(params.archr_genome) && !genome_ucsc_list.contains(params.archr_genome)) {
    exit 1, "Pls use --support_genome to show a list of supported genomes!"
  }
}

// Modules: local
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['csv':'']] )
include { GET_10XGENOMICS_FASTQ } from '../modules/local/get_10xgenomics_fastq'   addParams( options: modules['get_10xgenomics_fastq'] )

include { CELLRANGER_ATAC_COUNT } from '../modules/local/cellranger_atac_count'   addParams( options: modules['cellranger_atac_count'] )
include { CORRECT_BARCODE       } from '../modules/local/correct_barcode'         addParams( options: modules['correct_barcode'] )
include { CORRECT_BARCODE_PHENIQS } from '../modules/local/correct_barcode_pheniqs' addParams( options: modules['correct_barcode_pheniqs'] )
include { MATCH_READS           } from '../modules/local/match_reads'             addParams( options: modules['match_reads'] )
include { MATCH_READS_TRIMMED   } from '../modules/local/match_reads_trimmed'     addParams( options: modules['match_reads_trimmed'] )
include { FASTQC                } from '../modules/local/fastqc'                  addParams( options: modules['fastqc'] )

include { ADD_BARCODE_TO_READS       } from '../modules/local/add_barcode_to_reads'    addParams( options: modules['add_barcode_to_reads'] )
include { CUTADAPT         } from '../modules/local/cutadapt'    addParams( options: modules['cutadapt'] )

include { DOWNLOAD_FROM_UCSC; DOWNLOAD_FROM_UCSC as DOWNLOAD_FROM_UCSC2 } from '../modules/local/download_from_ucsc'    addParams( options: modules['download_from_ucsc'] )
include { DOWNLOAD_FROM_ENSEMBL } from '../modules/local/download_from_ensembl'    addParams( options: modules['download_from_ensembl'] )
// can be removed
include { GET_PRIMARY_GENOME        } from '../modules/local/get_primary_genome'    addParams( options: modules['get_primary_genome'] )
include { BWA_INDEX        } from '../modules/local/bwa_index'    addParams( options: modules['bwa_index'] )
include { BWA_MAP          } from '../modules/local/bwa_map'    addParams( options: modules['bwa_map'] )

include { BUILD_BSGENOME } from '../modules/local/build_bsgenome'
include { BUILD_TXDB } from '../modules/local/build_txdb'
include { PREP_GENOME } from '../modules/local/prep_genome'
include { PREP_GTF; PREP_GTF as PREP_GTF_ARCHR } from '../modules/local/prep_gtf'
include { BUILD_GENE_ANNOTATION } from '../modules/local/build_gene_annotation' addParams( options: modules['build_gene_annotation'] )
include { BUILD_GENOME_ANNOTATION } from '../modules/local/build_genome_annotation' addParams( options: modules['build_genome_annotation'] )
include { PREP_FRAGMENT } from '../modules/local/prep_fragment'

include { MINIMAP2_INDEX   } from '../modules/local/minimap2_index'    addParams( options: modules['minimap2_index'] )
include { MINIMAP2_MAP     } from '../modules/local/minimap2_map'    addParams( options: modules['minimap2_map'] )

include { BAM_FILTER       } from '../modules/local/bam_filter'    addParams( options: modules['bam_filter'] )
include { REMOVE_DUPLICATE } from '../modules/local/remove_duplicate'    addParams( options: modules['remove_duplicate'] )
include { QUALIMAP         } from '../modules/local/qualimap'    addParams( options: modules['qualimap'] )
include { GET_FRAGMENTS    } from '../modules/local/get_fragments'    addParams( options: modules['get_fragments'] )

include { DOWNLOAD_FROM_UCSC_GTF; DOWNLOAD_FROM_UCSC_GTF as DOWNLOAD_FROM_UCSC_GTF2 } from '../modules/local/download_from_ucsc_gtf'    addParams( options: modules['download_from_ucsc_gtf'] )
include { FIX_UCSC_GTF } from '../modules/local/fix_ucsc_gtf'    addParams( options: modules['fix_ucsc_gtf'] )
include { DOWNLOAD_FROM_ENSEMBL_GTF; DOWNLOAD_FROM_ENSEMBL_GTF as DOWNLOAD_FROM_ENSEMBL_GTF2 } from '../modules/local/download_from_ensembl_gtf'    addParams( options: modules['download_from_ensembl_gtf'] )
include { CELLRANGER_INDEX } from '../modules/local/cellranger_index'             addParams( options: modules['cellranger_index'] )

// For ArchR functions:
// include { ARCHR_GET_ANNOTATION } from '../modules/local/archr_get_annotation' addParams( options: modules['archr_get_annotation'] )
include { ARCHR_GET_ANNOTATION_BIOC } from '../modules/local/archr_get_annotation_bioc' addParams( options: modules['archr_get_annotation_bioc'] )
include { ARCHR_CREATE_ARROWFILES } from '../modules/local/archr_create_arrowfiles' addParams( options: modules['archr_create_arrowfiles'] )
include { ARCHR_CREATE_ARROWFILES_ANNOTATION } from '../modules/local/archr_create_arrowfiles_annotation' addParams( options: modules['archr_create_arrowfiles_annotation'] )
include { ARCHR_ADD_DOUBLETSCORES } from '../modules/local/archr_add_doubletscores' addParams( options: modules['archr_add_doubletscores'] )
include { ARCHR_ARCHRPROJECT } from '../modules/local/archr_archrproject' addParams( options: modules['archr_archrproject'] )
include { ARCHR_ARCHRPROJECT_ANNOTATION } from '../modules/local/archr_archrproject_annotation' addParams( options: modules['archr_archrproject_annotation'] )
include { ARCHR_FILTER_DOUBLETS } from '../modules/local/archr_filter_doublets' addParams( options: modules['archr_filter_doublets'] )
include { ARCHR_ARCHRPROJECT_QC } from '../modules/local/archr_archrproject_qc' addParams( options: modules['archr_archrproject_qc'] )
include { ARCHR_DIMENSION_REDUCTION } from '../modules/local/archr_dimension_reduction' addParams( options: modules['archr_dimension_reduction'] )
include { ARCHR_BATCH_CORRECTION } from '../modules/local/archr_batch_correction' addParams( options: modules['archr_batch_correction'] )
include { ARCHR_CLUSTERING } from '../modules/local/archr_clustering' addParams( options: modules['archr_clustering'] )
include { ARCHR_EMBEDDING } from '../modules/local/archr_embedding' addParams( options: modules['archr_embedding'] )
include { ARCHR_MARKER_GENE } from '../modules/local/archr_marker_gene' addParams( options: modules['archr_marker_gene'] )
include { ARCHR_SCRNASEQ_UNCONSTRAINED } from '../modules/local/archr_scrnaseq_unconstrained' addParams( options: modules['archr_scrnaseq_unconstrained'] )
include { ARCHR_SCRNASEQ_CONSTRAINED } from '../modules/local/archr_scrnaseq_constrained' addParams( options: modules['archr_scrnaseq_constrained'] )
include { ARCHR_PSEUDO_BULK_CLUSTERS } from '../modules/local/archr_pseudo_bulk_clusters' addParams( options: modules['archr_pseudo_bulk_clusters'] )
include { ARCHR_PSEUDO_BULK_CLUSTERS2 } from '../modules/local/archr_pseudo_bulk_clusters2' addParams( options: modules['archr_pseudo_bulk_clusters2'] )
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

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow PREPROCESS_10XGENOMICS {
  take:
    ch_samplesheet

  main:
    if (params.ref_cellranger) {
      // if cellranger index folder provided:
      log.info "Parameter --ref_cellranger supplied, will use it as index folder."
      GET_10XGENOMICS_FASTQ (ch_samplesheet)
      CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.fastq_folder, params.ref_cellranger)
    } else if (params.ref_fasta) {
      if (params.ref_gtf) {
        log.info "Parameter --ref_fasta/ref_gtf supplied, will build index."
        // Module: prep_genome
        PREP_GENOME (params.ref_fasta, "custom_genome")
        // Module: prep_gtf
        PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, params.ref_gtf)
        // Module: prepare cellranger index
        CELLRANGER_INDEX (PREP_GENOME.out.genome_fasta, PREP_GTF.out.gtf, PREP_GENOME.out.genome_name)
        // Module: prepare fastq folder
        GET_10XGENOMICS_FASTQ (ch_samplesheet)
        // Module: run cellranger-atac count
        CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.fastq_folder, CELLRANGER_INDEX.out.index_folder)
      } else {
        exit 1, "Pls supply --ref_gtf."
      }
    } else if (params.ref_cellranger_ensembl) {
      if (genome_ensembl_list.contains(params.ref_cellranger_ensembl)) {
        // if ensembl name supplied:
        // Module: download ensembl genome
        DOWNLOAD_FROM_ENSEMBL (params.ref_cellranger_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // Module: prep_genome
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // Module: download ensembl gtf
        DOWNLOAD_FROM_ENSEMBL_GTF (params.ref_cellranger_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // Module: prep_gtf
        PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_ENSEMBL_GTF.out.gtf)
        // Module: prepare cellranger index
        CELLRANGER_INDEX (PREP_GENOME.out.genome_fasta, PREP_GTF.out.gtf, PREP_GENOME.out.genome_name)
        // Module: prepare fastq folder
        GET_10XGENOMICS_FASTQ (ch_samplesheet)
        // Module: run cellranger-atac count
        CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.fastq_folder, CELLRANGER_INDEX.out.index_folder)
      } else {
        exit 1, "Pls use --support_genome to show a list of supported genomes!"
      }
    } else if (params.ref_cellranger_ucsc) {
      genome_ucsc_list = get_genome_ucsc()
      if (genome_ucsc_list.contains(params.ref_cellranger_ucsc)) {
        // Module: download ucsc genome
        DOWNLOAD_FROM_UCSC (params.ref_cellranger_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // Module: prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_gtf)
        // Module: extract primary genome
        // GET_PRIMARY_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta)
        // Module: download ucsc gtf
        DOWNLOAD_FROM_UCSC_GTF (params.ref_cellranger_ucsc)
        // Module: prep_gtf
        PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_gtf, DOWNLOAD_FROM_UCSC_GTF.out.gtf)
        // Module: fix gtf
        // FIX_UCSC_GTF (DOWNLOAD_FROM_UCSC_GTF.out.gtf, GET_PRIMARY_GENOME.out.genome_fasta)
        // Module: prepare cellranger index
        CELLRANGER_INDEX (PREP_GENOME.out.genome_fasta, PREP_GTF.out.gtf, PREP_GENOME.out.genome_name)
        // Module: prepare fastq folder
        GET_10XGENOMICS_FASTQ (ch_samplesheet)
        // Module: run cellranger-atac count
        CELLRANGER_ATAC_COUNT (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.fastq_folder, CELLRANGER_INDEX.out.index_folder)
      } else {
        exit 1, "Pls use --support_genome to show a list of supported genomes!"
      }
    } else {
      exit 1, "PREPROCESS_10XGENOMICS: --ref_cellranger_ucsc, or --ref_cellranger_ensembl, or --ref_fasta/ref_gtf must be specified!"
    }

    res_files = Channel.empty()

  emit:
    res_files // out[0]: res folders for MultiQC report
    CELLRANGER_ATAC_COUNT.out.fragments // out[1]: for split bed
    CELLRANGER_ATAC_COUNT.out.ch_fragment // out[2]: fragment ch for ArchR
    CELLRANGER_ATAC_COUNT.out.sample_name // out[3]: for split bam
    CELLRANGER_ATAC_COUNT.out.bam // out[4]: for split bam
}

workflow.onComplete {
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
