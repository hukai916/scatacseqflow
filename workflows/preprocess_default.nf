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
workflow PREPROCESS_DEFAULT {
  take:
    ch_samplesheet

  main:
    // log.info "INFO(2): --preprocess: default"
    GET_10XGENOMICS_FASTQ (ch_samplesheet)
    // module: fastQC
    FASTQC (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq,
    GET_10XGENOMICS_FASTQ.out.barcode_fastq)

    // module: barcode correction (optional) and add barcode: correct barcode fastq given whitelist and barcode fastq file
    if (!(params.barcode_whitelist)) {
      log.info "NOTICE: --barcode_whitelist: not supplied, skip barcode correction!"
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

    // module: mapping with bwa or minimap2: mark duplicate
    // bwa or minimap2
    if (params.mapper == 'bwa') {
      log.info "INFO: --mapper: bwa"
      if (params.ref_bwa_index) {
        BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, params.ref_bwa_index)
      } else if (params.ref_fasta) {
        log.info "INFO: --ref_fasta provided, use it for building bwa index."
        // module : prep_genome
        PREP_GENOME (params.ref_fasta, "custom_genome")
        // module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
        // module : bwa_map
        BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, BWA_INDEX.out.bwa_index_folder)
      } else if (params.ref_fasta_ensembl) {
        log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // module : download_from_ensembl
        DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // module: prep_genome
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
      } else if (params.ref_fasta_ucsc) {
        log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build bwa index, and map with bwa ..."
        // module : download_from_ucsc
        DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // module : prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_name)
        // module : extract primary sequence
        // GET_PRIMARY_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta)
        // module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
        // module : bwa_map
        BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, BWA_INDEX.out.bwa_index_folder)
      } else {
        exit 1, 'Parameter --ref_fasta_ensembl/--ref_fasta_ucsc: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
      }
    } else if (params.mapper == "minimap2") {
      log.info "INFO: --mapper: minimap2"
      if (params.ref_minimap2_index) {
        // use user provided bwa index for mapping
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, params.ref_minimap2_index)
      } else if (params.ref_fasta) {
        log.info "INFO: --ref_fasta provided, use it to build minimap2 index."
        // module : prep_genome
        PREP_GENOME (params.ref_fasta, "custom_genome")
        // module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, MINIMAP2_INDEX.out.minimap2_index)
      } else if (params.ref_fasta_ensembl) {
        log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // module : download_from_ensembl
        DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // module: PREP_GENOME
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, MINIMAP2_INDEX.out.minimap2_index)
      } else if (params.ref_fasta_ucsc) {
        log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // module : download_from_ucsc
        DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // module : prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_gtf)
        // module : get_primary_genome
        // GET_PRIMARY_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta)
        // module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, MINIMAP2_INDEX.out.minimap2_index)
      } else {
        exit 1, 'Parameter --ref_fasta_ucsc/--ref_fasta_ensembl: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
      }
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
    // REMOVE_DUPLICATE module:
    try {
      res_files = res_files.mix(REMOVE_DUPLICATE.out.remove_duplicate_summary.collect().ifEmpty([]))
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

  emit:
    res_files // out[0]: res folders for MultiQC report
    GET_FRAGMENTS.out.fragments // out[1]: for split bed
    GET_FRAGMENTS.out.ch_fragment // out[2]: fragment ch for ArchR
    REMOVE_DUPLICATE.out.sample_name // out[3]: for split bam
    REMOVE_DUPLICATE.out.bam // out[4]: for split bam
}

workflow.onComplete {
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
