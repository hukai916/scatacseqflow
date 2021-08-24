// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_FOOTPRINTING_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_footprinting_clusters', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc:0.5"

    // cache false

    input:
    path archr_project
    path archr_dir

    output:
    path "Plots/Footprints-*-Bias.pdf", emit: footprints
    path "save_archr_project/Plots/TSS-No-Normalization.pdf", emit: tss_no_normalization

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    # Footprinting of motif:
    motifPositions <- getPositions(proj)
    motifs <- c($options.motifs)
    markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

    seFoot <- getFootprints(
      ArchRProj = proj,
      positions = motifPositions[markerMotifs],
      groupBy = "Clusters"
      )
    plotName <- paste0("Footprints", "-", "$options.norm_method", "-Bias")
    plotFootprints(
      seFoot = seFoot,
      ArchRProj = proj,
      normMethod = "$options.norm_method",
      plotName = plotName,
      $options.args
      )

    # Footprinting of TSS (custom) Features
    seTSS <- getFootprints(
      ArchRProj = proj,
      positions = GRangesList(TSS = getTSS(proj)),
      groupBy = "Clusters",
      flank = $options.tss_flank
      )
    plotFootprints(
      seFoot = seTSS,
      ArchRProj = NULL,
      normMethod = "None",
      plotName = "TSS-No-Normalization",
      addDOC = FALSE,
      flank = $options.tss_flank,
      flankNorm = $options.flank_norm
      )

    ' > run.R

    Rscript run.R

    """
}
