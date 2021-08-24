// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_marker_peaks_in_tracks_clusters', publish_id:'') }

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
    path marker_peaks
    val gene_symbol
    val cluster_name

    output:
    path "Plots/Plot-Tracks-With-Features.pdf", emit: archr_tracks_with_features

    script:

    """
    echo '
    library(ArchR)
    markersPeaks <- readRDS("$marker_peaks")
    proj <- readRDS("$archr_project")

    p <- plotBrowserTrack(
      ArchRProj = proj,
      groupBy = "Clusters",
      geneSymbol = c("$gene_symbol"),
      features =  getMarkers(markersPeaks, cutOff = "$options.cutoff", returnGR = TRUE)["$cluster_name"],
      $options.args
    )

    plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    ' > run.R

    Rscript run.R

    """
}