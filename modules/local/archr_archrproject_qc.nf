// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_ARCHRPROJECT_QC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_archrproject_qc', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc:0.4"

    // cache false

    input:
    path archr_project
    // val archr_genome
    // val archr_thread

    output:
    path "Plots/p1.pdf", emit: pdf_p1
    // path quality_control, emit: quality_control // if using this syntax, the -resume won't work
    // path "QualityControl", emit: quality_control // using this, the -resume won't work either.
    // This is because the quality_control folder content gets updated after each run, and it will be used as input for itself, so each time, it rerun, the timestamp of this folder is newer.

    script:
    // for unknown reason, #!/usr/bin/R + direct R codes won't work
    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    # Create 
    p1 <- plotGroups(
      ArchRProj = proj,
      groupBy = "Sample",
      colorBy = "cellColData",
      name = "TSSEnrichment",
      plotAs = "ridges"
    )

    plotPDF(p1, name = "p1.pdf", ArchRProj = NULL, addDOC = FALSE, width = 4, height = 4)

    ' > run.R

    Rscript run.R

    """
}
