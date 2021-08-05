// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_CREATE_ARROWFILE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_create_arrowfile', publish_id:'') }

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
    val sample_name
    // val archr_genome
    // val archr_thread
    path fragment

    output:
    val sample_name, emit: sample_name
    // obj_arrowfile, emit: obj_arrowfile

    script:

    """
    #!/usr/bin/R

    library(ArchR)
    addArchRGenome("hg19")
    addArchRThreads(threads = 4)

    inputFiles <- $fragment
    names(inputFiles) <- $sample_name

    ArrowFiles <- createArrowFiles(
      inputFiles = inputFiles,
      sampleNames = names(inputFiles),
      minTSS = 4, #Dont set this too high because you can always increase later
      minFrags = 1000,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE
    )

    """
}
