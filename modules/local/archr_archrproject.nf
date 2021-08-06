// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_ARCHRPROJECT {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_archrproject', publish_id:'') }

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
    val archr_genome
    val archr_thread
    path arrowfile
    // path quality_control


    output:
    // val sample_name, emit: sample_name
    // path sample_name, emit: output_dir
    path "proj.rds", emit: archr_project

    script:
    // for unknown reason, #!/usr/bin/R + direct R codes won't work
    """
    echo '
    library(ArchR)

    addArchRGenome("$archr_genome")
    addArchRThreads(threads = $archr_thread)

    proj <- ArchRProject(
    ArrowFiles = "$arrowfile",
    outputDirectory = "$sample_name",
    $options.args)

    saveRDS(proj, file = "proj.rds")
    ' > run.R

    Rscript run.R

    """
}
