// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process DOWNLOAD_FROM_GENCODE_GTF {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_from_gencode_gtf', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/miniconda3_xenial:0.1"

    // cache false

    input:
    val genome_name

    output:
    path "*.gtf.gz", emit: gtf

    script:
    download_link = "https://hgdownload.soe.ucsc.edu/goldenPath/" + genome_name + "/bigZips/genes/" + genome_name + ".ncbiRefSeq.gtf.gz"

    """
    wget $download_link -o logfile.gtf.txt

    """
}
