// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Parse software version numbers
 */
process CORRECT_BARCODE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc_atac:0.1"

    // cache false

    input:
    path barcode_fastq
    path barcode_whitelist

    output:
    path "barcode_*", emit: corrected_barcode
    path "summary.txt", emit: corrected_barcode_summary

    script:

    """
    correct_barcode.R \
    --barcode_file=$barcode_fastq \
    --whitelist_file=$barcode_whitelist \
    --path_output_fq=./

    """
}
