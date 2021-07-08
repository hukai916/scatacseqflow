// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Parse software version numbers
 */
process CELLRANGER_ATAC_COUNT {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'cellranger_count', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/cellranger_atat_2.0.0:0.1"

    // cache false

    input:
    path fastq_folder
    path reference
    // val jobmode

    output:
    path "cellranger_atac_count_*", emit: cellranger_atac_count

    script:

    """
    cellranger-atac count \
    --id cellranger_atac_count_$fastq_folder \
    --fastqs $fastq_folder \
    --reference $reference \

    """
}
