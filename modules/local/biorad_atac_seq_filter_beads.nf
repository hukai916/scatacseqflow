// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Parse software version numbers
 */
process BIORAD_ATAC_SEQ_ALIGNMENT_QC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'biorad_bead_filtration', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/biorad_atac_seq_filter_beads:0.1"
    // cache false

    input:
    val sample_name
    path alignments
    val biorad_genome

    output:
    path "bead_filtration", emit: bead_filtration
    val sample_name, emit: sample_name

    script:

    """
    /filterBeads.sh -i $alignments -o bead_filtration -r $biorad_genome

    """
}
