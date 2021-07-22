// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Parse software version numbers
 */
process BWA_MAP {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'bwa_map', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/bwa_xenial:0.1"

    // cache false

    input:
    path read1_fastq
    path read2_fastq
    path bwa_index_folder

    output:
    path "*.bam", emit: bam

    script:
    genome_basename = genome_fasta.getName()

    """
    mkdir bwa_index
    ln $genome_fasta bwa_index/
    bwa index -a bwtsw bwa_index/$genome_basename

    """
}
