// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Parse software version numbers
 */
process GET_BIORAD_FASTQ {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'fastq_folder_biorad', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "ubuntu:xenial"

    // cache false

    input:
    tuple val(sample_name), val(path_fastq_1), val(path_fastq_2), val(path_barcode)

    output:
    path "fastq_*", emit: fastq_folder
    // path "fastq_*/*S1_L000_R1_001.fastq.gz", emit: read1_fastq
    // path "fastq_*/*S1_L000_R2_001.fastq.gz", emit: barcode_fastq
    // path "fastq_*/*S1_L000_R3_001.fastq.gz", emit: read2_fastq
    val sample_name, emit: sample_name

    script:
    """
    mkdir fastq_$sample_name

    fastq1=(`echo "$path_fastq_1" | tr ';' ' '`)
    fastq2=(`echo "$path_fastq_2" | tr ';' ' '`)

    for i in "\${fastq1[@]}"
    do
      ln -s "\$i" fastq_$sample_name/
    done

    for i in "\${fastq2[@]}"
    do
      ln -s "\$i" fastq_$sample_name/
    done

    """
}