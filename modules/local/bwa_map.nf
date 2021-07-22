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
    val sample_name
    path read1_fastq
    path read2_fastq
    path bwa_index_folder

    output:
    val sample_name, emit: sample_name
    path "*.sorted.bam", emit: bam

    script:

    """
    index_name=\$(basename $bwa_index_folder/*.bwt)
    index_name="\${filename%.*}"

    bwa mem -t $task.cpus $bwa_index_folder/\$index_name $read1_fastq $read2_fastq | samtools sort -@ $task.cpus -O bam -o ${sample_name}.sorted.bam

    """
}