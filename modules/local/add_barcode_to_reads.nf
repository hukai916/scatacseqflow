// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Parse software version numbers
 */
process ADD_BARCODE_TO_READS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'add_barcode_to_reads', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/add_barcode_to_reads:0.1"

    // cache false

    input:
    val sample_name
    path barcode_fastq
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name

    path "R1/*barcoded*", emit: read1_fastq
    path "R2/*barcoded*", emit: read2_fastq

    script:

    """
    echo "I am here!"
    mkdir -p R1/barcoded
    mkdir -p R2/barcoded

    # use the first read length from fastq file to determine the length since -b is required by sinto.
    filename=\$(basename -- "$barcode_fastq")
    extension="\${filename##*.}"

    if [[ "\$extension" == "gz" ]]
    then
      barcode_length=\$(zcat < $barcode_fastq | awk '{if(NR%4==2) print length(\$1)}' | head -n 1)
      echo "Using full length of the first record in barcode read fastq.gz file as -b to sinto."
    else
      # barcode_length=\$(cat < $barcode_fastq | awk '{if(NR%4==2)
      # print length(\$1)}' | head -n 1)
      echo "Using full length of the first record in barcode read fastq file as -b to sinto."
    fi
    """
}
