// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process CORRECT_BARCODE_PHENIQS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode_pheniqs', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/pheniqs_xenial:0.1"

    // cache false

    input:
    val sample_name
    path barcode_fastq
    path barcode_whitelist
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name
    path "barcode_corrected_*", emit: corrected_barcode
    path "summary_*.txt", emit: corrected_barcode_summary
    path "barcode_corrected*.R1.fastq.gz", emit: read1_fastq
    path "barcode_corrected*.R2.fastq.gz", emit: read2_fastq

    script:

    """
    # step1, interleave read and index files
    pheniqs mux -i $read1_fastq -i $barcode_fastq -i $read2_fastq --output ${sample_name}.cram

    # step2, prepare a json config file for pheniqs
    ## determine the index read length from index fastq file:
    barcode_length=\$(pheniqs mux -i $barcode_fastq | head -n 1000 | cut -f 10 | awk '{print length(\$0)}' | uniq -c | sort -r | head -n 1 | cut -f 3 -d " ")
    make_json.py $barcode_whitelist ${sample_name}.cram 3 0::,2:: 1::\$barcode_length ${sample_name}.json

    # step4, run pheniqs
    pheniqs mux --threads $task.cpus --decoding-threads $task.cpus --htslib-threads $task.cpus --config debarcode.json --output ${sample_name}.bam

    # step5, extract fastq from pheniqs output bam
    bam2fastq.py ${sample_name}.bam barcode_corrected_${sample_name}

    """
}
