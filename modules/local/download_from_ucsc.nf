// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process DOWNLOAD_FROM_UCSC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_from_ucsc', publish_id:'') }

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
    path "*.fa.gz", emit: genome_fasta
    path "md5sum.txt", emit: genome_md5
    val genome_name, emit: genome_name

    script:
    download_link = "https://hgdownload.soe.ucsc.edu/goldenPath/" + genome_name + "/bigZips/" + genome_name + ".fa.gz"
    md5_link = "https://hgdownload.soe.ucsc.edu/goldenPath/" + genome_name + "/bigZips/md5sum.txt"

    """
    wget $md5_link -o logfile.md5.txt
    wget $download_link -o logfile.genome.txt

    (cat \$(basename $md5_link) | grep \$( basename $download_link) || true) > md5_to_check.txt

    if [ -s md5_to_check.txt ]
    then
      md5sum -c md5_to_check.txt
    fi

    """
}
