// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Parse software version numbers
 */
process DOWNLOAD_FROM_ENSEMBL {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_from_ensembl', publish_id:'') }

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
    val ensembl_release

    output:
    path "*.fa.gz", emit: genome_fasta
    path "CHECKSUMS", emit: genome_md5

    script:
    def dict_genome_name = [homo_sapiens: "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", mus_musculus: "Mus_musculus.GRCm39.dna.primary_assembly.fa.gz", danio_rerio: "Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"]

    download_link = "http://ftp.ensembl.org/pub/release-" + ensembl_release + "/fasta/" + genome_name + "/dna/" + dict_genome_name[genome_name]

    md5_link = "http://ftp.ensembl.org/pub/release-" + ensembl_release + "/fasta/" + genome_name + "/dna/CHECKSUMS"

    """
    wget $md5_link -o logfile.md5.txt
    wget $download_link -o logfile.genome.txt

    (cat \$(basename $md5_link) | grep \$( basename $download_link) || true) > md5_to_check.txt

    if [ -s md5_to_check.txt ]
    then
      real=\$(cat md5_to_check.txt | cut -f 1,2 -d " ")
      measure=\$(sum \$( basename $download_link))

      if [ "\$real" == "\$measure" ]
      then
        echo "\$real,\$measure"
        exit 0
      else
        echo "\$real,\$measure"
        exit 1
      fi
    fi

    """
}
