// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_GENE_SCORE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_gene_score', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc:0.5"

    // cache false

    input:
    path archr_project

    output:
    path "proj_gene_score.rds", emit: archr_project

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- addUMAP(
      ArchRProj = proj,
      reducedDims = "IterativeLSI",
      name = "UMAP",
      $options.args
    )

    proj2 <- addTSNE(
      ArchRProj = proj2,
      reducedDims = "IterativeLSI",
      name = "TSNE",
      $options.args2
    )

    proj2 <- addUMAP(
      ArchRProj = proj,
      reducedDims = "Harmony",
      name = "UMAPHarmony",
      $options.args
    )

    ' > run.R

    Rscript run.R

    """
}
