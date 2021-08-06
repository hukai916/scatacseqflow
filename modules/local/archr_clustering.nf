// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_CLUSTERING {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_clustering', publish_id:'') }

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
    path "proj_clustering.rds", emit: archr_project
    path "Cluster-seurat-matrix.csv", emit: csv_cluster_seurat_matrix
    path "Cluster-scran-matrix.csv", emit: csv_cluster_scran_matrix
    path "Plots/Cluster-heatmap.pdf", emit: pdf_cluster_heatmap

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- addClusters(
      input = proj,
      reducedDims = "IterativeLSI",
      method = "Seurat",
      name = "Clusters",
      $options.args
    )
    proj2 <- addClusters(
      input = proj2,
      reducedDims = "IterativeLSI",
      method = "scran",
      name = "ScranClusters",
      $options.args2
    )

    saveRDS(proj2, file = "proj_clustering.rds")

    # Save text summary and heatmap summary
    cM <- confusionMatrix(paste0(proj2\$Clusters), paste0(proj2\$Sample))
    cM_scran <- confusionMatrix(paste0(proj2\$ScranClusters), paste0(proj2\$Sample))

    write.csv(cM, file="Cluster-seurat-matrix.csv")
    write.csv(cM_scran, file="Cluster-scran-matrix.csv")

    library(pheatmap)
    cM <- cM / Matrix::rowSums(cM)
    cM_scran <- cM_scran / Matrix::rowSums(cM)
    p1 <- pheatmap::pheatmap(
      mat = as.matrix(cM),
      color = paletteContinuous("whiteBlue"),
      border_color = "black"
    )
    p2 <- pheatmap::pheatmap(
      mat = as.matrix(cM_scran),
      color = paletteContinuous("whiteBlue"),
      border_color = "black"
    )
    plotPDF(p1, p2, name = "Cluster-heatmap.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)



    ' > run.R

    Rscript run.R

    """
}
