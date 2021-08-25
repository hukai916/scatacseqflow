// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_positive_tf_regulator_clusters', publish_id:'') }

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
    path "archr_coaccessibility_project.rds", emit: archr_project
    path "Plots/Plot-Tracks-With-Features.pdf", emit: plot_tracks_with_features

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Clusters")
    seZ <- seGroupMotif[rowData(seGroupMotif)\$seqnames=="z",]
    rowData(seZ)\$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
        rowMaxs(assay(seZ) - assay(seZ)[,x])
      }) %>% Reduce("cbind", .) %>% rowMaxs

    corGSM_MM <- correlateMatrices(
      ArchRProj = proj,
      useMatrix1 = "GeneScoreMatrix",
      useMatrix2 = "MotifMatrix",
      reducedDims = "IterativeLSI"
      )

    corGSM_MM\$maxDelta <- rowData(seZ)[match(corGSM_MM\$MotifMatrix_name, rowData(seZ)\$name), "maxDelta"]

    corGSM_MM <- corGSM_MM[order(abs(corGSM_MM\$cor), decreasing = TRUE), ]
    corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
    corGSM_MM\$TFRegulator <- "NO"
    corGSM_MM\$TFRegulator[which(corGSM_MM\$cor > 0.5 & corGSM_MM\$padj < 0.01 & corGSM_MM\$maxDelta > quantile(corGSM_MM\$maxDelta, 0.75))] <- "YES"
    sort(corGSM_MM[corGSM_MM\$TFRegulator=="YES",1])

    p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
      geom_point() +
      theme_ArchR() +
      geom_vline(xintercept = 0, lty = "dashed") +
      scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
      xlab("Correlation To Gene Score") +
      ylab("Max TF Motif Delta") +
      scale_y_continuous(
        expand = c(0,0),
        limits = c(0, max(corGSM_MM\$maxDelta)*1.05)
      )
    plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    ' > run.R

    Rscript run.R

    """
}
