#' Multiple comparisons analysis: Venn diagrams, Upset plots and Heatmaps
#'
#' Plot Venn diagramm from 2 to 5 comparisons or UpsetR plot if >5 comparisons. It also plots heatmaps between the comparisons determined, with the thresholds values (p-value and logFC) of interest.
#' The Venn and Upset plots function is based in previous work of Ferran Briansó and Miriam Mota ("https://raw.githubusercontent.com/rgonzalo/GitBackupScripts/master/VennDiagrams/CreateVennEuler.R").
#'
#' @inheritParams mc_venn_upset
#' @inheritParams mc_hm
#' @details
#' @import VennDiagram
#' @import writexl
#' @import limma
#' @import pheatmap
#' @import heatmaply
#' @import plotly
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com} and Esther Camacho
#' @examples
#' @return Plot Venn diagramm from 2 to 5 comparisons or UpsetR plot if >5 comparisons. Plot heatmap plots.
#' @keywords Venn diagrams Upset plots Heatmap
#' @references
#' @export

mc_all <- function(listofcsv, venn_comparNames, venn_compar, label, colFeat, colPval, pval, colFC, FC, pltR=TRUE,
                   pltPdf=TRUE, pltPng=FALSE, venn=TRUE, eul=FALSE, saveTables=TRUE, upsetPlot=FALSE, saveTables_extended=TRUE,
                   colors=rainbow(length(namescomp)), trans=0.5, cex1=1, rotation=0, position=rep(0,length(namescomp)), cex2=1, include, hm_dataset, targets, featureCol="Gene.Symbol", hm_comparNames, hm_compar, hm_groupexclude=rep(list(c()), length(hm_comparNames)), hm_pval, hm_pval.thr, hm_logFC="logFC", hm_logFC.thr, hm_palette=colorRampPalette(c("blue", "white", "red"))(n = 199), hm_clustCol=TRUE, hm_clustCol.dist="euclidian", hm_clustCol.method="complete", hm_clustRow.cor=TRUE, batcheffect=FALSE,
                   batchcolName=NULL, batchcolName2=NULL, batchNumcolName=NULL, hm_plots_interactive=FALSE, Venn_plots=TRUE, Heatmap_plots=TRUE,resultsSummFN, outputDir){
    if (Venn_plots){listsharedelems <- lapply(seq_along(venn_comparNames), function(v) {
        workingDir <- getwd()
        setwd(outputDir)
        namescomp <- names(listofcsv)[venn_compar[[v]]]
        mc_venn_upset(listofcsv=listofcsv, namescomp=namescomp, label = venn_comparNames[v], colFeat = colFeat[v],
                                  colPval = venn_pval[v], pval = venn_pval.thr[v], colFC=venn_FC_col[v], FC=venn_FC.thr[v], include=venn_include[v], pltR = FALSE,
                                  pltPdf = TRUE, pltPng=FALSE, venn = TRUE, eul = FALSE, saveTables = TRUE, upsetPlot=FALSE,
                                  saveTables_extended=Venn_sharedElems_extended, colors = rainbow(length(namescomp)), trans = 0.5,
                                  cex1 = 1, rotation=0, position=venn_pos[[v]], cex2 = 1, resultsSummFN=resultsSummFN, outputDir=resultsDir)
        setwd(workingDir)
    })}
    if (Heatmap_plots){mc_hm(hm_comparNames=hm_comparNames, hm_compar=hm_compar, hm_groupexclude=hm_groupexclude, hm_pval=hm_pval, hm_pval.thr=hm_pval.thr, hm_logFC=hm_logFC, hm_palette=hm_palette, hm_clustCol=hm_clustCol, hm_clustCol.dist=hm_clustCol.dist, hm_clustCol.method=hm_clustCol.method, hm_clustRow.cor=hm_clustRow.cor, batcheffect=batcheffect, batcheffect2=batcheffect2, batcheffect3=batcheffect3,batchcolName=batchcolName, batchcolName2=batchcolName2,batchcolName3=batchcolName3, resultsSummFN=resultsSummFN, outputDir=resultsDir)}
}
