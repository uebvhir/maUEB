#' GSEA analysis
#'
#' Performs GO and Reactome analysis for all the comparisons and store unfiltered results (no pvalue threshold) in lists. Creates dataframes of GO and Reactome results filtered by pvalue are obtained and performs a GO significance analysis, giving some plots with the most enriched GO terms or Reactome pathways.
#' @param annotPackage Annotation package (Clariom D/S, specific of species; e.g. "clariomdhumantranscriptcluster.db", "clariomshumantranscriptcluster.db", "mta10transcriptcluster.db", "clariomsmousetranscriptcluster.db")
#' @param organism_annot  <- Organism annotation package (e.g. "org.Hs.eg.db", "org.Mm.eg.db")
#' @param GO_comparisons <- Comparisons to be analysed, indicated by the corresponding position in 'listofcoef' or 'listofTopNames' (eg. 1:10 or c(1,2,8))
#' @param GO_pvalcutoff <- Thresholds of p-value for GO analysis (in this order: BP, CC, MF)
#' @param GO_pAdjustMethod <- Adjustment method for each GO category to analyse (CC, BP, MF) in each comparison
#' @param GO_minSetSize <- Minimum size of the pathways analysed in the analysis
#' @param GO_maxSetSize <- Maximum size of the pathways analysed in the analysis
#' @param GSEAresGO <-Object resulting from the function 'abs_gsea_GOunfilt'
#' @param Reac_comparisons <- Comparisons to be analysed, indicated by the corresponding position in 'listofcoef' or 'listofTopNames' (eg. 1:10 or c(1,2,8))
#' @param Reac_pvalcutoff <- Thresholds of p-value for GO analysis (in this order: BP, CC, MF)
#' @param Reac_pAdjustMethod <- Adjustment method for each GO category to analyse (CC, BP, MF) in each comparison
#' @param Reac_minSetSize <- Minimum size of the pathways analysed in the analysis
#' @param Reac_maxSetSize <- Maximum size of the pathways analysed in the analysis
#' @param GSEAresReac <-Object resulting from the function 'abs_gsea_ReactomePAunfilt'
#' @details
#' @import reactome.db AnnotationDbi
#' @import clusterProfiler
#' @import limma
#' @import DT
#' @import ggplot2
#' @import enrichplot
#' @import ReactomePA
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates lists of GO and Reactome unfiltered results and a table with number of terms found at different pvalue thresholds for the different categories analyzed.
#' @keywords GO Reactome GSEA ORA plots
#' @references
#' @export

abs_gsea_all <- function(listofTopNamed, namescomp, ranking_metric, geneColname,organism_annot, keyType, readable=TRUE, GO_categories, GO_minSetSize, GO_maxSetSize, resultsSummFN, saveRda=TRUE, saveGMT=TRUE, label, outputDir, rdaDir, gmtDir, GO_pvalcutoff, GO_pAdjustMethod, saveTables, GOPlots, GSEAresGO, GSEAresReac, Reac_pvalcutoff, Reac_pAdjustMethod, ReacPlots, col_entrez, Reac_minSetSize, Reac_maxSetSize){
if(GOAnalysis.unfilt){abs_gsea_GOunfilt(listofTopNamed=listofTopNamed, namescomp=namescomp, organism_annot= organism_annot, keyType="ENTREZID", GO_minSetSize=GO_minSetSize, GO_maxSetSize=GO_maxSetSize, resultsSummFN=resultsSummFN, saveRda=saveRda, saveGMT=saveGMT, label=label, outputDir=resultsDir, rdaDir=dataDir, gmtDir=resultsDir)}
if(GOAnalysis.filt){abs_gsea_GOfiltplot(GSEAresGO=GSEA_GOresults.unfilt, GO_pvalcutoff=GO_pvalcutoff, GO_pAdjustMethod=GO_pAdjustMethod, saveTables=saveTables, GOPlots=GOPlots, label=label, resultsSummFN=resultsSummFN, outputDir=resultsDir)}
if(ReactomeAnalysis.unfilt){abs_gsea_ReactomePAunfilt(listofTopNamed=listofTopNamed, namescomp=namescomp, organism=organism, organism_annot=organism_annot, Reac_minSetSize=Reac_minSetSize, Reac_maxSetSize=Reac_maxSetSize, resultsSummFN=resultsSummFN, saveRda=saveRda, saveGMT=saveGMT, label=label, outputDir=resultsDir, rdaDir=dataDir, gmtDir=resultsDir)}
if(ReactomeAnalysis.filt){abs_gsea_ReactomePAfiltplot(GSEAresReac=GSEA_ReactomePA.unfilt, Reac_pvalcutoff=Reac_pvalcutoff, Reac_pAdjustMethod=Reac_pAdjustMethod, saveTables=saveTables, ReacPlots=ReactomePlots, label=label, resultsSummFN=resultsSummFN, outputDir=resultsDir)}}
