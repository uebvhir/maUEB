#' ORA analysis
#'
#' Performs GO and Reactome overrepresentation analysis for all the comparisons and store unfiltered results (no pvalue threshold) in lists. Creates dataframes of GO and Reactome results filtered by pvalue and performs a GO and Reactome significance analysis, giving some plots with the most enriched GO terms and Reactome pathways.

#' @param annotPackage <- Annotation package (Clariom D/S, specific of species; e.g. "clariomdhumantranscriptcluster.db", "clariomshumantranscriptcluster.db", "mta10transcriptcluster.db", "clariomsmousetranscriptcluster.db")
#' @param organism_annot <- Organism annotation package (e.g. "org.Hs.eg.db", "org.Mm.eg.db")
#' @param organism <- Organism ("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly")
#' @param label <- Label to use in some file names
#' @param resultsSummFN <- Name of the file where a summary of the results will be saved
#' @param universe <- Universe
#' @param keyType <- ID used in the analysis for gene classification (e.g. "UNIPROT", "ENTREZID")
#' @param geneColname <- Name of the column used for the overrepresentation analysis (e.g. "EntrezID")
#' @param GO_categories <- List of categories that will be used in the analysis (e.g. "BP", "CC" and/or "MF")
#' @param GO_minSetSize <- Minimum set size for the analysis
#' @param GO_maxSetSize <- Maximum set size for the analysis
#' @param selected.pval.go <- "P.Value" or "adj.P.Val" (for each set of comparisons)
#' @param selected.pvalthr.go <- Pvalue threshold  (for each set of comparisons)
#' @param selected.logFC.go <- LogFC threshold
#' @param col_logFC.go <- Name of the column where the logFC is located
#' @param sign_logFC.go <- Label for deciding what kind of analysis perform: "up", "down" or "abs"
#' @param ORAresGO <- Object resulting from the function 'abs_ORA_GOunfilt'
#' @param sign <- vector of length namescomp composed of one number for comparison: +1 for up -1 for down or 0 for abs
#' @param GOPlots <- Whether to generate plots or not
#' @param GO_pvalcutoff <- Thresholds of p-value for GO analysis (in this order: BP, CC, MF)
#' @param GO_pAdjustMethod <- Adjustment method for each GO category to analyse (CC, BP, MF) in each comparison
#' @param ORAresReac <- Object resulting from the function 'abs_ORA_ReactomePAunfilt'
#' @param ReacPlots <- Whether to generate plots or not
#' @param Reac_pvalcutoff <- Thresholds of p-value for Reactome overrepresentation analysis (in this order: BP, CC, MF)
#' @param Reac_pAdjustMethod <- Adjustment method for each Reactome category to analyse (CC, BP, MF) in each comparison
#' @details
#' @import reactome.db AnnotationDbi
#' @import clusterProfiler
#' @import writexl
#' @import utils
#' @import limma
#' @import DT
#' @import ggplot2
#' @import enrichplot
#' @import ReactomePA
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates lists of GO unfiltered results and a table with number of terms found at different pvalue thresholds for the different categories analyzed.
#' @keywords GO Reactome ORA plots
#' @references
#' @export

abs_ora_all <- function(listofTopNamed, namescomp, universe, geneColname, readable=TRUE, selected.pval.go, selected.pvalthr.go, selected.logFC.go, col_logFC.go, sign_logFC.go, organism_annot, keyType, GO_categories, GO_minSetSize,
                        GO_maxSetSize, resultsSummFN, saveRda, saveGMT, label, outputDir, rdaDir, gmtDir, col_entrez, selected.pval.reac, selected.pvalthr.reac, selected.logFC.reac, col_logFC.reac,
                        sign_logFC.reac, organism, Reac_minSetSize, Reac_maxSetSize, ORAresGO, sign, GO_pvalcutoff, GO_pAdjustMethod, saveTables, GOPlots, ORAresReac, Reac_pvalcutoff, Reac_pAdjustMethod, ReacPlots){
            if(GOAnalysis.unfilt){ORA_GOresults.unfilt <- abs_ora_GOunfilt(listofTopNamed=listofTopNamed, namescomp=namescomp, universe=universe, geneColname=geneColname, selected.pval.go=selected.pval.go, selected.pvalthr.go=selected.pvalthr.go, selected.logFC.go=selected.logFC.go, col_logFC.go=col_logFC.go, sign_logFC.go=sign_logFC.go, organism_annot=organism_annot, keyType=keyType, GO_categories=GO_categories, GO_minSetSize=GO_minSetSize, GO_maxSetSize=GO_maxSetSize, resultsSummFN=resultsSummFN, saveRda=saveRda, saveGMT=saveGMT, label=label, outputDir=resultsDir, rdaDir=rdaDir, gmtDir=resultsDir, readable=readable)}
            if(GOAnalysis.filt){abs_ora_GOfiltplot(ORAresGO=ORA_GOresults.unfilt, sign=sign, GO_pvalcutoff=GO_pvalcutoff, GO_pAdjustMethod=GO_pAdjustMethod, saveTables=saveTables, GOPlots=GOPlots, label=label, resultsSummFN=resultsSummFN, outputDir=resultsDir)}
            if(ReactomeAnalysis.unfilt){ORA_ReactomePAresults.unfilt <- abs_ora_ReactomePAunfilt(listofTopNamed=listofTopNamed, namescomp=namescomp, universe=universe, organism=organism, organism_annot=organism_annot, col_entrez=col_entrez, readable=readable, selected.pval.reac=selected.pval.reac, selected.pvalthr.reac=selected.pvalthr.reac, selected.logFC.reac=selected.logFC.reac, col_logFC.reac=col_logFC.reac, sign_logFC.reac=sign_logFC.reac, keyType=keyType, Reac_minSetSize=Reac_minSetSize, Reac_maxSetSize=Reac_maxSetSize, resultsSummFN=resultsSummFN, saveRda=saveRda, saveGMT=saveGMT, label=label, outputDir=resultsDir, rdaDir=rdaDir, gmtDir=resultsDir)}
            if(ReactomeAnalysis.filt){abs_ora_ReactomePAfiltplot(ORAresReac=ORA_ReactomePAresults.unfilt, sign=sign, Reac_pvalcutoff=Reac_pvalcutoff, Reac_pAdjustMethod=Reac_pAdjustMethod, saveTables=saveTables, ReacPlots=ReacPlots, label=label, resultsSummFN=resultsSummFN, outputDir=resultsDir)}}
