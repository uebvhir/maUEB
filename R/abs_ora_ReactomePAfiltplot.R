#' Reactome filt ORA
#'
#' Creates dataframes of Reactome results filtered by pvalue and performs a Reactome overrepresentation analysis, giving some plots with the most enriched GO terms pathways.

#' @param ORAresReac <- Object resulting from the function 'abs_ORA_ReactomePAunfilt'
#' @param label <- Label to use in some file names
#' @param resultsSummFN <- Name of the file where a summary of the results will be saved
#' @param sign <- vector of length namescomp composed of one number for comparison: +1 for up -1 for down or 0 for abs
#' @param ReacPlots <- Whether to generate plots or not
#' @param Reac_pvalcutoff <- Thresholds of p-value for Reactome overrepresentation analysis (in this order: BP, CC, MF)
#' @param Reac_pAdjustMethod <- Adjustment method for each Reactome category to analyse (CC, BP, MF) in each comparison
#' @details
#' @import writexl
#' @import DT
#' @import ggplot2
#' @import enrichplot
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates dataframes of Reactome results filtered by pvalue and performs a Reactome overrespresentation analysis, giving some plots with the most overrepresented Reactome pathways.
#' @keywords Reactome ORA plots
#' @references
#' @export

abs_ora_ReactomePAfiltplot <- function(ORAresReac, sign, Reac_pvalcutoff, Reac_pAdjustMethod,
                                       saveTables, ReacPlots, label, resultsSummFN, outputDir) {
    #extreiem les llistes
    listselected <- ORAresReac$listselected
    listReactomePAresults <- ORAresReac$listReactomePAresults
    for (i in 1:length(listselected)){
        comparison <- names(listselected)[i]
        #get reactome result from list
        reactome.result <- listReactomePAresults[[comparison]]
        Reactab <- as.data.frame(reactome.result)
        #filter results by pvalue threshold set in parameters
        if (Reac_pAdjustMethod[i]=="none"){
            Reactab.f <- Reactab[Reactab$pvalue < Reac_pvalcutoff[i],]
        } else {Reactab.f <- Reactab[Reactab$p.adjust < Reac_pvalcutoff[i],]}
        #save summary of filtered results
        sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
        cat("\n-------------------------ORA-REACTOME----------------------------------")
        cat("\nComparison: ", comparison,"\n")
        cat("Reactome pval-cutoff: ", Reac_pvalcutoff[i],"\n")
        cat("Reactome pAdjustMethod: ", Reac_pAdjustMethod[i],"\n")
        cat("Number of terms found: ", nrow(Reactab.f))
        cat("\n")
        sink()
        #save filtered results as xls and html
        if (nrow(Reactab.f) > 0) {
            if (saveTables){
                writexl::write_xlsx(Reactab.f, file.path(outputDir, paste0("ReactomePAResultsORA.", comparison, label, ".xlsx")))
                #guardem llistes per Cytoscape-EnrichmentMap:
                nomscols <- c("ID", "Description", "pvalue", "p.adjust")
                Reactab.f1 <- Reactab.f[,nomscols]
                Reactab.f1$Sign <- sign[i]
                write.table(Reactab.f1, file = file.path(outputDir, paste0("ReactomePAResultsORA.", comparison, label, ".txt")), sep = "\t", dec = ".", row.names = FALSE)

                #Create HTML table
                ##create links for reactome pathways
                Reactab.f2 <- Reactab.f
                Reactab.f2$ID <- paste0("<a href=\'https://reactome.org/PathwayBrowser/#/", Reactab.f2$ID, "\'>", Reactab.f2$ID, "</a>")
                ##reformat dataframe to show
                Reactab.f2 <- Reactab.f2[,-c(7,9)] #delete columns that are not needed to show (qvalues,Counts)
                Reactab.f2$geneID <- gsub("/",", ",Reactab.f2$geneID) #change bars by commas
                Reactab.f2 <- Reactab.f2[order(Reactab.f2$pvalue),] #order results by pvalue
                Reactab.f2[,c("pvalue", "p.adjust")] <- format(Reactab.f2[,c("pvalue", "p.adjust")], digits=3, scientific=TRUE) #change pvalues to scientific format
                y<-datatable(Reactab.f2, escape=F, caption=paste0("ReactomePA Results for comparison: ", comparison),
                             options = list(autoWidth = TRUE,
                                            columnDefs = list(list(width = '15%', targets=0), list(width = '30%', targets=1), list(width="50%", targets=6)),
                                            pageLength = 50,
                                            lengthMenu = c(10, 50, 200, 500)),
                             rownames=F, width="100%")%>%
                    formatStyle(columns = 7, fontSize = '70%')
                saveWidget(y, file=paste0(outputDir, paste0("/ReactomePAResultsORA.", comparison, label, ".html")))
            }
            #Plots
            if (ReacPlots){
                pval <- "pvalue" #note: selecting the adjusted pvalue is not advised since usually is the same for the top terms selected and doesn't allow to visualize in different colors
                pdf(file = file.path(outputDir, paste0("ReactomePAplotsORA.",comparison, label, ".pdf")), width = 14, height = 14)
                if (as.numeric(R.Version()$major)>=4) {
                ##Dotplot
                ifelse (nrow(Reactab.f) < 15, ncateg <- nrow(Reactab.f), ncateg <- 15)
                print(dotplot(reactome.result, showCategory=ncateg, font.size = 15,
                              title = paste0("Reactome Pathway Analysis for ", comparison,". Dotplot"), color=pval))
                ##Enrichment map
                ifelse (nrow(Reactab.f) < 60, ncateg <- nrow(Reactab.f), ncateg <- 60)
                pt.reactome.result <- pairwise_termsim(reactome.result, method="JC", showCategory=ncateg)
                print(emapplot(pt.reactome.result, color=pval, showCategory=ncateg, min_edge=0.2, group_category=FALSE, group_legend=TRUE,
                               node_label="category", label_style="shadowtext", shadowtext=FALSE))
                ##Network
                print(cnetplot(reactome.result, categorySize = pval, colorEdge=TRUE,
                               showCategory = 5, circular = FALSE))
                } else {
                    ##Dotplot
                    ifelse (nrow(Reactab.f) < 15, ncateg <- nrow(Reactab.f), ncateg <- 15)
                    print(dotplot(reactome.result, showCategory=ncateg, font.size = 15,
                                  title = paste0("Reactome Pathway Analysis for ", comparison,". Dotplot"), color=pval))
                    ##Enrichment map
                    ifelse (nrow(Reactab.f) < 60, ncateg <- nrow(Reactab.f), ncateg <- 60)
                    print(emapplot(reactome.result, color=pval, showCategory=ncateg, min_edge=0.2, group_category=FALSE, group_legend=TRUE,
                                   node_label="category", label_style="shadowtext", shadowtext=FALSE))
                    ##Network
                    print(cnetplot(reactome.result, categorySize = pval, colorEdge=TRUE,
                                   showCategory = 5, circular = FALSE))
                }
                dev.off()
            }
        }
    }
}
