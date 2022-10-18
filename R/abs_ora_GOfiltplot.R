#' GO filt ORA
#'
#' Creates dataframes of GO results filtered by pvalue and performs a GO significance analysis, giving some plots with the most enriched GO terms pathways.

#' @param ORAresGO <- Object resulting from the function 'abs_ORA_GOunfilt'
#' @param label <- Label to use in some file names
#' @param resultsSummFN <- Name of the file where a summary of the results will be saved
#' @param sign <- vector of length namescomp composed of one number for comparison: +1 for up -1 for down or 0 for abs
#' @param GOPlots <- Whether to generate plots or not
#' @param GO_pvalcutoff <- Thresholds of p-value for GO analysis (in this order: BP, CC, MF)
#' @param GO_pAdjustMethod <- Adjustment method for each GO category to analyse (CC, BP, MF) in each comparison
#' @details
#' @import writexl
#' @import DT
#' @import ggplot2
#' @import enrichplot
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates dataframes of GO results filtered by pvalue and performs a GO overrespresentation analysis, giving some plots with the most overrepresented GO terms pathways.
#' @keywords GO ORA plots
#' @references
#' @export

abs_ora_GOfiltplot <- function(ORAresGO, sign, GO_pvalcutoff, GO_pAdjustMethod,
                               saveTables, GOPlots, label, resultsSummFN, outputDir) {
    #extreiem les llistes
    listselected <- ORAresGO$listselected
    listGOresults <- ORAresGO$listGOresults
    for (i in 1:length(listselected)){
        comparison <- names(listselected)[i]
        for (categ in names(listGOresults)){
            go.result <- listGOresults[[categ]][[i]]
            GOtab <- as.data.frame(go.result)
            #filter results by pvalue threshold set in parameters
            if (GO_pAdjustMethod[[i]][categ]=="none"){
                GOtab.f <- GOtab[GOtab$pvalue < GO_pvalcutoff[[i]][categ],]
            } else {GOtab.f <- GOtab[GOtab$p.adjust < GO_pvalcutoff[[i]][categ],]}
            #save summary of filtered results
            sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
            cat("\n-------------------------ORA-GO----------------------------------")
            cat("\nComparison: ", comparison,"\n")
            cat("GO category: ",categ,"\n")
            cat("GO pval-cutoff: ",GO_pvalcutoff[[i]][categ],"\n")
            cat("GO pAdjustMethod: ",GO_pAdjustMethod[[i]][categ],"\n")
            cat("Number of terms found: ", nrow(GOtab.f))
            cat("\n")
            sink()
            #save filtered results as xls and html
            if (nrow(GOtab.f) > 0) {
                if (saveTables){
                    writexl::write_xlsx(GOtab.f, file.path(outputDir, paste0("GOResultsORA.", categ, ".", comparison, label, ".xlsx")))
                    #guardem llistes per Cytoscape-EnrichmentMap:
                    nomscols <- c("ID", "Description", "pvalue", "p.adjust")
                    GOtab.f1 <- GOtab.f[,nomscols]
                    GOtab.f1$Sign <- sign[i]
                    write.table(GOtab.f1, file = file.path(outputDir, paste0("GOResultsORA.", categ, ".", comparison, label, ".txt")), sep = "\t", dec = ".", row.names = FALSE)
                    #Create HTML table
                    GOtab.f2 <- GOtab.f
                    ##include GO terms as links
                    GOtab.f2$ID <- paste0("<a href=\'http://amigo.geneontology.org/amigo/term/", GOtab.f2$ID, "\'>", GOtab.f2$ID, "</a>")
                    ##reformat dataframe
                    GOtab.f2 <- GOtab.f2[,-c(7,9)] #delete columns that are not needed to show (qvalues,Count)
                    GOtab.f2$geneID <- gsub("/",", ",GOtab.f2$geneID) #change bars to commas in 'geneID' column
                    GOtab.f2 <- GOtab.f2[order(GOtab.f2$pvalue),] #order results by pvalue
                    GOtab.f2[,c("pvalue", "p.adjust")] <- format(GOtab.f2[,c("pvalue", "p.adjust")], digits=3, scientific=TRUE) #convert pvalues to scientific format
                    y<-datatable(GOtab.f2, escape=F, caption=paste0("GO-", categ, " Results for comparison: ", comparison),
                                 options = list(autoWidth = TRUE,
                                                columnDefs = list(list(width = '30%', targets=1), list(width="50%", targets=6)),
                                                pageLength = 50,
                                                lengthMenu = c(10, 50, 200, 500)),
                                 rownames=F, width="100%")%>%
                        formatStyle(columns = 7, fontSize = '70%')
                    saveWidget(y, file=paste0(outputDir, paste0("/GOResultsORA.", categ, ".", comparison, label, ".html")))
                }
                if (GOPlots){
                    pval <- "pvalue" #note: selecting the adjusted pvalue is not advised since usually is the same for the top terms selected and doesn't allow to visualize in different colors
                    pdf(file=file.path(outputDir, paste0("GOplotsORA.",categ, ".", comparison, label, ".pdf")), width = 14, height = 14)
                    if (as.numeric(R.Version()$major)>=4) {
                    ##Dotplot
                    ifelse (nrow(GOtab.f) < 15, ncateg <- nrow(GOtab.f), ncateg <- 15)
                    print(dotplot(go.result, showCategory = ncateg, font.size = 15,
                                  title = paste0("GO Analysis for ",categ, ":", comparison), color=pval))
                    ##Enrichment map
                    ifelse (nrow(GOtab.f) < 60, ncateg <- nrow(GOtab.f), ncateg <- 60)
                    pt.go.result <- pairwise_termsim(go.result, method="JC", showCategory=ncateg)
                    print(emapplot(pt.go.result, color=pval, showCategory=ncateg, min_edge=0.2, group_category=FALSE, group_legend=TRUE,
                                   node_label="category", label_style="shadowtext", shadowtext=FALSE))
                    ##Network
                    print(cnetplot(go.result, categorySize = pval, colorEdge=TRUE,
                                   showCategory = 5, circular = FALSE))
                    } else {
                        ifelse (nrow(GOtab.f) < 15, ncateg <- nrow(GOtab.f), ncateg <- 15)
                        print(dotplot(go.result, showCategory = ncateg, font.size = 15,
                                      title = paste0("GO Analysis for ",categ, ":", comparison), color=pval))
                        ##Enrichment map
                        ifelse (nrow(GOtab.f) < 60, ncateg <- nrow(GOtab.f), ncateg <- 60)
                        print(emapplot(go.result, color=pval, showCategory=ncateg))
                        ##Network
                        print(cnetplot(go.result, categorySize = pval, colorEdge=TRUE,
                                       showCategory = 5, circular = FALSE))
                    }
                    dev.off()
                }
            }
        }
    }
}
