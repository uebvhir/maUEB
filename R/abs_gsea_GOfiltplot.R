#' GO results Filtering and Plots
#'
#' Filter GO results according to pvalue threshold and saves a pdf with different plots to visualize the results of enrichment analysis, using clusterProfiler package (dotplot, enrichment map and netplot)

#' @param GSEAresGO Object resulting from the function 'abs_gsea_GOunfilt'
#' @param GO_pvalcutoff List of p-value thresholds for filtering GO analysis results (in this order: BP, CC, MF) in each comparison
#' @param GO_pAdjustMethod List of pvalue adjustment method for each GO category to analyse (CC, BP, MF) in each comparison (eg. "BH", "bonferroni" or "none")
#' @param saveTables
#' @param GOPlots
#' @param label
#' @param resultsSummFN
#' @param outputDir
#' @details
#' @import DT ggplot2 enrichplot
#' @importFrom writexl write_xlsx
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Saves dataframes of GO results filtered by pvalue and plots to visualize the results.
#' @keywords GO GSEA plots
#' @references
#' @export

abs_gsea_GOfiltplot <- function(GSEAresGO, GO_pvalcutoff=rep(list(c(BP=0.05, CC=0.05, MF=0.05)),length(GSEAresGO)), GO_pAdjustMethod=rep(list(c(BP="BH", CC="BH", MF="BH")),length(GSEAresGO)), saveTables=TRUE, GOPlots=TRUE, label="", resultsSummFN="ResultsSummary_ABS-GSEA.txt", outputDir) {
    #extreiem les llistes
    listGOgeneFC <- GSEAresGO$listGOgeneFC
    listGOresults <- GSEAresGO$listGOresults
    for (i in 1:length(listGOgeneFC)){
        comparison <- names(listGOgeneFC)[i]
        for (categ in names(listGOresults)){
            gseago.result <- listGOresults[[categ]][[i]]
            GOtab <- as.data.frame(gseago.result)
            genes <- apply(GOtab, 1, function (x) {length(strsplit(x["core_enrichment"], "/")[[1]])})
            GOtab$NumGenes.core <- genes
            GOtab$GeneRatio <- round(GOtab$NumGenes.core/GOtab$setSize,3)
            #reordenem la taula
            GOtab <- GOtab[,c("ID", "Description", "setSize", "NumGenes.core", "GeneRatio", colnames(GOtab[,4:11]))]
            #filter results by pvalue threshold set in parameters
            if (GO_pAdjustMethod[[i]][categ]=="none"){
                GOtab.f <- GOtab[GOtab$pvalue < GO_pvalcutoff[[i]][categ],]
            } else {GOtab.f <- GOtab[GOtab$p.adjust < GO_pvalcutoff[[i]][categ],]}
            #save summary of filtered results
            sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
            cat("\n-------------------------GSEA-GO----------------------------------")
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
                    writexl::write_xlsx(GOtab.f, file.path(outputDir, paste0("GOResultsGSEA.", categ, ".", comparison, label, ".xlsx")))
                    #guardem llistes per Cytoscape-EnrichmentMap:
                    nomscols <- c("ID", "Description", "pvalue", "p.adjust", "NES", "enrichmentScore", "setSize", "core_enrichment")
                    GOtab.f1 <- GOtab.f[,nomscols]
                    GOtab.f1$core_enrichment <- gsub("/", ", ", GOtab.f1$core_enrichment)
                    write.table(GOtab.f1, file = file.path(outputDir, paste0("GOResultsGSEA.", categ, ".", comparison, label, ".txt")), sep = "\t", dec = ".", row.names = FALSE)

                    #Create HTML table
                    ##include GO terms as links
                    GOtab.f2 <- GOtab.f
                    GOtab.f2$ID <- paste0("<a href=\'http://amigo.geneontology.org/amigo/term/", GOtab.f2$ID, "\'>", GOtab.f2$ID, "</a>")
                    ##reformat dataframe
                    GOtab.f2 <- GOtab.f2[,-c(6,10:12)] #delete columns that are not needed to show (enrichmentScore,qvalues,rank,leading_edge)
                    GOtab.f2$core_enrichment <- gsub("/",", ",GOtab.f2$core_enrichment) #change bars to commas in 'core_enrichment' column
                    GOtab.f2$NES <- round(GOtab.f2$NES, 2)
                    GOtab.f2 <- GOtab.f2[order(GOtab.f2$pvalue),] #order results by pvalue
                    GOtab.f2[,7:8] <- format(GOtab.f2[,7:8], digits=3, scientific=TRUE) #convert pvalues to scientific format
                    y<-datatable(GOtab.f2, escape=F, caption=paste0("GO-", categ, " Results for comparison: ", comparison),
                                 options = list(autoWidth = TRUE,
                                                columnDefs = list(list(width = '15%', targets=0), list(width = '30%', targets=1), list(width="50%", targets=8)),
                                                pageLength = 50,
                                                lengthMenu = c(10, 50, 200, 500)),
                                 rownames=F, width="100%")%>%
                        formatStyle(columns = 9, fontSize = '70%')
                    saveWidget(y, file=paste0(outputDir, paste0("/GOResultsGSEA.", categ, ".", comparison, label, ".html")))
                }
                #GO plots
                if (GOPlots){
                    options(ggrepel.max.overlaps = 1000)
                    pval <- "pvalue" #note: selecting the adjusted pvalue is not advised since usually is the same for the top terms selected and doesn't allow to visualize in different colors
                    geneFC <- listGOgeneFC[[i]]
                    pdf(file=file.path(outputDir, paste0("GOplotsGSEA.",categ, ".", comparison, label, ".pdf")), width = 14, height = 14)
                    if (as.numeric(R.Version()$major)>=4) {
                    ##Dotplot
                    ifelse (nrow(GOtab.f) < 15, ncateg <- nrow(GOtab.f), ncateg <- 15)
                    dp <- dotplot(gseago.result, size="GeneRatio", showCategory = ncateg, font.size = 16,label_format=40,
                                  title = paste0("GO Analysis for ",categ, ":", comparison), color=pval,
                                  label_format=60, split=".sign") + facet_grid(.~.sign) +
                        theme(strip.text.x = element_text(size = 16), text=element_text(size=16))
                   print(dp)
                     ##Enrichment map
                    ifelse (nrow(GOtab.f) < 60, ncateg <- nrow(GOtab.f), ncateg <- 60)
                    pt.gseago.result <- pairwise_termsim(gseago.result, method="JC", showCategory=ncateg)
                    ep <- emapplot(pt.gseago.result, color=pval, showCategory=ncateg, min_edge=0.2, group_category=FALSE, group_legend=TRUE,
                                   node_label="category", label_style="shadowtext", shadowtext=FALSE)+ theme(legend.text=element_text(size=15), legend.title=element_text(size=15))
                    print(ep)
                    ##Network
                    cnp <- netplot(gseago.result, categorySize = pval, foldChange = geneFC, colorEdge=TRUE,
                                   showCategory = 5, circular = FALSE, cex_label_category=1.3)+
                        theme(legend.text=element_text(size=15), legend.title=element_text(size=15))
                    print(cnp)
                    } else {
                        ##Dotplot
                        ifelse (nrow(GOtab.f) < 15, ncateg <- nrow(GOtab.f), ncateg <- 15)
                        dp <- dotplot(gseago.result, size="GeneRatio", showCategory = ncateg, font.size = 16,
                                      title = paste0("GO Analysis for ",categ, ":", comparison), color=pval, split=".sign") + facet_grid(.~.sign)+
                            theme(strip.text.x = element_text(size = 16), text=element_text(size=16))
                        print(dp)
                        ##Enrichment map
                        ifelse (nrow(GOtab.f) < 60, ncateg <- nrow(GOtab.f), ncateg <- 60)
                        ep <-emapplot(gseago.result, color=pval, showCategory=ncateg)+ theme(legend.text=element_text(size=15), legend.title=element_text(size=15))
                        print(ep)
                        ##Network
                        cnp <- cnetplot(gseago.result, categorySize = pval, foldChange = geneFC, colorEdge=TRUE,
                                       showCategory = 5, circular = FALSE)+
                            theme(legend.text=element_text(size=15), legend.title=element_text(size=15))
                        print(cnp)
                    }
                    dev.off()
                }
            }
        }
    }
}

