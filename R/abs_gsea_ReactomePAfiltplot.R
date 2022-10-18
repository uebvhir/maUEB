#' Reactome filt
#'
#' Performs Creates dataframes of Reactome results filtered by pvalue are obtained and performs a Reactome significance analysis, giving some plots with the most enriched Reactome pathways.

#' @param annotPackage Annotation package (Clariom D/S, specific of species; e.g. "clariomdhumantranscriptcluster.db", "clariomshumantranscriptcluster.db", "mta10transcriptcluster.db", "clariomsmousetranscriptcluster.db")
#' @param organism_annot  <- Organism annotation package (e.g. "org.Hs.eg.db", "org.Mm.eg.db")
#' @param Reac_pvalcutoff <- Thresholds of p-value for GO analysis (in this order: BP, CC, MF)
#' @param Reac_pAdjustMethod <- Adjustment method for each GO category to analyse (CC, BP, MF) in each comparison
#' @param GSEAresReac <-Object resulting from the function 'abs_gsea_ReactomePAunfilt'
#' @details
#' @import DT
#' @import ggplot2
#' @import enrichplot
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates dataframes of Reactome results filtered by pvalue and performs a Reactome significance analysis, giving some plots with the most enriched Reactome pathways.
#' @keywords Reactome GSEA plots
#' @references
#' @export


abs_gsea_ReactomePAfiltplot <- function(GSEAresReac, Reac_pvalcutoff, Reac_pAdjustMethod, saveTables, ReacPlots, label, resultsSummFN, outputDir) {
    #extreiem les llistes
    listReactomePAgeneFC <- GSEAresReac$listReactomePAgeneFC
    listReactomePAresults <- GSEAresReac$listReactomePAresults
    #Filter Reactome results according to thresholds and save dataframes as html and xls
    for (i in 1:length(listReactomePAgeneFC)){
        comparison <- names(listReactomePAgeneFC)[i]
        #get reactome result from list
        gsea.result <- listReactomePAresults[[i]]
        Reactab <- as.data.frame(gsea.result)
        #afegim columnes "NumGenes.changed" i GeneRatio
        genes <- apply(Reactab, 1, function (x) {length(strsplit(x["core_enrichment"], "/")[[1]])})
        Reactab$NumGenes.core <- genes
        Reactab$GeneRatio <- round(Reactab$NumGenes.core/Reactab$setSize,3)
        #reordenem la taula
        Reactab <- Reactab[,c("ID", "Description", "setSize", "NumGenes.core", "GeneRatio", colnames(Reactab[,4:11]))]

        #filter results by pvalue threshold set in parameters
        if (Reac_pAdjustMethod[i]=="none"){
            Reactab.f <- Reactab[Reactab$pvalue < Reac_pvalcutoff[i],]
        } else {Reactab.f <- Reactab[Reactab$p.adjust < Reac_pvalcutoff[i],]}
        #save summary of filtered results
        sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
        cat("\n-------------------------GSEA-REACTOME----------------------------------")
        cat("\nComparison: ", comparison,"\n")
        cat("Reactome pval-cutoff: ", Reac_pvalcutoff[i],"\n")
        cat("Reactome pAdjustMethod: ", Reac_pAdjustMethod[i],"\n")
        cat("Number of terms found: ", nrow(Reactab.f))
        cat("\n")
        sink()
        #save filtered results as xls and html
        if (nrow(Reactab.f) != 0) {
            if (saveTables){
                writexl::write_xlsx(Reactab.f, file.path(outputDir, paste0("ReactomePAResultsGSEA.", comparison, label, ".xlsx")))
                #guardem llistes per Cytoscape-EnrichmentMap:
                nomscols <- c("ID", "Description", "pvalue", "p.adjust", "NES", "enrichmentScore", "setSize", "core_enrichment")
                Reactab.f1 <- Reactab.f[,nomscols]
                Reactab.f1$core_enrichment <- gsub("/", ", ", Reactab.f1$core_enrichment)
                write.table(Reactab.f1, file = file.path(outputDir, paste0("ReactomePAResultsGSEA.", comparison, label, ".txt")), sep = "\t", dec = ".", row.names = FALSE)

                #Create HTML table
                ##create links for reactome pathways
                Reactab.f2 <- Reactab.f
                Reactab.f2$ID <- paste0("<a href=\'https://reactome.org/PathwayBrowser/#/", Reactab.f2$ID, "\'>", Reactab.f2$ID, "</a>")
                ##reformat dataframe to show
                Reactab.f2 <- Reactab.f2[,-c(6,10:12)] #delete columns that are not needed to show (enrichmentScore,qvalues,rank,leading_edge)
                Reactab.f2$NES <- round(Reactab.f2$NES, 2)
                Reactab.f2$core_enrichment <- gsub("/",", ",Reactab.f2$core_enrichment) #change bars by commas
                Reactab.f2 <- Reactab.f2[order(Reactab.f2$pvalue),] #order results by pvalue
                Reactab.f2[,7:8] <- format(Reactab.f2[,7:8], digits=3, scientific=TRUE) #change pvalues to scientific format
                y<-datatable(Reactab.f2, escape=F, caption=paste0("ReactomePA Results for comparison: ", comparison),
                             options = list(autoWidth = TRUE,
                                            columnDefs = list(list(width = '15%', targets=0), list(width = '30%', targets=1), list(width="50%", targets=8)),
                                            pageLength = 50,
                                            lengthMenu = c(10, 50, 200, 500)),
                             rownames=F, width="100%")%>%
                    formatStyle(columns = 9, fontSize = '70%')
                saveWidget(y, file=paste0(outputDir, paste0("/ReactomePAResultsGSEA.", comparison, label, ".html")))
            }
            #Plots
            if (ReacPlots){
                options(ggrepel.max.overlaps = 1000)
                pval <- "pvalue" #note: selecting the adjusted pvalue is not advised since usually is the same for the top terms selected and doesn't allow to visualize in different colors
                geneFC <- listReactomePAgeneFC[[i]]
                pdf(file = file.path(outputDir, paste0("ReactomePAplotsGSEA.",comparison, label, ".pdf")), width = 14, height = 14)
                if (as.numeric(R.Version()$major)>=4) {
                ##Dotplot
                ifelse (nrow(Reactab.f) < 15, ncateg <- nrow(Reactab.f), ncateg <- 15)
                print(dotplot(gsea.result, showCategory=ncateg, font.size = 15, size="GeneRatio",
                              title = paste0("Reactome Pathway Analysis for ", comparison,". Dotplot"), color=pval,
                              label_format=60, split=".sign") + facet_grid(.~.sign))
                ##Enrichment map
                ifelse (nrow(Reactab.f) < 60, ncateg <- nrow(Reactab.f), ncateg <- 60)
                pt.gsea.result <- pairwise_termsim(gsea.result, method="JC", showCategory=ncateg)
                print(emapplot(pt.gsea.result, color=pval, showCategory=ncateg, min_edge=0.2, group_category=FALSE, group_legend=TRUE,
                               node_label="category", label_style="shadowtext", shadowtext=FALSE))

                ##Network
                print(cnetplot(gsea.result, categorySize = pval, foldChange = geneFC, colorEdge=TRUE,
                               showCategory = 5, circular = FALSE, cex_label_category=1.3))
                } else {
                    ##Dotplot
                    ifelse (nrow(Reactab.f) < 15, ncateg <- nrow(Reactab.f), ncateg <- 15)
                    print(dotplot(gsea.result, showCategory=ncateg, font.size = 15, size="GeneRatio",
                                  title = paste0("Reactome Pathway Analysis for ", comparison,". Dotplot"), color=pval,
                                  split=".sign") + facet_grid(.~.sign))
                    ##Enrichment map
                    ifelse (nrow(Reactab.f) < 60, ncateg <- nrow(Reactab.f), ncateg <- 60)
                    print(emapplot(gsea.result, color=pval, showCategory=ncateg))

                    ##Network
                    print(cnetplot(gsea.result, categorySize = pval, foldChange = geneFC, colorEdge=TRUE,
                                   showCategory = 5, circular = FALSE))
                }
                dev.off()
            }
        }
    }
}
