#' GO unfilt ORA
#'
#' Performs GO overrepresentation analysis for all the comparisons and store unfiltered results (no pvalue threshold) in lists.

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
#' @details
#' @import reactome.db AnnotationDbi
#' @import clusterProfiler
#' @import writexl
#' @import utils
#' @import limma
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates lists of GO unfiltered results and a table with number of terms found at different pvalue thresholds for the different categories analyzed.
#' @keywords GO ORA
#' @references
#' @export

abs_ora_GOunfilt <- function(listselgenesNamed=NULL, listofTopNamed, namescomp, universe, geneColname, readable, selected.pval.go, selected.pvalthr.go, selected.logFC.go, col_logFC.go, sign_logFC.go, organism_annot, keyType, GO_categories, GO_minSetSize, GO_maxSetSize, resultsSummFN, saveRda, saveGMT, label, outputDir, rdaDir, gmtDir){
    #Prepare lists for storing results
    if (!exists("listGOresults")) {listGOresults <- rep(list(list()), length(GO_categories)); names(listGOresults) <- GO_categories}
    if (!exists("listselected")) {listselected <- list()}
    #Prepare numTermsChanged table
    numTerms.df <- data.frame(pval_type=character(), pval_thr=numeric(), BP.NumTerms=integer(), CC.NumTerms=integer(), MF.NumTerms=integer(), stringsAsFactors=FALSE)
    numTerms.df[1:7,"pval_type"] <- c(rep("pvalue", 3), rep("p.adjust", 4))
    numTerms.df[,"pval_thr"] <- c(c(0.01,0.05,0.1), c(0.01,0.05,0.15,0.25))
    #Prepare selected genes for enrichment analysis
    ##If input is a list of Toptables
    if (length(listofTopNamed)>0){
        for (i in 1:length(namescomp)){
            comparison <- namescomp[i]
            #prepare data: select genes in top tab according to cutoffs (a millorar: q es pugui triar si logFC en absolut)
            topTab <- listofTopNamed[[comparison]]
            if (sign_logFC.go=="abs") topTab.sel <- topTab[topTab[,selected.pval.go[i]] < selected.pvalthr.go[i] & abs(topTab[,col_logFC.go]) > selected.logFC.go[i],]
            if (sign_logFC.go=="up") topTab.sel <- topTab[topTab[,selected.pval.go[i]] < selected.pvalthr.go[i] & topTab[,col_logFC.go] > selected.logFC.go[i],]
            if (sign_logFC.go=="down") topTab.sel <- topTab[topTab[,selected.pval.go[i]] < selected.pvalthr.go[i] & topTab[,col_logFC.go] < selected.logFC.go[i],]
            selected.genes <- as.character(topTab.sel[,geneColname])
            listselected[[comparison]] <- selected.genes
            #para guardar en el Rsession
            sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
            cat("\n-------------------------ORA-GO----------------------------------")
            cat("\nComparison: ", comparison)
            cat("\nUniverse: ", universe)
            cat("\nSelected genes-pval type: ", selected.pval.go[i])
            cat("\nSelected genes-pval threshold: ", selected.pvalthr.go[i])
            cat("\nSelected genes-logFC threshold: ", selected.logFC.go[i], "(", sign_logFC.go, ")")
            cat("\nNum. selected genes: ", length(selected.genes))
            cat("\nNum. unique selected genes: ", length(unique(selected.genes)))
            cat("\n")
            sink()
        }
    }
    ##If input is a list of genes
    if (length(listselgenesNamed)>0){
        listselected <- listselgenesNamed
    }
    #Enrichment GO ontology
    for (comparison in namescomp){
        selected.genes <- listselected[[comparison]]
        selected.genes <- unique(selected.genes)
        numTerms.dfi <- numTerms.df
        #Loop for the different GO categories
        for (categ in GO_categories){
            go.result <- enrichGO(gene = selected.genes,
                                  if (!is.null(universe)) eval(parse(text="universe=universe")),
                                  OrgDb = organism_annot,
                                  keyType = keyType,
                                  ont = categ,
                                  pAdjustMethod = "BH", # fix adjustment method so that it computes an adjusted pvalue, later on results will be filtered according to thresholds set in parameters
                                  pvalueCutoff  = 1, # fix to 1, later on results will be filtered according to thresholds set in parameters
                                  qvalueCutoff = 1, #important! sino per defecte agafa qvalueCutoff de 0.05 i no retorna gairebe res!!
                                  minGSSize = GO_minSetSize, # min set size
                                  maxGSSize = GO_maxSetSize, # max set size
                                  readable=readable) #whether to convert from entrez to GeneSymbol
            #store go result in list
            listGOresults[[categ]][[comparison]] <- go.result
            #save go result as dataframe (not filtered by pvalue)
            GOtab <- as.data.frame(go.result)
            #save unfiltered results
            writexl::write_xlsx(GOtab, file.path(outputDir, paste0("GOResultsORA.", categ, ".", comparison, label, ".unfilt.xlsx")))
            #fill numTermsChanged table with the number of terms found under different pvalue thresholds for GO significance
            for (j in 1:nrow(numTerms.dfi)){
                numTerms.dfi[j, paste0(categ,".NumTerms")] <- nrow(GOtab[(GOtab[, numTerms.dfi[j,"pval_type"]] < numTerms.dfi[j,"pval_thr"]),])
            }
        }
        #Save numTermsChanged table as xls
        writexl::write_xlsx(numTerms.dfi, file.path(outputDir, paste0("GO_numTermsChangedORA.", comparison, label, ".xlsx")))
    }
    #Save lists in Rda object
    if (saveRda) save(listGOresults, listselected, namescomp, file = file.path(rdaDir, paste0('ABS-ORA-GOresults-unfilt', label, '.Rda')), version=2)
    #save gmt files
    if (saveGMT) {
        lapply(GO_categories, function(categ) {
            ##tots els genesets independentment del tamany
            go.result <- listGOresults[[categ]][[1]]
            go.geneset <- go.result@geneSets
            go.ids <- names(go.geneset)
            go.descrip <- go2term(go.ids) #ojo q canvia l'ordre! important reordenar abans de fer la llista
            rownames(go.descrip) <- go.descrip$go_id
            go.descrip <- go.descrip[go.ids,]
            listgo.geneset <- list(genesets=go.geneset, geneset.names=go.ids, geneset.descriptions=go.descrip[go.ids,"Term"]) #important aquests noms
            write.gmt(listgo.geneset, file.path(gmtDir, paste0("GO-ORA",categ,".gmt")))
            #cal fer aixo, sino no ho reconeix be a enrichmentmap:
            gmt.f <- read.csv(file.path(gmtDir, paste0("GO-ORA",categ,".gmt")), header=FALSE, sep="\t")
            write.table(gmt.f, file.path(gmtDir, paste0("GO-ORA",categ,".gmt")), sep="\t", row.names=FALSE, col.names=FALSE)
        })
    }
    #return list of lists
    ORA_GOresults.unfilt <- list(listGOresults=listGOresults, listselected=listselected)
    return(ORA_GOresults.unfilt)
}

#function write.gmt from mw201608/msigdb
write.gmt=function(obj,filename){
    conn=file(filename,'w')
    for(i in 1:length(obj$genesets)){
        cat(obj$geneset.names[i],obj$geneset.descriptions[i],obj$genesets[[i]],file=conn,sep='\t')
        cat('\n',file=conn)
    }
    close(conn)
    return(invisible())
}
