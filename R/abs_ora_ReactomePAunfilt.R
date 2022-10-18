#' Reactome unfilt ORA
#'
#' Performs Reactome overrepresentation analysis for all the comparisons and store unfiltered results (no pvalue threshold) in lists.

#' @param annotPackage <- Annotation package (Clariom D/S, specific of species; e.g. "clariomdhumantranscriptcluster.db", "clariomshumantranscriptcluster.db", "mta10transcriptcluster.db", "clariomsmousetranscriptcluster.db")
#' @param organism_annot <- Organism annotation package (e.g. "org.Hs.eg.db", "org.Mm.eg.db")
#' @param organism <- Organism ("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly")
#' @param label <- Label to use in some file names
#' @param resultsSummFN <- Name of the file where a summary of the results will be saved
#' @param universe <- Universe
#' @param keyType <- ID used in the analysis for gene classification (e.g. "UNIPROT", "ENTREZID")
#' @param col_entrez <- Name of the column used for the overrepresentation analysis (e.g. "EntrezID")
#' @param Reac_minSetSize <- Minimum set size for the analysis
#' @param Reac_maxSetSize <- Maximum set size for the analysis
#' @param selected.pval.reac <- "P.Value" or "adj.P.Val" (for each set of comparisons)
#' @param selected.pvalthr.reac <- Pvalue threshold  (for each set of comparisons)
#' @param selected.logFC.reac <- LogFC threshold
#' @param col_logFC <- Name of the column where the logFC is located
#' @param sign_logFC <- Label for deciding what kind of analysis perform: "up", "down" or "abs"
#' @details
#' @import reactome.db AnnotationDbi
#' @import clusterProfiler
#' @import writexl
#' @import ReactomePA
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates lists of Reactome unfiltered results and a table with number of terms found at different pvalue thresholds for the different categories analyzed.
#' @keywords Reactome ORA
#' @references
#' @export

abs_ora_ReactomePAunfilt <- function(listselgenesNamed=NULL, listofTopNamed, namescomp, universe, col_entrez, readable, selected.pval.reac, selected.pvalthr.reac, selected.logFC.reac, col_logFC.reac, sign_logFC.reac, organism, organism_annot, keyType, Reac_minSetSize, Reac_maxSetSize, resultsSummFN, saveRda, saveGMT, label, outputDir, rdaDir, gmtDir){
    #Prepare lists for storing results
    if (!exists("listReactomePAresults")) {listReactomePAresults <- list()}
    if (!exists("listselected")) {listselected <- list()}
    #Prepare numTermsChanged table
    numTerms.df <- data.frame(pval_type=character(), pval_thr=numeric(), NumTerms=integer(), stringsAsFactors=FALSE)
    numTerms.df[1:7,"pval_type"] <- c(rep("pvalue", 3), rep("p.adjust", 4))
    numTerms.df[,"pval_thr"] <- c(c(0.01,0.05,0.1), c(0.01,0.05,0.15,0.25))
    #Prepare selected genes for enrichment analysis
    ##If input is a list of Toptables
    if (length(listofTopNamed)>0){
        for (i in 1:length(namescomp)){
            comparison = namescomp[i]
            #prepare data: select genes in top tab according to cutoffs (a millorar: q es pugui triar si logFC en absolut)
            topTab <- listofTopNamed[[comparison]]
            if (sign_logFC.reac=="abs") topTab.sel <- topTab[topTab[,selected.pval.reac[i]] < selected.pvalthr.reac[i] & abs(topTab[,col_logFC.reac]) > selected.logFC.reac[i],]
            if (sign_logFC.reac=="up") topTab.sel <- topTab[topTab[,selected.pval.reac[i]] < selected.pvalthr.reac[i] & topTab[,col_logFC.reac] > selected.logFC.reac[i],]
            if (sign_logFC.reac=="down") topTab.sel <- topTab[topTab[,selected.pval.reac[i]] < selected.pvalthr.reac[i] & topTab[,col_logFC.reac] < selected.logFC.reac[i],]
            selected.genes <- as.character(topTab.sel[,col_entrez])
            listselected[[comparison]] <- selected.genes
            #para guardar en el Rsession
            sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
            cat("\n-------------------------ORA-REACTOME----------------------------------")
            cat("\nComparison: ", comparison)
            cat("\nUniverse: ", universe)
            cat("\nSelected genes-pval type: ", selected.pval.reac[i])
            cat("\nSelected genes-pval threshold: ", selected.pvalthr.reac[i])
            cat("\nSelected genes-logFC threshold: ", selected.logFC.reac[i], "(", sign_logFC.reac, ")")
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
    #Reactome analysis for the comparisons selected
    for (comparison in namescomp){
        selected.genes <- listselected[[comparison]]
        selected.genes <- unique(selected.genes)
        numTerms.dfi <- numTerms.df
        #analysis
        reactome.result <- enrichPathway(gene = selected.genes,
                                         organism = organism,
                                         if (!is.null(universe)) eval(parse(text="universe=universe")),
                                         pAdjustMethod = "BH", # fix adjustment method so that it computes an adjusted pvalue, later on results will be filtered according to thresholds set in parameters
                                         pvalueCutoff  = 1, # fix to 1, later on results will be filtered according to thresholds set in parameters
                                         qvalueCutoff = 1, #important! sino per defecte agafa qvalueCutoff de 0.05 i no retorna gairebe res!!
                                         minGSSize = Reac_minSetSize,
                                         maxGSSize = Reac_maxSetSize,
                                         readable=readable) #convert from entrez to GeneSymbol
        #store reactome result in list
        listReactomePAresults[[comparison]] <- reactome.result
        #save reactome result as dataframe (not filtered by pvalue)
        Reactab <- as.data.frame(reactome.result)
        #save unfiltered results
        writexl::write_xlsx(Reactab, file.path(outputDir, paste0("ReactomePAResultsORA.", comparison, label, ".unfilt.xlsx")))
        #fill numTermsChanged table with the number of terms found under different pvalue thresholds for GO significance
        for (j in 1:nrow(numTerms.dfi)){
            numTerms.dfi[j, "NumTerms"] <- nrow(Reactab[(Reactab[, numTerms.dfi[j,"pval_type"]] < numTerms.dfi[j,"pval_thr"]),])
        }
        #Save numTermsChanged table as xls
        writexl::write_xlsx(numTerms.dfi, file.path(outputDir, paste0("Reactome_numTermsChangedORA.",comparison, label, ".xlsx")))
    }
    #Save lists in Rdata object
    if (saveRda) save(listReactomePAresults, listselected, file = file.path(rdaDir,paste0('ABS-ORA-ReactomePAresults-unfilt', label, '.Rda')), version=2) #return list of lists
    #save gmt files
    if (saveGMT) {
        ##tots els genesets independentment del tamany
        reac.result <- listReactomePAresults[[1]]
        reac.result.df <- as.data.frame(reac.result, stringsAsFactors=FALSE)
        reac.geneset <- reac.result@geneSets
        reac.ids <- names(reac.geneset)
        ##no hi ha la funcio go2term per reactome, utilitzo el paquet reactome.db per obtenir la description
        reac.descrip <- AnnotationDbi::select(reactome.db,
                                              keys=reac.ids,
                                              keytype="PATHID",
                                              column=c("PATHID","PATHNAME"),
                                              multivals="first")
        reac.descrip$PATHNAME <- gsub(".*: ","", reac.descrip$PATHNAME)
        rownames(reac.descrip) <- as.character(reac.descrip$PATHID)
        listreac.geneset <- list(genesets=reac.geneset, geneset.names=reac.ids, geneset.descriptions=reac.descrip[reac.ids,"PATHNAME"]) #important aquests noms
        write.gmt(listreac.geneset, file.path(gmtDir, paste0("ReacORA.gmt")))
        #cal fer aixo, sino no ho reconeix be a enrichmentmap:
        gmt.f <- read.csv(file.path(gmtDir, paste0("ReacORA.gmt")), header=FALSE, sep="\t")
        write.table(gmt.f, file.path(gmtDir, paste0("ReacORA.gmt")), sep="\t", row.names=FALSE, col.names = FALSE)
    }
    ORA_ReactomePAresults.unfilt <- list(listReactomePAresults=listReactomePAresults, listselected=listselected)
    return(ORA_ReactomePAresults.unfilt)
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
