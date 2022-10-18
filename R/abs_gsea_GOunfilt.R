#' GO unfilt
#'
#' Performs GO analysis for all the comparisons and store unfiltered results (no pvalue threshold) in lists.

#' @param listofTopNamed Named list of toptables
#' @param namescomp Name of the comparisons to be analysed, corresponding to names in list of toptables
#' @param ranking_metric Column from toptable by which features will be ranked
#' @param geneColname Column from toptable containing feature ID
#' @param keyType
#' @param readable
#' @param annotPackage Annotation package (Clariom D/S, specific of species; e.g. "clariomdhumantranscriptcluster.db", "clariomshumantranscriptcluster.db", "mta10transcriptcluster.db", "clariomsmousetranscriptcluster.db")
#' @param organism_annot  Organism annotation package (e.g. "org.Hs.eg.db", "org.Mm.eg.db")
#' @param GO_categories
#' @param GO_minSetSize Minimum size of the pathways analysed in the analysis
#' @param GO_maxSetSize Maximum size of the pathways analysed in the analysis
#' @param resultsSummFN
#' @param saveRda
#' @param saveGMT
#' @param label
#' @param outputDir
#' @param rdaDir
#' @param gmtDir
#' @details
#' @import clusterProfiler
#' @importFrom writexl write_xlsx
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates lists of GO unfiltered results and a table with number of terms found at different pvalue thresholds for the different categories analyzed.
#' @keywords GO GSEA
#' @references
#' @export

abs_gsea_GOunfilt <- function(listofTopNamed, namescomp, ranking_metric="logFC", geneColname="EntrezID", annotPackage, organism_annot, keyType="ENTREZID", readable=TRUE, GO_categories=c("BP","CC","MF"), GO_minSetSize=3, GO_maxSetSize=500, resultsSummFN="ResultsSummary_ABS-GSEA.txt", saveRda=TRUE, saveGMT=TRUE, label="", outputDir, rdaDir, gmtDir){
  #Prepare lists for storing results
  #Prepare lists for storing results
    require(annotPackage, character.only=TRUE) #required for GO/Reactome
    require(organism_annot, character.only=TRUE) #required for GO/Reactome
  listGOresults <- rep(list(list()), length(GO_categories)); names(listGOresults) <- GO_categories
  listGOgeneFC <- list() #required for plots
  #Prepare numTermsChanged table
  numTerms.df <- data.frame(pval_type=character(), pval_thr=numeric(), BP.NumTerms=integer(), CC.NumTerms=integer(), MF.NumTerms=integer(), stringsAsFactors=FALSE)
  numTerms.df[1:7,"pval_type"] <- c(rep("pvalue", 3), rep("p.adjust", 4))
  numTerms.df[,"pval_thr"] <- c(c(0.01,0.05,0.1), c(0.01,0.05,0.15,0.25))
  #GO analysis for the comparisons selected
  for (comparison in namescomp){
    numTerms.dfi <- numTerms.df
    #prepare data: reorder top tab by logFC and take the entrez
    topTab <- listofTopNamed[[comparison]]
    geneFC <- topTab[,ranking_metric]
    names(geneFC) <- as.character(topTab[,geneColname])
    #eliminem els duplicats, abans d'ordenar per logFC (es quedara amb els mes significatius ja que estava ordenada per pval)
    geneFC <- geneFC[!duplicated(names(geneFC))]
    geneFC = sort(geneFC, decreasing = TRUE)
    listGOgeneFC[[comparison]] <- geneFC #save ranked gene names in list
    #Loop for the different GO categories
    for (categ in GO_categories){
      set.seed(123) # since 'gseGO' performs a permutation test, setting the seed allows to obtain same results in different executions

        gseago.result <- gseGO(geneList = geneFC,
                             OrgDb = eval(parse(text=organism_annot)),
                             ont = categ,
                             keyType = keyType,
                             if (as.numeric(R.Version()$major)>=4) eps=0.0 else nPerm=1000,
                             exponent=1,
                             seed=FALSE,
                             pAdjustMethod = "BH", # fix adjustment method so that it computes an adjusted pvalue, later on results will be filtered according to thresholds set in parameters
                             pvalueCutoff  = 1, # fix to 1, later on results will be filtered according to thresholds set in parameters
                             minGSSize = GO_minSetSize, # min set size
                             maxGSSize = GO_maxSetSize, # max set size
                             verbose = FALSE)

      if (readable) gseago.result <- setReadable(gseago.result, organism_annot)
      #store gseago result in list
      listGOresults[[categ]][[comparison]] <- gseago.result
      #save gseago result as dataframe (not filtered by pvalue)
      GOtab <- as.data.frame(gseago.result)
      #save unfiltered results
      writexl::write_xlsx(GOtab, file.path(outputDir, paste0("GOResultsGSEA.", categ, ".", comparison, label, ".unfilt.xlsx")))
      #fill numTermsChanged table with the number of terms found under different pvalue thresholds for GO significance
      for (j in 1:nrow(numTerms.dfi)){
        numTerms.dfi[j, paste0(categ,".NumTerms")] <- nrow(GOtab[(GOtab[, numTerms.dfi[j,"pval_type"]] < numTerms.dfi[j,"pval_thr"]),])
      }
    }
    #Save numTermsChanged table as xls
    assign(paste0("GO_numTerms.", comparison), numTerms.dfi)
    writexl::write_xlsx(numTerms.dfi, file.path(outputDir, paste0("GO_numTermsChangedGSEA.",comparison, label, ".xlsx")))
  }
  #Save lists in Rdata object
  if (saveRda) save(listGOresults, listGOgeneFC, namescomp, file = file.path(rdaDir, paste0('ABS-GSEA-GOresults-unfilt',label,'.Rda')), version=2)
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
      write.gmt(listgo.geneset, file.path(gmtDir, paste0("GO",categ,"GSEA.gmt")))
      #cal fer aixo, sino no ho reconeix be a enrichmentmap:
      gmt.f <- read.csv(file.path(gmtDir, paste0("GO",categ,"GSEA.gmt")), header=FALSE, sep="\t")
      write.table(gmt.f, file.path(gmtDir, paste0("GO",categ,"GSEA.gmt")), sep="\t", row.names=FALSE, col.names=FALSE)
    })
  }
  #return list of lists
  GSEA_GOresults.unfilt <- list(listGOresults=listGOresults, listGOgeneFC=listGOgeneFC)
  return(GSEA_GOresults.unfilt)
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
