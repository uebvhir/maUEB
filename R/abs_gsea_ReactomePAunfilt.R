#' Gene Set Enrichment Analysis (GSEA) over Reactome Pathways database
#'
#' Performs GSEA analysis using clusterProfiler package over Reactome Pathways database for the results from differential expression analysis ranked by the specified metric. The unfiltered enrichment results are returned and saved as xls files.

#' @param listofTopNamed Named list of toptables
#' @param namescomp Name of the comparisons to be analysed, corresponding to names in list of toptables
#' @param ranking_metric Column from toptable by which features will be ranked
#' @param col_entrez Column from toptable containing EntrezID
#' @param readable Whether features should be converted to gene symbol. Defaults to TRUE.
#' @param organism Organism (eg. 'human', 'rat', 'mouse')
#' @param organism_annot Organism annotation package (e.g. "org.Hs.eg.db", "org.Mm.eg.db")
#' @param Reac_minSetSize Minimum size of the pathways analysed in the analysis
#' @param Reac_maxSetSize Maximum size of the pathways analysed in the analysis
#' @param resultsSummFN
#' @param saveRda
#' @param saveGMT
#' @param label
#' @param outputDir
#' @param rdaDir
#' @param gmtDir
#' @details
#' @import reactome.db AnnotationDbi ReactomePA clusterProfiler DOSE
#' @importFrom writexl write_xlsx
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Creates lists of Reactome unfiltered results and a table with number of terms found at different pvalue thresholds for the different categories analyzed.
#' @keywords Reactome GSEA
#' @references
#' @export

abs_gsea_ReactomePAunfilt <- function(listofTopNamed, namescomp, ranking_metric="logFC", col_entrez="EntrezID", readable=TRUE, organism, organism_annot, Reac_minSetSize=3, Reac_maxSetSize=500, resultsSummFN="ResultsSummary_ABS-GSEA.txt", saveRda=TRUE, saveGMT=TRUE, label="", outputDir, rdaDir, gmtDir){
    require(annotPackage, character.only=TRUE) #required for GO/Reactome
    require(organism_annot, character.only=TRUE) #required for GO/Reactome
    #Prepare lists for storing results
    listReactomePAresults <- list()
    listReactomePAgeneFC <- list()
    #Prepare numTermsChanged table
    numTerms.df <- data.frame(pval_type=character(), pval_thr=numeric(), NumTerms=integer(), stringsAsFactors=FALSE)
    numTerms.df[1:7,"pval_type"] <- c(rep("pvalue", 3), rep("p.adjust", 4))
    numTerms.df[,"pval_thr"] <- c(c(0.01,0.05,0.1), c(0.01,0.05,0.15,0.25))
    #Reactome analysis for the comparisons selected
    for (comparison in namescomp){
        numTerms.dfi <- numTerms.df
        #prepare ranked data. En aquest cas necessita el EntrezID
        topTab <- listofTopNamed[[comparison]]
        geneFC <- topTab[,ranking_metric]
        names(geneFC) <- as.character(topTab[,col_entrez])
        #eliminem els duplicats, abans d'ordenar per logFC (es quedara amb els mes significatius ja que estava ordenada per pval)
        geneFC <- geneFC[!duplicated(names(geneFC))]
        geneFC = sort(geneFC, decreasing = TRUE)
        listReactomePAgeneFC[[comparison]] <- geneFC #save ranked gene names in list
        #analysis
        set.seed(123)
        #Reactome necessita el entrezId com a input.
        gsea.result <- gsePathway(geneFC,
                                  if (as.numeric(R.Version()$major)>=4) eps=0.0 else nPerm=10000,
                                  exponent=1,
                                  seed=FALSE,
                                  pAdjustMethod = "BH", # fix adjustment method so that it computes an adjusted pvalue, later on results will be filtered according to thresholds set in parameters
                                  pvalueCutoff  = 1, # fix to 1, later on results will be filtered according to thresholds set in parameters
                                  minGSSize = Reac_minSetSize,
                                  maxGSSize = Reac_maxSetSize,
                                  verbose = TRUE,
                                  organism = organism)
        #convert from entrez to GeneSymbol
        if (readable) gsea.result <- setReadable(gsea.result, OrgDb = organism_annot)
        #store reactome result in list
        listReactomePAresults[[comparison]] <- gsea.result
        #save reactome result as dataframe (not filtered by pvalue)
        Reactab <- as.data.frame(gsea.result)
        #save unfiltered results
        writexl::write_xlsx(Reactab, file.path(outputDir, paste0("ReactomePAResultsGSEA.", comparison, label, ".unfilt.xlsx")))
        #fill numTermsChanged table with the number of terms found under different pvalue thresholds for GO significance
        for (j in 1:nrow(numTerms.dfi)){
            numTerms.dfi[j, "NumTerms"] <- nrow(Reactab[(Reactab[, numTerms.dfi[j,"pval_type"]] < numTerms.dfi[j,"pval_thr"]),])
        }
        #Save numTermsChanged table as xls
        assign(paste0("Reac_numTerms.", comparison), numTerms.dfi)
        writexl::write_xlsx(numTerms.dfi, file.path(outputDir, paste0("Reactome_numTermsChangedGSEA.", comparison, label, ".xlsx")))
    }
    #Save lists in Rdata object
    if (saveRda) save(listReactomePAresults, listReactomePAgeneFC, namescomp, file = file.path(rdaDir,paste0('ABS-GSEA-ReactomePAresults-unfilt', label, '.Rda')), version=2)
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
        write.gmt(listreac.geneset, file.path(gmtDir, paste0("ReacGSEA.gmt")))
        #cal fer aixo, sino no ho reconeix be a enrichmentmap:
        gmt.f <- read.csv(file.path(gmtDir, paste0("ReacGSEA.gmt")), header=FALSE, sep="\t")
        write.table(gmt.f, file.path(gmtDir, paste0("ReacGSEA.gmt")), sep="\t", row.names=FALSE, col.names = FALSE)
    }
    GSEA_ReactomePAresults.unfilt <- list(listReactomePAresults=listReactomePAresults, listReactomePAgeneFC=listReactomePAgeneFC)
    return(GSEA_ReactomePAresults.unfilt)
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
