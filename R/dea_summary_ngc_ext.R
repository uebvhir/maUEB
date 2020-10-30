#' Summarise results from differential expression analysis
#'
#' Create summary table of number of genes changed in the different contrasts under different thresholds
#' @param listofcsv List of toptables to summarize
#' @param listofcoef Name of comparisons, in the same order as listofcsv
#' @param B_thr.e Numeric vector of thresholds to be used for B statistic in the extended mode
#' @param Pval_thr.e Numeric vector of thresholds to be used for the p-value in the extended mode
#' @param adjPval_thr.e Numeric vector of thresholds to be used for the adjusted p-value in the extended mode
#' @param logFC_thr.e Numeric vector of thresholds to be used for the logFC in the extended mode
#' @param outputDir Name of the directory where the results summary will be saved
#' @details
#' @importFrom WriteXLS WriteXLS
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link{dea_summary_ngc}
#' @examples
#' @return A list with summaries for each toptable
#' @keywords numGenesChanged
#' @references
#' @export
dea_summary_ngc_ext <- function(listofcsv, listofcoef, B_thr.e=c(0), Pval_thr.e=c(0.01,0.05,0.1), adjPval_thr.e=c(0.01,0.05,0.15,0.25), logFC_thr.e=c(0,0.5,1,1.5,2), outputDir){
    stats_thr <- list(B=B_thr.e, adjPval=adjPval_thr.e, Pval=Pval_thr.e)
    #set names as defined in columns of top table
    stats_names <- c("B", "adj.P.Val", "P.Value") #names(stats_thr)
    stats_lengths <- sapply(stats_thr, length)
    #initialize dataframe
    numGenes.df <- data.frame(stat_type=character(), stat_thr=numeric(), logFC_thr=numeric(), NumGenes=integer(), NumGenesUp=integer(), NumGenesDown=integer(), stringsAsFactors=FALSE)
    numrows <- sum(stats_lengths)*length(logFC_thr.e)
    #fill parameters in dataframe
    numGenes.df[1:numrows,"logFC_thr"] <- rep(logFC_thr.e, sum(stats_lengths))
    numGenes.df[,"stat_type"] <- unlist(mapply(function(x,y) rep(x,y), stats_names, stats_lengths*length(logFC_thr.e)))
    numGenes.df[,"stat_thr"] <- unlist(mapply(function(x,y) rep(x,each=y), stats_thr, length(logFC_thr.e)))
    ngc_list <- list()
    for (i in 1:length(listofcsv)){
        top <- listofcsv[[i]]
        comp <- listofcoef[i]
        numGenes.dfi <- numGenes.df
        for (j in 1:nrow(numGenes.dfi)){
            if (numGenes.dfi[j,"stat_type"]=="B"){
                numGenes.dfi[j,"NumGenes"] <- nrow(top[(top[, numGenes.dfi[j,"stat_type"]] > numGenes.dfi[j,"stat_thr"] & abs(top$logFC) > numGenes.dfi[j,"logFC_thr"]),])
                numGenes.dfi[j,"NumGenesUp"] <- nrow(top[(top[, numGenes.dfi[j,"stat_type"]] > numGenes.dfi[j,"stat_thr"] & top$logFC > numGenes.dfi[j,"logFC_thr"]),])
                numGenes.dfi[j,"NumGenesDown"] <- nrow(top[(top[, numGenes.dfi[j,"stat_type"]] > numGenes.dfi[j,"stat_thr"] & top$logFC < -numGenes.dfi[j,"logFC_thr"]),])
            } else {
            numGenes.dfi[j,"NumGenes"] <- nrow(top[(top[, numGenes.dfi[j,"stat_type"]] < numGenes.dfi[j,"stat_thr"] & abs(top$logFC) > numGenes.dfi[j,"logFC_thr"]),])
            numGenes.dfi[j,"NumGenesUp"] <- nrow(top[(top[, numGenes.dfi[j,"stat_type"]] < numGenes.dfi[j,"stat_thr"] & top$logFC > numGenes.dfi[j,"logFC_thr"]),])
            numGenes.dfi[j,"NumGenesDown"] <- nrow(top[(top[, numGenes.dfi[j,"stat_type"]] < numGenes.dfi[j,"stat_thr"] & top$logFC < -numGenes.dfi[j,"logFC_thr"]),])
            }
        }
        ngc_list[[i]] <- numGenes.dfi
        WriteXLS(numGenes.dfi, ExcelFileName = file.path(outputDir, paste("numGenesChangedextended", i, comp, "xls", sep=".")),
                 row.names = FALSE)
    }
    names(ngc_list) <- listofcoef
    return(ngc_list)
}
