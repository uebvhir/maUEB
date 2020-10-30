#' Summarise results from differential expression analysis
#'
#' Create summary table of number of genes changed in the different contrasts under different thresholds
#' @param listofcsv List of toptables to summarize
#' @param listofcoef Name of comparisons, in the same order as listofcsv
#' @param B_thr Numeric vector of thresholds to be used for B statistic
#' @param Pval_thr Numeric vector of thresholds to be used for the p-value
#' @param adjPval_thr Numeric vector of thresholds to be used for the adjusted p-value
#' @param logFC_thr A number indicating the threshold to be used for logFC
#' @param reg A character vector indicating the direction of logFC change ("Up", "Down", and/or "All"). Default is c("Up", "Down"). If 'All' is indicated, it will take the absolute value of logFC.
#' @param outputDir Name of the directory where the results summary will be saved
#' @details
#' @importFrom WriteXLS WriteXLS
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link{dea_summary_ngc_ext}
#' @examples
#' @return A dataframe with summaries for each toptable
#' @keywords numGenesChanged
#' @references
#' @export
dea_summary_ngc <- function(listofcsv, listofcoef, B_thr=c(0), Pval_thr=c(0.01,0.05,0.1), adjPval_thr=c(0.01,0.05,0.15,0.25), logFC_thr=0, reg=c("Up", "Down"), outputDir){
    stats_thr <- list(B=B_thr, adjPval=adjPval_thr, Pval=Pval_thr)
    #set names as defined in columns of top table
    stats_names <- c("B", "adj.P.Val", "P.Value") #names(stats_thr)
    stats_lengths <- sapply(stats_thr, length)
    #initialize dataframe
    numGenes.df <- data.frame(stat_type=character(), stat_thr=numeric(), reg=character(), stringsAsFactors=FALSE)
    numrows <- sum(stats_lengths)*length(reg)
    #fill parameters in dataframe
    numGenes.df[1:numrows,"reg"] <- rep(reg, sum(stats_lengths))
    numGenes.df[,"stat_type"] <- unlist(mapply(function(x,y) rep(x,y), stats_names, stats_lengths*length(reg)))
    numGenes.df[,"stat_thr"] <- unlist(mapply(function(x,y) rep(x,each=y), stats_thr, length(reg)))
    numGenes.df$Stat <- paste0(numGenes.df$stat_type, numGenes.df$stat_thr, "_", numGenes.df$reg)
    for (i in 1:length(listofcsv)){
        top <- listofcsv[[i]]
        comp <- listofcoef[i]
        for (j in 1:nrow(numGenes.df)){
            if (numGenes.df[j,"stat_type"]=="B"){
                if (numGenes.df[j,"reg"]=="Up") {
                    numGenes.df[j,comp] <- nrow(top[(top[, numGenes.df[j,"stat_type"]] > numGenes.df[j,"stat_thr"] & top$logFC > logFC_thr),])
                }
                if (numGenes.df[j,"reg"]=="Down") {
                    numGenes.df[j,comp] <- nrow(top[(top[, numGenes.df[j,"stat_type"]] > numGenes.df[j,"stat_thr"] & top$logFC < -logFC_thr),])
                }
                if (numGenes.df[j,"reg"]=="All") {
                    numGenes.df[j,comp] <- nrow(top[(top[, numGenes.df[j,"stat_type"]] > numGenes.df[j,"stat_thr"] & abs(top$logFC) > logFC_thr),])
                }
                } else {
                    if (numGenes.df[j,"reg"]=="Up") {
                        numGenes.df[j,comp] <- nrow(top[(top[, numGenes.df[j,"stat_type"]] < numGenes.df[j,"stat_thr"] & top$logFC > logFC_thr),])
                    }
                    if (numGenes.df[j,"reg"]=="Down") {
                        numGenes.df[j,comp] <- nrow(top[(top[, numGenes.df[j,"stat_type"]] < numGenes.df[j,"stat_thr"] & top$logFC < -logFC_thr),])
                    }
                    if (numGenes.df[j,"reg"]=="All") {
                        numGenes.df[j,comp] <- nrow(top[(top[, numGenes.df[j,"stat_type"]] < numGenes.df[j,"stat_thr"] & abs(top$logFC) > logFC_thr),])
                    }
            }
        }
    }
    numGenes.df <- numGenes.df[,c("Stat", listofcoef)]
    WriteXLS(numGenes.df, ExcelFileName = file.path(outputDir, paste("numGenesChanged", "xls", sep=".")),
                 row.names = FALSE)
    return(numGenes.df)
}
