#' Summarise results from differential expression analysis
#'
#' Create summary table of number of genes changed in the different contrasts under different thresholds
#' @inheritParams dea_summary_ngc
#' @inheritParams dea_summary_ngc_ext
#' @param extended a logical value indicating whether extended summary tables should be created for each comparison. See \link{dea_summary_ngc_ext} for details.
#' @details
#' @importFrom WriteXLS WriteXLS
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link{dea_summary_ngc} \link{dea_summary_ngc_ext}
#' @examples
#' @return A dataframe with summaries for each toptable
#' @keywords numGenesChanged
#' @references
#' @export
dea_summary_ngc1 <- function(listofcsv, listofcoef, B_thr=c(0), Pval_thr=c(0.01,0.05,0.1), adjPval_thr=c(0.01,0.05,0.15,0.25), logFC_thr=0, reg=c("Up", "Down"),
                             extended=TRUE, B_thr.e=B_thr, Pval_thr.e=Pval_thr, adjPval_thr.e=adjPval_thr, logFC_thr.e=c(0,0.5,1,1.5,2), outputDir){
    ngc <- dea_summary_ngc(listofcsv=listofcsv, listofcoef=listofcoef, B_thr=B_thr, Pval_thr=Pval_thr, adjPval_thr=adjPval_thr, logFC_thr=logFC_thr, reg=reg, outputDir=outputDir)
    if (extended){
        ngc_list <- dea_summary_ngc_ext(listofcsv=listofcsv, listofcoef=listofcoef, B_thr.e=B_thr.e, Pval_thr.e=Pval_thr.e, adjPval_thr.e=adjPval_thr.e, logFC_thr.e=logFC_thr.e, outputDir=outputDir)
        ngc <- list(numGenesChanged=ngc, numGenesChanged_extended=ngc_list)
    }
    return(ngc)
}
