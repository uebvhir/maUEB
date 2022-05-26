#' Create design matrix, fit the linear model and compute contrasts and moderated t-tests
#'
#' Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
#' @inheritParams dea_lmdesign
#' @inheritParams dea_lmfit
#' @inheritParams dea_compare
#' @param ... Other arguments
#' @details
#' @importFrom limma is.fullrank lmFit duplicateCorrelation contrasts.fit makeContrasts eBayes
#' @export dea_all2
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[maUEB]{dea_lmdesign} \link[maUEB]{dea_lmfit} \link[maUEB]{dea_compare}
#' @examples
#' @return An MArrayLM object containing the result of the fits.
#' @keywords linear model fit contrasts moderated t-test Bayes
#' @references

dea_all2 <- function(data, targets, sampleNames, group, covariates, fmla=NULL, block.cor=NULL, contrastsv, moderated=TRUE, summaryFN, outputDir, ...) {
    design <- dea_lmdesign(targets, sampleNames, data, group, covariates, fmla, summaryFN, outputDir)
    fit <- dea_lmfit(data, targets, design, block.cor, summaryFN, outputDir)
    fit.main <- dea_compare(fit, contrastsv, design, moderated=TRUE, summaryFN, outputDir)
    listofcoef <- colnames(fit.main)
    listofcsv <- dea_toptab(listofcoef=listofcoef, fit.main=fit.main, eset=data, padjust.method="BH",
                            html_report=FALSE, #ReportHTMLTopTab,
                            html_ntop=500, html_group="Group",
                            outputDir=outputDir)
    numGenesChanged_all <- dea_summary_ngc1(listofcsv=listofcsv, listofcoef=listofcoef, B_thr=c(0), Pval_thr=c(0.01,0.05,0.1), adjPval_thr=c(0.01,0.05,0.15,0.25), logFC_thr=0, reg=c("Up", "Down"),
             extended=TRUE, B_thr.e=c(0), Pval_thr.e=c(0.01,0.05,0.1), adjPval_thr.e=c(0.01,0.05,0.15,0.25), logFC_thr.e=c(0,0.5,1,1.5,2), outputDir=outputDir)
    dea.reslist <- list(design=design, fit=fit, fit.main=fit.main, ngc=numGenesChanged_all)
    return(dea.reslist)
}
