#' Create design matrix, fit the linear model and compute contrasts and moderated t-tests
#'
#' Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
#' @inheritParams dea_lmdesign
#' @inheritParams dea_lmfit
#' @inheritParams dea_compare
#' @inheritParams dea_toptab
#' @inheritParams dea_summary_ngc1
#' @inheritParams dea_volcanoplot
#' @details
#' @importFrom limma is.fullrank lmFit duplicateCorrelation contrasts.fit makeContrasts eBayes
#' @import ggplot2 ggrepel
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[maUEB]{dea_lmdesign} \link[maUEB]{dea_lmfit} \link[maUEB]{dea_compare} \link[maUEB]{dea_toptab} \link[maUEB]{dea_summary_ngc1} \link[maUEB]{dea_volcanoplot}
#' @examples
#' @return A list containing the result of the different steps from the differential expression analysis (fit, fit.main, listofcsv, numGenesChanged)
#' @keywords linear model fit contrasts moderated t-test Bayes
#' @references
#' @export

dea_all1 <- function(data, targets, sampleNames, group, covariates, fmla=NULL, block.cor=NULL, contrastsv, moderated=TRUE, padjust.method="fdr",  html_report=FALSE, html_ntop=500, html_group="Group", B_thr=c(0), Pval_thr=c(0.01,0.05,0.1), adjPval_thr=c(0.01,0.05,0.15,0.25), logFC_thr=0, extended=TRUE, B_thr.e=c(0), Pval_thr.e=c(0.01,0.05,0.1), adjPval_thr.e=c(0.01,0.05,0.15,0.25), logFC_thr.e=c(0,0.5,1,1.5,2), dovolcanoplot=TRUE, volc_logFC=rep(1,length(contrastsv)), volc_pval=c(rep("adj.P.Val", length(contrastsv))), volc_pval.thr=rep(0.05, length(contrastsv)), volc_x0=rep(-3, length(contrastsv)), volc_x1=rep(+3, length(contrastsv)), volc_y0=rep(0, length(contrastsv)), volc_y1=rep(10, length(contrastsv)), n=6, cols=2, label="", summaryFN, outputDir){
    design <- dea_lmdesign(targets, sampleNames, data, group, covariates, fmla, summaryFN, outputDir)
    fit <- dea_lmfit(data, targets, design, block.cor, summaryFN, outputDir)
    fit.main <- dea_compare(fit, contrastsv, design, moderated=TRUE, summaryFN, outputDir)
    listofcoef <- colnames(fit.main)
    listofcsv <- dea_toptab(listofcoef=listofcoef, fit.main=fit.main, eset=data, padjust.method=padjust.method,
                            html_report=html_report, #ReportHTMLTopTab,
                            html_ntop=html_ntop, html_group=html_group,
                            outputDir=outputDir)
    numGenesChanged_all <- dea_summary_ngc1(listofcoef=listofcoef, listofcsv=listofcsv, B_thr=B_thr, Pval_thr=Pval_thr, adjPval_thr=adjPval_thr, logFC_thr=logFC_thr, extended=extended, B_thr.e=B_thr.e, Pval_thr.e=Pval_thr.e, adjPval_thr.e=adjPval_thr.e, logFC_thr.e=logFC_thr.e, outputDir=outputDir)
    if (dovolcanoplot) {dea_volcanoplot(listofcsv=listofcsv, listofcoef=listofcoef, volc_logFC=volc_logFC, volc_pval=volc_pval, volc_pval.thr=volc_pval.thr, volc_x0=volc_x0, volc_x1=volc_x1, volc_y0=volc_y0, volc_y1=volc_y1, n=n, cols=cols, outputDir=outputDir, label=label)}
    dea.results <- list(fit=fit, fit.main=fit.main, listofcsv=listofcsv, numGenesChanged=numGenesChanged_all)
    return(dea.results)
}
