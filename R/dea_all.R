#' Create design matrix, fit the linear model and compute contrasts and moderated t-tests
#'
#' Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
#' @inheritParams dea_lmdesign
#' @inheritParams dea_lmfit
#' @inheritParams dea_compare
#' @details
#' @importFrom limma is.fullrank lmFit duplicateCorrelation contrasts.fit makeContrasts eBayes
#' @export dea_all
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[maUEB]{dea_lmdesign} \link[maUEB]{dea_lmfit} \link[maUEB]{dea_compare}
#' @examples
#' @return An MArrayLM object containing the result of the fits.
#' @keywords linear model fit contrasts moderated t-test Bayes
#' @references

dea_all <- function(data, targets, sampleNames, group, covariates, fmla=NULL, block.cor=NULL, contrastsv, moderated=TRUE, summaryFN, outputDir){
    design <- dea_lmdesign(targets, sampleNames, data, group, covariates, fmla, summaryFN, outputDir)
    fit <- dea_lmfit(data, targets, design, block.cor, summaryFN, outputDir)
    fit.main <- dea_compare(fit, contrastsv, design, moderated=TRUE, summaryFN, outputDir)
    dea.reslist <- list(design=design, fit=fit, fit.main=fit.main)
    return(dea.reslist)
}
