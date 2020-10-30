#' Compute contrasts and moderated t-tests
#'
#' Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
#' @param fit An MArrayLM object created by lmFit
#' @param contrastsv A character vector with the contrasts to be performed (eg. c("Contrast1=Group2-Group1", "Contrast2=Group3vsGroup1"))
#' @param design the design matrix of the microarray experiment, with rows corresponding to arrays and columns to coefficients to be estimated.
#' @param moderated a logical value indicating whether empirical Bayes moderation of the standard errors towards a common value should be performed. Defaults to TRUE. See 'Details' for further information.
#' @param summaryFN Name of the file where the reporting summary of execution results will be saved
#' @param outputDir Name of the directory where the results summary will be saved
#' @details
#' @importFrom limma contrasts.fit makeContrasts eBayes
#' @export dea_compare
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[limma]{contrasts.fit} \link[limma]{eBayes}
#' @examples
#' @return An MArrayLM object containing the result of the fits.
#' @keywords linear model fit contrasts moderated t-test Bayes
#' @references

dea_compare <- function(fit, contrastsv, design, moderated=TRUE, summaryFN, outputDir){
    #Contrasts matrix
    contrastsMatrix <- eval(parse(text=paste("makeContrasts(", paste(contrastsv, collapse=", "), ", levels=design)")))
    #then compute the contrasts and moderated t-tests
    fit.main <- contrasts.fit(fit, contrastsMatrix)
    if (moderated) fit.main <- eBayes(fit.main)
    #Save design and contrast matrix and fit in Summary
    sink(file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    cat("\nCONTRASTS MATRIX:\n")
    sink()
    # print(as.data.frame(contrastsMatrix))
    capture.output(contrastsMatrix, file=file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    sink(file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    cat("\nEmpirical Bayes moderation of the errors performed?", moderated, "\n")
    cat("\nFIT.MAIN:\n")
    # fit.main
    # cat("\n")
    sink()
    capture.output(fit.main, file=file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    return(fit.main)
}
