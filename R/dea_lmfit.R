#' Fit multiple linear models to expression data
#'
#' This function uses 'lmFit' function from the 'limma' package to fit a linear model for each gene/probeset given a series of arrays. See 'Details' for extended information.
#' @param data An ExpressionSet instance or a matrix of log-expression values, with samples in columns and probes in rows, to which the linear model will be fitted. Column names should be the same and in the same order as the targets sample names.
#' @param targets Targets dataframe with the phenotypic information of the blocking factor (if provided) to be used to build the linear model. Can be extracted with pData(eset).
#' @param design the design matrix of the microarray experiment, with rows corresponding to arrays and columns to coefficients to be estimated.
#' @param block.cor Name of variable in targets containing the blocking factor accounting for random effects (if it exists) for which an inter-block correlation consensus will be computed and fit to the model. Default is NULL.
#' @param summaryFN Name of the file where the reporting summary of execution results will be saved
#' @param outputDir Name of the directory where the results summary will be saved
#' @details A linear model is fitted to the expression data for each probe. The expression data should be log-ratios for two-color array platforms or log-expression values for one-channel platforms.
#' The coefficients of the fitted models describe the differences between the RNA sources hybridized to the arrays.
#'
#' The function allows for missing values and accepts quantitative precision weights through the weights argument.
#'
#' If block is not NULL then different arrays are assumed to be correlated and a correlation consensus is computed and fit to the model. The value of the correlation consensus obtained is printed to the results summary and returned as field of the object 'fit' that is returned.
#' @importFrom limma lmFit
#' @importFrom limma duplicateCorrelation
#' @export dea_lmfit
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[limma]{lmFit}
#' @examples
#' @return An MArrayLM object containing the result of the fits.
#' @keywords linear model fit
#' @references

dea_lmfit <- function(data, targets, design, block.cor=NULL, summaryFN, outputDir){
   if (!is.null(block.cor)){
       #compute the correlation factor
        targets[, block.cor] <- factor(targets[,block.cor], levels=unique(targets[,block.cor]))
        corfit <- duplicateCorrelation(data, design, block=targets[,block.cor])
        fit <- lmFit(data, design, block=targets[, block.cor], correlation=corfit$consensus)
        sink(file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
        cat("\nFIT:\n")
        cat("Blocking correlation factor (consensus):", corfit$consensus, "\n")
        # fit
        # cat("\n")
        sink()
        capture.output(fit, file=file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    } else {
        fit <- lmFit(data, design)
        sink(file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
        cat("\nFIT:\n")
        # fit
        # cat("\n")
        sink()
        capture.output(fit, file=file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    }
    return(fit)
}
