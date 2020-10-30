#' Remove batch effects from expression data
#'
#' Removes batch/covariate effects from expression data using linear models (removeBatchEffect function from limma package)
#' @param data A matrix of expression values, with samples in columns and probes in rows.
#' @param targets Targets dataframe with the information of batch variables
#' @param batchFactors Vector with the name of variable(s) from targets dataframe indicating the batches/covariates to be adjusted for
#' @importFrom limma removeBatchEffect
#' @export removebatch
#' @return Returns an expression matrix with batch/covariate effects removed.
#' @note This function is not intended to be used prior to linear modelling. For linear modelling, it is better to include the batch factors in the linear model.
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[limma]{removeBatchEffect}
#' @examples
#' @keywords batch
#' @references

removebatch <- function(data, targets, batchFactors=NULL){
    fmla <- paste("~ ", paste(batchcolName, collapse=" + "), sep="+ ")
    batch.design <- model.matrix(as.formula(fmla), data=targets)
    batch.design <- batch.design[,-1] #treiem intercept
    data <- removeBatchEffect(data, covariates=batch.design) #function from limma package
    return(data)
}
