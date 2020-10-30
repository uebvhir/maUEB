#' Design matrix to be used in the linear model
#'
#' Constructs a design matrix from the main factor and covariates specified in parameters. A formula can be directly provided. The matrix is constructed by expanding factors to a set of dummy variables resulting in a matrix with rows corresponding to samples and columns to factors.
#' @param targets Targets dataframe with the phenotypic information of variables and covariates to build the linear model. Can be extracted with pData(eset)
#' @param sampleNames Name of variable in targets containing the sample names (the same identificators as the column names of the expression matrix)
#' @param data An ExpressionSet instance or a matrix of expression values, with samples in columns and probes in rows, to which the linear model will be fitted. Column names should be the same and in the same order as the targets sample names used to build the design matrix.
#' @param group Name of variable in targets containing the main factor of interest (usually: 'Group')
#' @param covariates Vector with the name of variable(s) to be used as covariates in the linear model. If not to consider any variable, set NULL.
#' @param fmla Formula to be used in the model to construct the design matrix (eg. ' ~ 0 + Group'). Default is NULL. If no formula is provided, it will be automatically built from the main factor and covariates specified in parameters such as: '0 + group + covariates'.
#' @param summaryFN Name of the file where the reporting summary of execution results will be saved
#' @param outputDir Name of the directory where the results summary will be saved
#' @importFrom limma is.fullrank
#' @export dea_lmdesign
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[stats]{model.matrix}
#' @examples
#' @return A design matrix
#' @keywords linear model design matrix
#' @references

dea_lmdesign <- function(targets, sampleNames, data, group, covariates, fmla=NULL, summaryFN, outputDir){
    #Before doing the design and contrast matrix, verify that order and number of samples in targets is the same as in the expression set
    #if not, you will need to arrange to make them coincide
    #(reason: it constructs the design matrix based only on targets sample name and then the fit is done based on ordering positions)
    if (sum(targets[,sampleNames]!=colnames(data))>0){
        stop("Attention: order of samples in targets and expression set is not the same. The model will not be constructed correctly.")
    }
    #Design matrix
    designfactors <- c(group, covariates)
    ##convert to factors
    for (i in designfactors){
        targets[,i] <- factor(targets[,i], levels=unique(targets[,i]))
    }
    ##Formula for linear model
    if (is.null(fmla)){
        fmla <- paste("~ 0 ", paste(designfactors, collapse=" + "), sep="+ ")
    }
    design <- model.matrix(as.formula(fmla), data=targets)
    ##Change names of "Group" columns of design matrix according to group levels
    colnames(design) <- gsub(group, "", colnames(design))
    rownames(design) <- targets[,sampleNames]
    if (!is.fullrank(as.matrix(design))) warning("The design matrix is not full ranked. You should remove correlated covariates.")
    #Save design and contrast matrix and fit in Summary
    sink(file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    cat("\nFormula of linear model: ", fmla)
    cat("\nDESIGN MATRIX:\n")
    cat("-Dimensions:", dim(design))
    cat("\n-Rank:", qr(as.matrix(design))$rank)
    cat("\n-Full rank:", is.fullrank(as.matrix(design)), "\n")
    # as.data.frame(design)
    sink()
    capture.output(design, file=file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    return(design)
}
