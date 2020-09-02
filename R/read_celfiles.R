#' A function to read cel files and create an ExpressionSet of raw data
#'
#' Reads the cel files specified in targets and creates an ExpressionSet
#' @param targets Targets object containing the filename of the celfiles to read and phenotypic data
#' @param inputDir Name of the directory containing the .CEL files to read.
#' @import Biobase oligo
#' @export read_celfiles
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso
#' @examples
#' @return rawData An ExpressionSet containing raw data
#' @keywords read cel ExpressionSet
#' @references

read_celfiles <- function(inputDir, targets){
    celfilesFN <- paste(inputDir, rownames(targets), sep="/")
    rawData <- read.celfiles(celfilesFN)
    pData(rawData) <- targets
    colnames(rawData) <- rownames(pData(rawData)) <- targets$ShortName
    return(rawData)
}
