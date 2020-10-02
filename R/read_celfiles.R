#' Read cel files and create an ExpressionSet of raw data
#'
#' Reads the cel files specified in FileName column of targets dataframe and creates an ExpressionSet
#' @param targets Targets object containing the filename of the celfiles to read and phenotypic data
#' @param inputDir Name of the directory containing the .CEL files to read.
#' @import Biobase oligo
#' @export read_celfiles
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso
#' @examples
#' @return Returns an ExpressionSet containing raw data
#' @keywords read cel ExpressionSet
#' @references

read_celfiles <- function(inputDir, targets){
    celfilesFN <- paste(inputDir, as.character(targets$FileName), sep="/")
    eset_raw <- read.celfiles(celfilesFN)
    pData(eset_raw) <- targets
    colnames(eset_raw) <- rownames(pData(eset_raw)) <- targets$ShortName
    return(eset_raw)
}
