#' A function to normalize microarray data
#'
#' Normalizes data using one of the available methods (see 'Details') and returns a FeatureSet of normalized data. The normalized expression values are saved into a csv/xls file.
#' @param data An ExpressionSet object with the raw data to be normalized.
#' @param norm.method Method used for normalization (currently: oligo::rma)
#' @param exonStudy If the analysis is to be performed at exon level (TRUE/FALSE). Defaults to FALSE.
#' @param norm.level Level of summarization (only for Exon/Gene arrays). Defaults to 'core' for summarization at transcript level. To summarize at probe level use target='probeset'.
#' @param bkgremove Whether to perform RMA background correction (TRUE/FALSE). Defaults to TRUE.
#' @param annotPkg Annotation package to be used to annotate the FeatureSet generated
#' @param outputFN File name of the datasheet saved containing normalized data (csv or xls file), without extension.
#' @param outputDir Name of the directory where the datasheet of normalized data will be saved
#' @details Currently using rma method from oligo package.
#' @importFrom oligo rma
#' @importFrom WriteXLS WriteXLS
#' @export normalization
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[oligo]{rma}
#' @examples
#' @return eset_norm An R object containing the normalized data (a FeatureSet)
#' @keywords normalization rma
#' @references

normalization <- function(data, norm.method="rma", exonStudy=FALSE, norm.level='core', bkgremove=TRUE, annotPkg=annotPackage, outputFN="Normalized.All", outputDir){
    if (norm.method=="rma"){
        if (exonStudy){
            eset_norm <- oligo::rma(data, target=norm.level, background=bkgremove)
        } else {eset_norm <- oligo::rma(data, background=bkgremove)}
    }
    #incloure aqui altres metodes possibles de normalitzacio
    ##Annotation of normalized expression set
    annotation(eset_norm) <- annotPkg
    ##Save normData as csv and xls if possible
    normData <- as.data.frame(exprs(eset_norm))
    write.csv2(normData, file.path(outputDir, paste0(outputFN, ".csv"))) #millor guardar csv que xls ja que xls nomes permet guardar 65535 files
    if ((nrow(normData)<=65535) & (ncol(normData)<=256)){
        WriteXLS(normData, ExcelFileName = file.path(outputDir, paste0(outputFN, ".xls")), row.names = TRUE)
    }
    return(eset_norm)
}


