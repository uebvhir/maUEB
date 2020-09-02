#' A function to save probeset annotations
#'
#' Uses aafTableAnn function from annaffy package to construct an aafTable object given a set of probe ids. Saves annotation for all summarized probes as csv and html files
#' @param data character vector of probe ids
#' @param annotPkg Annotation package to be used for annotations
#' @param outputFN File name of the datasheet saved containing annotation data, without extension.
#' @param saveHTML Whether annotation datasheet should be saved also as html (TRUE/FALSE). Defaults to TRUE.
#' @param title Title included in the html file.
#' @param outputDir Name of the directory where the files will be saved
#' @note In Exon studies, where probes are summarized at the probeset level but not transcript level, the probe annotations cannot be recovered using this function.
#' @return An aafTable object
#' @importFrom annaffy aafTableAnn
#' @export save_annotations
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[annaffy]{aafTableAnn}
#' @examples
#' @keywords annotation
#' @references

save_annotations <- function(data, annotPkg, outputFN="Annotations.AllGenes", saveHTML=TRUE, title="Annotations for all genes", outputDir){
    at <- aafTableAnn(data, annotPkg)
    saveHTML(at, file.path(outputDir, paste0(outputFN, ".html")), title)
    saveText(at, file.path(outputDir, paste0(outputFN, ".csv")), header=TRUE)
    return(at)
}


