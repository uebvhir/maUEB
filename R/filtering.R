#' Filtering of features in an ExpressionSet
#'
#' Uses 'nsFilter' function to remove features from an ExpressionSet. Data can be filtered by annotation and/or variance.
#' @param data An instance of ExpressionSet.
#' @inheritParams genefilter::nsFilter
#' @param outputFN File name of the datasheet saved containing normalized data (csv or xls file), without extension.
#' @param outputDir Name of the directory where the plot will be saved
#' @param summaryFN Name of the file where the reporting summary of execution results will be saved
#' @note For memory reasons, only half of the values are plotted.
#' @importFrom genefilter nsFilter
#' @importFrom WriteXLS WriteXLS
#' @export filtering
#' @return Returns an ExpressionSet with filtered features. The filtered expression values are saved into a csv/xls file.
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[genefilter]{nsFilter}
#' @examples
#' @keywords filter
#' @references

filtering <- function(data, outputFN, outputDir, summaryFN, feature.exclude="^AFFX",require.entrez=TRUE, remove.dupEntrez = TRUE,
                         var.filter=FALSE, var.func = IQR, var.cutoff = 0.5, filterByQuantile=TRUE,
                         require.GOBP = FALSE, require.GOCC = FALSE, require.GOMF = FALSE){
    filtered <- nsFilter(data, feature.exclude=feature.exclude,require.entrez=require.entrez, remove.dupEntrez = remove.dupEntrez,
                              var.filter=var.filter, var.func = eval(parse(text=var.func)), var.cutoff = var.cutoff, filterByQuantile=filterByQuantile,
                              require.GOBP = require.GOBP, require.GOCC = require.GOCC, require.GOMF = require.GOMF)
    eset_filtered <- filtered$eset
    filteredData <- as.data.frame(exprs(eset_filtered))
    #Save filtered data as csv and/or xls
    write.csv2(filteredData, file.path(outputDir,  paste0(outputFN, ".csv")))
    if ((nrow(filteredData)<=65535) & (ncol(filteredData)<=256)){
        WriteXLS(filteredData, ExcelFileName=file.path(outputDir, paste0(outputFN, ".xls")), row.names = TRUE)
    }
    #Save summary
    sink(file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    cat("\nNON-SPECIFIC FILTERING\n")
    cat("-Feature filtering:", feature.exclude, "\n")
    cat("-Filtering by Annotation:", require.entrez, "\n")
    cat("-Filtering by Variance:", var.filter, "\n")
    if (var.filter){
        cat("--Function used to filter by variance:", var.func, "\n")
        cat("--Percentage of features removed among total:", var.cutoff*100, "%\n")
    }
    cat("-Number of features before filtering:", nrow(exprs(data)), "\n")
    cat("-Number of features after filtering:", nrow(filteredData), "\n")
    print(as.data.frame(filtered$filter.log))
    cat("\n")
    sink()
    return(eset_filtered)
}
