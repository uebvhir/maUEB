#' Save toptables
#'
#' Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
#' @param listofcoef Coefficient names
#' @param fit.main MArrayLM object resulting from linear model fit
#' @param eset ExpressionSet object containing expression values used to fit the model
#' @param padjust.method Method for pvalue adjustment.
#' @param html_report true/false
#' @param html_ntop number of top genes
#' @param html_group group for boxplot
#' @param html_padjust_method pval to show in htmlreport
#' @param outputDir Name of the directory where the results summary will be saved
#' @details
#' @importFrom limma topTable
#' @import annotate
#' @import ReportingTools
#' @import WriteXLS
#' @export dea_toptab
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[limma]{topTable} \link[ReportingTools]{HTMLReport}
#' @examples
#' @return A list with csv2top
#' @keywords topTable
#' @references

dea_toptab <- function(listofcoef, fit.main, eset, padjust.method="fdr", html_report=FALSE, html_ntop=500, html_group="Group", html_padjust_method, outputDir){
    #Toptables
    listofTop <- lapply(listofcoef, function(x){topTable(fit.main, number = nrow(fit.main), coef = x, adjust = padjust.method)})
    names(listofTop) <- listofcoef
    #Annotate and add expression values as additional columns in a dataframe ('csvtopTab.XXX')
    selectedgenes <- as.data.frame(exprs(eset))
    listofcsv <- lapply(seq_along(listofTop), function(i){
        top <- listofTop[[i]]
        EntrezsA <- getEG (rownames(top), annotation(eset))
        SymbolsA <- getSYMBOL (rownames(top), annotation(eset))
        top.annot <- cbind.data.frame(SymbolsA, EntrezsA, top)
        top.annot <- merge(top.annot, selectedgenes, by=0)
        names(top.annot) <- c("ProbesetID","Gene.Symbol", "EntrezID", colnames(top),
                              colnames(selectedgenes))
        top.annot <- top.annot[order(-top.annot$B),]
        #posem el probeID com a rownames
        rownames(top.annot) <- top.annot$ProbesetID
        #Save annotated topTabs as xls or csv
        if ((nrow(top.annot)<=65535) & (ncol(top.annot)<=256)){
            WriteXLS(top.annot,
                     ExcelFileName = file.path(outputDir, paste("ExpressionAndTop_", listofcoef[i],".xls",sep="")), row.names = FALSE)
        } else {
            write.csv2(top.annot, file.path(outputDir, paste("ExpressionAndTop_", listofcoef[i],".csv",sep="")))
        }
        return(top.annot)})
    names(listofcsv) <- listofcoef
    #HTML report
    #Create HTML Report of TopTables (top 500 genes included)
    if (html_report){
        for(i in 1:length(listofcoef)){
            deReport <- HTMLReport(shortName = paste0("topTab.",listofcoef[i]) ,
                                   title = paste0("Analysis of Differential Expression for: TopTab.",
                                                  listofcoef[i]), reportDirectory = paste("./",basename(outputDir),sep=""))
            publish(fit.main, deReport, eSet = eset, factor = pData(eset)[,html_group], coef = listofcoef[i],
                    n = html_ntop, pvalueCutoff = 1, adjust.method=html_padjust_method[i]) #Gives the adjusted pvalue by default. For raw pvalue, use 'adjust.method="none"'. Then change column name in html file using a text editor.
            finish(deReport)
            #if adjust.method="none", the following 'if' modifies name of column from "adjusted pvalue" to "Pvalue" in html file
            if (html_padjust_method[i]=="none"){
                f <- file.path(resultsDir, paste0("topTab.",listofcoef[i], ".html"))
                x <- readLines(f)
                y <- gsub("Adjusted p-Value", "P.Value", x )
                cat(y, file=f, sep="\n")
            }
        }
    }
    return(listofcsv)
}
