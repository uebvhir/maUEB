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
#' @importFrom hwriter hwrite hwriteImage
#' @importFrom WriteXLS WriteXLS
#' @import annotate
#' @import ReportingTools
#' @import ggplot2
#' @export dea_toptab
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[limma]{topTable} \link[ReportingTools]{HTMLReport}
#' @examples
#' @return A list with csv2top
#' @keywords topTable
#' @references

dea_toptab <- function(listofcoef, fit.main, eset, padjust.method="fdr", html_report=FALSE, html_ntop=500, html_group="Group", outputDir){
    #Toptables
    listofTop <- lapply(listofcoef, function(x){topTable(fit.main, number = nrow(fit.main), coef = x, adjust = padjust.method)})
    names(listofTop) <- paste0("topTab.", listofcoef)
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
    names(listofcsv) <- paste0("csv2", names(listofTop))
    #HTML report
    #Create HTML Report of TopTables (top 500 genes included)
    if (html_report){
        #expression data for boxplots in html table
        exprdata <- listofcsv[[1]][,colnames(selectedgenes)]
        # if (batchRemove){
        #     exprdata <- removebatch(data=exprdata, targets=pData(eset), batchFactors="Batch")
        # }
        colstoshow <- c("ProbesetID","Gene.Symbol", "EntrezID", "logFC", "P.Value", "adj.P.Val", "B")
        for (j in 1:length(listofcoef)){
            figsDir <- file.path(outputDir, paste0("figurestopTab.", listofcoef[j]))
            if (dir.exists(figsDir)) unlink(figsDir, recursive = TRUE)
            dir.create(figsDir, showWarnings = FALSE)
            top.df <- listofcsv[[j]][1:html_ntop, colstoshow]
            #add link to entrezID
            top.df$EntrezID <- hwrite(as.character(top.df$EntrezID),
                                     link = paste0("http://www.ncbi.nlm.nih.gov/gene/", as.character(top.df$EntrezID)), table = FALSE)
            #add boxplots for each gene
            imagename <- c()
            for (g in 1:nrow(top.df)){
                imagename[g] <- paste0("boxplot.", top.df[g,"ProbesetID"], ".png")
                genexpr <- as.data.frame(t(exprdata[g,]))
                colnames(genexpr) <- "value"
                genexpr$Group <- pData(eset)[,html_group]
                bp <- ggplot(genexpr, aes(x=Group, y=value)) +
                    geom_boxplot()+
                    geom_point()+
                    theme_classic() +
                    labs(y="log2 normalized expression", title=top.df[g,"ProbesetID"]) +
                    # scale_y_continuous(limits=c(0,max(selectedgenes))) +
                    theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"),
                          axis.text=element_text(size=14), axis.title.y=element_text(vjust=6, size=14),
                          axis.title.x=element_text(vjust=-6, size=14),
                          plot.margin = unit(c(1,1,1,1), "cm"))
                png(file.path(figsDir, imagename[g]))
                # boxplot(genexpr[,1]~genexpr$Group, main=top.df[g,"ProbesetID"], ylab="log2-normalized intensity", xlab=html_group, cex.ylab=0.8)
                print(bp)
                dev.off()
            }
            top.df$Boxplot <- hwriteImage(file.path(figsDir, imagename), link=file.path(figsDir, imagename), table = FALSE, width=50, height=50)
            wd <- getwd()
            setwd(outputDir)
            deReport <- HTMLReport(shortName = paste0("topTab.",listofcoef[j]) ,
                                   title = paste("Analysis of Differential Expression for:",
                                                  listofcoef[j], "(top", html_ntop, "features)"), reportDirectory = ".")
            publish(top.df, deReport)
            finish(deReport)
            setwd(wd)
        }
    }
    return(listofcsv)
}
