#' A function to assess the source of batch effects in a microarray gene expression experiment
#'
#' Performs a Principal Variance Component Analysis (PVCA) using 'pvcaBatchAssess' function from 'pvca' package. A barplot is created with the percentage of total variance associated to each factor assessed. The effects associated to each factor are calculated by fitting all sources as random effects including two-way interaction terms in the Mixed Model to select principal components. See 'Details' for extended information.
#' @param data An instance of ExpressionSet.
#' @param factors A vector with the names of factor variables to be assessed
#' @param targetsPVCA If a different targets is to be used for PVCA, the name of the new targets csv file should be specified here. See 'Details'.
#' @param pct_threshold The percentile value of the minimum amount of the variabilities that the selected principal components need to explain
#' @param label A string to be included in the file name to be saved
#' @param inputDir Name of the directory from where the targetsPVCA csv file should be loaded, if specified.
#' @param outputDir Name of the directory where the plot will be saved
#' @param summaryFN Name of the file where the reporting summary of execution results will be saved
#' @details The approach leverages the strengths of two very popular data analysis methods: first, principal component analysis (PCA) is used to efficiently reduce data dimension while maintaining the majority of the variability in the data, and variance components analysis (VCA) fits a mixed linear model using factors of interest as random effects to estimate and partition the total variability.
#' @note At least two factors should be provided. In addition, for each factor, num of levels must be < num of samples/level, if not it gives an error. To overcome this issue you can create a simplified version of targets file with less levels just to create the pvca.
#' @import pvca
#' @export qc_pvca
#' @return Creates a barplot into the graphical device and returns a list with the loads of PCA.
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[pvca]{pvcaBatchAssess}
#' @examples
#' @keywords PVCA
#' @references

qc_pvca <- function(data, factors, targetsPVCA=NULL, pct_threshold=0.6, label=NULL, inputDir=NULL, outputDir, summaryFN){
    if (length(factors) < 2){stop("At least two factors are needed.")}
    if (!is.null(targetsPVCA) & (targets.pvcaFN != targetsFN)){
        targets.pvca <- read.table(file.path(inputDir, targetsPVCA), header = TRUE, row.names = 1, sep = ";", as.is=TRUE)
        pData(data.pvca) <- targets.pvca
    }
    #compute pvca and print output success in device and in ResultsSummary.txt
    (pvcaOutmessage <- try(pvcaObj <- pvcaBatchAssess (data, factors, pct_threshold), TRUE))
    sink(file.path(outputDir, summaryFN), split = FALSE, append = TRUE)
    cat("\nPVCA:\n")
    print(pvcaOutmessage)
    cat("\n")
    sink()
    #pvca barplot
    if (exists("pvcaObj")){
        #Create barplot and show plot on device
        bp <- barplot(pvcaObj$dat, xlab = "Effects",
                      ylab = "Weighted average proportion variance",
                      ylim= c(0,1.1),col = c("mediumorchid"), las=2,
                      main="PVCA estimation")
        axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.55, las=2)
        values = pvcaObj$dat
        new_values = round(values , 3)
        text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.5)
        #save plot in pdf
        pdf(file.path(outputDir, paste0("PVCA", label, ".pdf")))
        bp <- barplot(pvcaObj$dat, xlab = "Effects",
                      ylab = "Weighted average proportion variance",
                      ylim= c(0,1.1),col = c("mediumorchid"), las=2,
                      main="PVCA estimation")
        axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.55, las=2)
        values = pvcaObj$dat
        new_values = round(values , 3)
        text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.5)
        dev.off()
        return(pvcaObj)
    } else {return(pvcaOutmessage)}
}
