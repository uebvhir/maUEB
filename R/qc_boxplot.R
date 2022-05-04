#' Boxplot of raw or normalized expression data
#'
#' Creates a boxplot from raw or normalized expression data in log2 scale and saves it into a pdf file. See 'Details' for transformations performed on data.
#' @param data A FeatureSet-like object or ExpressionSet object.
#' @param group Name of the targets variable specifying the group to color samples. Defaults to "Group".
#' @param group.color A vector of group colors to be used to colour the bodies of the box plots. Defaults to "Colors".
#' @param samplenames Name of the targets variable specifying the sample names which will be printed under each boxplot. Defaults to "ShortName".
#' @param outputDir Name of the directory where the plot will be saved.
#' @param bp.main General title of the plot. Defaults to "Boxplot for array intensity".
#' @param label A label to be included in the title and/or file name to be saved.
#' @param cex.axis Numeric expansion factor for axis tick labels. Defaults to 0.5.
#' @param bp.ylab Title of y axis. Defaults to "log2 intensity".
#' @param cex.lab Numeric expansion factor for axis labels. Defaults to 0.7.
#' @param las Direction of axis tick labels. Defaults to 2.
#' @param bp.legend boolean indicating whether to insert legend into plot (TRUE/FALSE). Defaults to TRUE.
#' @param bp.legend.posx X coordinate of legend position. Defaults to 10.
#' @param bp.legend.posy Y coordinate of legend position. Defaults to 13.
#' @param bp.legend.cex Numeric expansion factor for legend text. Defaults to 0.6.
#' @param ... other arguments to be passed to boxplot function from oligo package, such as graphical parameters.
#' @details This function uses the 'boxplot' function from 'oligo' package, which takes as input an ExpressionSet. For raw data, expression values are transformed to log2 scale. For normalized data, no transformation is applied since they are already in log2 scale.
#' @note The boxplot methods for FeatureSet and Expression use a sample (via sample) of the probes/probesets to produce the plot. Therefore, the user interested in reproducibility is advised to use set.seed.
#' @importFrom oligo boxplot
#' @export qc_boxplot
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[oligo]{boxplot}
#' @examples
#' qc_boxplot(data=eset_raw, outputDir=resultsDir, label="RawData")
#' @keywords boxplot ExpressionSet FeatureSet
#' @references

qc_boxplot <- function(data, group="Group", group.color="Colors", samplenames="ShortName", outputDir, bp.main="Boxplot for array intensity", label, cex.axis=0.5, bp.ylab="log2 intensity", cex.lab=0.7, las=2, bp.legend=TRUE, bp.legend.posx=10, bp.legend.posy=13, bp.legend.cex=0.6, ...){
    ###Note: takes 'boxplot' function from package 'oligo'. For raw data (ExpressionSet), it transforms intensities to log2-scale (default: transfo=log2)
    targets <- pData(data)
    names <- as.character(targets[,samplenames])
    col <- targets[,group.color]
    ###Note: The boxplot methods for FeatureSet and Expression use a sample (via sample) of the probes/probesets to produce the plot. Therefore, the user interested in reproducibility is advised to use set.seed.
    set.seed(123)
    op <- par(mar = c(8, 5, 4, 2) + 0.1)
    oligo::boxplot(x=data, col=col, cex.axis=cex.axis, las=las, ylab=bp.ylab, cex.lab=cex.lab,
                   names=names, main=paste0(bp.main, ". ", label))
    if (bp.legend) {legend(x=bp.legend.posx, y=bp.legend.posy, legend=unique(targets[,group]), fill=unique(col), cex=bp.legend.cex, bg="white")}
    #ho he de fer dos vegades pq no mel deixa guardar com a objecte
    pdf(file.path(outputDir, paste0("Boxplot", label, ".pdf")))
    op <- par(mar = c(8, 5, 4, 2) + 0.1)
    oligo::boxplot(x=data, col=col, cex.axis=cex.axis, las=las, ylab=bp.ylab, cex.lab=cex.lab,
                   names=names, main=paste0(bp.main, ". ", label))
    if (bp.legend) {legend(x=bp.legend.posx, y=bp.legend.posy, legend=unique(targets[,group]), fill=unique(col), cex=bp.legend.cex, bg="white")}
    dev.off()
    par(op)
}
