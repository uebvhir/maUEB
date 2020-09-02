#' A function to create a boxplot from raw or normalized data
#'
#' Uses the boxplot function from package oligo to create a boxplot of raw or normalized data and saves it into a pdf file. For raw data, expression values are transformed to log2 scale. For normalized data, no transformation is applied since they are already in log2 scale.
#' @param data A FeatureSet-like object or ExpressionSet object.
#' @param group A vector or factor specifying the grouping variable to color samples
#' @param col A vector of colors to be used to colour the bodies of the box plots (for each group)
#' @param names A vector of sample names which will be printed under each boxplot
#' @param outputDir Name of the directory where the plot will be saved
#' @param label A string to be included in the file name to be saved
#' @param cex.axis Numeric expansion factor for axis tick labels
#' @param cex.lab Numeric expansion factor for axis labels
#' @param las Direction of axis tick labels
#' @param plotlegend Whether to insert legend into plot (TRUE/FALSE)
#' @param legend.posx X coordinate of legend position
#' @param legend.posy Y coordinate of legend position
#' @param legend.cex Numeric expansion factor for legend text
#' @param ... other arguments to be passed to boxplot function from oligo package, such as graphical parameters
#' @note The boxplot methods for FeatureSet and Expression use a sample (via sample) of the probes/probesets to produce the plot. Therefore, the user interested in reproducibility is advised to use set.seed.
#' @importFrom oligo boxplot
#' @export qc_boxplot
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[oligo]{boxplot}
#' @examples
#' createBoxplot(data=rawData, group=targets$Group, col=targets$Colors, names=targets$ShortName, outputDir, label="RawData")
#' @keywords boxplot ExpressionSet FEatureSet
#' @references

qc_boxplot <- function(data, group, col, names, outputDir, label, cex.axis=0.5, cex.lab=0.7, las=2, plotlegend=TRUE, legend.posx=10, legend.posy=13, legend.cex=0.6, ...){
    ###Note: takes 'boxplot' function from package 'oligo'. For raw data (ExpressionSet), it transforms intensities to log2-scale (default: transfo=log2)
    ###Note: The boxplot methods for FeatureSet and Expression use a sample (via sample) of the probes/probesets to produce the plot. Therefore, the user interested in reproducibility is advised to use set.seed.
    set.seed(123)
    op <- par(mar = c(8, 5, 4, 2) + 0.1)
    oligo::boxplot(data, col=col, cex.axis=cex.axis, las=las, ylab="log2 intensity", cex.lab=cex.lab,
                   names=names, main=paste0("Boxplot for array intensity: ", label))
    if (plotlegend) {legend(x=legend.posx, y=legend.posy, legend=unique(group), fill=unique(col), cex=legend.cex, bg="white")}
    #ho he de fer dos vegades pq no mel deixa guardar com a objecte
    pdf(file.path(outputDir, paste0("Boxplot", label, ".pdf")))
    oligo::boxplot(data, col=col, cex.axis=cex.axis, las=las, ylab="log2 intensity", cex.lab=cex.lab,
                   names=names, main=paste0("Boxplot for array intensity: ", label))
    if (plotlegend) {legend(x=legend.posx, y=legend.posy, legend=unique(group), fill=unique(col), cex=legend.cex, bg="white")}
    dev.off()
    par(op)
}
