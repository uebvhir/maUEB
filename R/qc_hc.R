#' Heatmap and a dendogram of sample distances
#'
#' Computes sample distances from a matrix of expression values and creates a heatmap of sample distances and a dendogram of the hierarchical clustering. The plots are saved into a pdf file.
#' @param data A matrix of expression values, with samples in columns and probes in rows.
#' @param dist.method The distance measure to be used, as used by 'dist()' function. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust.method The agglomeration method as used by 'hclust()' function. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param names A vector with sample names corresponding to columns of the data matrix (in the same order)
#' @param cexRow Size of labels in rows
#' @param cexCol Size of labels in columns
#' @param rect Whether a rectangle is to be plot on hierarchical tree to highlight main clusters
#' @param numClusters Number of main clusters to highlight by rect
#' @param outputDir Name of the directory where the plot will be saved
#' @param label A string to be included in the file name to be saved
#' @export qc_hc
#' @return Creates a heatmap of sample distances and a dendogram of the hierarchical clustering into the graphical device.
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[stats]{hclust}, \link[stats]{dist}, \link[stats]{heatmap}
#' @examples
#' @keywords hierarchical clustering heatmap samples
#' @references

qc_hc <- function(data, dist.method="euclidean", hclust.method="average", names, cexRow = 0.6, cexCol = 0.6, rect=TRUE, numclusters=2, outputDir, label){
    manDist <-  dist(t(data), method=dist.method)
    hc <- hclust(manDist, method=hclust.method)
    #fem els grafics
    heatmap(as.matrix(manDist),  col=heat.colors(16),labCol = names,
            labRow = names, cexRow = cexRow, cexCol = cexCol)
    plot(hc, labels=names,
         main="Hierarchical clustering of samples",  hang=-1, cex = 0.5)
    if (rect) {rect.hclust(hc, k = numclusters, border = "red")}
    #salvem els grafics en un pdf
    pdf(file.path(outputDir, paste0("QC", label, ".pdf")))
    heatmap(as.matrix(manDist),  col=heat.colors(16),labCol = names,
             labRow = names, cexRow = cexRow, cexCol = cexCol)
    plot(hc, labels=names,
         main="Hierarchical clustering of samples",  hang=-1, cex = 0.5)
    if (rect) {rect.hclust(hc, k = numclusters, border = "red")}
    dev.off()
}
