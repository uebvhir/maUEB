#' Quality Control plots and report generation for raw or normalized expression data
#'
#' Performs different Quality Control plots and report generation on the given data matrix.
#' @param data An ExpressionSet object containing expression data, and the targets as phenotypic data.
#' @param group Name of the targets variable specifying the group factor to color samples (eg. "Group")
#' @param group.color A vector of group colors to be used to colour the bodies of the box plots
#' @param samplenames Name of the targets variable specifying the sample names which will be printed under each boxplot
#' @param factors A vector with the names of factor variables contained in targets to color samples in pca plots, pvca and/or arrayQM report
#' @param colorlist Vector of RGB colors to be used for factor variables specified in 'factors' other than Group. If is NULL (default), a default colorlist will be used.
#' @param pca_scale a logical value indicating whether the variables in PCA should be scaled to have unit variance before the analysis takes place.
#' @param batchRemove a logical value indicating whether batch/covariate effects should be removed from the data prior to PCA
#' @param batchFactors Vector with the name of variable(s) from targets dataframe indicating the batches/covariates to be adjusted for
#' @param hc.method The agglomeration method as used by 'hclust()' function. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param hc.numclusters Number of main clusters to highlight by rect
#' @param outputDir Name of the directory where the plot will be saved
#' @param summaryFN Name of the file where the reporting summary of execution results will be saved
#' @param label A string to be included in the file name to be saved
#' @param doboxplot a logical value indicating whether boxplot should be created. Defaults to TRUE.
#' @param dohc a logical value indicating whether heatmap of sample distances should be created. Defaults to TRUE.
#' @param dopca a logical value indicating whether PCA plot should be created. Defaults to TRUE.
#' @param dopvca a logical value indicating whether PVCA plot should be created. Defaults to TRUE.
#' @param arrayQMreport a logical value indicating whether arrayQM report should be created. Defaults to TRUE.
#' @inheritParams qc_boxplot
#' @inheritParams qc_hc
#' @inheritParams qc_pca
#' @inheritParams qc_pvca
#' @export qc_all
#' @return Creates different QC plots into the folder specified.
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link{qc_boxplot}, \link{qc_hc}, \link{qc_pca}, \link{qc_pvca}
#' @examples
#' @keywords Quality control
#' @references

qc_all <- function(data, group, group.color, samplenames, factors, pca_scale, colorlist, batchRemove, batchFactors, hc.method, hc.numclusters=2, label, outputDir, summaryFN, doboxplot=TRUE, dopca=TRUE, dopvca=FALSE, dohc=TRUE, doarrayQMreport=FALSE) {
    targets <- pData(data)
    if (doboxplot) {qc_boxplot(data=data, group=group, group.color=group.color, samplenames=samplenames, outputDir=outputDir, bp.main="Boxplot for array intensity", label, cex.axis=0.5, bp.ylab="log2 intensity", cex.lab=0.7, las=2, bp.legend=TRUE, bp.legend.posx=10, bp.legend.posy=13, bp.legend.cex=0.6)}
    if (dohc) {qc_hc(data=exprs(data), hclust.method=hc.method, names=targets[,samplenames], cexRow = 0.6, cexCol = 0.6, rect=TRUE, numclusters=hc.numclusters, outputDir=outputDir, label=label)}
    if (doarrayQMreport) {
        require(arrayQualityMetrics)
        arrayQualityMetrics(data, outdir = file.path(outputDir, paste0("QCDir.", label)), force=TRUE, intgroup=factors)
        }
    if (dopca) {qc_pca1(data=exprs(data), scale=pca_scale, pca2D_factors=factors, targets=targets, col.group=group.color, colorlist=colorlist, names=targets[,samplenames], outputDir=outputDir, label=label)}
    if (dopvca) {qc_pvca(data=data, factors=factors,label=label, outputDir=resultsDir, summaryFN=summaryFN)}
}
