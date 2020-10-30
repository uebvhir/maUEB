#' Plot the standard deviations of probesets
#'
#' Computes the standard deviations of all probesets in the array and plots them in increasing order.
#' @param data A matrix of expression values, with samples in columns and probes in rows.
#' @param var_thr Percentile used for thresholding to be plotted as a red line
#' @param label A string to be included in the file name to be saved
#' @param outputDir Name of the directory where the plot will be saved
#' @note For memory reasons, only half of the values are plotted.
#' @export sdplot
#' @return Creates a plot into the graphical device and saves it as pdf in the output directory.
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[stats]{sd}
#' @examples
#' @keywords SD
#' @references

sdplot <- function(data, var_thr, label=NULL, outputDir){
    var <- apply(data,1,sd)
    var <- var[order(var)]
    is.odd <- function(x) x %% 2 != 0
    var2 <- var[is.odd(seq_along(var))]
    pdf(file.path(outputDir,paste0("SDplot", label, ".pdf")))
    plot(seq(1,length(var),2),var2, xlab="features",ylab="Standard deviation",las=2, main=paste0("Standard deviation distribution. ", label))
    abline(v=var_thr*length(var),lty=3, col="red")
    dev.off()
}
