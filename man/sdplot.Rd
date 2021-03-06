% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdplot.R
\name{sdplot}
\alias{sdplot}
\title{Plot the standard deviations of probesets}
\usage{
sdplot(data, var_thr, label = NULL, outputDir)
}
\arguments{
\item{data}{A matrix of expression values, with samples in columns and probes in rows.}

\item{var_thr}{Percentile used for thresholding to be plotted as a red line}

\item{label}{A string to be included in the file name to be saved}

\item{outputDir}{Name of the directory where the plot will be saved}
}
\value{
Creates a plot into the graphical device and saves it as pdf in the output directory.
}
\description{
Computes the standard deviations of all probesets in the array and plots them in increasing order.
}
\note{
For memory reasons, only half of the values are plotted.
}
\references{

}
\seealso{
\link[stats]{sd}
}
\author{
Mireia Ferrer Almirall \email{mireia.ferrer.vhir@gmail.com}
}
\keyword{SD}
