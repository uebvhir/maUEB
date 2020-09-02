#' A function to read targets file and convert categorical variables to factors
#'
#' Reads a targets file in csv2 format and converts categorical variables to factors
#' @param inputDir Name of the directory containing the targets file to read.
#' @param targetsFN Name of the targets file
#' @param targets.fact Covariables in targets file to be treated as factors
#' @export read_targets
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso read.table
#' @examples
#' @return targets An R object containing the targets file with categorical variables set as factors
#' @keywords read csv
#' @references

read_targets <- function(inputDir, targetsFN, targets.fact){
    targets <- read.table(file.path(inputDir, targetsFN), header = TRUE, row.names = 1, sep = ";", as.is=TRUE)
    for (i in targets.fact){
        #important to set 'levels=unique' so that order of levels is as it appears, and not alphabetic
        targets[,i] <- factor(targets[,i], levels=unique(targets[,i]))
    }
    return(targets)
}
