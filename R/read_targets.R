#' Read targets file and convert categorical variables to factors
#'
#' Reads a targets file in csv2 format and converts the categorical variables defined in parameter targets.fact to factors
#' @param inputDir Name of the directory containing the targets file to read.
#' @param targetsFN Name of the targets file to read (csv file with ; as field separator)
#' @param targets.fact Covariables in targets file to be treated as factors
#' @export read_targets
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso read.table
#' @examples
#' @return Returns a targets dataframe with the ShortName set as rownames and the categorical variables specified by targets.fact set as factors
#' @keywords read csv
#' @references

read_targets <- function(inputDir, targetsFN, targets.fact){
    targets <- read.table(file.path(inputDir, targetsFN), header = TRUE, row.names = NULL, sep = ";", as.is=TRUE)
    rownames(targets) <- as.character(targets$ShortName)
    for (i in targets.fact){
        #important to set 'levels=unique' so that order of levels is as it appears, and not alphabetic
        targets[,i] <- factor(targets[,i], levels=unique(targets[,i]))
    }
    return(targets)
}
