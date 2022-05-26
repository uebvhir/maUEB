#' Principal Component Analysis (PCA) from raw or normalized expression data
#'
#' Performs a principal components analysis on the given data matrix and creates a plot of the data projected onto the 2 or 3 first principal components. The plot is saved into a pdf file (2D) or as html (3D).
#' @param data A matrix of expression values, with samples in columns and probes in rows.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
#' @param dim Numeric vector with the number of the principal components to project the PCA. Defaults to first 2. For a 2D plot using PC2 and PC3 should be indicated as c(2,3)
#' @param pca3D Logical, whether a threedimensional plot is to be performed. Defaults to FALSE. If TRUE then dim should be 3 dimensions (default: the first 3 PC)
#' @param pca2D_factors A vector with the names of factor variables contained in targets to color samples in pca 2D plots. Default is "Group".
#' @param pca3D_factors A vector with the names of factor variables contained in targets to color samples in pca 3D plots. Default is "Group".
#' @param targets Targets dataframe with the information of factor variables to plot the pca
#' @param col.group Name of variable in targets containing the colors for Group
#' @param colorlist Vector of RGB colors to be used for factor variables specified in 'factors' other than Group. If is NULL (default), a default colorlist will be used.
#' @param names A vector with sample names corresponding to columns of the data matrix (in the same order)
#' @param batchRemove a logical value indicating whether batch/covariate effects should be removed from the data prior to PCA
#' @param batchFactors Vector with the name of variable(s) from targets dataframe indicating the batches/covariates to be adjusted for
#' @param outputDir Name of the directory where the plot will be saved
#' @param label A label to be included in the file name to be saved
#' @param size Size of labels
#' @param glineas Width of lines rising from labels
#' @import ggplot2 ggrepel
#' @importFrom plotly plot_ly layout as_widget
#' @importFrom magrittr %>%
#' @export qc_pca1
#' @return Creates a PCA plot into the graphical device, saves it as pdf and returns a list with the loads of PCA.
#' @author Mireia Ferrer Almirall \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[stats]{prcomp}, \link[ggplot2]{ggplot}
#' @examples
#' @keywords PCA
#' @references

qc_pca1 <- function(data, scale, dim=c(1,2,3), pca2D_factors=c("Group"), pca3D=FALSE, pca3D_factors=c("Group"), targets, col.group="Colors", colorlist=NULL, names, batchRemove=FALSE, batchFactors=NULL, outputDir, label, size = 1.5, glineas = 0.25){
    if (is.null(colorlist)){colorlist <- rainbow(ncol(data))}
    if (batchRemove){
        if (is.null(batchFactors)){stop("The variable names of batch factors should be provided.")}
        data <- removebatch(data=data, targets=targets, batchFactors=batchFactors)
        label <- paste0(label, "-corrbatch")
    }
    #compute PCA
    data.pca <- prcomp(t(data), scale=scale)
    loads <- round(data.pca$sdev^2/sum(data.pca$sdev^2)*100, 1)
    dims <- paste0("PC", dim)
    if (!pca3D){
        #creates pca plots colored by each factor specified and outputs into a list.
        ## here done using lapply. Can also be done using for loop but need to use local environment, otherwise it overwrites the plots in the list (eg. for (i in 1:length(factors)){listplots[[i]] <- local({i <- i; fact <- factors[i]; ... ; print(pcaplot)})}
        listplots <- lapply(pca2D_factors, function(fact){
            lev <- levels(targets[, fact])
            ifelse (fact=="Group", pcacolors <- unique(targets[, col.group]), pcacolors <- colorlist[1:length(lev)])
            pcatitle <- paste0("Principal Component Analysis for: ", label, ". ", fact)
            pcaplot <- ggplot(data.frame(data.pca$x), aes_string(x=dims[1], y=dims[2])) +
                theme_classic() +
                geom_hline(yintercept = 0, color = "gray70") +
                geom_vline(xintercept = 0, color = "gray70") +
                geom_point(aes(color = targets[,fact]), alpha = 0.55, size = 3) +
                coord_cartesian(xlim = c(min(data.pca$x[,dim[1]])-5,max(data.pca$x[,dim[1]])+5)) +
                # scale_fill_discrete(name = fact) +
                geom_text_repel(aes(y = eval(parse(text=dims[2])) + 0.25, label = names), segment.size = 0.25, size = size) +
                labs(x = c(paste(dims[1],loads[dim[1]],"%")), y=c(paste(dims[2],loads[dim[2]],"%"))) +
                ggtitle(pcatitle) +
                theme(plot.title = element_text(hjust = 0.5)) +
                scale_color_manual(name=fact, values=pcacolors)
            return(pcaplot)
        })
        #save as pdf
        pdf(file.path(outputDir, paste0("PCA", label, ".pdf")))
        for (j in listplots) {print(j)}
        dev.off()
        #print plot in device
        for (j in listplots) {print(j)}
    }
    if (pca3D){
        dataDf <- data.frame(data.pca$x)
        listplots3D <- lapply(pca3D_factors, function(fact3d){
            lev <- levels(targets[, fact3d])
            ifelse (fact3d=="Group", pcacolors <- unique(targets[, col.group]), pcacolors <- colorlist[1:length(lev)])
            pcatitle <- paste0("Principal Component Analysis for: ", label, ". ", fact3d)
            p <- plot_ly(dataDf, type = "scatter3d", x = dataDf[,dims[1]] , y = dataDf[,dims[2]], z= dataDf[,dims[3]],
                    text = rownames(dataDf), mode = "markers+text", showlegend =TRUE,
                    color = ~targets[,fact3d], colors = pcacolors) %>%
                layout(title=pcatitle,
                       scene = list(xaxis = list(title = paste(dims[1], loads[1], "%")),
                                    yaxis = list(title = paste(dims[2], loads[2], "%")),
                                    zaxis = list(title = paste(dims[3], loads[3], "%"))))
            #save plot as html file into outputdir
            htmlwidgets::saveWidget(as_widget(p), file.path(resultsDir, paste0("PCA3D", label, ".", fact3d, ".html")))
            return(p)
        })

    }
    return(loads)
}
