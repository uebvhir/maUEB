#' Venn diagrams and Upset plots
#'
#' Plot Venn diagramm from 2 to 5 comparisons or UpsetR plot if >5 comparisons. Function based in previous work of Ferran Briansó and Miriam Mota ("https://raw.githubusercontent.com/rgonzalo/GitBackupScripts/master/VennDiagrams/CreateVennEuler.R")

#' @param listofcsv list of the dataframes with the toptables
#' @param namescomp list with the names of each comparisons
#' @param label name that is going to appear in the file
#' @param colFeat column by which common genes are searched
#' @param colPval column to do the selection (normally Pvalue)
#' @param pval  p-value by which the selection is done
#' @param pltR to show or not the image in R
#' @param pltPdf to make or not the pdf
#' @param eul to draw or not euler plot
#' @param venn to perform or not the venn diagrams
#' @param upsetPlot to perform or not the upset plots. Normally when there are more than 5 comparisons or when this kind of plot is preferred
#' @param csv create the csv file with the common genes
#' @param colors colors to use for the plots (rainbow by defect)
#' @param trans transparency of colors
#' @param cex1 size of the labels
#' @param cex2 size of the numbers
#' @param rotation venn global rotation
#' @param position vector with the position in grades of each comparison label
#' @param FC logFC value to select genes
#' @param colFC column referring to logFC (eg. "logFC")
#' @param include direction in which genes are going to be selected (up / down / abs)
#' @param outputDir directory in which the results are going to be saved
#' @details
#' @import VennDiagram writexl grid
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com} and Esther Camacho
#' @examples
#' @return Plot Venn diagramm from 2 to 5 comparisons or UpsetR plot if >5 comparisons
#' @keywords Venn diagrams Upset plots
#' @references
#' @export

mc_venn_upset <- function(listofcsv, namescomp, label, colFeat, colPval, pval=0.05, colFC, FC=1, pltR=TRUE,
                          pltPdf=TRUE, pltPng=FALSE, venn=TRUE, eul=FALSE, saveTables=TRUE, upsetPlot=FALSE, saveTables_extended=TRUE,
                          colors=rainbow(length(namescomp)), trans=0.5, cex1=1, rotation=0, position=rep(0,length(namescomp)), cex2=1, include=rep("abs", length(namescomp)), margins_venn=c(5.1,4.1,4.1,2.1), resultsSummFN, outputDir){
    ## Initializing lists
    list_genes_sel <- list()
    topTabs <- listofcsv[namescomp]
    ## Reading input data
    for (i in 1:length(topTabs)) {
        if (include=="abs"){list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colPval] < pval & abs(topTabs[[i]][, colFC]) > FC])}
        if (include=="up"){list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colPval] < pval & topTabs[[i]][, colFC] > FC])}
        if (include=="down"){list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colPval] < pval & topTabs[[i]][, colFC] < FC])}
    }
    numgenes.sel <- sapply(list_genes_sel, length)
    namescomp_plot <- paste0(namescomp, paste0("\n(",numgenes.sel,")"))
    ## Creating Venn Diagram

    if (venn) {
        opt <- par(mar=margins_venn)
        if (include=="abs") {titulo <- paste0("Venn diagram for comparison: ", label, "\n(" , colPval, " < ", pval," & abs(logFC) > ", FC, ")")}
        if (include=="up") {titulo <- paste0("Venn diagram for comparison: ", label, "\nUp-regulated genes (" , colPval, " < ", pval," & ", colFC, " > ", FC, ")")}
        if (include=="down") {titulo <- paste0("Venn diagram for comparison: ", label, "\nDown-regulated genes (" , colPval, " < ", pval," & ", colFC, " < ", FC, ")")}
        venn.plot <- venn.diagram(list_genes_sel,
                                  category.names = namescomp_plot,
                                  fill = colors,
                                  alpha = trans,
                                  resolution = 800,
                                  cat.cex = cex1,
                                  cex = cex2,
                                  main.cex=cex2,
                                  main = titulo,
                                  filename = NULL,
                                  rotation.degree = rotation,
                                  cat.pos = position)
        if (pltPdf) {
            pdf(file.path(outputDir, paste0("VennDiagram.", include, ".",label, ".", colPval, pval, ".",colFC,FC,".pdf")))
            grid.draw(venn.plot)
            dev.off()
        }
        if (pltPng) {
            png(file.path(paste0("VennDiagram.", include, ".",label, ".", colPval, pval, ".", colFC,FC,".png")))
            grid.draw(venn.plot)
            dev.off()
        }
        if (pltR) {grid.draw(venn.plot)}

        #save summary of venn results
        sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
        cat("\n--------------------Venn settings:--------------------\n")
        cat("-Sets for comparison:", paste(venn_comparNames, collapse=", "))
        cat("\n-pvalue type:", paste(colPval, collapse=", "))
        cat("\n-pvalue thr:", paste(pval, collapse=", "))
        cat("\n-logFC:", paste(FC, collapse=", "))
        cat("\n")
        sink()
        par(opt)
    }

    ## Creating Euler Diagram
    if (eul) {
        set <- NULL
        for (i in 1:length(namescomp_plot)) {
            set <- c(set, rep(namescomp_plot[i],length(list_genes_sel[[i]])))
        }
        v <- venneuler(data.frame(elements = c(unlist(list_genes_sel)),
                                  sets = set))
        if (pltPdf) {
            pdf(file.path(outputDir, paste0("EulerDiagram", include, ".", label, ".", colPval, pval, ".pdf")))
            plot(v, main = paste0("Euler diagram (" , colPval, " < ", pval,")"))
            dev.off()
        }
        if (pltR) {plot(v, main = paste0("Euler diagram (" , colPval, " < ", pval,")"))}
    }

    ## Creating UpsetR plot
    if (upsetPlot) {
        if (include=="abs") {titulo <- paste0("UpSet plot for: ", label, "\n(" , colPval, " < ", pval," & abs(logFC) > ", FC, ")")}
        if (include=="up") {titulo <- paste0("UpSet plot for: ", label, "\nUp-regulated genes (" , colPval, " < ", pval," & logFC > ", FC, ")")}
        if (include=="down") {titulo <- paste0("UpSet plot for: ", label, "\nDown-regulated genes (" , colPval, " < ", pval," & logFC < -", FC, ")")}
        names(list_genes_sel) <- namescomp
        #això ens prepara l'objecte per fer el gràfic. la llista que surt no te rownames.
        data2upset <- fromList(list_genes_sel)
        plotupset <- upset(data2upset, sets=namescomp, keep.order=TRUE,order.by = "freq", nsets = length(whichcomp),  mainbar.y.label = "Genes in intersection",
                           sets.x.label = "Genes in comparison", main.bar.color = "#5536AB",
                           matrix.color = "#5536AB", sets.bar.color = "#5536AB")
        if (pltPdf) {
            pdf(file.path(outputDir, paste0("UpSetPlot", include, ".",label, ".", colPval, pval, ".", colFC,FC,".pdf")))
            print(plotupset)
            grid.text(titulo,x = 0.7, y=0.95, gp=gpar(fontsize=11))
            dev.off()
        }
        if (pltPng) {
            pdf(file.path(outputDir, paste0("UpSetPlot", include, ".",label, ".", colPval, pval, ".", colFC,FC,".png")))
            print(plotupset)
            grid.text(titulo,x = 0.7, y=0.95, gp=gpar(fontsize=11))
            dev.off()
        }
        if (pltR) {
            print(plotupset)
            grid.text(titulo,x = 0.7, y=0.95, gp=gpar(fontsize=11))
        }

    }
    ## Obtaining lists of combined elements and computing shared elements

    names(list_genes_sel) <- paste0(namescomp, "_")
    combs <-  unlist(lapply(1:length(list_genes_sel),
                            function(j) combn(names(list_genes_sel), j, simplify = FALSE)),
                     recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    #str(combs)

    elements <- lapply(combs, function(i) Setdiff(list_genes_sel[i], list_genes_sel[setdiff(names(list_genes_sel), i)]))
    n.elements <- sapply(elements, length)
    list_res <- list(elements = elements, n.elements = n.elements)

    seq.max <- seq_len(max(n.elements))
    mat <- sapply(elements, "[", i = seq.max)
    mat[is.na(mat)] <- ""
    sharedElements <- rbind(t(data.frame(n.elements)),data.frame(mat))

    ## Writing table of shared elements as csv
    if (saveTables) {
        writexl::write_xlsx(sharedElements, file.path(outputDir, paste0("sharedElements", include, ".", label, ".", colPval, pval,".",colFC,".",FC,".xlsx")))
    }
    if (saveTables_extended) {
      sharedElements.stack <- stack(sharedElements[-1,])
      sharedElements.stack <- sharedElements.stack[sharedElements.stack$values!="",]
      for (comparison in namescomp){
        listofcsv_comp <- topTabs[[comparison]]
        sharedElements.stack <- merge(sharedElements.stack, listofcsv_comp[,c(colFeat, colFC, colPval)], by=1, all.x=TRUE)
      }
      colnames(sharedElements.stack)[1:2] <- c(colFeat, "Venn.subset")
      colnames(sharedElements.stack)[grep(paste0("^",colFC), colnames(sharedElements.stack))] <- paste0(colFC,".", comparison)
      colnames(sharedElements.stack)[grep(paste0("^",colPval), colnames(sharedElements.stack))] <- paste0(colPval,".", comparison)
      if ((nrow(sharedElements.stack)<=65535) & (ncol(sharedElements.stack)<=256)){
        writexl::write_xlsx(sharedElements.stack,
                 file.path(outputDir, paste0("sharedElements.", include, ".", label, ".", colPval, pval,".",colFC,".",FC,".extended.xlsx")))
      } else {
        write.csv2(sharedElements.stack, file.path(outputDir, paste0("sharedElements.", include, ".", label, ".", colPval, pval,".",colFC,".",FC,".extended.csv")))
      }
    }
    listtables <- list(sharedElements=sharedElements, sharedElements.extended=sharedElements.stack)
    ## Returning table as data.frame
    return(listtables)
}
#####################################################################
###################################################################
###### Internal functions used to extract the lists of shared elements
###### both x and y are, in all cases, lists of character strings
###################################################################
Intersect <- function(x) {
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        intersect(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        intersect(x[[1]], Intersect(x[-1]))
    }
}

Union <- function(x) {
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        union(x[[1]], Union(x[-1]))
    }
}

Setdiff <- function(x, y) {
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
}
