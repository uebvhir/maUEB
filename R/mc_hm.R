#' Heatmap plots
#'
#' Plot heatmaps between the comparisons determined, with the thresholds values (p-value and logFC) of interest.
#'
#' @param listofcsv List of toptables
#' @param hm_dataset Dataset with expression values to plot on heatmap
#' @param targets Targets file
#' @param featureCol Column in targets file corresponding to features (eg. "Gene.Symbol")
#' @param hm_comparNames name of set of comparisons (without spaces/dots/_). Usually the same as used for Venn (eg. venn_comparNames)
#' @param hm_compar list with set of comparisons to select genes to do the Heatmap, positions (or names) according to 'listofcoef' or 'listofTopNames'. Usually the same as used for Venn (eg. venn_compar)
#' @param hm_groupexclude samples to exclude in each set of comparisons for Heatmap (length must be equal to number of comparisons). If no sample is to be excluded put c()
#' @param hm_pval "P.Value" or "adj.P.Val" (for each set of comparisons)
#' @param hm_pval.thr pvalue threshold (for each set of comparisons)
#' @param hm_logFC Name of column corresponding to logFC. Default is "logFC"
#' @param hm_logFC.thr logFC threshold (for each set of comparisons)
#' @param hm_palette colors used in the heatmap plots
#' @param hm_clustCol whether to perform hierarchical clustering of samples in  heatmap (column dendogram) (eg. TRUE/FALSE)
#' @param hm_clustCol.dist required if hm_clustCol=TRUE: method for distance calculation in hierarchical clustering of samples (default in heatmap.2 is "euclidian")
#' @param hm_clustCol.method required if hm_clustCol=TRUE: method for hierarchical clustering (default in heatmap.2 is "complete", other may be "average")
#' @param hm_clustRow.cor whether to cluster gene expression based on distance correlation (otherwise distance matrix is calculated based on euclidian distance (default in heatmap.2 function))
#' @param batcheffect whether to correct expression values for batch effect (eg. TRUE/FALSE)
#' @param batchcolName column name in 'targets' file containing Batch (eg. 'Batch')
#' @param batchcolName2 column name in 'targets' file containing categorical cofactor to remove from expression data (eg. 'Batch2')
#' @param batchNumcolName Vector of numerical variables contained in 'targets' file to be taken as covariates in removeBatcheffect function (eg. 'Day')
#' @param hm_fontsize
#' @param hm_fontsize_row
#' @param hm_fontsize_col
#' @param outputDir directory in which the results are going to be saved
#' @details
#' @import limma
#' @import pheatmap
#' @import heatmaply
#' @import plotly
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @examples
#' @return Heatmap plots between the comparisons determined, with the thresholds values (p-value and logFC) of interest. Plots are saved in pdf format in the outputDir.
#' @keywords Heatmap plots
#' @references
#' @export

mc_hm <- function(listofcsv, hm_dataset, targets, featureCol="Gene.Symbol", hm_comparNames, hm_compar, hm_groupexclude=rep(list(c()), length(hm_comparNames)), hm_pval, hm_pval.thr, hm_logFC="logFC", hm_logFC.thr, hm_palette=colorRampPalette(c("blue", "white", "red"))(n = 199), hm_clustCol=TRUE, hm_clustCol.dist="euclidian", hm_clustCol.method="complete", hm_clustRow.cor=TRUE, batcheffect=FALSE,
                  batchcolName=NULL, batchcolName2=NULL, batchNumcolName=NULL, hm_fontsize_row=rep(3, length(hm_comparNames)), hm_fontsize_col=rep(7, length(hm_comparNames)), hm_fontsize= rep(8, length(hm_comparNames)), hm_plots_interactive=FALSE, resultsSummFN, outputDir){
        if (batcheffect){
            if(is.null(batchcolName)) batch <- NULL else batch <- as.factor(targets[,batchcolName])
            if(is.null(batchcolName2)) batch2 <- NULL else batch2 <- as.factor(targets[,batchcolName2])
            if(is.null(batchNumcolName)) batch.num <- NULL else batch.num <- targets[,batchNumcolName] #batchNumcolName can be a vector (multiple columns)
            designmat <- model.matrix(~0+Group, targets)
            hm_dataset[,as.character(targets$ShortName)] <- removeBatchEffect(hm_dataset[,as.character(targets$ShortName)], batch=batch, batch2=batch2, covariates=batch.num, design=designmat) #function from limma package
        }
    list_hm <- list()
        for (h in 1:length(hm_compar)){
            whichcomp <- hm_compar[[h]]
            listofcsv_hm <- listofcsv[whichcomp]
            ##select subset of DEGs for heatmap according to logFC/pvalue thresholds (normally 100-200 DEGs for good visualization)
            hm_genes <- c()
            for (k in 1:length(listofcsv_hm)){
                hm_csvtop <- listofcsv_hm[[k]]
                hm_genesel <- hm_csvtop[(hm_csvtop[,hm_pval[h]] < hm_pval.thr[h] & abs(hm_csvtop[, hm_logFC]) > hm_logFC.thr[h]), featureCol]
                ##set a list with all genes from the subset of comparisons
                hm_genes <- c(hm_genes, as.vector(hm_genesel))
            }
            length(hm_genes)
            ##filter out duplicate genes
            hm_genes.unique <- unique(hm_genes)
            length(hm_genes.unique)
            ##take expression values of selected genes
            Expr.genes2hm <- hm_dataset[which(hm_dataset[,featureCol] %in% hm_genes.unique),]
            colnames(Expr.genes2hm)
            dim(Expr.genes2hm)
            ##select samples to include in heatmap
            if(length(hm_groupexclude) > 0 & length(hm_groupexclude[[h]]) > 0) {
                samples.selected <- targets[!(targets$Group %in% hm_groupexclude[[h]]),"ShortName"]
            } else {
                samples.selected <- targets[, "ShortName"]
            }
            expr.selected <- Expr.genes2hm[,c(featureCol, samples.selected)] #AQUí, aquest objecte no es crea bé. penso que el problema pot estar en la X de la design matrix davant del nom.
            rownames(expr.selected) <- expr.selected[, featureCol] #set gene symbols as rownames
            expr.selected <- as.matrix(expr.selected[,-1]) #set dataframe as matrix
            dim(expr.selected) #must be equal to the number of samples
            #calculem el dendograma per les columnes (samples) si així s'ha especificat a paràmetres
            if (hm_clustCol) {
                clustCol <- hclust(dist(t(expr.selected), method=hm_clustCol.dist), method = hm_clustCol.method)
                filename.clustcol <-""
                dendroCol <- as.dendrogram(clustCol)
            } else {
                clustCol <- FALSE
                filename.clustcol <- ".noCol"
                dendroCol <- FALSE
            }
            #colors pels heatmaps corresponent a les mostres incloses
            targets1 <- targets
            rownames(targets1) <- targets1$ShortName
            targets1 <- targets1[colnames(expr.selected),]
            hm_samplecolors <- unique(targets1$Colors)
            names(hm_samplecolors) <- unique(targets1$Group)
            hm_samplecolors <- list(Group=hm_samplecolors)
            #calculem el dendograma per les files (gens) en base a distancia de correlacio (en el basicpipe ho fa aixi) si així s'ha indicat a parametres
            if (hm_clustRow.cor) {
                clustRow <- hclust(as.dist(1 - cor(t(expr.selected))),  method = "average")
                filename.cor <- ".corRow"
                dendroRow <- as.dendrogram(clustRow)
            } else {clustRow <- TRUE;dendroRow <- TRUE; filename.cor <- ""} #sino deixa dendrorow=TRUE pq utilitzi el clustering com ve per defecte a la funcio heatmap.2
            #save summary of heatmap results
            sink(file.path(outputDir, resultsSummFN), split = FALSE, append = TRUE)
            cat("\n-----------------------Heatmap settings---------------------------\n")
            cat("-Sets for comparison:", paste(hm_comparNames, collapse=", "))
            cat("\n-pvalue type:", paste(hm_pval, collapse=", "))
            cat("\n-pvalue thr:", paste(hm_pval.thr, collapse=", "))
            cat("\n-logFC:", paste(hm_logFC.thr, collapse=", "))
            cat("\n-Column clustering:", hm_clustCol, "; dist=", hm_clustCol.dist, "; method=", hm_clustCol.method)
            cat("\n-Row clustering by correlation distance:", hm_clustRow.cor)
            cat("\n")
            sink()
            ##CREATE HEATMAPS with pheatmap
            ###plot
            hm <- pheatmap(expr.selected,
                        cluster_rows = clustRow,
                        cluster_cols = clustCol, #determines if samples are reordered according to hierarchical clustering
                        main = paste0("Heatmap for ", hm_comparNames[h], ":\n", hm_pval[h], " < ", hm_pval.thr[h], " & abs(logFC) > ", hm_logFC.thr[h],
                                     " (", nrow(expr.selected), " genes)"),
                        scale = "row",
                        col = hm_palette,
                        border_color=NA,
                        fontsize=hm_fontsize[h],
                        fontsize_row = hm_fontsize_row[h],
                        fontsize_col = hm_fontsize_col[h],
                        annotation_col = targets1["Group"],
                        angle_col="45",
                        annotation_colors = hm_samplecolors)
            pdf(file.path(outputDir, paste0("Heatmap", hm_comparNames[h], filename.clustcol, ".", hm_pval[h], hm_pval.thr[h], ".", hm_logFC, hm_logFC.thr[h], ".pdf")))
            print(hm)
            dev.off()
            list_hm[[h]] <- hm

            ##CREATE INTERACTIVE HEATMAPS
            if (hm_plots_interactive) {
                heatmaply(expr.selected, main = paste0("HeatMap for ", hm_comparNames[h], ":\n", hm_pval[h], " < ", hm_pval.thr[h], " & abs(logFC) > ", hm_logFC.thr[h],
                                                   " (", nrow(expr.selected), " genes)"),
                        margins=c(50,50,80,50),
                        Colv = dendroCol,
                        Rowv = dendroRow,
                        scale = "row",
                        col = hm_palette,
                        sepcolor = "white",
                        sepwidth = c(0.05,0.05),
                        cexRow = 0.5,
                        cexCol = 0.5,
                        key = TRUE,
                        keysize = 2, #tamaño de la "key" si es pequeña no salen los colores
                        density.info = "histogram",
                        tracecol = NULL,
                        srtCol = 30 ,
                        ColSideColors = targets1$Colors,
                        revC=TRUE,
                        width=600,
                        height=800,
                        branches_lwd=0.3,
                        file=paste0(outputDir, "/Heatmap", hm_comparNames[h], filename.clustcol, ".", hm_pval[h], hm_pval.thr[h], ".", hm_logFC, hm_logFC.thr[h], ".html"))
            }
        }
    names(list_hm) <- hm_comparNames
    return(list_hm)
    }


