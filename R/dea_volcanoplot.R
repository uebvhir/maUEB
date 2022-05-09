#' Volcano plot of results from differential expression analysis
#'
#' Creates a volcano plot of the results from differential expression analysis. Dots are colored according to threshold.
#' @param listofcsv List of toptables to summarize
#' @param listofcoef Name of comparisons, in the same order as listofcsv
#' @param volc_logFC Numeric vector of thresholds to be used for logFC (absolute logFC) for each comparison. Genes above this threshold will be colored on the plot.
#' @param volc_pval Character vector indicating the name of the variable in the toptable containing the p-values to be used for thresholding (eg. "P.Value" or "adj.P.Val").
#' @param volc_pval.thr Numeric vector of thresholds to be used for the p-value for each comparison. Genes above this threshold will be colored on the plot.
#' @param volc_x0 Numeric vector of edge x inferior limits for each comparison
#' @param volc_x1 Numeric vector of edge x superior limits for each comparison
#' @param volc_y0 Numeric vector of edge y inferior limits for each comparison
#' @param volc_y1 Numeric vector of edge y superior limits for each comparison
#' @param n Number of plots per sheet. For multiple plots, n=6 is recommended. Default is 1.
#' @param cols Number of columns per sheet. For multiple plots, cols=2 is recommended. Default is 1.
#' @param outputDir Name of the directory where the plot will be saved
#' @param label A string to be included in the file name to be saved
#' @details
#' @import ggplot2 ggrepel
#' @author Mireia Ferrer \email{mireia.ferrer.vhir@@gmail.com}
#' @seealso \link[ggplot2]{ggplot}
#' @examples
#' @return Creates a volcano plot into the graphical device, saves it as pdf.
#' @keywords volcano
#' @references
#' @export

dea_volcanoplot <- function(listofcsv, listofcoef, volc_logFC=rep(1,length(listofcoef)), volc_pval=c(rep("adj.P.Val", length(listofcoef))), volc_pval.thr=rep(0.05, length(listofcoef)), volc_x0=rep(-3, length(listofcoef)), volc_x1=rep(+3, length(listofcoef)), volc_y0=rep(0, length(listofcoef)), volc_y1=rep(10, length(listofcoef)), n=1, cols=1, outputDir, label=""){
    #aquí es miren els punts que sortiran amb nom
    plotes <- list()
    for(i in 1:length(listofcsv)){
        top <- listofcsv[[i]]
        coef <- listofcoef[i]
        #s'afegeix una columna per dir quins gens es pintaran després al volcano
        top$volcThreshold <- as.factor(top[,volc_pval[i]]<volc_pval.thr[i] & abs(top[,"logFC"])>volc_logFC[i])
        g <- ggplot(data = top,
                    aes(x = logFC, y = -log10(eval(parse(text=volc_pval[i]))), colour = volcThreshold)) +
            scale_colour_manual(values = c("TRUE"="purple","FALSE"="lightblue"), breaks="TRUE", labels=paste0(volc_pval[i], "<",volc_pval.thr[i], " & abs(logFC)>", volc_logFC[i]), name="Threshold") +
            geom_point(alpha = 0.4, size = 0.75) +
            xlim(c(volc_x0[i], volc_x1[i])) +
            xlab("log2 fold change") + ylab(paste("-log10 ", volc_pval[i]))
        pl <- g + theme_bw() +
            ggtitle(paste0("Volcano for: ", listofcoef[i])) +
            theme(plot.title = element_text(lineheight = 0.8, face = "bold", hjust = 0.5, size = 10),
                  axis.title.x = element_text(size=8), axis.title.y = element_text(size=8), legend.title = element_blank(),
                  legend.key.size = unit(0.5, "cm"), legend.key = element_rect(fill="grey90", size=0.1), legend.text = element_text(size=8),
                  legend.position="top", legend.margin=margin(0,-2,-2,0), legend.box.margin=margin(-2,-8,-8,-6)) +
            geom_vline(xintercept = c(-volc_logFC[i], volc_logFC[i]), linetype = "dotted") +
            geom_text_repel(data = subset(top, volcThreshold == "TRUE")[1:20,], aes(label = Gene.Symbol),
                            size = 2, vjust = -0.25, hjust = 1.1, segment.size=0.1, col="gray20") +
            scale_y_continuous(trans = "log1p", limits = c(volc_y0[i], volc_y1[i]))
        plotes[[i]] <- pl
    }
    ###Saves volcano plots in pdf file in groups of n plots/sheet
    q <- length(plotes)%/%n #integer part of division: how many groups of 6
    r <- length(plotes)%%n #remainder of division
    ifelse (r>0, r1<-1, r1<-0)
    i1<-seq(1,length(plotes),by=n)
    i2<-c((1:q)*n,q*n+r)
    pdf(file.path(outputDir, paste0("Volcanos", label, ".pdf")))
    for (j in 1:(q+r1)){
        multiplot(plotlist = plotes[i1[j]:i2[j]], cols = cols)
    }
    dev.off()
    #to print plot in device
    for (j in 1:(q+r1)){
        multiplot(plotlist = plotes[i1[j]:i2[j]], cols = cols)
    }
}

#funció per fer disposar gràfics per files i columnes dintre d'un "for" amb ggplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

