#' @name QCKReportTools
#' @title Report tools for quality check.
#' @description
#' Visualization tools for quality check reports.
#' 
#' @concepts Quality_check
#' 
#' @section Installation of latticeExtra.
#' Download latticeExtra
#' sudo R CMD INSTALL latticeExtra_0.6-24.tar.gz
#' 
#' @import latticeExtra
#' @import ggplot2
library(latticeExtra)
library("ggplot2")

## -----------------------------------------------------------------------------
## Tests / plots used in the report

#' @title Heatmap
#' 
#' @description Heatmap
#' 
#' @details 
#' Display a heatmap including a dendrogram in the legend.
#' 
#' @param dataset A table of intensities.
#' @param distMethod One of the dist() methods to compute distances (base 
#'  package).
#' @param title Title for the heatmap.
displayHeatmap <- function(dataset, distMethod, title="Heatmap") {
    
    ## NB: dist computes distances between rows, so we have to transpose data 
    ## to compute distances between samples
    distances <- as.matrix(dist(t(dataset), method=distMethod))
    
    ## Create dendrogram
    distancesRow <- as.dendrogram(hclust(as.dist(distances)))
    rowOrdered <- order.dendrogram(distancesRow)
    legend <- list(
        top=list(
            fun=dendrogramGrob,
            args=list(x=distancesRow, side="top"))
        )
    
    ## Create heatmap including dendrogram
    p <- levelplot(distances[rowOrdered, rowOrdered],
              scales=list(x=list(rot=90)), 
              xlab='',
              ylab='',
              legend=legend,
              main=title)
    
    return(p)
}

#' @title Linear regression
#' 
#' @description Linear regression
#' 
#' @details 
#' Display a linear regression between two columns of a dataset.
#' 
#' @param dataset A data.frame of intensities.
#' @param indexX Column index for X values.
#' @param indexY Column index for Y values.
#' @param annX Legend position on X.
#' @param annY Legend position on Y.
#' @param show00 Should axis be shown? 
displayLinearReg <- function(dataset, 
                             indexX=1, indexY=2, 
                             annX=10, annY=19, show00=FALSE) {
    
    model <- lm(dataset[,indexY]~dataset[,indexX])
    
    eq <- paste(
        c(
            'y = ', 
            format(coef(model)[1], digits=2),
            ifelse((coef(model)[2]>=0), ' + ', ' - '),
            format(abs(coef(model)[2]), digits=3),
            ' . x'
        ), collapse=''
    )    
    r2val <- paste( c('r2 = ', format(summary(model)$r.squared, digits=3) ), 
                    collapse='')    
    p <- ggplot(data=dataset, 
                aes_string(x=names(dataset)[indexX], y=names(dataset)[indexY]))
    if(isTRUE(show00)) {
        p <- p + geom_point() + 
            geom_hline(yintercept=0, colour="grey") + 
            geom_vline(xintercept=0, colour="grey")
    }
    p <- p + geom_point() + 
        geom_smooth(method=lm, formula=y~x, se=FALSE, color="blue") +
        annotate("text", x=annX, y=annY+0.2, hjust=0, vjust=0, label=eq) + 
        annotate("text", x=annX, y=annY-0.2, hjust=0, vjust=1, label=r2val) 
    
    return(p)
}

