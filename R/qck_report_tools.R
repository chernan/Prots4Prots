#' Download latticeExtra
#' sudo R CMD INSTALL latticeExtra_0.6-24.tar.gz
library(latticeExtra)

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


