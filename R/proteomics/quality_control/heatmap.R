#' Download latticeExtra
#' sudo R CMD INSTALL latticeExtra_0.6-24.tar.gz
library(latticeExtra)

#' method : one of the dist methods to compute distances
display_heatmap <- function(data.norm, dist.method, title="Heatmap") {
    
    #NB: dist computes distances between rows, so we have to transpose data to compute distances between samples
    distances <- as.matrix(dist(t(data.norm), method=dist.method))
    
    distances.row <- as.dendrogram(hclust(as.dist(distances)))
    row.ordered <- order.dendrogram(distances.row)
    legend <- list(
        top=list(
            fun=dendrogramGrob,
            args=list(x=distances.row, side="top"))
        )
    levelplot(distances[row.ordered, row.ordered],
              scales=list(x=list(rot=90)), 
              xlab='',
              ylab='',
              legend=legend,
              main=title)
}


#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
report_heatmap <- function(dataset, outfolder, output_file, dist.method="euclidean") {
    
    tempoutput <- paste(c(outfolder, '/report_heatmap_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(dataset, tempoutput, sep="\t") 
    
    library("ggplot2")
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Quality check
-------------------------

',
        paste(c('Checking quality of the samples by computing distances (',dist.method,') between them, and applying a hierarchical clustering.'), collapse='')
    )

cat(
'
```{r, echo=FALSE, fig.width=8, fig.height=8}

',
    'display_heatmap <- ')
    print(display_heatmap)
cat(
    paste(c('matrixdata <- as.matrix(read.table("',tempoutput,'", stringsAsFactors=FALSE))'), collapse=''),
    paste(c('display_heatmap(matrixdata, dist.method="',dist.method,'")'), collapse=''),
    '```',
sep = "\n")
    
cat(
'

---------------------------------------------------------------------------
'
)
    sink(file=NULL)
    
    return(dataset)
}
