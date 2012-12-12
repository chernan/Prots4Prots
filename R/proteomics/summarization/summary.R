library(plyr)


#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_NA <- function(dataset, dataset.norm, outfolder, output_file) {
    #no summary, just to test pipeline
    data_proteins <- new("ExpressionSet", 
                         exprs = dataset.norm, 
                         phenoData = phenoData(dataset),
                         featureData = featureData(dataset)
    )
    
    tempoutput <- paste(c(outfolder, '/report_summarypep_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.exprs(data_proteins, tempoutput, sep="\t") 
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Summary
-------------------------

No summary method was applied.
')
    
    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
')
    
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('cols.exp <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Experiment")]),collapse='", "'), '")'), collapse=''),
        paste(c('cols.ctl <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Control")]),collapse='", "'), '")'), collapse=''),
        'data.to.test1 <- apply(matrixdata[,cols.exp], 1, mean)',
        'data.to.test2 <- apply(matrixdata[,cols.ctl], 1, mean)',
        'display_normalization_MAplot(data.to.test1,data.to.test2)',
        sep = "\n"
    )
    
    cat('
```
')
    
cat('

---------------------------------------------------------------------------
'
    )
    sink(file=NULL)
    
    return(data_proteins)
}


#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_mean <- function(dataset, dataset.norm, outfolder, output_file) {
    #mean by prot
    aggr.dataset <- aggregate(dataset.norm, 
              by = list(Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
              mean)
    
    data_proteins <- new("ExpressionSet", 
                         exprs = aggr.dataset[,names(dataset.norm)], 
                         phenoData = phenoData(dataset),
                         featureData = new("AnnotatedDataFrame",data.frame(aggr.dataset[,1]))
    )
    
    tempoutput <- paste(c(outfolder, '/report_summarypep_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(aggr.dataset, tempoutput, sep="\t") 
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Summary
-------------------------

Mean by protein group.
')
    
    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
')
    
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('cols.exp <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Experiment")]),collapse='", "'), '")'), collapse=''),
        paste(c('cols.ctl <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Control")]),collapse='", "'), '")'), collapse=''),
        'data.to.test1 <- apply(matrixdata[,cols.exp], 1, mean)',
        'data.to.test2 <- apply(matrixdata[,cols.ctl], 1, mean)',
        'display_normalization_MAplot(data.to.test1,data.to.test2)',
        sep = "\n"
    )
    
    cat('
```
')
    
cat('

---------------------------------------------------------------------------
'
)
    sink(file=NULL)
    
    return(data_proteins)
}



#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_median <- function(dataset, dataset.norm, outfolder, output_file) {
    #median by prot
    aggr.dataset <- aggregate(dataset.norm, 
                              by = list(Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
                              median)
    
    data_proteins <- new("ExpressionSet", 
                         exprs = aggr.dataset[,names(dataset.norm)], 
                         phenoData = phenoData(dataset),
                         featureData = new("AnnotatedDataFrame",data.frame(aggr.dataset[,1]))
    )
    
    tempoutput <- paste(c(outfolder, '/report_summarypep_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(aggr.dataset, tempoutput, sep="\t") 
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Summary
-------------------------

Median by protein group.
')
    
    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
')
    
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('cols.exp <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Experiment")]),collapse='", "'), '")'), collapse=''),
        paste(c('cols.ctl <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Control")]),collapse='", "'), '")'), collapse=''),
        'data.to.test1 <- apply(matrixdata[,cols.exp], 1, mean)',
        'data.to.test2 <- apply(matrixdata[,cols.ctl], 1, mean)',
        'display_normalization_MAplot(data.to.test1,data.to.test2)',
        sep = "\n"
    )
    
    cat('
```
')
    
cat('
    
---------------------------------------------------------------------------
'
    )
    sink(file=NULL)
    
    return(data_proteins)
}



#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_trimmed20 <- function(dataset, dataset.norm, outfolder, output_file) {
    #trimmed mean by prot
    aggr.dataset <- aggregate(dataset.norm, 
                              by = list(Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
                              FUN=function(x){mean(x,trim=0.1)}
    )
    
    data_proteins <- new("ExpressionSet", 
                         exprs = aggr.dataset[,names(dataset.norm)], 
                         phenoData = phenoData(dataset),
                         featureData = new("AnnotatedDataFrame",data.frame(aggr.dataset[,1]))
    )
    
    tempoutput <- paste(c(outfolder, '/report_summarypep_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(aggr.dataset, tempoutput, sep="\t") 
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Summary
-------------------------

Trimmed mean by protein group (20%).
')
    
    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
')
    
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('cols.exp <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Experiment")]),collapse='", "'), '")'), collapse=''),
        paste(c('cols.ctl <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Control")]),collapse='", "'), '")'), collapse=''),
        'data.to.test1 <- apply(matrixdata[,cols.exp], 1, mean)',
        'data.to.test2 <- apply(matrixdata[,cols.ctl], 1, mean)',
        'display_normalization_MAplot(data.to.test1,data.to.test2)',
        sep = "\n"
    )
    
    cat('
 ```
')
    
cat('
    
---------------------------------------------------------------------------
'
    )
    sink(file=NULL)
    
    return(data_proteins)
}



#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_mostvariablepep <- function(dataset, dataset.norm, outfolder, output_file) {
    #Look 
    temp.df <- data.frame(dataset.norm, "Majority.protein.IDs"=featureData(dataset)[["Majority.protein.IDs"]])
    aggr.dataset <- ddply(.data=temp.df, 
                          .variables="Majority.protein.IDs", 
                          .fun=function(df) {
                                variances <- apply(df[,names(dataset.norm)], 1, var)
                                max.var <- max(variances)
                                return(df[(which(variances == max.var))[1],])
                          }
    )

    data_proteins <- new("ExpressionSet", 
                         exprs = data.frame(aggr.dataset[,names(dataset.norm)]), 
                         phenoData = phenoData(dataset),
                         featureData = new("AnnotatedDataFrame",data.frame("Majority.protein.IDs"=aggr.dataset[["Majority.protein.IDs"]]))
    )
    
#     aggr.dataset <- by(dataset.norm, 
#                        INDICES=featureData(dataset)[["Majority.protein.IDs"]],
#                        FUN=function(df) {
#                            variances <- apply(df, 1, var)
#                            max.var <- max(variances)
#                            return(df[(which(variances == max.var))[1],])
#                        })
    
#     data_proteins <- new("ExpressionSet", 
#                          exprs = data.frame(as.matrix(t(sapply(names(aggr.dataset), FUN=function(name){ return(aggr.dataset[[name]]) } )))), 
#                          phenoData = phenoData(dataset),
#                          featureData = new("AnnotatedDataFrame",data.frame("Majority.protein.IDs"=names(aggr.dataset),row.names=names(aggr.dataset)))
#     )
    
    tempoutput <- paste(c(outfolder, '/report_summarypep_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(aggr.dataset, tempoutput, sep="\t") 
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Summary
-------------------------

Most variable peptide for each protein.
')
    
    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
')
    
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('cols.exp <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Experiment")]),collapse='", "'), '")'), collapse=''),
        paste(c('cols.ctl <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Control")]),collapse='", "'), '")'), collapse=''),
        'data.to.test1 <- apply(matrixdata[,cols.exp], 1, mean)',
        'data.to.test2 <- apply(matrixdata[,cols.ctl], 1, mean)',
        'display_normalization_MAplot(data.to.test1,data.to.test2)',
        sep = "\n"
    )
    
    cat('
```
')
    
cat('
    
---------------------------------------------------------------------------
'
    )
    sink(file=NULL)
    
    return(data_proteins)
}

#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_mostintensepep <- function(dataset, dataset.norm, outfolder, output_file) {

    temp.df <- data.frame(dataset.norm, "Majority.protein.IDs"=featureData(dataset)[["Majority.protein.IDs"]])
    aggr.dataset <- ddply(.data=temp.df, 
                          .variables="Majority.protein.IDs", 
                          .fun=function(df) {
                              sums <- apply(df[,names(dataset.norm)], 1, sum)
                              max.sum <- max(sums)
                              return(df[(which(sums == max.sum))[1],])
                          }
    )
    
    data_proteins <- new("ExpressionSet", 
                         exprs = data.frame(aggr.dataset[,names(dataset.norm)]), 
                         phenoData = phenoData(dataset),
                         featureData = new("AnnotatedDataFrame",data.frame("Majority.protein.IDs"=aggr.dataset[["Majority.protein.IDs"]]))
    )

    tempoutput <- paste(c(outfolder, '/report_summarypep_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(aggr.dataset, tempoutput, sep="\t") 
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Summary
-------------------------

Most insense peptide for each protein.
')

    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
')
    
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('cols.exp <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Experiment")]),collapse='", "'), '")'), collapse=''),
        paste(c('cols.ctl <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Control")]),collapse='", "'), '")'), collapse=''),
        'data.to.test1 <- apply(matrixdata[,cols.exp], 1, mean)',
        'data.to.test2 <- apply(matrixdata[,cols.ctl], 1, mean)',
        'display_normalization_MAplot(data.to.test1,data.to.test2)',
        sep = "\n"
    )
    
    cat('
```
        ')

cat('
---------------------------------------------------------------------------
'
    )
    sink(file=NULL)
    
    return(data_proteins)
}

#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_sumintens <- function(dataset, dataset.norm, outfolder, output_file) {
    #mean by prot
    aggr.dataset <- aggregate(dataset.norm, 
                              by = list(Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
                              sum)
    
    data_proteins <- new("ExpressionSet", 
                         exprs = aggr.dataset[,names(dataset.norm)], 
                         phenoData = phenoData(dataset),
                         featureData = new("AnnotatedDataFrame",data.frame(aggr.dataset[,1]))
    )
    
    tempoutput <- paste(c(outfolder, '/report_summarypep_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(aggr.dataset, tempoutput, sep="\t") 
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Summary
        -------------------------
        
        Sum intensities of peptides by protein group.
        ')
    
    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
        ')
    
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('cols.exp <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Experiment")]),collapse='", "'), '")'), collapse=''),
        paste(c('cols.ctl <- c("',paste(sampleNames(data_proteins[,(data_proteins$Design == "Control")]),collapse='", "'), '")'), collapse=''),
        'data.to.test1 <- apply(matrixdata[,cols.exp], 1, mean)',
        'data.to.test2 <- apply(matrixdata[,cols.ctl], 1, mean)',
        'display_normalization_MAplot(data.to.test1,data.to.test2)',
        sep = "\n"
    )
    
    cat('
```
        ')
    
    cat('

        ---------------------------------------------------------------------------
        '
)
    sink(file=NULL)
    
    return(data_proteins)
}
