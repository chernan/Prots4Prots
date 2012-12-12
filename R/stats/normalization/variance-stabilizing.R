#' VSN Variance Stabilizing Method
#' 
#' Default settings is lts.quantile=0.9 but lts.quantile=0.5 is more robust. 
#' The reason why lts.quantile=0.5 is not the default is that the estimator 
#' with lts.quantile=0.9 is more efficient (more precise with less data) 
#' if the fraction of differentially expressed genes is not that large.
#' 
#' 
#' R package
#' Article

# For reproducible results using VSN package
set.seed(0)
# Variance-stabilizing normalization
library(vsn)


#' call just VSN ...
#' output normalized data are glogged
apply_vsn09 <- function(data_nn) {

    # as.matrix(cbind(experiment_nn, control_nn))
    
    # normalized data are also "logged" transformed 
    data_normalized <- data.frame(justvsn(as.matrix(data_nn), verbose=FALSE)) #, backgroundsubtract=FALSE

    # or:
    #fit <- vsn2(as.matrix(data_nn), subsample=100000)
    #x <- predict(fit, newdata=as.matrix(data_nn))
    
    return(data_normalized)
}

#' VSN with a lts parameter more robust to big quantities of outliers
apply_vsn05 <- function(data_nn) {
    
    # normalized data are also "logged" transformed 
    data_normalized <- data.frame(justvsn(as.matrix(data_nn), verbose=FALSE, lts.quantile=0.5))
    
    return(data_normalized)
}

#' Apply VSN on provided data.
#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_vsn <- function(dataset, output_folder_temp, output_file, use.robust.fit=FALSE) {
    
    dataset.norm <- matrix(c(''))
    # Apply VSN
    if(use.robust.fit==FALSE) {
        dataset.norm <- apply_vsn09(dataset)
    }
    else {
        dataset.norm <- apply_vsn05(dataset)
    }
    
    # Generate report
    sink(output_file, append=TRUE)

cat('
Normalization
-------------------------

```{r citationVSN, echo=FALSE, warning=FALSE}
vsn_citation <- citation("vsn")
vsn_description <- packageDescription("vsn")
```
Normalisation was achieved using the Variance Stabilizing method (R package version `r vsn_description$Version`).
')

    if(use.robust.fit==FALSE) {
        cat("Default settings (lts.quantile=0.9) were used.\n")
    }
    else {
        cat("Settings more robust to high number of outliers were used (lts.quantile=0.5).\n")
    }

    all_tests_normalization_Rmd(dataset.norm, title="VS Normalization", outfolder=output_folder_temp)
    
    
    tempoutput <- paste(c(output_folder_temp, '/vsn_tests_normalization_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(dataset.norm, tempoutput, sep="\t") 
    cat('',
        'Specific test for VSN.',
        '',
        '```{r, echo=FALSE, fig.width=10, fig.height=6}',
        '',
        sep="\n")
    cat('vsn_test_normalization <- ')
    print(vsn_test_normalization)
    cat( 
        paste(c('matrixdata <- as.matrix(read.table("',tempoutput,'", stringsAsFactors=FALSE))'), collapse=''),
        'vsn_test_normalization(matrixdata)',
        '',
        '```',
        '',
        '',
        sep = "\n"
    )
        
cat('
>
> R package for VSN : published in `r vsn_citation$title` by `r vsn_citation$author`
> 

---------------------------------------------------------------------------
')
    sink(file=NULL)
    
    return(dataset.norm)
}



vsn_test_normalization <- function(matrixdata) {
    meanSdPlot(matrixdata)
}


#' Apply VSN on provided data.
#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
# apply_byDesign_and_report_vsn <- function(dataset, output_folder_images, output_file) {
#     
#     # Apply VSN
#     dataset.norm <- data.frame( apply_vsn09(dataset[,1:4]), apply_vsn09(dataset[,5:8]) )
#     
#     # Generate report
#     sink(output_file, append=TRUE)
#     
#     cat('
# Normalization
# -------------------------
# 
# ```{r citationVSN, echo=FALSE, warning=FALSE}
# vsn_citation <- citation("vsn")
# vsn_description <- packageDescription("vsn")
# ```
# Normalisation was achieved using the Variance Stabilizing method (R package version `r vsn_description$Version`), and was applied by design.
# 
# ')
# 
#     all_tests_normalization_Rmd(dataset.norm, title="VS Normalization", outfolder=output_folder_images)
# 
#     cat('
# >
# > R package for VSN : published in `r vsn_citation$title` by `r vsn_citation$author`
# > 
# 
# ---------------------------------------------------------------------------
# ')
#     sink(file=NULL)
#     
#     return(dataset.norm)
# }
