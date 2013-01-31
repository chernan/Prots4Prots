#' Limma: linear models for microarray data. Normalization functions.
#' (WithinArrays, BetweenArrays)
#' 
#' 
#' R package
#' Article

# For reproducible results
set.seed(0)
# Limma library
library(limma)

# Print-tip loess normalization by default, but we don't have print tip.
# Then loess or robust spline.
# Does a background correction.
# Oshlack et al. show that loess normalization can tolerate up to about 30% asymmetric 
# differential expression while still giving good results.
#
# But object of class list, RGList or MAList containing red and green intensities 
# constituting two-color microarray data.
# Ref: limma help
apply_limma_withinArrays <- function(data_nn) {

    data_normalized <- normalizeWithinArrays(as.matrix(data_nn), method="loess")

    return(data_normalized)
}

# Loess normalization doesn’t affect the A-values.
# Applying quantile normalization to the A-values makes the distributions 
# essentially the same across arrays as well as channels.
# Quantile normalization is noisier.
# method="vsn" gives an interface to the variance-stabilizing norm. methods of the vsn package
# Ref: limma help
apply_limma_betweenArrays <- function(data_nn) {
    
    data_normalized <- normalizeBetweenArrays(as.matrix(data_nn), method="quantile")
    
    return(data_normalized)
}


#' Apply Normalization on provided data.
#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_limmaNorm <- function(dataset, output_folder_temp, output_file, betweenArrays=FALSE) {
    
    dataset.norm <- matrix(c(''))
    # Apply Normalization
    if(betweenArrays==FALSE) {
        dataset.norm <- apply_limma_withinArrays(dataset)
    }
    else {
        dataset.norm <- apply_limma_betweenArrays(dataset)
    }
    
    # Generate report
    sink(output_file, append=TRUE)

cat('
Normalization
-------------------------

```{r citationLIMMAnorm, echo=FALSE, warning=FALSE}
limma_citation <- citation("limma")
limma_description <- packageDescription("limma")
```
')

    if(betweenArrays==FALSE) {
        cat("Normalization was achieved using the normalizeWithinArrays(\"loess\") function of Limma (linear models for microarray data) package (R package version `r limma_description$Version`).
\n")
    }
    else {
        cat("Normalization was achieved using the normalizeBetweenArrays(\"quantile\") function of Limma (linear models for microarray data) package (R package version `r limma_description$Version`).
\n")
    }
    
    
    all_tests_normalization_Rmd(dataset.norm, title="LIMMA Normalization", outfolder=output_folder_temp)
    
cat('
>
> R package for LIMMA : published in `r limma_citation$title` by `r limma_citation$author`
> 
> Methodology paper of normalization functions :
> Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. Methods 31, 265-273
> 

---------------------------------------------------------------------------
')
    sink(file=NULL)
    
    return(dataset.norm)
}


