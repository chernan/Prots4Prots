#' Inter-Quartile-Normalization
#' 
#' R package
#' Article

# For reproducible results using LPE package
set.seed(0)
# Local-Pooled-Error method package
library(LPE)

#' call IQR from LPE package
#' normalized data are also transformed
apply_iqr_from_lpe1 <- function(data_nn) {
    
    data_normalized <- preprocess(
        # non normalized data
        data_nn, 
        # is MAS5 ok for proteomics data?
        data.type="MAS5", 
        # thresholding value' below which all data would be thresholded
        threshold=1,
        # lowess normalization needs to be performed?
        LOWESS=FALSE
    )
    
    # normalized data are also transformed
    return(data_normalized)
}


apply_and_report_iqr <- function(dataset, output_folder_images, output_file) {
    dataset.norm <- apply_iqr_from_lpe1(dataset)
    
    sink(output_file, append=TRUE)
    cat('
Normalization
-------------------------

```{r citationLPE, echo=FALSE, warning=FALSE}
lpe_citation <- citation("LPE")
```
Normalisation was achieved by using the Inter-Quartile-Range method, as available in the LPE R package (`r lpe_citation$note`).

')

    all_tests_normalization_Rmd(dataset.norm, title="IQR normalization", outfolder=output_folder_images)

cat('
> 
> R package for LPE : published in `r lpe_citation$title` by `r lpe_citation$author`
>

---------------------------------------------------------------------------
        ')
    sink(file=NULL)
    
    return(dataset.norm)
}

# apply_byDesign_and_report_iqr <- function(dataset, output_folder_images, output_file) {
#     dataset.norm <- data.frame( apply_iqr_from_lpe1(dataset[,1:4]), apply_iqr_from_lpe1(dataset[,5:8]) )
#     
#     sink(output_file, append=TRUE)
#     cat('
# Normalization
# -------------------------
# 
# ```{r citationLPE, echo=FALSE, warning=FALSE}
# lpe_citation <- citation("LPE")
# ```
# Normalisation was achieved by using the Inter-Quartile-Range method, as available in the LPE R package (`r lpe_citation$note`), and applying by design.
# 
# ')
#     
#     all_tests_normalization_Rmd(dataset.norm, title="IQR normalization", outfolder=output_folder_images)
#     
#     cat('
# > 
# > R package for LPE : published in `r lpe_citation$title` by `r lpe_citation$author`
# >
# 
# ---------------------------------------------------------------------------
# ')
#     sink(file=NULL)
#     
#     return(dataset.norm)
# }
