#' LOESS
#' 
#' Oshlack et al. show that loess normalization can tolerate up to about 30% asymmetric 
#' differential expression while still giving good results.
#' A. Oshlack, D. Emslie, L. Corcoran, and G. K. Smyth. Normalization of boutique two-
#' color microarrays with a high proportion of differentially expressed probes. Genome
#' Biology, 8:R2, 2007.
#' 
#' R package
#' Article

# For reproducible results
set.seed(0)

# Loess normalization doesn’t affect the A-values.
# Ref: limma help
apply_loess <- function(data_nn) {
    
    data_normalized <- loess(as.matrix(data_nn))
    
    return(data_normalized)
}


#' Apply Normalization on provided data.
#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_limmaNorm <- function(dataset, output_folder_temp, output_file, betweenArrays=FALSE) {
    
    dataset.norm <- apply_loess(dataset)
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Normalization
-------------------------

')

    cat("Normalization was achieved using \"loess\" method available in stats package of base R distribution (`r version$version.string`).")
    
    all_tests_normalization_Rmd(dataset.norm, title="LOESS_Normalization", outfolder=output_folder_temp)

    cat('
>
> Base stats R package (loess function).
> 
> Reference:
> W. S. Cleveland, E. Grosse and W. M. Shyu (1992) 
> Local regression models. Chapter 8 of Statistical Models in S
> eds J.M. Chambers and T.J. Hastie, Wadsworth & Brooks/Cole.

---------------------------------------------------------------------------
')
    sink(file=NULL)
    
    return(dataset.norm)
}

