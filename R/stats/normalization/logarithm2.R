#' Log2 Normalization
#' 
#' R package
#' Article


#' Simple normalisation by applying log2
apply_log2 <- function(data_nn) {
    
    data_normalized <- log(data_nn, base=2)
    
    return(data_normalized)
}

apply_and_report_log2 <- function(dataset, output_folder_temp, output_file) {
    dataset.norm <- apply_log2(dataset)
    
    sink(output_file, append=TRUE)
    cat('
Normalization
-------------------------

Base 2 logarithm was applied on the raw intensities. Note that this transformation is only for comparison with other normalization methods.

')

    all_tests_normalization_Rmd(dataset.norm, title="Log2 normalization", outfolder=output_folder_temp)

cat('

---------------------------------------------------------------------------
')        
    sink(file=NULL)
    
    return(dataset.norm)
}

