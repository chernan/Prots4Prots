#' RMA
#' 
#' 
#' R package
#' Article

# Installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")

# For reproducible results
set.seed(0)
# Library
library("affy")

# Data need to be an ExpressionSet
apply_rma <- function(data_nn) {
    
#     data_normalized <- rma(new("AffyBatch",exprs=as.matrix(data_nn)))
    data_normalized <- vsnrma(data_nn, verbose=FALSE)
    
    return(data_normalized)
}


#' Apply Normalization on provided data (ExpressionSet).
#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
apply_and_report_rma <- function(dataset, output_folder_temp, output_file) {
    
    dataset.norm <- apply_rma(dataset)
    
    # Generate report
    sink(output_file, append=TRUE)
    
    cat('
Normalization
-------------------------

')

    cat("Normalization was achieved using \"RMA\" method.\n")
    
    
    all_tests_normalization_Rmd(dataset.norm, title="RMA Normalization", outfolder=output_folder_temp)
    
    cat('
>

---------------------------------------------------------------------------
')
    sink(file=NULL)
    
    return(dataset.norm)
}

