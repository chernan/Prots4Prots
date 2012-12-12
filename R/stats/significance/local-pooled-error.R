#' Local-pooled-error test
#' 
#' !!! CAUTION !!!
#' LPE library is primarily used for analyzing data between two conditions. To use it for
#' paired data, see PLPE library. For using LPE in multiple conditions, use HEM library.
#' 
#' R package
#' http://www.bioconductor.org/packages/release/bioc/html/LPE.html
#' See also
#' http://www.bioconductor.org/packages/2.10/bioc/html/PLPE.html

# For reproducible results using LPE package
set.seed(0)
# Local-pooled-error library
library(LPE)


#' Function apply_lpe_method
#' Apply log2 on values inside the function
apply_lpe_method <- function(experiment, control, bin_val=0.01, rownames) {
    control <- log2(control)
    var.Naive <- baseOlig.error(control, q=bin_val)
    
    experiment <- log2(experiment)
    var.Activated <- baseOlig.error(experiment, q=bin_val)
    
    #compute z-values using LPE method
    lpe.val <- data.frame(lpe(experiment, control, 
                              var.Activated, var.Naive,
                              probe.set.name=rownames)
    )
    
    # round z-values
    lpe.val <- round(lpe.val, digits=2)
    
    # compute p-value from z-values
    p.vals <- 2*pnorm(-abs( lpe.val[,"z.stats"] ))
    lpe.val <- cbind(lpe.val, p.vals)
    
    return(lpe.val)
}

#' Function apply_lpe_method
#' Suppose that data is already log2'ed
apply_lpe_method_sslog2 <- function(experiment, control, bin_val=0.01, rownames) {
    var.Naive <- baseOlig.error(control, q=bin_val)
    var.Activated <- baseOlig.error(experiment, q=bin_val)
    
    #compute z-values using LPE method
    lpe.val <- data.frame(lpe(experiment, control, 
                              var.Activated, var.Naive,
                              probe.set.name=rownames)
    )
    
    # round z-values
#     lpe.val <- round(lpe.val, digits=2)
    
    # compute p-value from z-values
    p.vals <- 2*pnorm(-abs( lpe.val[,"z.stats"] ))
    lpe.val <- cbind(lpe.val, p.vals)
    
    return(lpe.val)
}

#' Needs an lpe.val output of LPE library
order_by_pvalue <- function(lpe.val) {
    order_pval <- order(lpe.val[,"p.vals"], decreasing=FALSE)
    return(lpe.val[order_pval, ])
}

get_lpepval_header <- function() {
    return("p.vals")
}




apply_and_report_lpe <- function(experiment, control, output_folder_temp, output_file, threshold_pval) {
    lpe.val <- apply_lpe_method_sslog2(experiment, control, rownames=rownames(experiment))
    
    matrixdata <- as.matrix(data.frame(
        p.values = lpe.val[,"p.vals"],
        fold.change = lpe.val[,"median.1"] - lpe.val[,"median.2"]
    ))
    
    sink(output_file, append=TRUE)
    cat('
Significantly different proteins
--------------------------------------------------------------------------------------------------------
        
Data analysis was performed using LPE.

')

    all_tests_significance_Rmd(matrixdata, threshold_pval, output_folder_temp)
    
    cat('        
> 
> LPE citation.
>

---------------------------------------------------------------------------
')        
    sink(file=NULL)
    
    return(matrixdata)
}
