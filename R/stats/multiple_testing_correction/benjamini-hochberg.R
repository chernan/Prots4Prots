
#' Benjamini-Hochberg correction
#'
#'
#' The "BH" (aka "fdr") and "BY" method of Benjamini, Hochberg, and Yekutieli control 
#' the false discovery rate, the expected proportion of false discoveries amongst the 
#' rejected hypotheses. The false discovery rate is a less stringent condition than 
#' the family-wise error rate, so these methods are more powerful than the others. 




#' Manual correction, but with a pre-set alpha threshold
manual_BH <- function(nb_p_values, threshold_pval) {
    fdr.BH <- ((1:nb_p_values)*threshold_pval/nb_p_values)
    return(fdr.BH)
}

#' Using p.adjust() function available in R
#' p_values a vector/data.frame slice containing the sorted p-values
#' round if superior to -1, rounds up p-values
#' 
#' Need data to be sorted out before?
R_adjust_BH <- function(p_values, round=-1) {
    fdr.BH <- p.adjust(p_values, "BH")
    if(round > -1) {
        fdr.BH <- round(fdr.BH, round)
    }
    return(fdr.BH)
}

#' Local-pooled-error adjustment
#' Needs an lpe object, output of a Local-Pooled-error method
lpe_adjust_BH <- function(lpe.val) {
    fdr.BH <- fdr.adjust(lpe.val, adjp="BH")
    return(fdr.BH)
}

# # compute FDR adjusted p-values (Benjamini-Hochberg correction)
# order_pval <- order(lpe.val[,"p.vals"], decreasing=TRUE)
# lpe.val <- lpe.val[order_pval, ]
# fdr.BH <- p.adjust(lpe.val[,"p.vals"], "BH")                         # R
# #fdr.BH <- (1:nrow(lpe.val))*threshold_pval/nrow(lpe.val) # manual BH correction, but with a pre-set alpha
# #fdr.BH <- fdr.adjust(lpe.val, adjp="BH")                # LPE package adjustment
# #fdr.BH <- round(fdr.BH, 4)
# 
# lpe.val <- cbind(lpe.val, fdr.BH)
# 
# # fdr.resamp <- fdr.adjust(lpe.val, adjp="resamp", iterations=2)
# # fdr.resamp

#install.packages(c("wordcloud","tm"),repos="http://cran.r-project.org")
# library(wordcloud)
# library(tm)

apply_and_report_BH <- function(matrixdata, output_folder_temp, output_file, threshold_pval) {
    order_pval <- order(matrixdata[,"p.values"], decreasing=FALSE)
    
    pvalues.corr <- R_adjust_BH(c(matrixdata[order_pval,"p.values"]))
    
    #ordered by increasing p-value
    matrixdata2 <- data.frame(
        p.values = matrixdata[order_pval,"p.values"],
        fold.change = matrixdata[order_pval,"fold.change"],
        p.values.corrected = pvalues.corr,
        row.names= c(1:nrow(matrixdata))[order_pval]
    )
    if(any(grepl("significant", names(matrixdata)))) {
        matrixdata2[["significant"]] <- matrixdata[order_pval,"significant"]
    }
    #re-order to have same as input data
    matrixdata3 <- matrixdata2[ order(as.numeric(row.names(matrixdata2))), ]
    matrixdata3 <- data.frame(matrixdata3, row.names=row.names(matrixdata))
    
    sink(output_file, append=TRUE)
    cat('
Significantly different proteins (with multiple testing correction)
--------------------------------------------------------------------------------------------------------

Benjamini-Hochberg method was chosen for multiple testing correction.
')

    all_tests_corrected_Rmd(matrixdata3, threshold_pval, output_folder_temp)
    
    cat('        
>
> Benjamini-Hochberg correction : published in "Controlling the false discovery rate: a practical and powerful approach to multiple testing" by Yoav Benjamini & Yosef Hochberg (1995)
>
>

---------------------------------------------------------------------------
')        
    sink(file=NULL)
    
    return(matrixdata3)
}
