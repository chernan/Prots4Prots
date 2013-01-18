#' Multiple testing correction
#' 
#' When testing multiple times a null hypothesis, incorrect rejection is more 
#' likely to occur when considering all tests together. Several methods have 
#' been developped in order to correct this problem and maintain FDR to an 
#' acceptable level.
#' 
#' @seealso http://en.wikipedia.org/wiki/Multiple_comparisons
#' 
sapply(
    list.files(file.path(getwd(), 'R'), 
               pattern="mtc_[^a][^l]{2}.+\\.R", 
               full.names=TRUE),
    source
)

#' @title Multiple testing correction
#' 
#' @description Multiple testing correction
#' 
#' @details
#' Apply multiple testing correction on the outcome of statistical tests.
#' 
#' @param mtcMethod A label indicating which method should be applied.
#'   Available methods include:
#'     bh : Benjamini-Hochberg FDR
#' @param matrixData Object containg the p-values to correct
#'   Any data.frame containing a column entitled "p.values" containing values 
#'   to be adjusted. Other columns should include "fold.change", can include 
#'   "significant".
#' @param outputFolderTemp Where the temporary files will be created
#' @param outputFileNameRmd Where report paragraphs will be written
#' @param thresholdPVal FDR threshold
#' @return A data frame with data in matrixData (p.values, fold.change, 
#'   significant) and the corrected p-values.
applyMultipleTestingCorrection <- function(mtcMethod, 
                                           matrixData, 
                                           outputFolderTemp, 
                                           outputFileNameRmd, 
                                           thresholdPVal) {
    
    mtcPVal <- switch(
        mtcMethod,
        bh = applyAndReportBH(matrixData, outputFolderTemp, 
                              outputFileNameRmd, thresholdPVal),
        stop(
            paste0(
                c(
                    "Unknown method in applyMultipleTestingCorrection call : \'", 
                    mtcMethod,
                    "\'."
                ), 
                sep=''), 
            call=TRUE)
    )
    
    return(mtcPVal)
}