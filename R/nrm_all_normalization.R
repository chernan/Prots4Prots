#' @name Normalization
#' @title Normalization methods
#' @description 
#' List of normalization methods
#' 
#' @keywords Normalization
#' 
#' @section Introduction:
#' Normalization is the process of making different samples comparable by 
#' adjusting them to a common scale.
#' 
#' @section References:
#' http://en.wikipedia.org/wiki/Normalization_%28statistics%29
#' 
sapply(
    list.files(file.path(getwd(), 'R'), 
               pattern="nrm_[^a][^l]{2}.+\\.R", 
               full.names=TRUE),
    source
)


#' @title Normalization methods.
#' 
#' @description Normalization methods.
#' 
#' @details
#' Apply normalization methods on intensity data. Transform data in log-scale. 
#' Stop if given an incorrect method.
#' 
#' @param normMethod A label indicating which method should be applied.
#'   Available methods include:
#'    log2 : log2 transformation (for comparison)
#'    iqr : inter-quartile range normalization
#'    vsn : variance stabilizing normalization
#'    vsn05 : variance stabilizing normalization (robust to outliers)
#' @param dataset A ExpressionSet with intensities to normalize.
#' @param outputFolderTemp Where the temporary files will be created
#' @param outputFileNameRmd Where report paragraphs will be written
#' @return A data frame with normalized data (same format as input)
#' @export
applyNormalization <- function(normMethod, dataset, 
                               outputFolderTemp, outputFileNameRmd) {
    
    outPVal <- switch(
        normMethod,
        log2 = applyAndReportLog2(exprs(dataset), 
                                  outputFolderTemp, outputFileNameRmd),
        iqr = applyAndReportIQR(exprs(dataset), 
                                outputFolderTemp, outputFileNameRmd),
        vsn = applyAndReportVSN(exprs(dataset), 
                                outputFolderTemp, outputFileNameRmd),
        vsn05 = applyAndReportVSN(exprs(dataset), 
                                  outputFolderTemp, outputFileNameRmd, 
                                  useRobustFit=TRUE),
        stop(
            paste0(
                c(
                    "Unknown method in applyNormalization call : \'", 
                    normMethod,
                    "\'."
                ), 
                sep=''), 
            call=TRUE)
    )
    
    return(outPVal)
}
