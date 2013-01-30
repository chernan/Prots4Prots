#' @name Summarization
#' @title Summarization methods
#' @description
#' List of Summarization methods.
#' 
#' @keywords Summarization
#' 
#' @section Introduction:
#' Summarization is the process of summarizing peptide information up to the 
#' protein level.
#' The protein groups have already been computed by identification tools, 
#' like MaxQuant/Andromeda.
#' 
#' 
sapply(
    list.files(file.path(getwd(), 'R'), 
               pattern="sum_[^a][^l]{2}.+\\.R", 
               full.names=TRUE),
    source
)

#' @title Summarization methods.
#' 
#' @description Summarization methods.
#' 
#' @details
#' Apply summarization methods on intensity data. 
#' Stop if given an incorrect method.
#' 
#' @param summMethod A label indicating which method should be applied.
#'   Available methods include:
#'    NonA : "non applicable", doesn't summarize. Ex. of protein-level data.
#'    meanPep : mean of peptides
#'    medianPep : median of peptides
#'    trimmed20Pep : trimmed mean (20%)
#'    variablePep : most variable peptide
#'    intensePep : most intense peptide
#'    sumIntens : sum of all peptide intensities for each reporter
#' @param dataset A ExpressionSet with raw intensities.
#' @param outVal A data.frame with normalized intensities.
#' @param outputFolderTemp Where the temporary files will be created
#' @param outputFileNameRmd Where report paragraphs will be written
#' @return A ExpressionSet with normalized data and additional info from raw 
#'  ExpressionSet.
applySummarization <- function(summMethod, dataset, outVal, 
                               outputFolderTemp, outputFileNameRmd) {
    pgESet <- switch(
        summMethod,
        NonA = applyAndReportSumNA(dataset, outVal, 
                                   outputFolderTemp, outputFileNameRmd),
        meanPep = applyAndReportSumMean(dataset, outVal, 
                                        outputFolderTemp, outputFileNameRmd),
        medianPep = applyAndReportSumMedian(dataset, outVal, 
                                            outputFolderTemp, outputFileNameRmd),
        trimmed20Pep = applyAndReportSumTrimmedMean20(dataset, outVal, 
                                                      outputFolderTemp, 
                                                      outputFileNameRmd),
        variablePep = applyAndReportSumMostVariablePep(dataset, outVal, 
                                                       outputFolderTemp, 
                                                       outputFileNameRmd),
        intensePep = applyAndReportSumMostintensePep(dataset, outVal, 
                                                     outputFolderTemp, 
                                                     outputFileNameRmd),
        sumIntens = applyAndReportSumIntens(dataset, outVal, 
                                            outputFolderTemp, outputFileNameRmd),
        stop(
            paste0(
                c(
                    "Unknown method in applySummarization call : \'", 
                    summMethod,
                    "\'."
                ), 
                sep=''), 
            call=TRUE)
    )
    
    return(pgESet)
}