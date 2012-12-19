#' Benjamini-Hochberg procedure
#' 
#' Benjamini-Hochberg procedure to control FDR.
#' 
#' The Benjamini-Hochberg procedure controls the false discovery rate, the 
#' expected proportion of false discoveries amongst the positive discoveries 
#' (rejected null hypotheses). This procedure, as well as the Benjamini-
#' Yekutieli, are less stringent than the family-wise error rate, making them 
#' more powerful than other multiple testing correction approaches. 
#' 
#' @references http://en.wikipedia.org/wiki/False_discovery_rate#Controlling_procedures
#' @references "Controlling the false discovery rate: a practical and powerful 
#' approach to multiple testing" by Yoav Benjamini & Yosef Hochberg (1995)


library("LPE")

#' @title Wrapper to R Base p.adjust()
#' 
#' @description Wrapper to R Base p.adjust()
#' 
#' @details 
#' Using \code{p.adjust} function available in R Base. p-values really 
#' need to be sorted out beforehand?
#' 
#' @param pValues a vector/data.frame slice containing the sorted p-values
#' @return a vector of corrected p-values
rAdjustBH <- function(pValues) {
    
    fdrBH <- p.adjust(pValues, "BH")
    
    return(fdrBH)
}

#' @title Wrapper to implementation in package LPE
#' 
#' @description Wrapper to implementation in package LPE
#' 
#' @details 
#' Implementation of the Benjamini-Hochberg procedure as available in the 
#' Local-Pooled Error (LPE) package
#' 
#' @param lpeVal A lpe object, output of a Local-Pooled-error method
#' @return a table of thresholds for a set of pre-defined FDR values
lpeAdjustBH <- function(lpeVal) {

    
    fdrBH <- fdr.adjust(lpeVal, adjp="BH")
    
    return(fdrBH)
}

#' @title Apply Benjamini-Hochberg and report
#' 
#' @description Apply Benjamini-Hochberg and report
#' 
#' @details 
#' Calls \code{R_adjust_BH} and write a report in a Rmd file. Data generated in 
#' this step are also saved in an intermediary table (CSV file).
#' 
#' @param matrixData Any data.frame containing a column entitled "p.values".  
#'   This column should contain the p-values to be adjusted.
#'   Other columns should include "fold.change", can include "significant".
#' @param outputFolderTemp Temporary folder to store data generated at this step of the analysis.
#' @param outputFile Report of current analysis step will be appended to this Rmd file.
#' @param thresholdPVal Chosen False Discovery rate value.
#' @return a data frame with data in matrixData (p.values, fold.change, 
#'   significant) and the corrected p-values.
applyAndReportBH <- function(matrixData, outputFolderTemp, outputFile, 
                             thresholdPVal) {
    
    orderPVal <- order(matrixData[, "p.values"], decreasing=FALSE)
    
    pValuesCorr <- rAdjustBH( c(matrixData[orderPVal, "p.values"]))
    
    ## ordered by increasing p-value
    matrixData2 <- data.frame(
        p.values = matrixData[orderPVal, "p.values"],
        fold.change = matrixData[orderPVal, "fold.change"],
        p.values.corrected = pValuesCorr,
        row.names = c(1:nrow(matrixData))[orderPVal]
    )
    if (any(grepl("significant", names(matrixData)))) {
        matrixData2[["significant"]] <- matrixData[orderPVal,"significant"]
    }
    ## re-order to have same as input data
    matrixData3 <- matrixData2[order(as.numeric(row.names(matrixData2))), ]
    matrixData3 <- data.frame(matrixData3, row.names=row.names(matrixData))
    
    sink(outputFile, append=TRUE)
    cat('',
        'Significantly different proteins (with multiple testing correction)',
        '---------------------------------------------------------------------',
        '',
        'Benjamini-Hochberg method was chosen for multiple testing correction.',
        '',
        sep="\n")

    allTestsCorrectedRmd(matrixData3, thresholdPVal, outputFolderTemp)
    
    cat('',        
        '>',
        '> Citation :',
        '> "Controlling the false discovery rate: a practical and powerful approach to multiple testing"',
        '> Yoav Benjamini & Yosef Hochberg (1995)',
        '>',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n")        
    sink(file=NULL)
    
    return(matrixData3)
}
