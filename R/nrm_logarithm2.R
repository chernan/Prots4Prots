#' @name Log2
#' @title Log2 transformation
#' @description
#' Log2 transformation.
#' 
#' @concepts Normalization
#' 
#' @section Introduction:
#' Transforms data to convert them to a linear scale.
#' 
NULL

#' @title Apply log2.
#' 
#' @description Apply log2.
#' 
#' @details 
#' Apply log 2 to transforms data and convert them to a linear scale.
#' Data should not contain 0 values.
#'  
#' @param dataNotNorm Data to be normalized.
#' @return A dataframe containing normlized data.
applyLog2 <- function(dataNotNorm) {
    
    dataNormalized <- log(dataNotNorm, base=2)
    
    return(dataNormalized)
}

#' @title Apply log2 transformation and report result
#' 
#' @description Apply log2 transformation and report result
#' 
#' @details 
#' Apply log2 transformation on a dataset. 
#' TODO Use ExpressionSet object as input
#' 
#' @param dataNotNorm Dataset to be normalized, should not contain 0 values.
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A data frame with normalized data.
applyAndReportLog2 <- function(dataNotNorm, outputFolderTemp, outputFile) {
    
    ## Transform data
    dataNormalized <- applyLog2(dataNotNorm)
    
    ## Generate report
    
    cat('',
        'Normalization (log2 transformation)',
        '---------------------------------------------------------------------',
        '',
        'Base 2 logarithm was applied on the raw intensities to convert them to a linear scale. Note that this transformation is only for comparison with other normalization methods. It\'s not a normalization method.',
        '',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsNormalizationRmd(dataNormalized, outputFile=outputFile,
                             outFolder=outputFolderTemp)
    
    cat('',
        '> Log2 transformation : require no citation.',
        '>',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)        
    
    return(dataNormalized)
}

