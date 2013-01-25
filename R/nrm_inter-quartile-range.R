#' @name InterQuartileRange
#' @title InterQuartile Range Normalization
#' @description 
#' InterQuartile Range Normalization.
#' 
#' @concepts Normalization
#' 
#' @section Introduction:
#' "The median and interquartile range are computed for the set of all gene 
#' expression values for each array. The median of the medians and the median 
#' of the IQRs are calculated. For each array, values are initially adjusted 
#' by subtracting the median, and then they are multiplied by a constant 
#' specific to each array chosen such that the scaled values have an IQR equal 
#' to the median of the IQRs for all arrays (i.e., ratio of the median of all 
#' IQR values to the array IQR). Finally, the median of the medians is added 
#' to all arrays. This produces expression values in which all arrays have the 
#' same median and same IQR." From ???.
#' 
#' This method is not sensitive to outliers as it normalizes data relatively 
#' to the inter-quartile range only (between 25th and 75th percentiles).
#' 
#' NB: there are many methods to compute a percentile (other than the median) 
#' and they may not give the same result...
#' 
#' @section References:
#' 
#' Hyndman, R. J., & Fan, Y. (1996). 
#' Sample quantiles in statistical packages. 
#' The American Statistician, 50(4), 361-365.
#' Retrieved from http://www.tandfonline.com/doi/abs/10.1080/00031305.1996.10473566
#' 
#' @import LPE
# set.seed(0)
library(LPE)

#' @title Apply interquartile range normalization.
#' 
#' @description Apply interquartile range normalization.
#' 
#' @details 
#' Apply interquartile range normalization and log2-transform data. This 
#' function uses the implementation available in the LPE R package.
#'  
#' @param dataNotNorm Data to be normalized.
#' @return A dataframe containing normlized data.
applyIqrFromLpe <- function(dataNotNorm) {
    
    dataNormalized <- preprocess(
        # non normalized data
        dataNotNorm, 
        # is MAS5 ok for proteomics data?
        data.type="MAS5", 
        # thresholding value' below which all data would be thresholded
        threshold=1,
        # lowess normalization needs to be performed?
        LOWESS=FALSE
    )
    
    return(dataNormalized)
}


#' @title Apply IQR normalization and report result
#' 
#' @description Apply IQR normalization and report result
#' 
#' @details 
#' Apply IQR normalization on a dataset. 
#' TODO Use ExpressionSet object as input
#' 
#' @param dataNotNorm Dataset to be normalized.
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A data frame with normalized data.
applyAndReportIQR <- function(dataNotNorm, outputFolderTemp, outputFile) {
    
    dataNormalized <- applyIqrFromLpe(dataNotNorm)
    
    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc( 
            runif(1) * 10000)), 
        collapse='')
    
    cat(' ',
        'Normalization (IQR)',
        '---------------------------------------------------------------------',
        '',
        paste(
            c('```{r applyAndReportIQR', 
              execLabel, ', echo=FALSE, warning=FALSE}'),
            collapse=''),
        'lpe_citation <- citation("LPE")',
        '```',
        'Normalisation was achieved by using the Inter-Quartile-Range method, as available in the LPE R package (`r lpe_citation$note`).',
        '',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsNormalizationRmd(dataNormalized, outputFile=outputFile,
                             outFolder=outputFolderTemp)
    
    cat('',
        '> ',
        '> Please cite these articles whenever using results from this software in a publication :',
        '> ',
        '> Software article (R package) :',
        '> `r lpe_citation$title` by `r lpe_citation$author`',
        '> ',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataNormalized)
}
