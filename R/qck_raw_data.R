#' @name QualityCheck
#' @title Quality check
#' @description
#' Quality check of raw data.
#' 
#' @concepts Quality_check
#' 
#' @section Introduction:
#' Data to be treated need to be checked for integrity.
#' 
NULL

#' @title Report QC on raw data
#' 
#' @description Report QC on raw data
#' 
#' @details 
#' Report QC on raw data. 
#' TODO Use ExpressionSet object as input
#' 
#' @param dataset Dataset to check.
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @param distMethod Distance method, for the heatmap.
reportQualityCheckRawData <- function(dataset, outputFolderTemp, outputFile, 
                                      distMethod="euclidean") {
    
    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc( 
            runif(1) * 10000)), 
        collapse='')

    tempOutput <- paste(c(outputFolderTemp, '/report_heatmap_Rmd_', 
                          execLabel,'.txt'), collapse='')
    write.table(dataset, tempOutput, sep="\t") 
    
    # Generate report    
    cat('',
        'Quality check',
        '---------------------------------------------------------------------',
        '',
        '',
        paste(
            c('Checking quality of the samples by computing distances (',
              distMethod,
              ') between them, and applying a hierarchical clustering.'), 
            collapse=''),
        '',
        paste(
            c('```{r reportQualityCheckRawData', execLabel, 
              ', echo=FALSE, fig.width=8, fig.height=8}'),
            collapse=''),
        '',
        '', 
        sep = "\n", file=outputFile, append=TRUE)
    
    cat('displayHeatmap <- ', file=outputFile, append=TRUE)
#     print(displayHeatmap)
    cat(paste(deparse(displayHeatmap), collapse="\n"),
        file=outputFile, append=TRUE)

    cat(' ',
        paste(
            c('matrixdata <- as.matrix(read.table("',
              tempOutput,
              '", stringsAsFactors=FALSE))'), 
            collapse=''),
        paste(
            c('displayHeatmap(matrixdata, distMethod="',
              distMethod, '")'), 
            collapse=''),
        '```',
        '',
        '',
        '---------------------------------------------------------------------',
        '',
        sep = "\n", file=outputFile, append=TRUE)
    
    return(dataset)
}
