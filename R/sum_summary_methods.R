library(plyr)

#' @title Report summarization
#' 
#' @description Report summarization (knitr)
#' 
#' @details 
#' Write a summary report of a summarization step, i.e. summarization of 
#' peptide information to protein level.
#' Data at this step will be saved in a temporary file.
#' 
#' @param dataProteins A ExpressionSet containing the summarized intensities.
#' @param outputFile Report output file (R markdown format)
#' @param outFolder Where the temporary file will be saved
#' 
allTestsSummarizationRmd <- function(dataProteins, 
                                     outputFile, outFolder) {
    
    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc(runif(1) * 10000)), 
        collapse='')
    
    tempOutput <- paste(
        c(outFolder, '/all_tests_summarization_Rmd_data_', execLabel, '.txt'), 
        collapse='')
    write.exprs(dataProteins, tempOutput, sep="\t") 
    
    cat('',
        paste(
            c('```{r allTestsSummarizationRmd', 
              execLabel, ', echo=FALSE, fig.width=12, fig.height=6}'),
            collapse=''),
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    cat('displayNormalizationMAplot <- ', file=outputFile, append=TRUE)
    #     print(displayNormalizationMAplot)
    cat(paste(deparse(displayNormalizationMAplot), collapse="\n"),
        file=outputFile, append=TRUE)
    
    cat(' ',
        paste(
            c('matrixData <- read.table("', tempOutput, '", stringsAsFactors=FALSE)'), 
            collapse=''),
        paste(c('colsExp <- c("', 
                paste(
                    sampleNames(
                        dataProteins[, (dataProteins$Design == "Experiment")]), 
                    collapse='", "'), 
                '")'), 
              collapse=''),
        paste(c('colsCtl <- c("', 
                paste(
                    sampleNames(
                        dataProteins[, (dataProteins$Design == "Control")]), 
                    collapse='", "'), 
                '")'), 
              collapse=''),
        'dataToTest1 <- apply(matrixData[, colsExp], 1, mean)',
        'dataToTest2 <- apply(matrixData[, colsCtl], 1, mean)',
        'displayNormalizationMAplot(dataToTest1, dataToTest2)',
        ' ',
        '```',
        ' ',
        sep="\n", file=outputFile, append=TRUE)
    
}

#' @title Report summarization NA
#' 
#' @description Report summarization NA
#' 
#' @details 
#' Report no summarization, just to test pipeline, or in case of data already 
#' at the protein level. 
#' 
#' @param dataset A ExpressionSet containing raw data
#' @param datasetNorm A dataframe containing normalized data
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A ExpressionSet containing (non-)summarized data as expression data.
applyAndReportSumNA <- function(dataset, datasetNorm, 
                                outputFolderTemp, outputFile) {
    
    #no summary applied, just to test pipeline
    
    dataProteins <- new("ExpressionSet", 
                        exprs = datasetNorm, 
                        phenoData = phenoData(dataset),
                        featureData = featureData(dataset)
    )
    
    ## Generate report    
    cat('',
        'Summary',
        '-------------------------',
        '',
        'No summary method was applied.',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsSummarizationRmd(dataProteins, outputFile, outputFolderTemp)
    
    cat('',
        '',
        '---------------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataProteins)
}


#' @title Report summarization Mean
#' 
#' @description Report summarization Mean
#' 
#' @details 
#' Report summarization by computing the mean of intensities for each channel, 
#' for all peptides of a protein.
#' 
#' @param dataset A ExpressionSet containing raw data
#' @param datasetNorm A dataframe containing normalized data
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A ExpressionSet containing summarized data as expression data.
applyAndReportSumMean <- function(dataset, datasetNorm, 
                                  outputFolderTemp, outputFile) {
    
    #mean by prot
    aggrDataset <- aggregate(
        datasetNorm, 
        by = list(
            Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
        mean)
    
    dataProteins <- new(
        "ExpressionSet", 
        exprs = aggrDataset[, names(datasetNorm)], 
        phenoData = phenoData(dataset),
        featureData = new("AnnotatedDataFrame", data.frame(aggrDataset[, 1]))
    )
    
    ## Generate report
    cat('',
        'Summary',
        '-------------------------',
        '',
        'Mean by protein group.',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsSummarizationRmd(dataProteins, outputFile, outputFolderTemp)
    
    cat('',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataProteins)
}



#' @title Report summarization Median
#' 
#' @description Report summarization Median
#' 
#' @details 
#' Report summarization by computing the median of intensities for each channel, 
#' for all peptides of a protein.
#' 
#' @param dataset A ExpressionSet containing raw data
#' @param datasetNorm A dataframe containing normalized data
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A ExpressionSet containing summarized data as expression data.
applyAndReportSumMedian <- function(dataset, datasetNorm, 
                                    outputFolderTemp, outputFile) {
    
    #median by prot
    aggrDataset <- aggregate(
        datasetNorm, 
        by = list(Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
        median)
    
    dataProteins <- new(
        "ExpressionSet", 
        exprs = aggrDataset[, names(datasetNorm)], 
        phenoData = phenoData(dataset),
        featureData = new("AnnotatedDataFrame", data.frame(aggrDataset[, 1]))
    )
    
    ## Generate report    
    cat('',
        'Summary',
        '-------------------------',
        '',
        'Median by protein group.',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsSummarizationRmd(dataProteins, outputFile, outputFolderTemp)
    
    cat('',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataProteins)
}



#' @title Report summarization Trimmed Mean
#' 
#' @description Report summarization Trimmed Mean
#' 
#' @details 
#' Report summarization by computing the trimmed mean (20%) of intensities 
#' for each channel, for all peptides of a protein.
#' 
#' @param dataset A ExpressionSet containing raw data
#' @param datasetNorm A dataframe containing normalized data
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A ExpressionSet containing summarized data as expression data.
applyAndReportSumTrimmedMean20 <- function(dataset, datasetNorm, 
                                           outputFolderTemp, outputFile) {
    
    #trimmed mean by prot
    aggrDataset <- aggregate(
        datasetNorm, 
        by = list(Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
        FUN=function(x){mean(x,trim=0.1)}
    )
    
    dataProteins <- new(
        "ExpressionSet", 
        exprs = aggrDataset[, names(datasetNorm)], 
        phenoData = phenoData(dataset),
        featureData = new("AnnotatedDataFrame", data.frame(aggrDataset[, 1]))
    )
    
    ## Generate report    
    cat('',
        'Summary',
        '-------------------------',
        '',
        'Trimmed mean by protein group (20%).',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsSummarizationRmd(dataProteins, outputFile, outputFolderTemp)
    
    cat('',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataProteins)
}



#' @title Report summarization Most variable peptide
#' 
#' @description Report summarization Most variable peptide
#' 
#' @details 
#' Report summarization by selecting the most variable peptide for each protein.
#' 
#' @param dataset A ExpressionSet containing raw data
#' @param datasetNorm A dataframe containing normalized data
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A ExpressionSet containing summarized data as expression data.
applyAndReportSumMostVariablePep <- function(dataset, datasetNorm, 
                                             outputFolderTemp, outputFile) {
    #Look 
    tempDF <- data.frame(
        datasetNorm, 
        "Majority.protein.IDs"=featureData(dataset)[["Majority.protein.IDs"]])
    aggrDataset <- ddply(.data=tempDF, 
                         .variables="Majority.protein.IDs", 
                         .fun=function(df) {
                             variances <- apply(df[, names(datasetNorm)], 1, var)
                             max.var <- max(variances)
                             return(df[(which(variances == max.var))[1], ])
                         }
    )
    
    dataProteins <- new(
        "ExpressionSet", 
        exprs = data.frame(aggrDataset[, names(datasetNorm)]), 
        phenoData = phenoData(dataset),
        featureData = new(
            "AnnotatedDataFrame", 
            data.frame(
                "Majority.protein.IDs"=aggrDataset[["Majority.protein.IDs"]]))
    )
    
    ## Generate report    
    cat('',
        'Summary',
        '-------------------------',
        '',
        'Most variable peptide for each protein.',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsSummarizationRmd(dataProteins, outputFile, outputFolderTemp)
    
    cat('',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataProteins)
}

#' @title Report summarization Most intense peptide
#' 
#' @description Report summarization Most intense peptide
#' 
#' @details 
#' Report summarization by selecting the most intense peptide of a protein.
#' 
#' @param dataset A ExpressionSet containing raw data
#' @param datasetNorm A dataframe containing normalized data
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A ExpressionSet containing summarized data as expression data.
applyAndReportSumMostintensePep <- function(dataset, datasetNorm, 
                                            outputFolderTemp, outputFile) {
    
    tempDF <- data.frame(
        datasetNorm, 
        "Majority.protein.IDs"=featureData(dataset)[["Majority.protein.IDs"]])
    aggrDataset <- ddply(.data=tempDF, 
                         .variables="Majority.protein.IDs", 
                         .fun=function(df) {
                             sums <- apply(df[, names(datasetNorm)], 1, sum)
                             max.sum <- max(sums)
                             return(df[(which(sums == max.sum))[1],])
                         }
    )
    
    dataProteins <- new(
        "ExpressionSet", 
        exprs = data.frame(aggrDataset[, names(datasetNorm)]), 
        phenoData = phenoData(dataset),
        featureData = new(
            "AnnotatedDataFrame",
            data.frame(
                "Majority.protein.IDs"=aggrDataset[["Majority.protein.IDs"]]))
    )
    
    ## Generate report    
    cat('',
        'Summary',
        '-------------------------',
        '',
        'Most insense peptide for each protein.',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsSummarizationRmd(dataProteins, outputFile, outputFolderTemp)
    
    cat('',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataProteins)
}

#' @title Report summarization Sum Intensities
#' 
#' @description Report summarization Sum intensities
#' 
#' @details 
#' Report summarization by computing the sum of intensities for each channel, 
#' for all peptides of a protein.
#' 
#' @param dataset A ExpressionSet containing raw data
#' @param datasetNorm A dataframe containing normalized data
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @return A ExpressionSet containing summarized data as expression data.
#' Output a new paragraph in an Rmd file describing the outcome of applying this method, including code to generate graphs.
applyAndReportSumIntens <- function(dataset, datasetNorm, 
                                    outputFolderTemp, outputFile) {
    
    #mean by prot
    aggrDataset <- aggregate(
        datasetNorm, 
        by = list(
            Majority.protein.IDs = featureData(dataset)[["Majority.protein.IDs"]]), 
        sum)
    
    dataProteins <- new(
        "ExpressionSet", 
        exprs = aggrDataset[, names(datasetNorm)], 
        phenoData = phenoData(dataset),
        featureData = new("AnnotatedDataFrame", data.frame(aggrDataset[, 1]))
    )
    
    ## Generate report    
    cat('',
        'Summary',
        '-------------------------',
        '',
        'Sum intensities of peptides by protein group.',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsSummarizationRmd(dataProteins, outputFile, outputFolderTemp)
    
    cat('',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataProteins)
}
