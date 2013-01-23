
test_rep_allTestsSummarizationRmd <- function() {
    
    ## Setup data
    normalData <- replicate(8, rnorm(n=500, mean=11.5, sd=2.2))
    normalData <- 2^normalData
    
    pData <- data.frame(Design=c(rep("Experiment",4), rep("Control",4)))
    rownames(pData) <- paste("V", 1:8, sep='')
    
    dataSet <- new("ExpressionSet", exprs = normalData, 
                   phenoData = new("AnnotatedDataFrame", pData))
    
    ## Run function... as test?
    allTestsSummarizationRmd(dataSet,
                             stdout(),
                             file.path(getwd(), 'temp')
    )
    
}

test_rep_applyAndReportSumNA <- function() {

    ## Setup data
    normalData <- replicate(8, rnorm(n=500, mean=11.5, sd=2.2))
    normalData2 <- 2^normalData
    
    pData <- data.frame(Design=c(rep("Experiment",4), rep("Control",4)))
    rownames(pData) <- paste("V", 1:8, sep='')
    
    dataSet <- new("ExpressionSet", exprs = normalData2, 
                   phenoData = new("AnnotatedDataFrame", pData))
    
    ## Run function... as test?
    applyAndReportSumNA(dataSet, normalData, 
                        file.path(getwd(), 'temp'), stdout())    
    
}

test_rep_applyAndReportSumMean <- function() {
    
    ## Setup data
    normalData <- data.frame(replicate(8, rnorm(n=500, mean=11.5, sd=2.2)))
    names(normalData) <- paste("V", 1:8, sep='')
    normalData2 <- 2^normalData
    
    pData <- data.frame(
        Design=c(rep("Experiment", 4), rep("Control", 4)),
        row.names = paste("V", 1:8, sep='')
    )


    fData <- data.frame(
        Majority.protein.IDs = rep(paste("Protein", 1:100, sep=''), each=5)
    )
    
    dataSet <- new("ExpressionSet", exprs = normalData2, 
                   phenoData = new("AnnotatedDataFrame", pData),
                   featureData = new("AnnotatedDataFrame", fData))
    
    ## Run function... as test?
    applyAndReportSumMean(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}

test_rep_applyAndReportMedian <- function() {
    
    ## Setup data
    normalData <- data.frame(replicate(8, rnorm(n=500, mean=11.5, sd=2.2)))
    names(normalData) <- paste("V", 1:8, sep='')
    normalData2 <- 2^normalData
    
    pData <- data.frame(
        Design=c(rep("Experiment", 4), rep("Control", 4)),
        row.names = paste("V", 1:8, sep='')
    )
    
    
    fData <- data.frame(
        Majority.protein.IDs = rep(paste("Protein", 1:100, sep=''), each=5)
    )
    
    dataSet <- new("ExpressionSet", exprs = normalData2, 
                   phenoData = new("AnnotatedDataFrame", pData),
                   featureData = new("AnnotatedDataFrame", fData))
    
    ## Run function... as test?
    applyAndReportMedian(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportTrimmedMean20 <- function() {
    
    ## Setup data
    normalData <- data.frame(replicate(8, rnorm(n=500, mean=11.5, sd=2.2)))
    names(normalData) <- paste("V", 1:8, sep='')
    normalData2 <- 2^normalData
    
    pData <- data.frame(
        Design=c(rep("Experiment", 4), rep("Control", 4)),
        row.names = paste("V", 1:8, sep='')
    )
    
    
    fData <- data.frame(
        Majority.protein.IDs = rep(paste("Protein", 1:100, sep=''), each=5)
    )
    
    dataSet <- new("ExpressionSet", exprs = normalData2, 
                   phenoData = new("AnnotatedDataFrame", pData),
                   featureData = new("AnnotatedDataFrame", fData))
    
    ## Run function... as test?
    applyAndReportTrimmedMean20(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportMostVariablePep <- function() {
    
    ## Setup data
    normalData <- data.frame(replicate(8, rnorm(n=500, mean=11.5, sd=2.2)))
    names(normalData) <- paste("V", 1:8, sep='')
    normalData2 <- 2^normalData
    
    pData <- data.frame(
        Design=c(rep("Experiment", 4), rep("Control", 4)),
        row.names = paste("V", 1:8, sep='')
    )
    
    
    fData <- data.frame(
        Majority.protein.IDs = rep(paste("Protein", 1:100, sep=''), each=5)
    )
    
    dataSet <- new("ExpressionSet", exprs = normalData2, 
                   phenoData = new("AnnotatedDataFrame", pData),
                   featureData = new("AnnotatedDataFrame", fData))
    
    ## Run function... as test?
    applyAndReportMostVariablePep(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportMostintensePep <- function() {
    
    ## Setup data
    normalData <- data.frame(replicate(8, rnorm(n=500, mean=11.5, sd=2.2)))
    names(normalData) <- paste("V", 1:8, sep='')
    normalData2 <- 2^normalData
    
    pData <- data.frame(
        Design=c(rep("Experiment", 4), rep("Control", 4)),
        row.names = paste("V", 1:8, sep='')
    )
    
    
    fData <- data.frame(
        Majority.protein.IDs = rep(paste("Protein", 1:100, sep=''), each=5)
    )
    
    dataSet <- new("ExpressionSet", exprs = normalData2, 
                   phenoData = new("AnnotatedDataFrame", pData),
                   featureData = new("AnnotatedDataFrame", fData))
    
    ## Run function... as test?
    applyAndReportMostintensePep(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportSumintens <- function() {
    
    ## Setup data
    normalData <- data.frame(replicate(8, rnorm(n=500, mean=11.5, sd=2.2)))
    names(normalData) <- paste("V", 1:8, sep='')
    normalData2 <- 2^normalData
    
    pData <- data.frame(
        Design=c(rep("Experiment", 4), rep("Control", 4)),
        row.names = paste("V", 1:8, sep='')
    )
    
    
    fData <- data.frame(
        Majority.protein.IDs = rep(paste("Protein", 1:100, sep=''), each=5)
    )
    
    dataSet <- new("ExpressionSet", exprs = normalData2, 
                   phenoData = new("AnnotatedDataFrame", pData),
                   featureData = new("AnnotatedDataFrame", fData))
    
    ## Run function... as test?
    applyAndReportSumintens(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}
