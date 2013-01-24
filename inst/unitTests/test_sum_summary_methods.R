
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

test_rep_applyAndReportSumMedian <- function() {
    
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
    applyAndReportSumMedian(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportSumTrimmedMean20 <- function() {
    
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
    applyAndReportSumTrimmedMean20(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportSumMostVariablePep <- function() {
    
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
    applyAndReportSumMostVariablePep(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportSumMostintensePep <- function() {
    
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
    applyAndReportSumMostintensePep(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}


test_rep_applyAndReportSumIntens <- function() {
    
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
    applyAndReportSumIntens(dataSet, normalData, 
                          file.path(getwd(), 'temp'), stdout())    
    
}
