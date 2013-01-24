
test_rep_applySummarization <- function() {

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
    
    
    # Test
    applySummarization(
        'NonA',
        dataSet, normalData,
        file.path(getwd(), 'temp'), 
        stdout()
    )
    
    applySummarization(
        'meanPep',
        dataSet, normalData,
        file.path(getwd(), 'temp'), 
        stdout()
    )

    applySummarization(
        'medianPep',
        dataSet, normalData,
        file.path(getwd(), 'temp'), 
        stdout()
    )
    
    applySummarization(
        'trimmed20Pep',
        dataSet, normalData,
        file.path(getwd(), 'temp'), 
        stdout() 
    )
    
    applySummarization(
        'variablePep',
        dataSet, normalData,
        file.path(getwd(), 'temp'), 
        stdout()
    )

    applySummarization(
        'intensePep',
        dataSet, normalData,
        file.path(getwd(), 'temp'), 
        stdout()
    )
    
    applySummarization(
        'sumIntens',
        dataSet, normalData,
        file.path(getwd(), 'temp'), 
        stdout() 
    )
    
    checkException(
        applySummarization(
            '',
            dataSet, normalData,
            file.path(getwd(), 'temp'), 
            stdout()
            ),
        silent = TRUE
    )
    
}
