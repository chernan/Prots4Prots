
test_rep_applyNormalization <- function() {

    ## Setup data
    normalData <- replicate(8, rnorm(n=3000, mean=11.5, sd=2.2))
    normalData <- 2^normalData
    pData <- data.frame(Design=c(rep("Experiment",4), rep("Control",4)))
    rownames(pData) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8" )
    
    dataSet <- new("ExpressionSet", exprs = normalData, 
                   phenoData = new("AnnotatedDataFrame",pData))
    
    
    # Test
    applyNormalization(
        'log2',
        dataSet, 
        file.path(getwd(), 'temp'), 
        stdout() #file.path(getwd(), 'temp', 'out_test_nrm_applyNormalization.Rmd')
    )
    
    applyNormalization(
        'iqr',
        dataSet, 
        file.path(getwd(), 'temp'), 
        stdout() #file.path(getwd(), 'temp', 'out_test_nrm_applyNormalization.Rmd')
    )

    applyNormalization(
        'vsn',
        dataSet, 
        file.path(getwd(), 'temp'), 
        stdout() #file.path(getwd(), 'temp', 'out_test_nrm_applyNormalization.Rmd')
    )
    
    applyNormalization(
        'vsn05',
        dataSet, 
        file.path(getwd(), 'temp'), 
        stdout() #file.path(getwd(), 'temp', 'out_test_nrm_applyNormalization.Rmd')
    )
    
    checkException(
        applyNormalization(
            '',
            dataSet, 
            file.path(getwd(), 'temp'), 
            stdout() #file.path(getwd(), 'temp', 'out_test_nrm_applyNormalization.Rmd')
            ),
        silent = TRUE
    )
    
}
