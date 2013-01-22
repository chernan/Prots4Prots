
test_rep_allTestsNormalizationRmd <- function() {
    
    ## Setup data
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    ## Run function... as test?
    allTestsNormalizationRmd(normalData, 
                         stdout(),
                         file.path(getwd(), 'temp')
    )
}


test_gra_displayNormalizationViolin <- function() {

    ## Setup data
    normalData <- matrix(
        replicate(8, rnorm(n=500, mean=11, sd=2.5)),
        ncol=8, nrow=500
    )
    rownames(normalData) <- as.character(1:500)
    colnames(normalData) <- as.character(1:8)
    
    ## Run function... as test?
    print(displayNormalizationViolin(normalData))
    
}

test_gra_displayNormalizationMAplot <- function() {
    
    ## Setup data
    normalData1 <- replicate(4, rnorm(n=500, mean=11, sd=2.5))
    normalData2 <- replicate(4, rnorm(n=500, mean=11, sd=2.5))
    dataToTest1 <- apply(normalData1, 1, mean)
    dataToTest2 <- apply(normalData2, 1, mean)
    
    ## Run function... as test?
    print(displayNormalizationMAplot(dataToTest1, dataToTest2))
    
}

test_gra_displayNormalizationQQplot <- function() {
    
    ## Setup data
    normalData <- rnorm(n=500, mean=11, sd=2.5)
    
    ## Run function... as test?
    displayNormalizationQQplot(normalData)
    
}
