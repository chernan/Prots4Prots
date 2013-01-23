
test_rep_reportQualityCheckRawData <- function() {
    
    normalData <- replicate(8, rnorm(n=3000, mean=11.5, sd=2.2))
    normalData <- 2^normalData
    
    ## Test
    reportQualityCheckRawData(
        normalData, 
        file.path(getwd(), 'temp'), 
        stdout(), #file.path(getwd(), 'temp', 'out_test_rep_reportQualityCheckRawData.Rmd'),
        "euclidean"
    )
    
}

