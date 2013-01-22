
test_fun_applyIqrFromLpe <- function() {
    
    ## Setup data
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    normalData <- 2^normalData
    
    ## Setup data
    applyIqrFromLpe(normalData)
    
    ## TODO check names of elements before / after
    
}

test_rep_applyAndReportIQR <- function() {

    ## Setup data
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    normalData <- 2^normalData
    
    ## Test
    applyAndReportIQR(
        normalData, 
        file.path(getwd(), 'temp'),
        stdout()
    )
    
}

