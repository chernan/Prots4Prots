
test_fun_applyIqrFromLpe <- function() {
    
    set.seed(0)
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    ## Setup data
    applyIqrFromLpe(normalData)
    
    ## TODO check names of elements before / after
    
}

test_rep_applyAndReportIQR <- function() {

    ## Setup data
    set.seed(0)
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    ## Test
    applyAndReportIQR(
        normalData, 
        file.path(getwd(), 'temp'),
        stdout()
    )
    
}

