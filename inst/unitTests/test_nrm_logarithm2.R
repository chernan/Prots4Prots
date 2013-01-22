
test_fun_applyLog2 <- function() {
    
    ## Test on normal data
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    applyLog2(normalData)
    
    ## Test all zero data
    zeroData <- matrix(c(0, 0, 0, 
                         0, 0, 0, 
                         0, 0, 0), nrow=3, ncol=3)
    checkEquals(
        applyLog2(zeroData),
        matrix(c(-Inf, -Inf, -Inf, 
                 -Inf, -Inf, -Inf, 
                 -Inf, -Inf, -Inf), nrow=3, ncol=3)
    )
    
    ## TODO check names of elements before / after
    
}

test_rep_applyAndReportLog2 <- function() {

    ## Setup data
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    ## Test
    applyAndReportLog2(
        normalData, 
        file.path(getwd(), 'temp'),
        stdout()
    )
    
}

