
test_fun_applyVSN09 <- function() {
    
    ## Setup data
    set.seed(0)
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    applyVSN09(normalData)
    
    ## TODO check names of elements before / after
    
}

test_fun_applyVSN05 <- function() {
    
    ## Setup data
    set.seed(0)
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    applyVSN05(normalData)
    
    ## TODO check names of elements before / after
    
}

test_rep_applyAndReportVSN <- function() {

    ## Setup data
    set.seed(0)
    normalData <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    ## Test
    applyAndReportVSN(
        normalData, 
        file.path(getwd(), 'temp'),
        stdout(),
        useRobustFit=FALSE)
    applyAndReportVSN(
        normalData, 
        file.path(getwd(), 'temp'),
        stdout(),
        useRobustFit=TRUE)
    
}

