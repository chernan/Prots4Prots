
test_rep_allTestsCorrectedRmd <- function() {
    
    ## Setup dataset
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1),
        p.values.corrected = p.adjust(p, "BH")
    )
    
    ## Run function... as test?
    allTestsCorrectedRmd(matrixData, 
                         0.05,
                         stdout(),
                         file.path(getwd(), 'temp')
    )
}


test_gra_displayCorrectedVolcanoPlot <- function() {

    ## Setup data
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1),
        p.values.corrected = p.adjust(p, "BH")
    )
    
    ## Run function... as test?
    displayCorrectedVolcanoPlot(matrixData, 0.05, 
                                title="test_displayCorrectedVolcanoPlot")
    displayCorrectedVolcanoPlot(matrixData, 0.05)
    
}
