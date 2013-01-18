
test_rep_allTestsSignificanceRmd <- function() {
    
    ## Setup dataset
    set.seed(0)
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1)
    )
    
    ## Run function... as test?
    allTestsSignificanceRmd(matrixData, 
                            0.05,
                            stdout(),
                            file.path(getwd(), 'temp')
    )
    
}


test_gra_displaySignificanceVolcanoPlot <- function() {
    
    ## Setup data
    set.seed(0)
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1)
    )
    
    ## Run function... as test?
    displaySignificanceVolcanoPlot(matrixData, 0.05, 
                                   title="test_displaySignificanceVolcanoPlot")
    displaySignificanceVolcanoPlot(matrixData, 0.05)
    
}
