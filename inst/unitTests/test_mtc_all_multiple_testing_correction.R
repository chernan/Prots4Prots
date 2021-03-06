
test_rep_applyMultipleTestingCorrection <- function() {

    ## Setup data
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1)
    )
    
    # Test
    checkIdentical(
        applyMultipleTestingCorrection(
            'bh',
            matrixData, 
            file.path(getwd(), 'temp'), 
            stdout(), #file.path(getwd(), 'temp', 'out_test_applyMultipleTestingCorrection.Rmd'), 
            0.05)[["p.values.corrected"]],
        p.adjust(p, "BH")
    )
    
    checkException(
        applyMultipleTestingCorrection(
            '',
            matrixData, 
            file.path(getwd(), 'temp'), 
            stdout(), #file.path(getwd(), 'temp', 'out_test_applyMultipleTestingCorrection.Rmd'), 
            0.05),
        silent = TRUE
    )
}
