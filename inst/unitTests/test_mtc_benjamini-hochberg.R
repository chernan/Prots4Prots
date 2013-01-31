
test_fun_rAdjustBH <- function() {
    
    ## Bsaic tests
    checkEquals(rAdjustBH(c()), numeric(0))
    checkEquals(rAdjustBH(c(0, 0)), c(0, 0))
    
    ## Test with real values
    x <- rnorm(50, mean=c(rep(0,25), rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    checkEqualsNumeric(
        rAdjustBH(p), 
        p.adjust(p, "BH")
    )
}

test_fun_lpeAdjustBH <- function() {
    library(LPE)
    
    ## Setup data (based on a LPE example)
    data(Ley)
    Ley[,2:7] <- preprocess(Ley[,2:7],data.type="MAS5")
    subsetLey <- Ley[1:1000,]
    var1 <- baseOlig.error(subsetLey[, 2:4], q=0.01)
    var2 <- baseOlig.error(subsetLey[, 5:7], q=0.01)
    lpeResult <- lpe(subsetLey[, 2:4], subsetLey[, 5:7], var1, var2,
                      probe.set.name=subsetLey[, 1])
    ## Test
    checkEquals(
        lpeAdjustBH(lpeResult),
        fdr.adjust(lpeResult, adjp="BH")
    )
}

test_rep_applyAndReportBH <- function() {
    
    ## Setup data
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1)
    )
    ## Test
    checkIdentical(
        applyAndReportBH(matrixData, 
                         file.path(getwd(), 'temp'), 
                         stdout(), #file.path(getwd(), 'temp', 'out_test_applyAndReportBH.Rmd'), 
                         0.05)[["p.values.corrected"]],
        p.adjust(p, "BH")
    )
    
}
