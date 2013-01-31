
test_fun_applyLpeMethod <- function() {
    library(LPE)
    
    ## This test doesn't work when run from test suite only! 
    # Error in smooth.spline(base.var.adap[, 1], base.var.adap[, 2], df = df) : 
    # missing or infinite values in inputs are not allowed 
    ## Ok on command line... :-/
#     set.seed(0)
#     normalData1 <- log2(replicate(3, rnorm(n=500, mean=11, sd=2.5)))
#     normalData2 <- log2(replicate(3, rnorm(n=500, mean=11, sd=2)))
#     
#     varActivated <- baseOlig.error(normalData2)
#     varNaive <- baseOlig.error(normalData1)
#     lpeVal <- lpe(normalData2, normalData1, varActivated, varNaive)
#     
#     ## Computation of p-values
#     checkEqualsNumeric(
#         applyLpeMethod(normalData2, normalData1, rowNames=1:500)$z.stats,
#         lpeVal$z.stats
#     )
    
    
    ## Setup data (based on a LPE example)
    data(Ley)
    Ley[,2:7] <- preprocess(Ley[,2:7],data.type="MAS5")
    subsetLey <- Ley[1:1000,]
    var1 <- baseOlig.error(subsetLey[,2:4], q=0.01)
    var2 <- baseOlig.error(subsetLey[,5:7], q=0.01)
    lpeResult <- lpe(subsetLey[,2:4], subsetLey[,5:7], var1, var2,
                     probe.set.name=subsetLey[,1])
    
    ## Computation of p-values
    checkEquals(
        applyLpeMethod(subsetLey[,2:4], subsetLey[,5:7], 
                       rowNames=subsetLey[,1])$z.stats,
        lpeResult$z.stats
    )
    
    ## TODO check names of elements before / after
    
}

test_rep_applyAndReportLPE <- function() {

    ## Setup data
    data(Ley)
    Ley[,2:7] <- preprocess(Ley[,2:7],data.type="MAS5")
    subsetLey <- Ley[1:1000,]
    
    ## Test
    applyAndReportLPE(
        subsetLey[,2:4], subsetLey[,5:7], 
        file.path(getwd(), 'temp'),
        stdout(),  # file.path(getwd(), 'temp', 'out_test_applyAndReportLPE.Rmd'), 
        0.05)[["p.values"]]
    
}

