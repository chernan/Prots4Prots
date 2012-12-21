
test_fun_apply2SampleTTest <- function() {
    
    ## Inspired from tests in R base documentation for t.test()
    checkEquals(
        apply2SampleTTest(1:10, 7:20)$p.value, 
        apply2SampleTTest(1:10, 7:20, 0.95, FALSE, FALSE)$p.value
    )
    checkEquals(
        apply2SampleTTest(1:10, 7:20)$p.value, 
        t.test(1:10, y=c(7:20))$p.value
    )
    checkEquals(
        apply2SampleTTest(1:10, 7:20, 0.95, FALSE, FALSE)$p.value, 
        t.test(1:10, y=c(7:20))$p.value
    )

    checkEquals(
        apply2SampleTTest(1:10, c(7:20, 200))$p.value, 
        t.test(1:10, y=c(7:20, 200))$p.value
    )
    checkEquals(
        apply2SampleTTest(1:10, c(7:20, 200), 0.95, FALSE, FALSE)$p.value, 
        t.test(1:10, y=c(7:20, 200))$p.value
    )
    
}

test_fun_applyTTests <- function() {
    
    set.seed(0)
    normalData1 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    normalData2 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    
    ## Computation of p-values
    checkEqualsNumeric(
        applyTTests(normalData1, normalData2, 1:500, 0.05, FALSE)$pval,
        sapply(1:500, 
               FUN = function(x) { 
                   return(t.test(normalData1[x, ], normalData2[x, ])$p.value)
               }
        )
    )
    
    ## Computation of fold changes - non paired
    checkEqualsNumeric(
        applyTTests(normalData1, normalData2, 1:500, 0.05, FALSE)$foldchange,
        sapply(1:500, 
               FUN = function(x) { 
                   tTest <- t.test(normalData1[x, ], normalData2[x, ], 
                                   paired=FALSE)
                   return(tTest$estimate[["mean of x"]]-tTest$estimate[["mean of y"]])
               }
        )
    )
    
    ## Computation of fold changes - paired
    checkEqualsNumeric(
        applyTTests(normalData1, normalData2, 1:500, 0.05, TRUE)$foldchange,
        sapply(1:500, 
               FUN = function(x) { 
                   tTest <- t.test(normalData1[x, ], normalData2[x, ], 
                                   paired=TRUE)
                   return(tTest$estimate[["mean of the differences"]])
               }
        )
    )
    
    ## TODO check order of elements before / after
    
}

test_fun_applyAndReportTTests <- function() {

    ## Setup data
    set.seed(0)
    normalData1 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    normalData2 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
        
    ## Test
    checkEqualsNumeric(
        applyAndReportTTests(
            normalData1, normalData2, 
            file.path(getwd(), 'temp'),
            file.path(getwd(), 'temp', 'out_test_applyAndReportTTests.Rmd'), 
            0.05)[["p.values"]],
        sapply(1:500, 
               FUN = function(x) { 
                   return(t.test(normalData1[x, ], 
                                 normalData2[x, ], paired=FALSE)$p.value)
               }
        )
    )
    
}

test_rep_reportTTest <- function() {
    
    set.seed(0)
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1)
    )

    reportTTest(
        matrixData,
        0.05,
        stdout(),  ##file.path(getwd(), 'temp', 'out_test_reportTTest.Rmd'),
        file.path(getwd(), 'temp'),
        FALSE, NULL, FALSE
    )
}