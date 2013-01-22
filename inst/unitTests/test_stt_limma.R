
test_fun_applyLimma <- function() {
    
    normalData1 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    normalData2 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    
    fit <- lmFit(
        cbind(normalData1, normalData2),
        method = "ls",
        design = cbind(Grp1=1, 
                       Grp2vs1=c(rep(0,ncol(normalData1)), 
                                 rep(1,ncol(normalData2)))
        )
    )
    fit <- eBayes(fit)
    limmaVals <- topTable(fit, number=nrow(normalData2), coef=2)
    limmaVals <- limmaVals[order(as.numeric(row.names(limmaVals))),]
    
    checkEquals(
        applyLimma(normalData2, normalData1)$P.Value,
        limmaVals$P.Value
    )
    
    
}

test_fun_applyPairedLimma <- function() {
    
    normalData1 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    normalData2 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    
    Replicates <- factor(rep(c(1:3), times=2))
    Treat <- factor(rep(c("C", "T"), each=3), levels=c("C","T"))

    fit <- lmFit(
        cbind(normalData1, normalData2),
        method = "ls",
        design = model.matrix(~Replicates+Treat)
    )
    fit <- eBayes(fit)
    limmaVals <- topTable(fit, number=500, coef="TreatT")
    limmaVals <- limmaVals[order(as.numeric(row.names(limmaVals))),]
    
    checkEquals(
        applyPairedLimma(normalData2, normalData1)$P.Value,
        limmaVals$P.Value
    )
    
    
}

test_rep_applyAndReportLimma <- function() {
    
    ## Setup data
    normalData1 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    normalData2 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    
    ## Test
    applyAndReportLimma(
        normalData2, normalData1, 
        file.path(getwd(), 'temp'),
        stdout(),
        0.05)[["p.values"]]
    
}

test_rep_reportLimma <- function() {
    
    ## Setup data
    x <- rnorm(50, mean=c(rep(0,25), rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1)
    )
    
    ## Test
    reportLimma(
        matrixData,
        0.05,
        stdout(),  ##file.path(getwd(), 'temp', 'out_test_reportTTest.Rmd'),
        file.path(getwd(), 'temp')
    )
    
}
