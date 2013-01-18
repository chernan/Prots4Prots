
test_fun_applySam <- function() {
    
    set.seed(0)
    control <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    experiment <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    
    data <- as.matrix(data.frame(control, experiment))
    rowNames <- rownames(experiment)
    row.names(data) <- rowNames
    
    design <- c(rep(1, ncol(control)), rep(2, ncol(experiment)))
    
    samFit <- SAM(x=data, y=design, 
                  logged2=TRUE, 
                  resp.type="Two class unpaired", fdr.output=0.05,
                  genenames=rowNames, geneid=rowNames)
    
    pValues <- samr.pvalues.from.perms(samFit$samr.obj$tt, 
                                       samFit$samr.obj$ttstar)
    
    ## As SAM applies randomosation, p-values are never exactly the same...
    checkEqualsNumeric(
        applySam(experiment, control, 0.05)$p.values,
        pValues,
        tolerance=0.005
    )
    
}

test_fun_applyPairedSam <- function() {
    
    set.seed(0)
    control <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    experiment <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    
    data <- as.matrix(data.frame(control, experiment))
    rowNames <- rownames(experiment)
    row.names(data) <- rowNames
    
    design <- rep(c(1:ncol(control)), times=2)
    design[1:ncol(control)] <- design[1:ncol(control)] * -1
    
    samFit <- SAM(x=data, y=design, 
                  logged2=TRUE, 
                  resp.type="Two class paired", fdr.output=0.05,
                  genenames=rowNames, geneid=rowNames,
    )
        
    pValues <- samr.pvalues.from.perms(samFit$samr.obj$tt, 
                                       samFit$samr.obj$ttstar)
    
    ## As SAM applies randomosation, p-values are never exactly the same...
    checkEqualsNumeric(
        applyPairedSam(experiment, control, 0.05)$p.values,
        pValues,
        tolerance=0.005
    )
    
}


test_rep_applyAndReportSam <- function() {
    
    ## Setup data
    set.seed(0)
    normalData1 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    rownames(normalData1) <- as.character(1:500)
    normalData2 <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    rownames(normalData2) <- as.character(1:500)
    
    ## Test
    applyAndReportSam(
        normalData2, normalData1, 
        file.path(getwd(), 'temp'),
        stdout(),
        0.05)[["p.values"]]
    
}

test_rep_reportSam <- function() {
    
    ## Setup data
    set.seed(0)
    x <- rnorm(50, mean=c(rep(0,25), rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))

    significant <- matrix(rep('STABLE', 50))
    significant[1:10] <- 'LO'
    significant[11:20] <- 'UP'
    
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1),
        significant = significant
    )
    
    ## Test
    reportSam(
        matrixData,
        0.05,
        stdout(),  ##file.path(getwd(), 'temp', 'out_test_reportTTest.Rmd'),
        file.path(getwd(), 'temp')
    )
    
}

test_gra_displaySamSignificanceVolcanoPlot <- function() {
    
    ## Setup data
    set.seed(0)
    x <- rnorm(50, mean=c(rep(0,25),rep(3,25)))
    p <- 2*pnorm( sort(-abs(x)))

    significant <- matrix(rep('STABLE', 50))
    significant[1:10] <- 'LO'
    significant[11:20] <- 'UP'
    
    matrixData <- data.frame(
        p.values = p,
        fold.change = rnorm(50, mean=0, sd=1),
        significant = significant
    )
    
    ## Run function... as test?
    displaySamSignificanceVolcanoPlot(matrixData, 0.05, 
                                      title="test_displaySamCorrectedVolcanoPlot")
    displaySamSignificanceVolcanoPlot(matrixData, 0.05)
    
}
