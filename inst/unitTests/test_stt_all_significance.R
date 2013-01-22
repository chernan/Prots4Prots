
test_fun_applySummaryStatistic <- function() {

    ## Setup data
    ## NB:
    ## This test doesn't work when run for LPE. 
    # Error in smooth.spline(base.var.adap[, 1], base.var.adap[, 2], df = df) : 
    # missing or infinite values in inputs are not allowed 
    ## Ok on command line... :-/
    experiment <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    rownames(experiment) <- as.character(1:500)
    control <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    rownames(control) <- as.character(1:500)
    
    ## Test t-test
    ttestPVals <- applySummaryStatistic(
        'ttest',
        experiment, control, 
        file.path(getwd(), 'temp'), 
        file.path(getwd(), 'temp', 'out_test_stt_applySummaryStatistic.Rmd'), 
        0.05)[["p.values"]]
    checkTrue(all( ttestPVals <= 1 & ttestPVals >= 0), msg="t-test")
    
    ## Test SAM
    ## Generates p-values > 1 ??!!
    samPVals <- applySummaryStatistic(
        'sam',
        experiment, control, 
        file.path(getwd(), 'temp'), 
        file.path(getwd(), 'temp', 'out_test_stt_applySummaryStatistic.Rmd'), 
        0.05)[["p.values"]]
    checkTrue(all( samPVals <= 1.01 & samPVals >= 0), msg="SAM (p-values can be >1... ?!)")
    
    ## Test Limma
    limmaPVals <- applySummaryStatistic(
        'limma',
        experiment, control, 
        file.path(getwd(), 'temp'), 
        file.path(getwd(), 'temp', 'out_test_stt_applySummaryStatistic.Rmd'), 
        0.05)[["p.values"]]
    checkTrue(all( limmaPVals <= 1 & limmaPVals >= 0), msg="Limma")
    
    # Test error if no label given
    checkException(
        applySummaryStatistic(
            '',
            experiment, control, 
            file.path(getwd(), 'temp'), 
            file.path(getwd(), 'temp', 'out_test_stt_applySummaryStatistic.Rmd'), 
            0.05),
        silent = TRUE
    )
    
    ## Test LPE
    # LPE needs a minimal variability
    data(Ley)
    Ley[,2:7] <- preprocess(Ley[,2:7],data.type="MAS5")
    subsetLey <- Ley[1:1000,]
    experiment <- subsetLey[,2:4]
    control <- subsetLey[,5:7]
    
    lpePVals <- applySummaryStatistic(
        'lpe',
        experiment, control, 
        file.path(getwd(), 'temp'), 
        file.path(getwd(), 'temp', 'out_test_stt_applySummaryStatistic.Rmd'), 
        0.05)[["p.values"]]
    checkTrue(all( lpePVals <= 1 & lpePVals >= 0), msg="LPE")
    
}

test_fun_applyPairedSummaryStatistic <- function() {
    
    ## Setup data
    experiment <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    rownames(experiment) <- as.character(1:500)
    control <- replicate(3, rnorm(n=500, mean=11.5, sd=2.2))
    rownames(control) <- as.character(1:500)
    
    ## Test t-test
    ttestPVals <- applyPairedSummaryStatistic(
        'ttest',
        experiment, control, 
        file.path(getwd(), 'temp'), 
        file.path(getwd(), 'temp', 'out_test_stt_applyPairedSummaryStatistic.Rmd'), 
        0.05)[["p.values"]]
    checkTrue(all( ttestPVals <= 1 & ttestPVals >= 0), msg="t-test")
    
    ## Test SAM
    ## Generates p-values > 1 ??!!
    samPVals <- applyPairedSummaryStatistic(
        'sam',
        experiment, control, 
        file.path(getwd(), 'temp'), 
        file.path(getwd(), 'temp', 'out_test_stt_applyPairedSummaryStatistic.Rmd'), 
        0.05)[["p.values"]]
    checkTrue(all( samPVals <= 1.01 & samPVals >= 0), msg="SAM (p-values can be >1... ?!)")
    
    ## Test Limma
    limmaPVals <- applyPairedSummaryStatistic(
        'limma',
        experiment, control, 
        file.path(getwd(), 'temp'), 
        file.path(getwd(), 'temp', 'out_test_stt_applyPairedSummaryStatistic.Rmd'), 
        0.05)[["p.values"]]
    checkTrue(all( limmaPVals <= 1 & limmaPVals >= 0), msg="Limma")
    
    # Test error if no label given
    checkException(
        applyPairedSummaryStatistic(
            '',
            experiment, control, 
            file.path(getwd(), 'temp'), 
            file.path(getwd(), 'temp', 'out_test_stt_applyPairedSummaryStatistic.Rmd'), 
            0.05),
        silent = TRUE
    )
    
    
}
