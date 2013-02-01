#' @name SignificanceTests
#' @title Significance tests
#' @description 
#' List of Significance tests
#' 
#' @keywords Significance
#' 
#' @section Introduction:
#' Testing a null hypothesis which is "each protein is not differencially 
#' expressed". 
#' The obtained p-values are probabilities of an observation (or more extreme 
#' observations), supposing that the null hypothesis is true under a set of 
#' assumptions (like for instance normality, homoscedasticity, and random 
#' sampling for t test).
#' NB: these tests should be followed by multiple testing correction.
#' 
#' Note also that the use of classical statistical tests has been prone to high 
#' criticism. Modern statistics propose good alternatives including Bayesian 
#' statistics (see Limma package).
#' 
#' @section References:
#' 
#' Johnson, D. (1999). 
#' The Insignificance of Statistical Significance Testing. 
#' The journal of wildlife management, 63(3), 763. 
#' doi:10.2307/3802789
#' 
#' Erceg-Hurn, D. M., & Mirosevich, V. M. (2008). 
#' Modern robust statistical methods: an easy way to maximize the accuracy and 
#' power of your research. 
#' The American psychologist, 63(7), 591-601. 
#' doi:10.1037/0003-066X.63.7.591
#' 
#' Jeanmougin, M., De Reynies, A., Marisa, L., Paccard, C., Nuel, G., & 
#' Guedj, M. (2010). 
#' Should we abandon the t-test in the analysis of gene expression microarray 
#' data: a comparison of variance modeling strategies. 
#' PloS one, 5(9), e12336. 
#' doi:10.1371/journal.pone.0012336
#' 
#' http://en.wikipedia.org/wiki/Statistical_hypothesis_testing
#' http://en.wikipedia.org/wiki/Statistical_hypothesis_testing#Controversy
#' http://cseweb.ucsd.edu/users/goguen/courses/275f00/stat.html
#' 
sapply(
    list.files(file.path(getwd(), 'R'), 
               pattern="stt_[^a][^l]{2}.+\\.R", 
               full.names=TRUE),
    source
)

#' @title Test for non-paired data.
#' 
#' @description Test for differentially expressed proteins on non-paired data.
#' 
#' @details
#' Apply statistical tests on non-paired data to find differentially 
#' expressed proteins. Stop if given an incorrect method.
#' 
#' @param signMethod A label indicating which method should be applied.
#'   Available methods include:
#'    ttest : t-test
#'    lpe : local-pooled-error method
#'    sam : significance analysis of microarrays
#'    limma : linear models for microarray analysis
#' @param experiment A dataframe with normalized intensities.
#' @param control A dataframe with normalized intensities.
#' @param outputFolderTemp Where the temporary files will be created
#' @param outputFileNameRmd Where report paragraphs will be written
#' @param thresholdPVal FDR threshold
#' @return A data frame with data (p.values, fold.change, 
#'   significant).
applySummaryStatistic <- function(signMethod, 
                                  experiment, control, 
#' @export
                                  outputFolderTemp, outputFileNameRmd, 
                                  thresholdPVal) {
    outPVal <- switch(
        signMethod,
        ttest = applyAndReportTTests(experiment, control, 
                                     outputFolderTemp, outputFileNameRmd, 
                                     thresholdPVal, isPaired=FALSE),
        lpe = applyAndReportLPE(experiment, control, 
                                outputFolderTemp, outputFileNameRmd, 
                                thresholdPVal),
        sam = applyAndReportSam(experiment, control, 
                                outputFolderTemp, outputFileNameRmd, 
                                thresholdPVal, isPaired=FALSE),
        limma = applyAndReportLimma(experiment, control, 
                                    outputFolderTemp, outputFileNameRmd, 
                                    thresholdPVal, isPaired=FALSE),
        stop(
            paste0(
                c(
                    "Unknown method in applySummaryStatistic call : \'", 
                    signMethod,
                    "\'."
                ), 
                sep=''), 
            call=TRUE)
    )
    return(outPVal)
}

#' @title Test for paired data.
#' 
#' @description Test for differentially expressed proteins on paired data.
#' 
#' @details
#' Apply statistical tests on paired data to find differentially 
#' expressed proteins. Stop if given an incorrect method.
#' Note that LPE method is not included: use instead PLPE package (to be 
#' implemented).
#' 
#' @param signMethod A label indicating which method should be applied.
#'   Available methods include:
#'    ttest : t-test
#'    sam : significance analysis of microarrays
#'    limma : linear models for microarray analysis
#' @param experiment A dataframe with normalized intensities.
#' @param control A dataframe with normalized intensities.
#' @param outputFolderTemp Where the temporary files will be created
#' @param outputFileNameRmd Where report paragraphs will be written
#' @param thresholdPVal FDR threshold
#' @return A data frame with data (p.values, fold.change, 
#'   significant).
applyPairedSummaryStatistic <- function(signMethod, 
                                        experiment, control, 
                                        outputFolderTemp, outputFileNameRmd, 
                                        thresholdPVal) {
    outPVal <- switch(
        signMethod,
        ttest = applyAndReportTTests(experiment, control, 
                                     outputFolderTemp, outputFileNameRmd, 
                                     thresholdPVal, isPaired=TRUE),
        sam = applyAndReportSam(experiment, control, 
                                outputFolderTemp, outputFileNameRmd, 
                                thresholdPVal, isPaired=TRUE),
        limma = applyAndReportLimma(experiment, control, 
                                    outputFolderTemp, outputFileNameRmd, 
                                    thresholdPVal, isPaired=TRUE),
        stop(
            paste0(
                c(
                    "Unknown method in applySummaryStatistic call : \'", 
                    signMethod,
                    "\'."
                ), 
                sep=''), 
            call=TRUE)
    )
    
    return(outPVal)
}
