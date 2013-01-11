#' t-test
#' 
#' A t statistic is the distance between two data sets in units of standard 
#' deviation.
#' This test has low statistical power when applied to small sample size.
#' 
#' Student's t-test relies on certain assumptions. Welch's t-test can be used 
#' as an alternative (see paragraph : 2/ 'Homoscedasticity, or equality of 
#' variance between samples')
#' 
#' A two-sample t-test can be unpaired or paired.
#' A comparison between paired samples has greater power, because of a similar 
#' noise factor in paired samples.
#' 
#' Assumptions
#' 
#' 1/ Normality
#' 
#' The two samples should follow a normal distribution under the null 
#' hypothesis (e.g if there is no difference betwween them). 
#' This can be tested by Shapiro-Wilk (parametric) or Kolmogorov–Smirnov 
#' (non-parametric) tests, or assessed by a normal quantile plot.
#' 
#' 2/ Homoscedasticity, or equality of variance between samples
#' 
#' This assumption is often neglected.
#' In case of difference in variance between samples (heteroscedasticity), the 
#' probability of obtaining a "significant" result may be greater than the 
#' desired alpha level (even if the null hypothesis is true). 
#' In other words, using the t-test when the variances are not equal will 
#' result in a rate of false positive superior to the desired value.
#' 
#' To avoid that problem, Welch's t-test can be used instead of Student's t-test. 
#' In this case the variance is calculated separately for each sample, and a 
#' modified degree of freedom is used. But sample sizes have to be the same.
#' 
#' Some tests can evaluate if the difference in variance between samples is 
#' significant or not.
#' - Bartlett : If the distributions are nearly normal, otherwise it may just be 
#'   testing for non-normality (being sensitive to departures from normality).
#' - Levene : If samples are small, or data are not normal (or you don’t know). 
#'   Relies on the mean.
#' - Brown–Forsythe : Non-parametric alternative to the Levene test, relies on 
#'   either the median or the trimmed mean. To be used in case of strong 
#'   evidence of departure from normality (after a Q-Q plot, for instance).
#' - F test : Same as Bartlett's test. Extremely sensitive to departures from 
#'   normality.
#' - Q-Q plot : For a graphical assessment, but preferably if samples are nor 
#'   small.
#' 
#' 3/ Independence
#' 
#' Data should be sampled independently from the two populations being compared.
#' This is in general not testable, and need to understand the performed 
#' experimental design.
#' 
#' NB: use of the t-test requires no citation. But it should be mentionned 
#' whether Student's or Welch's test was used.
#' 
#' @seealso http://en.wikipedia.org/wiki/Student%27s_t-test


#' @title Wrapper to R Base t.test()
#' 
#' @description Wrapper to R Base t.test()
#' 
#' @details 
#' Using \code{t.test} function available in R Base to perform Student's or 
#' Welch's t test, depending on the value of isHomoscedastic.
#' 
#' @param xExp Dataset for experiment (condition 1).
#' @param xControl Dataset for control (condition 2).
#' @param confidence Desired confidence level (default to 0.95, as t.test)
#' @param isHomoscedastic If set to FALSE, Welch's test will be performed, 
#'  otherwise Student's test will be used (Welch's is the default in t.test). 
#' @param isPaired If the samples are paired (ex: sample 1,2,3 in condition 1 vs 
#'  sample 1,2,3 in condition 2)
#' @return The test result object as described in R documentation
apply2SampleTTest <- function(xExp, xControl, confidence=0.95, 
                              isHomoscedastic=FALSE, isPaired=FALSE) {
    
    testResult <- t.test(x=xExp, y=xControl, na.action=na.omit, 
                         var.equal=isHomoscedastic, conf.level=confidence, 
                         paired=isPaired)
    
    return (testResult)
}


#' @title Apply multiple Welch's t-tests
#' 
#' @description Apply Welch's t-tests between rows of two datasets
#' 
#' @details 
#' Apply Welch's t-tests between rows of two datasets. 
#' As the fold change will also be computed, data must be log2-transformed.
#' TODO Change Bartlett to Levene variance test
#' TODO Test variance between *rows* + *MTC*
#' TODO Use ExpressionSet object as input
#' TODO Return p-values of normality and variance tests!
#' 
#' @param experiment Dataset for experiment (condition 1).
#' @param control Dataset for control (condition 2).
#' @param significance Maximum threshold for the p-values (e.g 0.05)
#' @param isHomoscedastic If set to FALSE, Welch's test will be performed, 
#'  otherwise Student's test will be used. 
#' @param isPaired If the samples are paired (ex: sample 1,2,3 in condition 1  
#'  vs sample 1,2,3 in condition 2)
#' @return A data frame with two columns "pval" p-value of the t-test, 
#'  and "foldchange" computed as either the difference between the means 
#'  (non-paired) or the mean of the differences (paired))
applyTTests <- function(experiment, control, rowNames, significance, isPaired) {
    
    nbRow <- nrow(experiment)
    confidence <- (1 - significance)
    
    #     stackedData <- stack(data.frame(cbind(experiment,control), row.names=rowNames))
    #     isHomoscedastic <- (bartlett.test(stackedData$values, 
    #                                        g=stackedData$ind)$p.value < significance)
    
    resultT <- t(
        apply(matrix(c(1:nbRow), ncol=1), 1, 
              FUN = function(x) { 
                  tTest <- apply2SampleTTest(
                      as.matrix(experiment[x, ]), 
                      as.matrix(control[x, ]), 
                      confidence, FALSE, isPaired)
                  if(!isPaired) {
                      return(c(
                          pval=tTest$p.value, 
                          foldchange=tTest$estimate[["mean of x"]]-tTest$estimate[["mean of y"]]
                      ))
                  } else {
                      return(c(
                          pval=tTest$p.value, 
                          foldchange=tTest$estimate[["mean of the differences"]]
                      ))
                  }
              }
        )
    )
    resultT <- data.frame(resultT, row.names=rowNames)
    return(resultT)
}

#' @title Apply multiple Welch's t-tests and report result
#' 
#' @description Apply Welch's t-tests between rows of two datasets
#' 
#' @details 
#' Apply Welch's t-tests between rows of two datasets. 
#' TODO Change Bartlett to Levene variance test
#' TODO Test variance between *rows* + *MTC*
#' TODO Use ExpressionSet object as input
#' TODO Return p-values of normality and variance tests!
#' 
#' @param experiment Dataset for experiment (condition 1).
#' @param control Dataset for control (condition 2).
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @param thresholdPVal Maximum threshold for the p-values (e.g 0.05)
#' @param isPaired If the samples are paired (ex: sample 1,2,3 in condition 1  
#'  vs sample 1,2,3 in condition 2). Default to FALSE
#' @return A data frame with two columns "p.values" "fold.change", in the 
#'  same order as the input data
applyAndReportTTests <- function(experiment, control, outputFolderTemp, 
                                 outputFile, thresholdPVal, isPaired=FALSE) {
    
    datasetTTest <- applyTTests(experiment, control, 
                                rownames(experiment), 
                                thresholdPVal, isPaired)
    
    matrixData <- data.frame(
        p.values = datasetTTest$pval,
        fold.change = datasetTTest$foldchange,
        row.names = rownames(experiment)
    )
    
    #     stacked.data <- stack(data.frame(cbind(experiment,control), row.names=rownames(experiment)))
    #     is.homoscedastic <- (bartlett.test(stacked.data$values, g=stacked.data$ind)$p.value < significance)
    
    #Erreur dans shapiro.test(vector1d) : sample size must be between 3 and 5000
    #     if(nrow(experiment) > 3 & nrow(experiment) < 5000) {
    #         is.norm.by.sample <- apply(cbind(experiment,control), 2, FUN=function(vector1d){return(shapiro.test(vector1d)$p.value)} )
    #     }else {
    #         is.norm.by.sample <- NULL
    #     }
    
    reportTTest(matrixData, thresholdPVal, outputFile, outputFolderTemp, 
                FALSE, NULL, isPaired)
    
    return(matrixData)
}

#' @title Apply multiple Welch's t-tests and report result
#' 
#' @description Apply Welch's t-tests between rows of two datasets
#' 
#' @details 
#' Apply Welch's t-tests between rows of two datasets. 
#' 
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @param matrixData A data frame with two columns "p.values" "fold.change"
#' @param thresholdPVal Maximum threshold for the p-values (e.g 0.05)
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param isHomoscedastic deprecated FALSE, Welch's test was performed
#' @param isNormBySample deprecated NULL
#' @param isPaired If the samples were paired (ex: sample 1,2,3 in condition 1  
#'  vs sample 1,2,3 in condition 2).
reportTTest <- function(matrixData, thresholdPVal, outputFile, outputFolderTemp, 
                          isHomoscedastic, isNormBySample, isPaired) {
    
    cat('',
        'Significantly different proteins',
        '---------------------------------------------------------------------',
        '',
        '',
        sep="\n", file=outputFile, append=TRUE)
    cat('Data analysis was performed using t-test ',
        ifelse(isPaired, 
               'for paired samples ', 'for non-paired samples '),
        ifelse(isHomoscedastic == TRUE, '
               (equal variance - Student)', '(unequal variance - Welch)'),
        '.',
        "\n",
        "\n",
        sep='', file=outputFile, append=TRUE)
    
    #     if(!all(sapply(is.norm.by.sample, is.null))) {
    #         cat('Check normality of input data using Shapiro\'s test, for each sample :',
    #             '',
    #             'Experiment | Control',
    #             '------------- | -------------',
    #             apply( matrix(format.pval(is.norm.by.sample, digits=4),
    #                           ncol=2, byrow=FALSE) ,
    #                    1,
    #                    FUN=function(x) {paste(x,collapse=' | ',sep="\n")}),
    #             '',
    #             sep="\n")
    #     }
    
    allTestsSignificanceRmd(matrixData, thresholdPVal, 
                            outputFile, outputFolderTemp)
    cat('',
        '> t-test : require no citation.',
        '>',
        '',
        '---------------------------------------------------------------------',
        '',
        sep="\n", file=outputFile, append=TRUE)        
    
}





