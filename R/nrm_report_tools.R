
#' @title Report results of normalization
#' 
#' @description Report results of normalization (knitr)
#' 
#' @details 
#' Write a summary report of normalization on a dataset. 
#' Suppose a 2*nb experimental design, as level for displayNormalizationMAplot 
#' call...
#' Data at this step will be saved in a temporary file.
#' TODO Use ExpressionSet object as input?
#' 
#' @param matrixData Dataset to check for normality.
#'  expDesign Experimental design for the dataset, categories.
#' @param outputFile Report to be extended.
#' @param outFolder Where the temporary file will be saved.
allTestsNormalizationRmd <- function(matrixData, 
                                     outputFile, outFolder) {
    
    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc( 
            runif(1) * 10000)), 
        collapse='')
    
    ## Write data in a file usable later when the report will be transformed 
    ## into HTML
    tempOutput <- paste(
        c(outFolder, '/all_tests_normalization_Rmd_data_', execLabel, '.txt'), 
        collapse='')
    write.table(matrixData, tempOutput, sep="\t")
    
    cat(' ',
        paste(
            c('```{r allTestsNormalizationRmd', 
              execLabel, ', echo=FALSE, fig.width=12, fig.height=6}'),
            collapse=''),
        ' ',
        sep="\n", file=outputFile, append=TRUE)
    
    cat('displayNormalizationViolin <- ', file=outputFile, append=TRUE)
#     print(displayNormalizationViolin)
    cat(paste(deparse(displayNormalizationViolin), collapse="\n"), '',
        sep="\n", file=outputFile, append=TRUE)
    cat('displayNormalizationQQplot <- ', file=outputFile, append=TRUE)
#     print(displayNormalizationQQplot)
    cat(paste(deparse(displayNormalizationQQplot), collapse="\n"), '',
        sep="\n", file=outputFile, append=TRUE)
    cat('displayNormalizationMAplot <- ', file=outputFile, append=TRUE)
#     print(displayNormalizationMAplot)
    cat(paste(deparse(displayNormalizationMAplot), collapse="\n"), ' ',
        sep="\n", file=outputFile, append=TRUE)
    
    cat('',
        paste(
            c('layout(matrix(c(', paste(c(1:ncol(matrixData)), collapse=', '), 
              '), 1,', ncol(matrixData), ', byrow=TRUE))'), 
            collapse=''),
        paste(c('matrixData <- as.matrix(read.table("',
                tempOutput,
                '", stringsAsFactors=FALSE))'), 
              collapse=''),
        'displayNormalizationViolin(matrixData)',
        'tempQQ <- apply(matrixData, 2, FUN=displayNormalizationQQplot)',
        'nbColumnsPerExp <- ncol(matrixData)/2',
        'dataToTest1 <- apply(matrixData[, 1:nbColumnsPerExp], 1, mean)',
        'dataToTest2 <- apply(matrixData[, (nbColumnsPerExp + 1):(2 * nbColumnsPerExp)], 1, mean)',
        'displayNormalizationMAplot(dataToTest2, dataToTest1)',
        '',
        '```',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
}

## -----------------------------------------------------------------------------
## Tests / plots used in the report

#' @title Violin plot
#' 
#' @description Violin plots using ggplot2
#' 
#' @details 
#' Display violin plots using the ggplot2 library.
#' TODO Provide exp design instead of supposing it
#' 
#' @param matrixData A data.frame containing the p-values (3 columns necessary: 
#'  "p.values", "p.values.corrected" and "fold.change")
displayNormalizationViolin <- function(matrixData) {
    library("ggplot2")
    
    # Stack data by pasting column after column
    dataStacked <- c(matrixData)
    
    rowNames <- row.names(matrixData)
    colNames <- colnames(matrixData)
    numProts <- length(rowNames)
    numCond <- 2
    numReplicates <- ncol(matrixData) / numCond
    numAllExp <- numCond * numReplicates
    
    intensitiesDF <- data.frame(
        Intensities = dataStacked,
        Protein.groups = factor(rep(rowNames, numAllExp)),
        Reporters = factor( rep(colNames, rep(c(numProts), numAllExp)) ),
        Conditions = factor( 
            rep(c("Control", "Exp."), 
                c(numProts * numReplicates, numProts * numReplicates)) )
        )
    
    p <- ggplot(data=intensitiesDF, aes(x=Reporters, y=Intensities)) 
    p + geom_violin(scale = "count", mapping=aes(fill = Conditions), 
                    trim=TRUE, alpha=.50) +
        stat_summary(aes(group=Reporters), 
                     fun.y=mean, fun.ymin=min, fun.ymax=max, 
                     fill="red", shape=21, size=1) +
        ylab('glog2( Intensities )') +
        xlab('') +
        theme(legend.position="bottom")
}

#' Both vectors should be log'ed before!
displayNormalizationMAplot <- function(red1d, green1d) {
    
    aVals <- (red1d + green1d)/2
    mVals <- red1d - green1d
    dat <- data.frame(aVals, mVals)
    
    sdFoldchange <- sd(mVals)
    
#     plot(aVals, mVals, pch=16)
    ggplot(dat, aes(x=aVals, y=mVals)) + 
        geom_point(shape=16) + 
        geom_abline(slope = 0, intercept = 0)
    
}

displayNormalizationQQplot <- function(vector1d) {
    qqnorm(vector1d, main="", ylab='')
    qqline(vector1d)
}

###############################################################################
## Non used :

# #' Computes the mean for each row in order to display histogram and QQplot
# allTestsNormalization <- function(matrixData, title, outFolder) {
#     
#     #Graphics outputs
# #     pdf(paste( c(outFolder, '/all_tests_normalization.pdf'), collapse='') )
#     jpeg(paste(c(outFolder, '/all_tests_normalization.jpg'), collapse=''))
#     layout(
#         matrix(c(1, 2, 3), 1, 3, byrow=TRUE)
#     )
#     
#     displayNormalizationBoxplot(matrixData, title=title)
#     
#     dataToTest <- apply(matrixData, 1, mean)
#     displayNormalizationHistogram(dataToTest)
#     displayNormalizationQQplot(dataToTest)
#     
#     dev.off()
#     
#     #Non-graphics output
#     testNormalizationShapiro(dataToTest, title)
#     
# }

# displayNormalizationBoxplot <- function(matrixData, title) {
#     boxplot(matrixData, main=title)
# }

# displayNormalizationDensities <- function(matrixData) {
#     
#     for(name in names(matrixData)) {
#         vector1d <- dataset[,name]
#         plot(density(vector1d), main=name)
# 
#         # Print mean and standard deviation
#         meanv <- mean(vector1d)
#         sdv <- sd(vector1d)
#         abline(v=c(meanv-sdv,meanv,meanv+sdv), col="green")
#     }
# 
# }

# displayNormalizationHistogram <- function(vector1d, x_label = "Dataset") {
# 
#     #Display histogram of values
#     histogram <- hist(vector1d, breaks=100, col="red", xlab=x_label, main="Histogram with Normal Curve")
#     
#     # Compute Normal Curve (Thanks to Peter Dalgaard)
#     xfit <- seq( min(vector1d, na.rm = TRUE), max(vector1d, na.rm = TRUE), length=3000)
#     yfit <- dnorm( xfit, mean=mean(vector1d, na.rm = TRUE), sd=sd(vector1d, na.rm = TRUE) )
#     yfit <- yfit * diff( histogram$mids[1:2]) * length(vector1d)
#     
#     # Print Normal Curve
#     lines(xfit, yfit, type="l", col="blue", lwd=2)
#     
#     # Print mean and standard deviation
#     meanv <- mean(vector1d)
#     sdv <- sd(vector1d)
#     abline(v=c(meanv-sdv,meanv,meanv+sdv), col="green")
# 
#     
# }



##########################################################
# Statistical tests
# Normality
#
#There are 5 major tests used:
#    Shapiro-Wilk W test
#Anderson-Darling test
#Martinez-Iglewicz test
#Kolmogorov-Smirnov test
#D’Agostino Omnibus test
#NB: Power of all is weak if N < 1


#' Shapiro-Wilk W test
#'Developed by Shapiro and Wilk (1965).
#'One of the most powerful overall tests.
#'It is the ratio of two estimates of variance (actual
#'calculation is cumbersome; adversely affected by ties).
#'Test statistic is W; roughly a measure of the straightness
#'of the quantile-quantile plot.
#'The closer W is to 1, the more normal the sample is.
# test_normalization_shapiro <- function(vector1d) {
#     stest <- shapiro.test(vector1d)
# }

#' Kolmogorov-Smirnov Test
#'Calculates expected normal distribution and compares
#'it with the observed distribution.
#'Uses cumulative distribution functions.
#'Based on the max difference between two distributions.
#'Poor for discrimination below N = 30.
#'Power to detect differences is low.
#'Historically popular.
# test_normality_kolmogorov_smirnov <- function(vector1d) {
#     #ks.test(F500, "pnorm", mean=a[1], sd=a[2])
# }