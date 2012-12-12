


#' Computes the mean for each row in order to display histogram and QQplot
all_tests_normalization <- function(matrixdata, title, outfolder) {
    
    #Graphics outputs
#     pdf( paste( c(outfolder,'/all_tests_normalization.pdf'), collapse='') )
    jpeg(paste(c(outfolder,'/all_tests_normalization.jpg'), collapse=''))
    layout(
        matrix(c(1,2,3),1,3,byrow=TRUE)
    )
    
    display_normalization_boxplot(matrixdata, title=title)
    
    data.to.test <- apply(matrixdata, 1, mean)
    display_normalization_histogram(data.to.test)
    display_normalization_qqplot(data.to.test)
    
    dev.off()
    
    #Non-graphics output
    test_normalization_shapiro(data.to.test, title)
    
}

#' Suppose a 2* nb experimental design, level of display_normalization_MAplot call
all_tests_normalization_Rmd <- function(matrixdata, title, outfolder) {
    tempoutput <- paste(c(outfolder, '/all_tests_normalization_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(matrixdata, tempoutput, sep="\t") 
    
    cat('
```{r, echo=FALSE, fig.width=12, fig.height=6}
')
    
    cat('display_normalization_violin <- ')
    print(display_normalization_violin)
    cat('display_normalization_qqplot <- ')
    print(display_normalization_qqplot)
    cat('display_normalization_MAplot <- ')
    print(display_normalization_MAplot)
    
    cat( 
paste(c('layout(matrix(c(',paste(c(1:ncol(matrixdata)),collapse=','),'),1,',ncol(matrixdata),',byrow=TRUE))'), collapse=''),
paste(c('matrixdata <- as.matrix(read.table("',tempoutput,'", stringsAsFactors=FALSE))'), collapse=''),
'display_normalization_violin(matrixdata)',
'temp.qq <- apply(matrixdata, 2, FUN=display_normalization_qqplot)',
'nb.columns.per.exp <- ncol(matrixdata)/2',
'data.to.test1 <- apply(matrixdata[,1:nb.columns.per.exp], 1, mean)',
'data.to.test2 <- apply(matrixdata[,(nb.columns.per.exp+1):(2*nb.columns.per.exp)], 1, mean)',
'display_normalization_MAplot(data.to.test2,data.to.test1)',
sep = "\n"
    )

    cat('
```
')
    
}

####################################################################
#' On matrix-like data

display_normalization_boxplot <- function(matrixdata, title) {
    boxplot(matrixdata, main=title)
}

# display_normalization_densities <- function(matrixdata) {
#     
#     for(name in names(matrixdata)) {
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

display_normalization_violin <- function(matrixdata) {
    library("ggplot2")
    
    # Stack data by pasting column after column
    data.stacked <- c(matrixdata)
    
    row.names <- row.names(matrixdata)
    col.names <- colnames(matrixdata)
    numprots <- length(row.names)
    num.cond <- 2
    num.replicates <- ncol(matrixdata)/num.cond
    numallexp <- num.cond * num.replicates
    
    intensitiesdf <- data.frame(
        Intensities = data.stacked,
        Protein.groups = factor(rep(row.names, numallexp)),
        Reporters = factor( rep(col.names, rep(c(numprots), numallexp)) ),
        Conditions = factor( rep(c("Control","Exp."), c(numprots*num.replicates, numprots*num.replicates)) )
        )
    
    p <- ggplot(data=intensitiesdf, aes(x=Reporters, y=Intensities)) 
    p + geom_violin(scale = "count", mapping=aes(fill = Conditions), trim=TRUE, alpha=.50) +
        stat_summary(aes(group=Reporters), fun.y=mean, fun.ymin=min, fun.ymax=max, fill="red", shape=21, size=1) +
        ylab('glog2( Intensities )') +
        xlab('') +
        theme(legend.position="bottom")
}

####################################################################
#' Only on 1D vector

#' Both vectors should be log'ed before!
display_normalization_MAplot <- function(red1d, green1d) {
    
    
    a.vals <- (red1d + green1d)/2
    m.vals <- red1d - green1d
    dat <- data.frame(a.vals,m.vals)
    
    sd_foldchange <- sd(m.vals)
    
#     plot(a.vals, m.vals, pch=16)
    ggplot(dat, aes(x=a.vals, y=m.vals)) + geom_point(shape=16) + geom_abline(slope = 0, intercept = 0)
    
}


display_normalization_histogram <- function(vector1d, x_label = "Dataset") {

    #Display histogram of values
    histogram <- hist(vector1d, breaks=100, col="red", xlab=x_label, main="Histogram with Normal Curve")
    
    # Compute Normal Curve (Thanks to Peter Dalgaard)
    xfit <- seq( min(vector1d, na.rm = TRUE), max(vector1d, na.rm = TRUE), length=3000)
    yfit <- dnorm( xfit, mean=mean(vector1d, na.rm = TRUE), sd=sd(vector1d, na.rm = TRUE) )
    yfit <- yfit * diff( histogram$mids[1:2]) * length(vector1d)
    
    # Print Normal Curve
    lines(xfit, yfit, type="l", col="blue", lwd=2)
    
    # Print mean and standard deviation
    meanv <- mean(vector1d)
    sdv <- sd(vector1d)
    abline(v=c(meanv-sdv,meanv,meanv+sdv), col="green")

    
}


display_normalization_qqplot <- function(vector1d) {
    qqnorm(vector1d, main="", ylab='')
    qqline(vector1d)
}


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
test_normalization_shapiro <- function(vector1d) {
    stest <- shapiro.test(vector1d)
}

#' Kolmogorov-Smirnov Test
#'Calculates expected normal distribution and compares
#'it with the observed distribution.
#'Uses cumulative distribution functions.
#'Based on the max difference between two distributions.
#'Poor for discrimination below N = 30.
#'Power to detect differences is low.
#'Historically popular.
test_normality_kolmogorov_smirnov <- function(vector1d) {
    #ks.test(F500, "pnorm", mean=a[1], sd=a[2])
}