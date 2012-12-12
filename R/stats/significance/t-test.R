#! Student's t-test
#! 
#! R package
#! Article

## Suppose sample gaussian and with same variance!
t_test_pval <- function(xexp, xcontrol, confidence=0.95, is.homoscedastic, is.paired) {
    test_res <- t.test(x=xexp, y=xcontrol, na.action=na.omit, var.equal=is.homoscedastic, conf.level=confidence, paired=is.paired)
    return (test_res)
}


# t-test
apply_t_test <- function(experiment, control, rownames, threshold, is.paired) {
    nbrow <- nrow(experiment)
    confidence <- (1 - threshold)
    
    stacked.data <- stack(data.frame(cbind(experiment,control), row.names=rownames))
    is.homoscedastic <- (bartlett.test(stacked.data$values, g=stacked.data$ind)$p.value < 0.05)
    
    resultt <- t(
        apply(matrix(c(1:nbrow),ncol=1), 1, 
                   FUN = function(x) { 
                          ttest <- t_test_pval(as.matrix(experiment[x,]), as.matrix(control[x,]), confidence, is.homoscedastic, is.paired)
                          if(!is.paired) {
                              return(c(
                                  pval=ttest$p.value, 
                                  foldchange=ttest$estimate[["mean of x"]]-ttest$estimate[["mean of y"]]
                              ))
                          } else {
                              return(c(
                                  pval=ttest$p.value, 
                                  foldchange=ttest$estimate[["mean of the differences"]]
                              ))
                          }
                      }
              )
    )
    resultt <- data.frame(resultt, row.names=rownames)
    return(resultt)
}

apply_and_report_t_test <- function(experiment, control, output_folder_temp, output_file, threshold_pval, is.paired=FALSE) {
    dataset.ttest <- apply_t_test(experiment, control, rownames(experiment), threshold_pval, is.paired)
    
    matrixdata <- data.frame(
        p.values = dataset.ttest$pval,
        fold.change = dataset.ttest$foldchange,
        row.names = rownames(experiment)
    )
    
    stacked.data <- stack(data.frame(cbind(experiment,control), row.names=rownames(experiment)))
    is.homoscedastic <- (bartlett.test(stacked.data$values, g=stacked.data$ind)$p.value < 0.05)

    #Erreur dans shapiro.test(vector1d) : sample size must be between 3 and 5000
    if(nrow(experiment) > 3 & nrow(experiment) < 5000) {
        is.norm.by.sample <- apply(cbind(experiment,control), 2, FUN=function(vector1d){return(shapiro.test(vector1d)$p.value)} )
    }else {
        is.norm.by.sample <- NULL
    }

    report_t_test(output_file, matrixdata, threshold_pval, output_folder_temp, is.homoscedastic, is.norm.by.sample, is.paired)
        
    return(matrixdata)
}

report_t_test <- function(output_file, matrixdata, threshold_pval, output_folder_temp, is.homoscedastic, is.norm.by.sample, is.paired) {
    sink(output_file, append=TRUE)
    cat('
Significantly different proteins
--------------------------------------------------------------------------------------------------------

')
    cat('Data analysis was performed using Student\'s t-test ',
        ifelse(is.paired, 'for paired samples ', 'for non-paired samples '),
        ifelse(is.homoscedastic==TRUE,'(equal variance)', '(unequal variance)'),
        '.',
        "\n",
        "\n",
        sep='')
        
    if(!all(sapply(is.norm.by.sample, is.null))) {
        cat('Check normality of input data using Shapiro\'s test, for each sample :',
            '',
            'Experiment | Control',
            '------------- | -------------',
            apply( matrix(format.pval(is.norm.by.sample, digits=4), ncol=2, byrow=FALSE) , 1, FUN=function(x) {paste(x,collapse=' | ',sep="\n")}),
            '',
            sep="\n")
    }
    
    all_tests_significance_Rmd(matrixdata, threshold_pval, output_folder_temp)
    
    cat('        
> 
> Student\'s t-Test : require no citation.
>

---------------------------------------------------------------------------
')        
    sink(file=NULL)
    
}



#Bartlett's or Levene's test 
# Generally speaking, the Levene test works well in the ANOVA framework, providing there are small to moderate deviations from the normality. 
# In this case, it outperfoms the Bartlett test. If the distribution are nearly normal, however, the Bartlett test is better. 
# I've also heard of the Brown–Forsythe test as a non-parametric alternative to the Levene test. Basically, it relies on either the median or 
# the trimmed mean (as compared to the mean in the Levene test). According to Brown and Forsythe (1974), a test based on the mean provided 
# the best power for symmetric distributions with moderate tails.
# In conclusion, I would say that if there is strong evidence of departure from the normality (as seen e.g., with the help of a Q-Q plot), 
# then use a non-parametric test (FK or BF test); otherwise, use Levene or Bartlett test.
# Use Bartlett’s test if your data follow a normal, bell-shaped distribution. 
# If your samples are small, or your data are not normal (or you don’t know whether they’re normal), use Levene’s test.
# 
#Test homoscedasticity (Bartlett\'s test)
#If the variances are different from each other (heteroscedasticity), the probability of obtaining a "significant" result may be greater than the desired alpha level (even if the null hypothesis is true).


