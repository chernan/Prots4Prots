#' Limma
#' 
#' Article

# #Print-tip loess normalization:
# MA <- normalizeWithinArrays(RG)
# #Estimate the fold changes and standard errors by fitting a linear model for each gene. The
# #design matrix indicates which arrays are dye-swaps.
# fit <- lmFit(MA, design=c(-1,1,-1,1))
# #Apply empirical Bayes smoothing to the standard errors.
# fit <- eBayes(fit)

# For reproducible results
set.seed(0)
# Limma library
library("limma")

# Values must be in log2
apply_limma <- function(experiment, control) {

    # Estimate the fold changes and standard errors by fitting a linear model for each gene. 
    # The design matrix indicates which columns are from a different condition.
    fit <- lmFit(
        cbind(control,experiment),
        method = "ls",
        design = cbind(Grp1=1,Grp2vs1=c(rep(0,ncol(control)),rep(1,ncol(experiment))))
        )
    
    # Apply empirical Bayes smoothing to the standard errors.
    #Need a proportion of supposedly differentially expressed genes?
    fit <- eBayes(fit) # , proportion=0.01
    
    # Compute fold change and p-values
    limma.vals <- topTable(fit, number=nrow(experiment), coef=2)

    # Order by id
    limma.vals <- limma.vals[order(as.numeric(row.names(limma.vals))),]
    
    return(limma.vals)
}


apply_paired_limma <- function(experiment, control) {
    
    Replicates <- factor(rep(c(1:ncol(control)), times=2))
    Treat <- factor(rep(c("C","T"), each=ncol(control)), levels=c("C","T"))

#     dupcor <- duplicateCorrelation(
#         cbind(control,experiment), 
#         block=c(1,2,3,1,2,3),
#         design=cbind(Grp1=1,Grp2vs1=c(rep(0,ncol(control)),rep(1,ncol(experiment)))))

    # Estimate the fold changes and standard errors by fitting a linear model for each gene. 
    fit <- lmFit(
        cbind(control,experiment), 
        method = "ls",
        design = model.matrix(~Replicates+Treat)
    )
    
    # Apply empirical Bayes smoothing to the standard errors.
    fit <- eBayes(fit) # , proportion=0.01
    # Compute fold change and p-values
    limma.vals <- topTable(fit, number=nrow(experiment), coef="TreatT")
    
    limma.vals <- limma.vals[order(as.numeric(row.names(limma.vals))),]
    
    return(limma.vals)
}


apply_and_report_limma <- function(experiment, control, output_folder_temp, output_file, threshold_pval, is.paired=FALSE) {
    if(!is.paired) {
        limma.vals <- apply_limma(experiment, control)
    } else {
        limma.vals <- apply_paired_limma(experiment, control)
    }
    
    matrixdata <- data.frame(
        p.values = limma.vals$P.Value,
        fold.change = limma.vals$logFC,
        row.names = row.names(experiment)
    )
    
    report_limma(output_file, matrixdata, threshold_pval, output_folder_temp)
        
    return(matrixdata)
}



report_limma <- function(output_file, matrixdata, threshold_pval, output_folder_temp) {
    sink(output_file, append=TRUE)
    cat('
Significantly different proteins
--------------------------------------------------------------------------------------------------------

Data analysis was performed using Limma modeling.

')
    
    all_tests_significance_Rmd(matrixdata, threshold_pval, output_folder_temp)
    
    cat('        
> 
> Limma citation.
>

---------------------------------------------------------------------------
')        
    sink(file=NULL)
    
}
