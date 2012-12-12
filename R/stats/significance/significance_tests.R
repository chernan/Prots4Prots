
all_tests_significance_Rmd <- function(matrixdata, threshold_pval, outfolder) {
    tempoutput <- paste(c(outfolder, '/all_tests_significance_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(matrixdata, tempoutput, sep="\t") 
    
    cat('
```{r, echo=FALSE, fig.width=14, fig.height=10}
')
    
    cat('display_significance_volcanoplot <- ')
    print(display_significance_volcanoplot)
    
    cat( 
        paste(c('matrixdata <- as.matrix(read.table("',tempoutput,'", stringsAsFactors=FALSE))'), collapse=''),
        paste(c('display_significance_volcanoplot(matrixdata, threshold_pval=',threshold_pval,')'), collapse=''),
        sep = "\n"
    )
    
    cat('
```

')

    cat(
        paste(c("Without multiple testing correction, <b>`r length(which(as.numeric(matrixdata[,\"p.values\"]) < ",threshold_pval,"))`</b> significant protein groups were found, with a p-value inferior to ",threshold_pval,"."))
    )
    
}

display_significance_volcanoplot <- function(matrixdata, threshold_pval, title="Volcano plot") {
    
    pvals <- matrixdata[,"p.values"]
    foldchange <- matrixdata[,"fold.change"]
    
    plot( foldchange, -log2(pvals),
          main=title, sub=paste("(threshold p-value = ",threshold_pval,")"),
          xlab="log2( Fold change )", 
          ylab="-log2( p-value )",
          col="gray", 
          pch=16, cex.lab = 1, cex.axis = 1, cex.main = 1)
    grid(col="lightgray")
    abline(h=-log2(threshold_pval), col="orange", lty=2)
    sd_foldchange <- sd(foldchange)
    abline(v=c(c(-2,-1,0,1,2)*sd_foldchange), col="gray", lty=2)
    
    is_test_sig <- as.numeric(pvals) < threshold_pval
    points( foldchange[is_test_sig], -log2(pvals[is_test_sig]), pch=16, col="darkorange")
    #text( -log2(pvals[is_test_sig]), foldchange[is_test_sig], rownames(matrixdata)[is_test_sig], cex=0.7, pos=4, col="darkorange")
    
    legend("bottomright", 
           legend=c("Non significant", "Significant", paste("p-value threshold (",threshold_pval,")")), 
           col=c("gray", "darkorange", "orange"),
           fill=c("gray", "darkorange", 0),
           border=c("gray", "gray", 0),
           lty = c(0, 0, 2),
           merge = TRUE,
           cex = 0.8)
}


###########################################################################################################################
# For corrected p-values


all_tests_corrected_Rmd <- function(matrixdata, threshold_pval, outfolder) {
    tempoutput <- paste(c(outfolder, '/all_tests_BH_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(matrixdata, tempoutput, sep="\t") 
    
    cat('
```{r, echo=FALSE, fig.width=14, fig.height=10}
')
    
    cat('display_corrected_volcanoplot <- ')
    print(display_corrected_volcanoplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('display_corrected_volcanoplot(matrixdata, threshold_pval=',threshold_pval,')'), collapse=''),
        sep = "\n"
    )
    
    cat('
```

')

    cat(
        paste(c("Before multiple testing correction, <b>`r length(which(as.numeric(matrixdata[,\"p.values\"]) < ",threshold_pval,"))`</b> protein groups had a p-value inferior to ",threshold_pval,"."), collapse=''),
        paste(c("After Benjamini-Hochberg correction of the p-values for multiple-testing, <b>`r length(which(as.numeric(matrixdata[,\"p.values.corrected\"]) < ",threshold_pval,"))`</b> protein groups were found significant, on a total of `r length(as.numeric(matrixdata[,\"p.values\"]))`, with an FDR of ",threshold_pval,"."), collapse=''),
        sep="\n"
    )
    
}

display_corrected_volcanoplot <- function(matrixdata, threshold_pval, title="Volcano plot") {
    
    pvals <- matrixdata[,"p.values"]
    pvals.correct <- as.numeric(matrixdata[,"p.values.corrected"])
    foldchange <- matrixdata[,"fold.change"]
    
    plot( foldchange, -log2(pvals),
          main="Volcano plot", sub=paste("(threshold p-value = ",threshold_pval,")"),
          xlab="log2( Fold change )", 
          ylab="-log2( p-value )",
          col="gray", 
          pch=16, cex.lab = 1, cex.axis = 1, cex.main = 1)
    grid(col="lightgray")
    sd_foldchange <- sd(foldchange)
    abline(v=c(c(-2,-1,0,1,2)*sd_foldchange), col="gray", lty=2)
    
    is_test_sig <- as.numeric(pvals) < threshold_pval
    abline( h=-log2(threshold_pval), col="orange", lty=2)
    points( foldchange[is_test_sig], -log2(pvals[is_test_sig]), pch=16, col="darkorange")
    #text( foldchange[is_test_sig], -log2(pvals[is_test_sig]), rownames(matrixdata)[is_test_sig], cex=0.7, pos=4, col="darkorange")
    
    is_corrected_sig <- pvals.correct < threshold_pval
    if(length( which(is_corrected_sig == TRUE) ) >0) {
        order.pvals <- order(pvals, decreasing=FALSE)
        reordered.foldchange <- foldchange[order.pvals]
        reordered.pvals <- pvals[order.pvals]
        reordered.pvals.correct <- pvals.correct[order.pvals]
        is_corrected_sig <- reordered.pvals.correct < threshold_pval
        
        max.pvalcorrect <- max(reordered.pvals.correct[is_corrected_sig])
        ind.maxpval <- which(reordered.pvals.correct[is_corrected_sig]==max.pvalcorrect,arr.ind=TRUE)
        abline( h=-log2( max(reordered.pvals[ ind.maxpval ]) ), col="red", lty=2)
        points( reordered.foldchange[is_corrected_sig], -log2(reordered.pvals[is_corrected_sig]), pch=16, col="darkred")
    }
    
    legend("bottomright", 
           legend=c("Non significant", "Non corrected", "Significant after correction", paste("Non-corrected threshold"), paste("Corrected threshold (FDR=",threshold_pval,")")), 
           col=c("gray", "darkorange", "darkred", "orange", "red"),
           fill=c("gray", "darkorange", "darkred", 0, 0),
           border=c("gray", "gray", "gray", 0, 0),
           lty = c(0, 0, 0, 2, 2),
           merge = TRUE,
           cex = 0.8)
}