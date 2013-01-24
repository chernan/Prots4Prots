#' @title Report p-values of a statistical test
#' 
#' @description Report p-values of a statistical test (knitr)
#' 
#' @details 
#' Write a summary report of a significance test step. 
#' Display a volcano plot enhancing significant p-values before multiple 
#' testing correction, as well as the number of significant proteins.
#' Data at this step will be saved in a temporary file.
#' 
#' @param matrixData A data.frame containing the p-values (2 columns necessary: 
#'  "p.values" and "fold.change")
#' @param thresholdPVal Threshold to be used for significant p-values
#' @param outputFile Report output file (R markdown format)
#' @param outFolder Where the temporary file will be saved
#' 
allTestsSignificanceRmd <- function(matrixData, thresholdPVal, 
                                    outputFile, outFolder) {
    
    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc(runif(1) * 10000)), 
        collapse='')
    
    tempOutput <- paste(
        c(outFolder, '/all_tests_significance_Rmd_data_', execLabel, '.txt'), 
        collapse='')
    write.table(matrixData, tempOutput, sep="\t") 
    
    cat('',
        paste(
            c('```{r allTestsSignificanceRmd', 
              execLabel, ', echo=FALSE, fig.width=14, fig.height=10}'),
            collapse=''),
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    cat('displaySignificanceVolcanoPlot <- ', file=outputFile, append=TRUE)
#     print(displaySignificanceVolcanoPlot)
    cat(paste(deparse(displaySignificanceVolcanoPlot), collapse="\n"),
        file=outputFile, append=TRUE)
    
    cat('',
        paste(
            c('matrixData <- as.matrix(read.table("',
              tempOutput,
              '", stringsAsFactors=FALSE))'), collapse=''),
        paste(
            c('displaySignificanceVolcanoPlot(matrixData, thresholdPVal=',
              thresholdPVal,
              ')'), collapse=''),
        '```',
        '',
        '',
        paste(
            c('Without multiple testing correction, ',
              '<b>`r length(which(as.numeric(matrixData[,\"p.values\"]) < ',
              thresholdPVal,
              '))`</b> significant protein groups were found ',
              'with a p-value inferior to ',
              thresholdPVal,
              ' on a total of `r length(as.numeric(matrixData[,\"p.values\"]))`',
              '.'), 
            collapse=''),
        sep="\n", file=outputFile, append=TRUE
    )
    
}

## -----------------------------------------------------------------------------
## Tests / plots used in the report

#' @title Volcano plot
#' 
#' @description Volcano plot for non-corrected p-values
#' 
#' @details 
#' Display a volcano plot enhancing significant p-values before multiple 
#' testing correction.
#' 
#' @param matrixData A data.frame containing the p-values (2 columns necessary: 
#'  "p.values" and "fold.change")
#' @param thresholdPVal Threshold to be used for significant p-values
#' @param title (optional) A title for the plot
displaySignificanceVolcanoPlot <- function(matrixData, thresholdPVal, 
                                           title="Volcano plot") {
    
    pVals <- matrixData[, "p.values"]
    foldChange <- matrixData[, "fold.change"]
    
    plot(foldChange, -log2(pVals),
         main=title, sub=paste("(threshold p-value = ", thresholdPVal, ")"),
         xlab="log2( Fold change )", 
         ylab="-log2( p-value )",
         col="gray", 
         pch=16, cex.lab=1, cex.axis=1, cex.main=1)
    grid(col="lightgray")
    abline(h=-log2(thresholdPVal), col="orange", lty=2)
    abline(v=c(c(-2,-1,0,1,2)*sd(foldChange)), col="gray", lty=2)
    
    is_test_sig <- as.numeric(pVals) < thresholdPVal
    points(foldChange[is_test_sig], -log2(pVals[is_test_sig]), 
           pch=16, col="darkorange")
    #text(-log2(pVals[is_test_sig]), foldChange[is_test_sig], 
    #     rownames(matrixData)[is_test_sig], cex=0.7, pos=4, col="darkorange")
    
    legend("bottomright", 
           legend = c("Non significant", "Significant", 
                      paste("p-value threshold (", thresholdPVal, ")")), 
           col = c("gray", "darkorange", "orange"),
           fill = c("gray", "darkorange", 0),
           border = c("gray", "gray", 0),
           lty = c(0, 0, 2),
           merge = TRUE,
           cex = 0.8)
}


