
#' @title Report p-values corrected for multiple-testing
#' 
#' @description Report p-values corrected for multiple-testing (knitr)
#' 
#' @details 
#' Write a summary report of a multiple testing correction step. 
#' Display a volcano plot enhancing p-values still significant after multiple 
#' testing correction, as well as the number of significant proteins before 
#' and after correction.
#' Data at this step will be saved in a temporary file.
#' 
#' @param matrixData A data.frame containing the p-values (3 columns necessary: 
#'  "p.values", "p.values.corrected" and "fold.change")
#' @param thresholdPVal Threshold to be used for significant p-values
#' @param outFolder Where the temporary file will be saved
allTestsCorrectedRmd <- function(matrixData, thresholdPVal, 
                                 outputFile, outFolder) {
    
    ## Write data in a file usable later when the report will be transformed 
    ## into HTML
    tempOutput <- paste(
        c(outFolder, '/all_tests_BH_Rmd_', 
          format(Sys.time(), "%Y%m%d%H%M%S"), 
          trunc(runif(1)*10000),'.txt'), 
        collapse='')
    write.table(matrixData, tempOutput, sep="\t") 
    
    
    ## Write report text
    cat('',
        '```{r, echo=FALSE, fig.width=14, fig.height=10}',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    cat('displayCorrectedVolcanoPlot <- ', file=outputFile, append=TRUE)
#     print(displayCorrectedVolcanoPlot)
    cat(paste(deparse(displayCorrectedVolcanoPlot), collapse="\n"),
        file=outputFile, append=TRUE)
    
    cat('',
        paste(c('matrixData <- read.table("', 
                tempOutput, 
                '", stringsAsFactors=FALSE)'), 
              collapse=''),
        paste(c('displayCorrectedVolcanoPlot(matrixData, thresholdPVal=', 
                thresholdPVal, 
                ')'), 
              collapse=''),
        '```',
        '',
        paste(c("Before correction of the p-values for multiple-testing, ", 
                "<b>`r length(which(as.numeric(matrixData[,\"p.values\"]) < ", 
                thresholdPVal, "))`</b> ", 
                "protein groups had a p-value inferior to ", 
                thresholdPVal, "."), 
              collapse=''),
        '',
        paste(c("After correction of the p-values for multiple-testing, ",
                "<b>`r length(which(as.numeric(matrixData[,\"p.values.corrected\"]) < ", 
                thresholdPVal, "))`</b> ",
                "protein groups were found significant, ",
                "on a total of `r length(as.numeric(matrixData[,\"p.values\"]))`, ",
                "with an FDR of ", 
                thresholdPVal, "."), 
              collapse=''),
        sep="\n", file=outputFile, append=TRUE
    )
    
}

## -----------------------------------------------------------------------------
## Tests / plots used in the report

#' @title Volcano plot (corrected)
#' 
#' @description Volcano plot for corrected p-values
#' 
#' @details 
#' Display a volcano plot enhancing p-values still significant after multiple 
#' testing correction
#' 
#' @param matrixData A data.frame containing the p-values (3 columns necessary: 
#'  "p.values", "p.values.corrected" and "fold.change")
#' @param thresholdPVal Threshold to be used for significant p-values
#' @param title (optional) A title for the plot
displayCorrectedVolcanoPlot <- function(matrixData, thresholdPVal, 
                                        title="Volcano plot") {
    
    pVals <- matrixData[, "p.values"]
    pValsCorrected <- as.numeric(matrixData[, "p.values.corrected"])
    foldChange <- matrixData[, "fold.change"]
    
    plot(foldChange, -log2(pVals),
         main=title, 
         sub=paste("(FDR = ", thresholdPVal, ")"),
         xlab="log2( Fold change )", 
         ylab="-log2( p-value )",
         col="gray", 
         pch=16, cex.lab=1, cex.axis=1, cex.main=1)
    grid(col="lightgray")
    sdFoldChange <- sd(foldChange)
    abline(v=c(c(-2,-1,0,1,2)*sdFoldChange), col="gray", lty=2)
    
    isTestSig <- as.numeric(pVals) < thresholdPVal
    abline(h=-log2(thresholdPVal), col="orange", lty=2)
    points(foldChange[isTestSig], -log2(pVals[isTestSig]), 
           pch=16, col="darkorange")
    # text(foldChange[isTestSig], -log2(pVals[isTestSig]), 
    #      rownames(matrixData)[isTestSig], cex=0.7, pos=4, col="darkorange")
    
    isCorrectedSig <- pValsCorrected < thresholdPVal
    if (length(which(isCorrectedSig == TRUE)) > 0) {
        
        orderPVals <- order(pVals, decreasing=FALSE)
        reorderedFoldChange <- foldChange[orderPVals]
        reorderedPVals <- pVals[orderPVals]
        reorderedPValsCorrected <- pValsCorrected[orderPVals]
                                                  
        isCorrectedSig <- reorderedPValsCorrected < thresholdPVal
        maxPValsCorrected <- max(reorderedPValsCorrected[isCorrectedSig])
        indexMaxPVal <- which(
            reorderedPValsCorrected[isCorrectedSig] == maxPValsCorrected, 
            arr.ind=TRUE
        )
        abline( h=-log2(max(reorderedPVals[indexMaxPVal])), 
                col="red", lty=2)
        points( reorderedFoldChange[isCorrectedSig], 
                -log2(reorderedPVals[isCorrectedSig]), 
                pch=16, col="darkred")
    }
    
    legend("bottomright", 
           legend = c("Non significant", 
                      "Non corrected", 
                      "Significant after correction", 
                      "Non-corrected threshold", 
                      paste("Corrected threshold (FDR=", thresholdPVal, ")")
                      ), 
           col = c("gray", "darkorange", "darkred", "orange", "red"),
           fill = c("gray", "darkorange", "darkred", 0, 0),
           border = c("gray", "gray", "gray", 0, 0),
           lty = c(0, 0, 0, 2, 2),
           merge = TRUE,
           cex = 0.8)
}