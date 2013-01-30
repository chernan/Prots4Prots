#' @name SAM
#' @title Significance Analysis of Microarray data
#' @description
#' Significance Analysis of Microarray data (SAM)
#' 
#' @keywords Significance
#' 
#' @section Introduction:
#' p-values are computed by bootstraping.
#' 
#' @section References:
#' Tusher, V., Tibshirani, R., & Chu, G. (2001). 
#' Significance analysis of microarrays applied to the ionizing radiation response. 
#' Proceedings of the National Academy of Sciences of the United States of America, 98(9), 5116-21. 
#' doi:10.1073/pnas.091062498
#' 
#' @section Installation:
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("impute")
#' install.packages("samr")
#' 
#' @import samr
# library("samr")
NULL

#' @title Apply SAM method
#' 
#' @description Apply SAM method on non-paired data.
#' 
#' @details 
#' Apply SAM method on two-class non-paired data, supposing that data has been 
#' log2-transformed before. 
#'  
#' @param experiment Dataset for experiment (condition 1).
#' @param control Dataset for control (condition 2).
#' @param thresholdPVal Desired FDR.
#' @return A list containing p-values, and the output object of SAM().
#' @export
applySam <- function(experiment, control, thresholdPVal) {
    
    data <- data.frame(control, experiment)
    rowNames <- rownames(experiment)
    row.names(data) <- rowNames
    
    design <- c(rep(1, ncol(control)), rep(2, ncol(experiment)))
    
    capture.output(
        samFit <- SAM(x=data, y=design, 
                      logged2=TRUE, 
                      resp.type="Two class unpaired", fdr.output=thresholdPVal,
                      genenames=rowNames, geneid=rowNames)
    )
    capture.output(
        pValues <- samr.pvalues.from.perms(samFit$samr.obj$tt, 
                                           samFit$samr.obj$ttstar)
    )
    returnValues <- list(sam.fit=samFit,
                         p.values=pValues
    )
    
    return(returnValues)
}

#' @title Apply SAM method on paired data
#' 
#' @description Apply SAM method on paired data.
#' 
#' @details 
#' Apply SAM method, supposing that data has been log2-transformed before. 
#' Data should be paired, meaning that first column of experiment is paired 
#' with the first column of control, etc.
#'  
#' @param experiment Dataset for experiment (condition 1).
#' @param control Dataset for control (condition 2).
#' @param thresholdPVal Desired FDR.
#' @return A list containing p-values, and the output object of SAM().
#' @export
applyPairedSam <- function(experiment, control, thresholdPVal) {
    
    data <- as.matrix(data.frame(control, experiment))
    rowNames <- rownames(experiment)
    row.names(data) <- rowNames
    
    design <- rep(c(1:ncol(control)), times=2)
    design[1:ncol(control)] <- design[1:ncol(control)] * -1
    
    capture.output(
        samFit <- SAM(x=data, y=design, 
                      logged2=TRUE, 
                      resp.type="Two class paired", fdr.output=thresholdPVal,
                      genenames=rowNames, geneid=rowNames,
        )
    )
    capture.output(
        pValues <- samr.pvalues.from.perms(samFit$samr.obj$tt, 
                                           samFit$samr.obj$ttstar)
    )
    returnValues <- list(sam.fit = samFit,
                         p.values = pValues
    )
    
    return(returnValues)
}

#' @title Apply SAM and report result
#' 
#' @description Apply SAM and report result
#' 
#' @details 
#' Apply SAM between two datasets. 
#' Note that SAM dosen't give any way to know which genes/proteins are 
#' considered up/down regulated, excepted by reading which ones are listed 
#' in siggenes.table. There is no way to know which limit was applied to 
#' achieve significance.
#' TODO Use ExpressionSet object as input
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
#' @return A data frame with three columns "p.values" "fold.change" 
#'  "significant", in the same order as the input data
applyAndReportSam <- function(experiment, control, 
                              outputFolderTemp, outputFile, 
                              thresholdPVal, isPaired=FALSE) {
    
    ## Run corresponding function if data are paired or not
    if(!isPaired) {
        samrVals <- applySam(experiment, control, thresholdPVal)
    } else {
        samrVals <- applyPairedSam(experiment, control, thresholdPVal)
    }
    
    
    ## Re-shape data
    ## Note that SAM dosen't give any way to know which genes/proteins are 
    ## considered up/down regulated, excepted by reading which ones are listed 
    ## in siggenes.table. There is no way to know which limit was applied to 
    ## achieve significance.
    samFit <- samrVals$sam.fit
    pValues <- samrVals$p.values
    significant <- matrix(rep('STABLE', nrow(experiment)))
    names(significant) <- as.character(row.names(experiment))
    
    ## Header:
    ## Gene ID  Gene Name  Score(d)  Numerator(r)  Denominator(s+s0)  Fold Change  q-value(%)
    ## Get 2nd column, but no colnames?? but a print shows the col.names!
    
    ## 2* Special cases for '1' necessary, otherwise get an error:
    ##  Erreur dans sam.fit$siggenes.table$genes.lo[1, 2] : 
    ##      nombre de dimensions incorrect
    
    if(samFit$siggenes.table$ngenes.lo == 1) {
        significant[samFit$siggenes.table$genes.lo[2]] <- 'LO'
    }
    if(samFit$siggenes.table$ngenes.lo > 1) {
        pos <- samFit$siggenes.table$ngenes.lo
        significant[samFit$siggenes.table$genes.lo[1:pos,2]] <- 'LO'
    }
    
    if(samFit$siggenes.table$ngenes.up == 1 ) {
        significant[samFit$siggenes.table$genes.up[2]] <- 'UP'
    }
    if(samFit$siggenes.table$ngenes.up > 1 ) {
        pos <- samFit$siggenes.table$ngenes.up
        significant[samFit$siggenes.table$genes.up[1:pos,2]] <- 'UP'
    }
    
    
    ## Final returned value
    matrixData <- data.frame(
        p.values = pValues,
        fold.change = log2(samFit$samr.obj$foldchange),
        significant = significant
    )
    
    ## Report results
    reportSam(matrixData, thresholdPVal, outputFile, outputFolderTemp)
    
    return(matrixData)
}

#' @title Report result of Limma
#' 
#' @description Report result of Limma
#' 
#' @details 
#' Report Limma on two datasets. 
#' 
#' @param matrixData A data frame with two columns "p.values" "fold.change"
#' @param thresholdPVal Maximum threshold for the p-values (e.g 0.05)
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
reportSam <- function(matrixData, thresholdPVal, 
                      outputFile, outputFolderTemp) {
    
    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc( 
            runif(1) * 10000)), 
        collapse='')

    cat('',
        'Significantly different proteins',
        '---------------------------------------------------------------------',
        '',
        '',
        'Data analysis was performed using Significance Analysis of Microarray data (R package version `r packageDescription("samr")$Version`).',
        '',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    tempOutput <- paste(c(outputFolderTemp, '/samr_tests_significance_Rmd_data_', 
                          execLabel, '.txt'), 
                        collapse='')
    write.table(matrixData, tempOutput, sep="\t") 
    
    cat('',
        paste(
            c('```{r reportSam2', execLabel, 
              ', echo=FALSE, fig.width=14, fig.height=10}'),
            collapse=''),
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    cat('displaySamSignificanceVolcanoPlot <- ', file=outputFile, append=TRUE)
    #     print(displaySamSignificanceVolcanoPlot)
    cat(paste(deparse(displaySamSignificanceVolcanoPlot), collapse="\n"),
        file=outputFile, append=TRUE)
    
    
    cat('',
        paste(c('matrixData <- read.table("',
                tempOutput,
                '", stringsAsFactors=FALSE)'), 
              collapse=''),
        paste(c(
            'displaySamSignificanceVolcanoPlot(matrixData, thresholdPVal=',
            thresholdPVal,
            ')'),
              collapse=''),
        '',
        '```',
        '', 
        paste(
            c("SAM found <b>`r length(which(matrixData[[\"significant\"]]=='UP'))`</b> significant up-regulated and <b>`r length(which(matrixData[[\"significant\"]]=='LO'))`</b> down-regulated protein groups, with an estimated FDR of ",
              thresholdPVal,
              "."),
            collapse=''), 
        ' ',       
        '> ',
        '> Please cite these articles whenever using results from this software in a publication :',
        '> ',
        '> Method article :',
        '> Significance analysis of microarrays applied to the ionizing radiation response. ',
        '> Tusher, V., Tibshirani, R., & Chu, G. (2001). ',
        '> Proceedings of the National Academy of Sciences of the United States of America, 98(9), 5116-21. ',
        '> doi:10.1073/pnas.091062498',
        '> ',
        '> Software article (R package) :',
        '> `r citation("samr")$textVersion`',
        '> ',
        '',
        '---------------------------------------------------------------------',
        '', 
        sep="\n", file=outputFile, append=TRUE)        
    
}

displaySamSignificanceVolcanoPlot <- function(matrixData, thresholdPVal, 
                                              title="Volcano plot") {
    
    pVals <- matrixData[, "p.values"]
    foldChange <- matrixData[, "fold.change"]
    samSignificant <- matrixData[, "significant"]
    
    plot( foldChange, -log2(pVals),
          main=title, sub=paste("(FDR = ", thresholdPVal, ")"),
          xlab="log2( Fold change )", 
          ylab="-log2( p-value )",
          col="gray", 
          pch=16, cex.lab = 1, cex.axis = 1, cex.main = 1)
    grid(col="lightgray")
    sdFoldChange <- sd(foldChange)
    abline(v=c(c(-2,-1,0,1,2) * sdFoldChange), col="gray", lty=2)
    
    if(any(samSignificant == 'UP') | any(samSignificant == 'LO')) {
        points( foldChange[samSignificant == 'UP'], 
                -log2(pVals[samSignificant == 'UP']), 
                pch=16, col="darkred")
        points( foldChange[samSignificant == 'LO'], 
                -log2(pVals[samSignificant == 'LO']), 
                pch=16, col="darkgreen")
    }
    
    legend("bottomright", 
           legend=c("Non significant", "Significant up-regulated", 
                    "Significant down-regulated", 
                    paste("Significance threshold (FDR=", thresholdPVal, ")")), 
           col=c("gray", "darkred", "darkgreen", "orange"),
           fill=c("gray", "darkred", "darkgreen", 0),
           border=c("gray", "gray", "gray", 0),
           lty=c(0, 0, 0, 2),
           merge=TRUE)
}

