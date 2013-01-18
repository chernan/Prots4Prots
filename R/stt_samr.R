#' Significance Analysis of Microarray data
#' 
#' p-values are computed by bootstraping.
#' 
#' Tusher, V., Tibshirani, R., & Chu, G. (2001). 
#' Significance analysis of microarrays applied to the ionizing radiation response. 
#' Proceedings of the National Academy of Sciences of the United States of America, 98(9), 5116–21. 
#' doi:10.1073/pnas.091062498
#' 
#' How to install:
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("impute")
#' install.packages("samr")
#' 


# set.seed(0)
library("samr")



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
applySam <- function(experiment, control, thresholdPVal) {
    
    data <- data.frame(control, experiment)
    rowNames <- rownames(experiment)
    row.names(data) <- rowNames
    
    design <- c(rep(1, ncol(control)), rep(2, ncol(experiment)))
    
    samFit <- SAM(x=data, y=design, 
                  logged2=TRUE, 
                  resp.type="Two class unpaired", fdr.output=thresholdPVal,
                  genenames=rowNames, geneid=rowNames)
    
    pValues <- samr.pvalues.from.perms(samFit$samr.obj$tt, 
                                       samFit$samr.obj$ttstar)
    
    returnValues <- list(sam.fit=samFit,
                         p.values=pValues
    )
    
    #     significant <- matrix(rep('STABLE', nrow(experiment)))
    #     names(significant) <- rowNames
    # 
    #     #Get 2nd column, but no colnames?? but a print shows the col.names!
    #     #Header:
    #     #Gene ID  Gene Name  Score(d)  Numerator(r)  Denominator(s+s0)  Fold Change  q-value(%)
    #     
    #     # 2* Special cases for '1' necessary, otherwise get an error:
    #     #  Erreur dans sam.fit$siggenes.table$genes.lo[1, 2] : 
    #     #      nombre de dimensions incorrect
    #     
    #     if(sam.fit$siggenes.table$ngenes.lo == 1) {
    #         significant[sam.fit$siggenes.table$genes.lo[2]] <- 'LO'
    #     }
    #     if(sam.fit$siggenes.table$ngenes.lo > 1) {
    #         pos <- sam.fit$siggenes.table$ngenes.lo
    #         significant[sam.fit$siggenes.table$genes.lo[1:pos,2]] <- 'LO'
    #     }
    #     
    #     if(sam.fit$siggenes.table$ngenes.up == 1 ) {
    #         significant[sam.fit$siggenes.table$genes.up[2]] <- 'UP'
    #     }
    #     if(sam.fit$siggenes.table$ngenes.up > 1 ) {
    #         pos <- sam.fit$siggenes.table$ngenes.up
    #         significant[sam.fit$siggenes.table$genes.up[1:pos,2]] <- 'UP'
    #     }
    #     
    #     return.values <- data.frame(
    #         p.values = p.values,
    #         fold.change = log2(sam.fit$samr.obj$foldchange),
    #         significant = significant
    #         )
    
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
applyPairedSam <- function(experiment, control, thresholdPVal) {
    
    data <- as.matrix(data.frame(control, experiment))
    rowNames <- rownames(experiment)
    row.names(data) <- rowNames
    
    design <- rep(c(1:ncol(control)), times=2)
    design[1:ncol(control)] <- design[1:ncol(control)] * -1
    
    samFit <- SAM(x=data, y=design, 
                  logged2=TRUE, 
                  resp.type="Two class paired", fdr.output=thresholdPVal,
                  genenames=rowNames, geneid=rowNames,
    )
    
    pValues <- samr.pvalues.from.perms(samFit$samr.obj$tt, 
                                       samFit$samr.obj$ttstar)
    
    returnValues <- list(sam.fit = samFit,
                         p.values = pValues
    )
    
    #     significant <- matrix(rep('STABLE',nrow(experiment)))
    #     names(significant) <- rowNames
    #     
    #     #Get 2nd column, but no colnames?? but a print shows the col.names!
    #     #Header:
    #     #Gene ID  Gene Name  Score(d)  Numerator(r)  Denominator(s+s0)  Fold Change  q-value(%)
    #     
    #     # 2* Special cases for '1' necessary, otherwise get an error:
    #     #  Erreur dans sam.fit$siggenes.table$genes.lo[1, 2] : 
    #     #      nombre de dimensions incorrect
    #     
    #     if(sam.fit$siggenes.table$ngenes.lo == 1) {
    #         significant[sam.fit$siggenes.table$genes.lo[2]] <- 'LO'
    #     }
    #     if(sam.fit$siggenes.table$ngenes.lo > 1) {
    #         pos <- sam.fit$siggenes.table$ngenes.lo
    #         significant[sam.fit$siggenes.table$genes.lo[1:pos,2]] <- 'LO'
    #     }
    #     
    #     if(sam.fit$siggenes.table$ngenes.up == 1 ) {
    #         significant[sam.fit$siggenes.table$genes.up[2]] <- 'UP'
    #     }
    #     if(sam.fit$siggenes.table$ngenes.up > 1 ) {
    #         pos <- sam.fit$siggenes.table$ngenes.up
    #         significant[sam.fit$siggenes.table$genes.up[1:pos,2]] <- 'UP'
    #     }
    #     
    #     return.values <- data.frame(
    #         p.values = p.values,
    #         fold.change = log2(sam.fit$samr.obj$foldchange),
    #         significant = significant
    #     )
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
    
    cat('',
        'Significantly different proteins',
        '---------------------------------------------------------------------',
        '',
        '```{r echo=FALSE, warning=FALSE}',
        'samCitation <- citation("samr")',
        'samDescription <- packageDescription("samr")',
        '```',
        '',
        'Data analysis was performed using Significance Analysis of Microarray data (R package version `r samDescription$Version`).',
        '',
        '',
        sep="\n", file=outputFile, append=TRUE)
    
    tempOutput <- paste(c(outputFolderTemp, '/samr_tests_significance_Rmd_', 
                          format(Sys.time(), "%Y%m%d%H%M%S"), 
                          trunc(runif(1) * 10000), '.txt'), 
                        collapse='')
    write.table(matrixData, tempOutput, sep="\t") 
    
    cat('',
        '```{r, echo=FALSE, fig.width=14, fig.height=10}',
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
        '> Proceedings of the National Academy of Sciences of the United States of America, 98(9), 5116–21. ',
        '> doi:10.1073/pnas.091062498',
        '> ',
        '> Software article (R package) :',
        '> `r samCitation$textVersion`',
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

###########################################################################################################

# uppercent <- function(normal_data, percent, datafolder) {
#     copydata <- normal_data
#     selected_rows <- runif(nrow(normal_data)) < percent/100
#     delta <- 0.4*abs(rowMeans(normal_data[selected_rows,]))
#     copydata[selected_rows,1:2] <- normal_data[selected_rows,1:2] + delta
#     copydata[selected_rows,3:4] <- normal_data[selected_rows,3:4] - delta
#     
#     copydata <- 2^copydata
#     dataSet <- new("ExpressionSet", exprs = copydata)
#     
#     tempoutput <- paste(c(datafolder, '/normal8up',percent,'_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
#     write.table(cbind(copydata,selected_rows), tempoutput, sep="\t")
#     
#     return(list(dataset=dataSet,type="prot",file=tempoutput))
# }
# 
# 
# 
# 
# test_samr <- function() {
# 
#     source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/01_quality_control/heatmap.R")
#     source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/00_normalization/normalization_test.R")
#     source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/00_normalization/variance-stabilizing.R")
#     source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/00_normalization/samr.R")
#     
#     
#     #Set output directory as working directory. All output files/images/etc. will be created there
#     oldwd <- getwd()
#     basedir <- '/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/out_summary_4plex'
#     setwd(basedir)
#     
#     datafolder <- "/home/chernan/Workspace/GitHub/Prots4Prots/test_data/normal8ups_4plex"
#     dir.create(datafolder)
# 
#     normal_data <- replicate(4, rnorm(n=3000, mean=11.5, sd=2.2))
#     dataset.name <- "normal8up0"
#     dataset.obj <- uppercent(normal_data, 0, datafolder)
#     dataset <- dataset.obj[['dataset']]
#     all.methods.label <- paste(c('vsn05','samr'), collapse='_')
#     
#     ## Setup environment ####
#     
#     #Create working directory
#     working.dir <- paste(c(getwd(),'/', dataset.name), collapse='')
#     dir.create(working.dir)
#     output_folder <- paste(c(working.dir,'/', all.methods.label), collapse='')
#     dir.create(output_folder)
#     
#     #Create temp folder for temp files used in Rmd report
#     output_folder_temp <- paste(c(output_folder,'/temp'), collapse='')
#     dir.create(output_folder_temp)
#     #Create Rmd report file
#     output_file_name_Rmd <- paste(c(output_folder,'/output_',all.methods.label,'.Rmd'), collapse='')
#     write(c(''),output_file_name_Rmd)
#     
#     # Analysis
#     
#     #First quality checks before normalization
#     report_heatmap(exprs(dataset), output_folder_temp, output_file_name_Rmd, dist.method="manhattan")
#     
#     # Normalization
#     out.val <- apply_and_report_vsn(exprs(dataset), output_folder_temp, output_file_name_Rmd, use.robust.fit=TRUE)
#     
#     #Second quality check after normalization
#     report_heatmap(out.val, output_folder_temp, output_file_name_Rmd, dist.method="manhattan")
#     
#     #Get data to compare, depending on experimental design
#     #Temporary! >_<; 
#     experiment <- out.val[,1:2]
#     control <- out.val[,3:4]
#     threshold_pval <- 0.001
#     
#     #Test significance
# #     out.pval <- apply_and_report_samr(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval)
#     rownames<- rownames(experiment)
#     data <- data.frame(cbind(control, experiment))
#     row.names(data) <- rownames
#     
#     design <- c(rep(1,ncol(control)), rep(2,ncol(experiment)))
#     
#     sam.fit <- SAM(x=data, y=design, logged2=FALSE, 
#                    resp.type="Two class unpaired", fdr.output=threshold_pval)
#     
#     p.values <- samr.pvalues.from.perms(sam.fit$samr.obj$tt, sam.fit$samr.obj$ttstar)
#     
#     significant <- matrix(rep('STABLE',nrow(experiment)))
#     names(significant) <- rownames
# 
#     #Get 2nd column, but no colnames?? but a print shows the col.names!
#     #Header:
#     #Gene ID  Gene Name  Score(d)  Numerator(r)  Denominator(s+s0)  Fold Change  q-value(%)
#     
#     # 2* Special cases for '1' necessary, otherwise get an error:
#     #  Erreur dans sam.fit$siggenes.table$genes.lo[1, 2] : 
#     #      nombre de dimensions incorrect
# 
#     if(sam.fit$siggenes.table$ngenes.lo == 1) {
#         significant[sam.fit$siggenes.table$genes.lo[2]] <- 'LO'
#     }
#     if(sam.fit$siggenes.table$ngenes.lo > 1) {
#         pos <- sam.fit$siggenes.table$ngenes.lo
#         significant[sam.fit$siggenes.table$genes.lo[1:pos,2]] <- 'LO'
#     }
#     
#     if(sam.fit$siggenes.table$ngenes.up == 1 ) {
#         significant[sam.fit$siggenes.table$genes.up[2]] <- 'UP'
#     }
#     if(sam.fit$siggenes.table$ngenes.up > 1 ) {
#         pos <- sam.fit$siggenes.table$ngenes.up
#         significant[sam.fit$siggenes.table$genes.up[1:pos,2]] <- 'UP'
#     }
#     
#     return.values <- data.frame(
#         p.values = p.values,
#         fold.change = log2(sam.fit$samr.obj$foldchange),
#         significant = significant
#     )
#     
#     #Multiple testing correction
#     apply_and_report_BH(out.pval, output_folder_temp, output_file_name_Rmd, threshold_pval)
#     
#     #Knit report
#     output_file_name_md <- paste(c(output_folder,'/output_',all.methods.label,'.md'), collapse='')
#     output_file_name_html <- paste(c(output_folder,'/output_',all.methods.label,'.html'), collapse='')    
#     knit(output_file_name_Rmd, out=output_file_name_md)
#     knitr::knit2html(output_file_name_md, output=output_file_name_html)
# 
#     #Back to precedent working directory
#     setwd(oldwd)
#     
# }