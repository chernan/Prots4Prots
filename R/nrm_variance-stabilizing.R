#' @name VSN
#' @title Variance Stabilizing Normalization
#' @description
#' Variance Stabilizing Normalization (VSN) is based on a generalized logarithm 
#' transformation and aims at stabilizing variance inside microarray data.
#' 
#' @keywords Normalization
#' 
#' @section Introduction:
#' Usage of VSN was suggested by Karp (2010) but its usage is not common as the 
#' scale change can make the fold change ratio difficult to interpret, and 
#' applying it can lead to negative values. This is well summarized in 
#' Mahoney (2011) :
#' "Karp et al proposed a variance stabilizing transformation to stabilize the 
#' variance across all proteins[15]. This has been applied to other data types 
#' such as mRNA and miRNA[33-35]. We prefer use of WLS for the following reasons: 
#' 1) constant variance is generally not expected across all proteins in a study; 
#' 2) downstream statistical analyses assume equal variances between groups 
#' within a protein but do not assume constant variance across all proteins; 
#' and 3) the scale is constant across all proteins with WLS. The interpretation 
#' is difficult after the proposed transformation, with fold changes on the 
#' transformed scale interpreted on the raw scale at low intensities, log scale 
#' at high intensities, and a sliding hybrid of these two scales over the 
#' continuum of middle range abundance values."
#' 
#' #' Usage note:
#' Default settings is lts.quantile=0.9 but lts.quantile=0.5 is more robust. 
#' The reason why lts.quantile=0.5 is not the default is that the estimator 
#' with lts.quantile=0.9 is more efficient (more precise with less data) 
#' if the fraction of differentially expressed genes is not that large.
#' From Bioconductor case studies.
#' 
#' @section References:
#' 
#' Huber, W., Von Heydebreck, A., Sultmann, H., Poustka, A., & Vingron, M. (2002). 
#' Variance stabilization applied to microarray data calibration and to the 
#' quantification of differential expression. 
#' Bioinformatics, 18(Suppl 1), S96-S104. 
#' doi:10.1093/bioinformatics/18.suppl_1.S96
#' 
#' Karp, N. A., Huber, W., Sadowski, P. G., Charles, P. D., Hester, S. V, 
#' Lilley, K. S. (2010). 
#' Addressing accuracy and precision issues in iTRAQ quantitation. 
#' Molecular & cellular proteomics : MCP, 9(9), 1885-97. 
#' doi:10.1074/mcp.M900628-MCP200
#' 
#' Mahoney, D. W., Therneau, T. M., Heppelmann, C. J., Higgins, L., 
#' Benson, L. M., Zenka, R. M., Jagtap, P., et al. (2011). 
#' Relative quantification: characterization of bias, variability and fold 
#' changes in mass spectrometry data from iTRAQ-labeled peptides. 
#' Journal of proteome research, 10(9), 4325-33. 
#' doi:10.1021/pr2001308
#' 
#' http://bioinfo.cnio.es/files/training/Microarray_Course/3_UBio_Normalization_Theory.pdf
#' 
#' @import vsn
# library(vsn)
NULL

#' @title Apply variance-stabilizing normalization.
#' 
#' @description Apply variance-stabilizing normalization.
#' 
#' @details 
#' Apply variance-stabilizing normalization using generalized logarithm. Use 
#' default parameter: lts.quantile = 0.9
#'  
#' @param dataNotNorm Data to be normalized.
#' @return A dataframe containing normlized data.
applyVSN09 <- function(dataNotNorm) {

    dataNormalized <- data.frame(
        justvsn(as.matrix(dataNotNorm), verbose=FALSE)
    ) #, backgroundsubtract=FALSE ?

    # equivalent to :
    #fit <- vsn2(as.matrix(dataNotNorm), subsample=100000)
    #x <- predict(fit, newdata=as.matrix(dataNotNorm))
    
    return(dataNormalized)
}

#' @title Apply robust variance-stabilizing normalization.
#' 
#' @description Apply robust variance-stabilizing normalization.
#' 
#' @details 
#' Apply variance-stabilizing normalization using generalized logarithm. Use 
#' parameter: lts.quantile = 0.5, more robust to high number of outliers.
#'  
#' @param dataNotNorm Data to be normalized.
#' @return A dataframe containing normlized data.
applyVSN05 <- function(dataNotNorm) {
    
    dataNormalized <- data.frame(
        justvsn(as.matrix(dataNotNorm), verbose=FALSE, lts.quantile=0.5)
    )
    
    return(dataNormalized)
}

#' @title Apply VSN and report result
#' 
#' @description Apply VSN and report result
#' 
#' @details 
#' Apply VSN on a dataset. 
#' TODO Use ExpressionSet object as input
#' 
#' @param dataNotNorm Dataset to be normalized.
#' @param outputFolderTemp Temporary folder to store data generated at this 
#'  step of the analysis.
#' @param outputFile Report of current analysis step will be appended to 
#'  this Rmd file.
#' @param useRobustFit If the more robust version of the method should be used.
#' @return A data frame with normalized data.
applyAndReportVSN <- function(dataNotNorm, 
                                 outputFolderTemp, outputFile, 
                                 useRobustFit=FALSE) {
    
    # Apply VSN
    dataNormalized <- matrix(c(''))
    if(useRobustFit==FALSE) {
        dataNormalized <- applyVSN09(dataNotNorm)
    }
    else {
        dataNormalized <- applyVSN05(dataNotNorm)
    }
    
    # Generate report
    
    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc( 
            runif(1) * 10000)), 
        collapse='')
    
    cat('',
        'Normalization (VSN)',
        '---------------------------------------------------------------------',
        '',
        'Normalisation was achieved using the Variance Stabilizing method (R package version `r packageDescription("vsn")$Version`).',
        '',
        ifelse(
            useRobustFit==FALSE,
            "Default settings (lts.quantile=0.9) were used.",
            "Settings more robust to high number of outliers were used (lts.quantile=0.5)."
        ),
        sep="\n", file=outputFile, append=TRUE)
    
    allTestsNormalizationRmd(dataNormalized, outputFile=outputFile,
                             outFolder=outputFolderTemp)
    
    tempOutput <- paste(
        c(outputFolderTemp, '/vsn_tests_normalization_Rmd_data_', 
          execLabel, '.txt'), 
        collapse='')
    write.table(dataNormalized, tempOutput, sep="\t")
    cat('',
        'Specific test for VSN.',
        '',
        paste(
            c('```{r applyAndReportVSN2', execLabel, 
              ', echo=FALSE, fig.width=10, fig.height=6}'),
            collapse=''),
        '',
        sep="\n", file=outputFile, append=TRUE
    )
    cat( 
        paste(
            c('matrixdata <- as.matrix(read.table("', 
              tempOutput, 
              '", stringsAsFactors=FALSE))'), 
            collapse=''
        ),
        'meanSdPlot(matrixdata)',
        '',
        '```',
        '',
        sep="\n", file=outputFile, append=TRUE
    )
        
    cat(' ',
        '> ',
        '> Please cite these articles whenever using results from this software in a publication :',
        '> ',
        '> Method article :',
        '> Variance stabilization applied to microarray data calibration and to the quantification of differential expression. ',
        '> Huber, W., Von Heydebreck, A., Sultmann, H., Poustka, A., & Vingron, M. (2002).',
        '> Bioinformatics, 18(Suppl 1), S96-S104. ',
        '> doi:10.1093/bioinformatics/18.suppl_1.S96',
        '> ',
        '> Software article (R package) :',
        '> `r citation("vsn")$textVersion`',
        '> ',
        '> Usage for quantitative proteomics data suggested by :',
        '> Karp, N. A., Huber, W., Sadowski, P. G., Charles, P. D., Hester, S. V, Lilley, K. S. (2010). ',
        '> Addressing accuracy and precision issues in iTRAQ quantitation. ',
        '> Molecular & cellular proteomics : MCP, 9(9), 1885-97. ',
        '> doi:10.1074/mcp.M900628-MCP200',
        ' ',
        '---------------------------------------------------------------------',
        ' ',
        sep="\n", file=outputFile, append=TRUE)
    
    return(dataNormalized)
}



