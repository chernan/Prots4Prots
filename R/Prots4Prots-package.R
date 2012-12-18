#' Protocols for Proteomics
#'
#' A package comparing different statistical tools in order to find 
#' differentially expressed proteins. 
#' 
#' The R folder contains the toolbox, tools brought together by topic. 
#' The inst folder contains use cases of the toolbox.
#' 
#' Quantitative Proteomics
#' 
#' This package divides adata analysis pipeline into multiple steps:
#' \enumerate{
#'   \item Loading data.
#'   \item Pre-processing
#'     \itemize {
#'       \item Normalization
#'       \item Treating missing values
#'       \item (Summarization)
#'     }
#'   \item Statistics
#'   \item Multiple testing correction
#'   \item Selection of Differentially Expressed proteins
#'   \item Inclusion of other sources of information.
#' }
#'
#' Each step of the pipeline should be able to generate a report paragraph, 
#' itself transformed in the end into a HTML report through the use of 
#' \code{\link{knitr}}.
#' 
#' Loading data (LDD)
#' A function should be provided to load data into an ExpressionSet object.
#' The exprs part of the object sould contain the intensities for each channel.
#' The experimental design (which channels are Control, and which are 
#' Experiments)
#' 
#'  (Treating missing values) (TMS)
#'  Only option available for now: remove incomplete information.
#'  Imputation should be done after normalization!
#'  TODO: check that normalization mehtods can accept missing values.
#'  
#' Pre-processing
#' 
#'  Normalization (NRM)
#'  This step includes the usage of Variance Stabilizing Normalization, VSN 
#'  with a more robust parameter, and the Inter Quartile Range normalization.
#'  As a comparison, transformation by logarithm, base 2, is also provided.
#'  
#'  (Summarization) (SUM)
#'  This step is optional if the peptide intensities have already been 
#'  summarized at the protein group level.
#'  Possible summarizations include mean, trimmed mean (20%), median, most 
#'  intense peptide, most variable peptide.
#' 
#' Statistics (STT)
#' Available statistics include:
#' \enumerate{
#'   \item t-test (Welsh)
#'   \item Local-pooled error
#'   \item LIMMA
#'   \item Significance Analysis of Microarrays
#' }
#' 
#' Multiple testing correction (MTC)
#' Only option available for now: Benjamini-Hochberg FDR.
#' 
#' Selection of Differentially Expressed proteins (SDE)
#' All methods use a static threshold for the p-value.
#' SAM method uses bootstraping.
#' 
#' Other categories for utilities functions:
#' Generation of reports (RPT)
#'
#' @docType package
#' @name Prots4Prots 
NULL