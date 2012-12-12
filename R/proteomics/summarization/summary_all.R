#' Summary
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/02_summarization/summary.R")

#' Summarization
#' 
#' NonA : non applicable. Doesn't summarize. Ex. of protein-level data.
#' meanPep : mean of peptides
#' medianPep : median of peptides
#' trimmed20Pep : trimmed mean (20%)
#' variablePep : most variable peptide
#' intensePep : most intense peptide
#' sumIntens : sum of all intensities for each reporter
apply_summarization <- function(summ.method, dataset, out.val, output_folder_temp, output_file_name_Rmd) {
    pg.eset <- switch(summ.method,
                      "NonA" = apply_and_report_NA(dataset, out.val, output_folder_temp, output_file_name_Rmd),
                      meanPep = apply_and_report_mean(dataset, out.val, output_folder_temp, output_file_name_Rmd),
                      medianPep = apply_and_report_median(dataset, out.val, output_folder_temp, output_file_name_Rmd),
                      trimmed20Pep = apply_and_report_trimmed20(dataset, out.val, output_folder_temp, output_file_name_Rmd),
                      variablePep = apply_and_report_mostvariablepep(dataset, out.val, output_folder_temp, output_file_name_Rmd),
                      intensePep = apply_and_report_mostintensepep(dataset, out.val, output_folder_temp, output_file_name_Rmd),
                      sumIntens = apply_and_report_sumintens(dataset, out.val, output_folder_temp, output_file_name_Rmd)
    )
    
}