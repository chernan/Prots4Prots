#' Significance tests
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/03_significance/significance_tests.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/03_significance/t-test.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/03_significance/local-pooled-error.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/03_significance/samr.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/03_significance/limma.R")

#' Summary statistic
#' 
#' ttest : t-test
#' lpe : local-pooled-error
#' samr : significance analysis of microarrays
#' limma : linear models for microarray analysis
apply_summary_statistic <- function(sign.method, experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval) {
    out.pval <- switch(sign.method,
                       ttest = apply_and_report_t_test(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval, is.paired=FALSE),
                       lpe = apply_and_report_lpe(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval),
                       samr = apply_and_report_samr(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval, is.paired=FALSE),
                       limma = apply_and_report_limma(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval, is.paired=FALSE)
    )
    return(out.pval)
}

#' Paired statistics
#' 
#' ttest : paired t-test
#' lpe : local-pooled-error
#' samr : significance analysis of microarrays
#' limma : linear models for microarray analysis
apply_paired_summary_statistic <- function(sign.method, experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval) {
    out.pval <- switch(sign.method,
                       ttest = apply_and_report_t_test(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval, is.paired=TRUE),
#                        lpe = apply_and_report_lpe(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval),
                       samr = apply_and_report_samr(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval, is.paired=TRUE),
                       limma = apply_and_report_limma(experiment, control, output_folder_temp, output_file_name_Rmd, threshold_pval, is.paired=TRUE)
    )
    return(out.pval)
}