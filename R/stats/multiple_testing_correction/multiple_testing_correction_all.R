#' Multiple testing correction
#' Display methods are also in significance_tests.R
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/03_significance/significance_tests.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/04_multiple_testing_correction/benjamini-hochberg.R")


#' Multiple testing correction
#' 
#' bh : Banjamini-Hochberg
apply_multiple_testing_correction <- function(mtc.method, out.pval, output_folder_temp, output_file_name_Rmd, threshold_pval) {
    mtc.pval <- switch(mtc.method,
                       bh = apply_and_report_BH(out.pval, output_folder_temp, output_file_name_Rmd, threshold_pval)
    )
    return(mtc.pval)
}