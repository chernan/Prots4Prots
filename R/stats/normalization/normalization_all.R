#' Normalizations
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/01_normalization/normalization_test.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/01_normalization/logarithm2.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/01_normalization/inter-quartile-range.R")
source("/home/chernan/Workspace/GitHub/Prots4Prots/Prots4Prots/R/01_normalization/variance-stabilizing.R")


#' Normalization methods
#' 
#' log2 : log2 transformation (for comparison)
#' iqr : inter-quartile range normalization
#' vsn : variance stabilizing normalization
#' vsn05 : variance stabilizing normalization (robust to outliers)
apply_normalization <- function(norm.method, dataset, output_folder_temp, output_file_name_Rmd) {
    
    out.val <- switch(norm.method,
                      log2 = apply_and_report_log2(exprs(dataset), output_folder_temp, output_file_name_Rmd),
                      iqr = apply_and_report_iqr(exprs(dataset), output_folder_temp, output_file_name_Rmd),
                      vsn = apply_and_report_vsn(exprs(dataset), output_folder_temp, output_file_name_Rmd),
                      vsn05 = apply_and_report_vsn(exprs(dataset), output_folder_temp, output_file_name_Rmd, use.robust.fit=TRUE)
#                       limmaBetween = apply_and_report_limmaNorm(exprs(dataset), output_folder_temp, output_file_name_Rmd, betweenArrays=TRUE),
#                       limmaWithin = apply_and_report_limmaNorm(exprs(dataset), output_folder_temp, output_file_name_Rmd, betweenArrays=FALSE)
#                       vsnrma = apply_and_report_rma(dataset, output_folder_temp, output_file_name_Rmd)
    )
    
    return(out.val)
}
