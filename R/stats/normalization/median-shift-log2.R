#' Median Normalization
#' 
#' R package
#' Article


#' 8-plex only
apply_median_norm <- function(data_nn) {
    
    experiment <- exprs(data_nn[,(data_nn$Design == "Experiment")])
    control <- exprs(data_nn[,(data_nn$Design == "Control")])

    median.exp21 <- median(experiment[,2]/experiment[,1])
    median.exp31 <- median(experiment[,3]/experiment[,1])
    median.exp41 <- median(experiment[,4]/experiment[,1])
    median.ctr21 <- median(control[,2]/control[,1])
    median.ctr31 <- median(control[,3]/control[,1])
    median.ctr41 <- median(control[,4]/control[,1])
    
    
    data_normalized <- data.frame(
        V1 = log(experiment[,1], base=2),
        V2 = log(experiment[,2]/median.exp21, base=2),
        V3 = log(experiment[,3]/median.exp31, base=2),
        V4 = log(experiment[,4]/median.exp41, base=2),
        
        V5 = log(control[,1], base=2),
        V6 = log(control[,2]/median.ctr21, base=2),
        V7 = log(control[,3]/median.ctr31, base=2),
        V8 = log(control[,4]/median.ctr41, base=2)
        )
       
    return(data_normalized)
}

apply_and_report_median_norm <- function(dataset, output_folder_temp, output_file) {
    dataset.norm <- apply_median_norm(dataset)
    
    sink(output_file, append=TRUE)
    cat('
Normalization
-------------------------

Intensities were corrected for median of log ratios compared to first reported.
Base 2 logarithm was applied on the intensities.

')

    all_tests_normalization_Rmd(dataset.norm, title="Median normalization", outfolder=output_folder_temp)

cat('

---------------------------------------------------------------------------
')        
    sink(file=NULL)
    
    return(dataset.norm)
}

