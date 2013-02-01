## Setup environment ###########################################################

## To avoid transformation of small numbers into scientific notation causing 
## incompatibilities with further sorting methods (googleVis by ex.)
oldScipen <- options(scipen=500)

## Set output directory as working directory. All output files/images/etc. 
## will be created there
setwd('/home/chernan/Workspace/testRpackage2/Prots4Prots/inst/protocols/quantitation/Label-free')
oldwd <- getwd()

################################################################################

## List datasets to be analyzed ################################################
source("./load_data.R")
dataset_list <- getDatasetDpulex246daysDMSOM4pg()


## Apply analysis methods ######################################################

basedir <- '/home/chernan/Workspace/testRpackage2/Prots4Prots/inst/protocols/quantitation/Label-free/outputs/'
dir.create(basedir)
setwd(basedir)

all.methods <- list(
    Normalization = c("vsn05"),
    Significance = c("lpe")
)

#Combine all (normalization, test) methods
combined_methods_prot <- as.matrix(expand.grid(all.methods))

# Reproducible analysis
set.seed(0)
threshold_pval <- 0.05

library(Prots4Prots)

#For each dataset, apply all possible combination of methods
output <- lapply(
    names(dataset_list), 
    FUN=function(dataset.name) {
        working.dir <- paste(c(getwd(),'/', threshold_pval, '/',dataset.name), 
                             collapse='')
        dir.create(working.dir, recursive=TRUE)        
        currentdata <- dataset_list[[dataset.name]]
        
        out.val <- apply(
            combined_methods_prot, 1, 
            FUN = function(methods) {
                dataset <- currentdata[['dataset']]
                all.methods.label <- paste(c(methods["Normalization"], 
                                             methods["Significance"]), 
                                           collapse='_')
                
                ## Setup environment ####
                print("Setup")
                #Create working directory
                output_folder <- paste(c(working.dir,'/', all.methods.label), 
                                       collapse='')
                dir.create(output_folder)
                #Create temp folder for temp files used in Rmd report
                output_folder_temp <- paste(c(output_folder,'/temp'), 
                                            collapse='')
                dir.create(output_folder_temp)
                #Create Rmd report file
                output_file_name_Rmd <- paste(c(output_folder, '/output_', 
                                                all.methods.label, '.Rmd'), 
                                              collapse='')
                write(c(''),output_file_name_Rmd)
                
                print("Quality check 1")
                reportQualityCheckRawData(
                    dataset=exprs(dataset), 
                    outputFolderTemp=output_folder_temp, 
                    outputFile=output_file_name_Rmd, 
                    distMethod="euclidean")
                
                print("Treating missing values")
                dataset <- filterDataset(dataset)
                
                print("Normalization")
                out.val <- applyNormalization(
                    normMethod=methods["Normalization"], 
                    dataset=dataset, 
                    outputFolderTemp=output_folder_temp, 
                    outputFileNameRmd=output_file_name_Rmd)
                
                print("Quality check 2")
                reportQualityCheckRawData(
                    dataset=out.val, 
                    outputFolderTemp=output_folder_temp, 
                    outputFile=output_file_name_Rmd, 
                    distMethod="euclidean")
                
                print("Summarization")
                pg.eset <- applySummarization(
                    summMethod="NonA", 
                    dataset=dataset, 
                    outVal=out.val, 
                    outputFolderTemp=output_folder_temp, 
                    outputFileNameRmd=output_file_name_Rmd)
                
                #Get data to compare, depending on experimental design
                experiment <- exprs(pg.eset[, (pg.eset$Design == "Experiment")])
                control <- exprs(pg.eset[, (pg.eset$Design == "Control")])
                
                print("Summary statistic")
                out.pval <- applySummaryStatistic(
                    signMethod=methods["Significance"], 
                    experiment=experiment, control=control, 
                    outputFolderTemp=output_folder_temp, 
                    outputFileNameRmd=output_file_name_Rmd, 
                    thresholdPVal=threshold_pval)
                
                print("Multiple testing correction")
                out.pval.corr <- applyMultipleTestingCorrection(
                    mtcMethod="bh", matrixData=out.pval, 
                    outputFolderTemp=output_folder_temp, 
                    outputFileNameRmd=output_file_name_Rmd, 
                    thresholdPVal=threshold_pval)
                
                #Append copyright to report 
                
                #Knit report
                print("Report")
                output_file_name_html <- generateKnitrReport(
                    outputLabel=all.methods.label, 
                    outputFileNameRmd=output_file_name_Rmd, 
                    outputFolder=output_folder)
                
                #Export summarized data as tsv for clients/us
            }
        )
    }
)

#Back to precedent working directory
setwd(oldwd)

