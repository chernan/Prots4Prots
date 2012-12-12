#' SAM library for R
#' 
#' 
#' 

# source("http://bioconductor.org/biocLite.R")
# biocLite("impute")
# install.packages("samr")

# For reproducible results
set.seed(0)
# SAM library for R
library("samr")



#' Function apply_samr_method
#' two class unpaired comparison
apply_samr_method <- function(experiment, control, threshold_pval, rownames) {
    
    data <- data.frame(control, experiment)
    row.names(data) <- rownames
    
    design <- c(rep(1,ncol(control)), rep(2,ncol(experiment)))
    
#     dataset <- list(x=data, y=design, 
#                  geneid=rownames,
#                  genenames=rownames, 
#                  logged2=TRUE)
#     
#     samr.obj<-samr(dataset, 
#                    resp.type="Two class unpaired", nperms=100)

    sam.fit <- SAM(x=data, y=design, logged2=TRUE, 
                 resp.type="Two class unpaired", fdr.output=threshold_pval,
                   genenames = rownames, geneid = rownames)
    
    p.values <- samr.pvalues.from.perms(sam.fit$samr.obj$tt, sam.fit$samr.obj$ttstar)
    
    significant <- matrix(rep('STABLE',nrow(experiment)))
    names(significant) <- rownames

    #Get 2nd column, but no colnames?? but a print shows the col.names!
    #Header:
    #Gene ID  Gene Name  Score(d)  Numerator(r)  Denominator(s+s0)  Fold Change  q-value(%)
    
    # 2* Special cases for '1' necessary, otherwise get an error:
    #  Erreur dans sam.fit$siggenes.table$genes.lo[1, 2] : 
    #      nombre de dimensions incorrect
    
    if(sam.fit$siggenes.table$ngenes.lo == 1) {
        significant[sam.fit$siggenes.table$genes.lo[2]] <- 'LO'
    }
    if(sam.fit$siggenes.table$ngenes.lo > 1) {
        pos <- sam.fit$siggenes.table$ngenes.lo
        significant[sam.fit$siggenes.table$genes.lo[1:pos,2]] <- 'LO'
    }
    
    if(sam.fit$siggenes.table$ngenes.up == 1 ) {
        significant[sam.fit$siggenes.table$genes.up[2]] <- 'UP'
    }
    if(sam.fit$siggenes.table$ngenes.up > 1 ) {
        pos <- sam.fit$siggenes.table$ngenes.up
        significant[sam.fit$siggenes.table$genes.up[1:pos,2]] <- 'UP'
    }
    
    return.values <- data.frame(
        p.values = p.values,
        fold.change = log2(sam.fit$samr.obj$foldchange),
        significant = significant
        )
    return(return.values)
}

#' Function apply_paired_samr_method
#' two class paired comparison
apply_paired_samr_method <- function(experiment, control, threshold_pval, rownames) {
    
    data <- as.matrix(data.frame(control, experiment))
    row.names(data) <- rownames
    
    design <- rep(c(1:ncol(control)), times=2)
    design[1:ncol(control)] <- design[1:ncol(control)]*-1
    
    sam.fit <- SAM(x=data, y=design, logged2=TRUE, 
                   resp.type="Two class paired", fdr.output=threshold_pval,
                   genenames = rownames, geneid = rownames,
                   )
    
    p.values <- samr.pvalues.from.perms(sam.fit$samr.obj$tt, sam.fit$samr.obj$ttstar)
    
    significant <- matrix(rep('STABLE',nrow(experiment)))
    names(significant) <- rownames
    
    #Get 2nd column, but no colnames?? but a print shows the col.names!
    #Header:
    #Gene ID  Gene Name  Score(d)  Numerator(r)  Denominator(s+s0)  Fold Change  q-value(%)
    
    # 2* Special cases for '1' necessary, otherwise get an error:
    #  Erreur dans sam.fit$siggenes.table$genes.lo[1, 2] : 
    #      nombre de dimensions incorrect
    
    if(sam.fit$siggenes.table$ngenes.lo == 1) {
        significant[sam.fit$siggenes.table$genes.lo[2]] <- 'LO'
    }
    if(sam.fit$siggenes.table$ngenes.lo > 1) {
        pos <- sam.fit$siggenes.table$ngenes.lo
        significant[sam.fit$siggenes.table$genes.lo[1:pos,2]] <- 'LO'
    }
    
    if(sam.fit$siggenes.table$ngenes.up == 1 ) {
        significant[sam.fit$siggenes.table$genes.up[2]] <- 'UP'
    }
    if(sam.fit$siggenes.table$ngenes.up > 1 ) {
        pos <- sam.fit$siggenes.table$ngenes.up
        significant[sam.fit$siggenes.table$genes.up[1:pos,2]] <- 'UP'
    }
    
    return.values <- data.frame(
        p.values = p.values,
        fold.change = log2(sam.fit$samr.obj$foldchange),
        significant = significant
    )
    return(return.values)
}


apply_and_report_samr <- function(experiment, control, output_folder_temp, output_file, threshold_pval, is.paired=FALSE) {
    if(!is.paired) {
        samr.vals <- apply_samr_method(experiment, control, threshold_pval, rownames=rownames(experiment))
    } else {
        samr.vals <- apply_paired_samr_method(experiment, control, threshold_pval, rownames=rownames(experiment))
    }
        
    sink(output_file, append=TRUE)
    cat('
Significantly different proteins
--------------------------------------------------------------------------------------------------------

Data analysis was performed using SAM.

')

    tempoutput <- paste(c(output_folder_temp, '/samr_tests_significance_Rmd_', format(Sys.time(), "%Y%m%d%H%M%S"),trunc(runif(1)*10000),'.txt'), collapse='')
    write.table(samr.vals, tempoutput, sep="\t") 
    
    cat('
```{r, echo=FALSE, fig.width=14, fig.height=10}
')
    
    cat('display_SAM_significance_volcanoplot <- ')
    print(display_SAM_significance_volcanoplot)
    
    cat( 
        paste(c('matrixdata <- read.table("',tempoutput,'", stringsAsFactors=FALSE)'), collapse=''),
        paste(c('display_SAM_significance_volcanoplot(matrixdata, threshold_pval=',threshold_pval,')'), collapse=''),
        sep = "\n"
    )
    
    cat('
```
')

    cat(
        paste(c("SAM found <b>`r length(which(matrixdata[[\"significant\"]]=='UP'))`</b> significant up-regulated and <b>`r length(which(matrixdata[[\"significant\"]]=='LO'))`</b> down-regulated protein groups, with an estimated FDR of ",threshold_pval,"."))
    )
    
    cat('        
> 
> SAM citation.
>

---------------------------------------------------------------------------
')        
    sink(file=NULL)
    
    return(samr.vals)
}

display_SAM_significance_volcanoplot <- function(matrixdata, threshold_pval, title="Volcano plot") {
    
    pvals <- matrixdata[,"p.values"]
    foldchange <- matrixdata[,"fold.change"]
    samr.significant <- matrixdata[,"significant"]
    
    plot( foldchange, -log2(pvals),
          main=title, sub=paste("(FDR = ",threshold_pval,")"),
          xlab="log2( Fold change )", 
          ylab="-log2( p-value )",
          col="gray", 
          pch=16, cex.lab = 1, cex.axis = 1, cex.main = 1)
    grid(col="lightgray")
    sd_foldchange <- sd(foldchange)
    abline(v=c(c(-2,-1,0,1,2)*sd_foldchange), col="gray", lty=2)
        
    if(any(samr.significant=='UP') | any(samr.significant=='LO')) {
        points( foldchange[samr.significant=='UP'], -log2(pvals[samr.significant=='UP']), pch=16, col="darkred")
        points( foldchange[samr.significant=='LO'], -log2(pvals[samr.significant=='LO']), pch=16, col="darkgreen")
    }
    
    legend("bottomright", 
           legend=c("Non significant", "Significant up-regulated", "Significant down-regulated", paste("Significance threshold (FDR=",threshold_pval,")")), 
           col=c("gray", "darkred", "darkgreen", "orange"),
           fill=c("gray", "darkred", "darkgreen", 0),
           border=c("gray", "gray", "gray", 0),
           lty = c(0, 0, 0, 2),
           merge = TRUE)
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