library(Biobase)


#'#' Load data
getDatasetDpulex246daysDMSOM4pg <- function() {
    
    data_file <- c('/home/chernan/Workspace/testRpackage2/Prots4Prots/inst/protocols/quantitation/Label-free/data/proteinGroups.txt')
    
    ## Load data from MQ iTRAQ
    temp_proteins <- read.delim(data_file, stringsAsFactors=FALSE, 
                                row.names=NULL, header=TRUE)
    
    ## Select intensities
    reporterIntHeaders2d <- c(
        "LFQ.intensity.DMSO.2d.G", "LFQ.intensity.DMSO.2d.H", 
        "LFQ.intensity.M4.2d.A", "LFQ.intensity.M4.2d.B"
    )
    exprs2d <- temp_proteins[, reporterIntHeaders2d]
    reporterIntHeaders4d <- c(
        "LFQ.intensity.DMSO.4d.I", "LFQ.intensity.DMSO.4d.J",
        "LFQ.intensity.M4.4d.C", "LFQ.intensity.M4.4d.D"
    )
    exprs4d <- temp_proteins[, reporterIntHeaders4d]
    reporterIntHeaders6d <- c(
        "LFQ.intensity.DMSO.6d.K", "LFQ.intensity.DMSO.7d.L", 
        "LFQ.intensity.M4.7d.E", "LFQ.intensity.M4.7d.F"
    )
    exprs6d <- temp_proteins[, reporterIntHeaders6d]
    
    ## Phenotypic data
    pData2d <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData2d) <- reporterIntHeaders2d
    pData4d <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData4d) <- reporterIntHeaders4d
    pData6d <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData6d) <- reporterIntHeaders6d
    
    ## Feature data
    infoHeaders <- c("Majority.protein.IDs", "Fasta.headers", "PEP", 
                     "Peptides", "Only.identified.by.site", "Reverse", 
                     "Contaminant")
    featuredata <- temp_proteins[, infoHeaders]
    names(featuredata) <- c("Majority.protein.IDs", "Fasta.headers", "PEP", 
                            "Peptides", "Only.identified.by.site", "Reverse", 
                            "Contaminant")
    
    ## ExpressionSet construction
    data_proteins2d <- new("ExpressionSet", 
                           exprs = exprs2d, 
                           phenoData = new("AnnotatedDataFrame", pData2d), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    data_proteins4d <- new("ExpressionSet", 
                           exprs = exprs4d, 
                           phenoData = new("AnnotatedDataFrame", pData4d), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    data_proteins6d <- new("ExpressionSet", 
                           exprs = exprs6d, 
                           phenoData = new("AnnotatedDataFrame", pData6d), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    
    return(
        list(
            dpulex2dDMSOM4 = list(dataset=data_proteins2d, type="prot"),
            dpulex4dDMSOM4 = list(dataset=data_proteins4d, type="prot"),
            dpulex6dDMSOM4 = list(dataset=data_proteins6d, type="prot")
        )
    )
    
}

filterDataset <- function(expSet, threshold=0) {
    
    ## Remove contaminants
    whichOkProteins <- 
        pData(featureData(expSet))[["Only.identified.by.site"]] != '+' & 
        pData(featureData(expSet))[["Reverse"]] != '+' &
        pData(featureData(expSet))[["Contaminant"]] != '+'
    
    ## Remove non-complete observations
    whichOkIntensitites <- 
        apply(exprs(expSet), 1, FUN=function(x){all(!is.na(x))}) &
        apply(exprs(expSet), 1, FUN=function(x){all(x>threshold)})   
    
    return(expSet[(whichOkProteins & whichOkIntensitites), ])
    
}

## Save as tsv
exportDatasetDpulexDMSOM4pg <- function(expSet) {
    
    orderprot <- order( row.names(all_intensities_PEP_sorted), decreasing=FALSE)
    ratioorder <- order( row.names(fold_change), decreasing=FALSE)
    lpe.val.order <- order( row.names(lpe.val), decreasing=FALSE)
    
    results <- cbind(
        all_intensities_PEP_sorted[orderprot,] ,
        fold_change[ratioorder,],
        lpe.val[lpe.val.order,]
    )
    
    #Modify format of UniProt IDs for further displays
    results[,"Majority.protein.IDs"] <- gsub(";", "; ", results[,"Majority.protein.IDs"])
    # round z-values
    # lpe.val <- round(lpe.val, digits=2)
    
    #Improve headers
    final_headers <- c(
        "Majority.protein.IDs", "Peptides", "PEP", "Log2(fold_change)",
        "flag.outlier.x", "flag.outlier.y",
        "p.vals", "lpe.fdr.BH"
    )
    client_headers <- c(
        "UniProt IDs", "Peptides", "MaxQuant PEP", "Log2(fold_change)", 
        "Outlier DMSO", "Outlier M4",
#         "Outlier DMSO", "Outlier C1",
#         "Outlier DMSO", "Outlier C2",
        "p-value", "p-value BH corrected"
    )
    client_results <- data.frame(results[,final_headers])
    names(client_results) <- client_headers
    
    #For us
    write.csv(results, "MaxQuant_LFQ_6-plex_non-linked_vsnlpe_full_DMSOC1.csv")
    #For clients
    write.csv(client_results, "MaxQuant_LFQ_6-plex_non-linked_vsnlpe_final_DMSOC1.csv")
    
}
