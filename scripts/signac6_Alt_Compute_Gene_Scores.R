#Computing gene scores 
#With adaptations from Granja et al. (2019) 

#Multimodal profiling of the transcriptional regulatory landscape
#of developing mouse cortex identifies Neurog2 as a key epigenome remodeler
#Cite Noack et al. (2022) 

library(cicero)
library(Signac)
library(Seurat)
library(monocle3)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
set.seed(1)

source('scripts/figures/config.R')
source('scripts/aux_functions.R')

####################################################
#Functions
####################################################

grToFeature <- function(gr){
    peakinfo <- data.frame(
        row.names = paste(seqnames(gr),start(gr),end(gr),sep="_"),
        site_name = paste(seqnames(gr),start(gr),end(gr),sep="_"),
        chr = gsub("chr","",as.character(seqnames(gr))),
        bp1 = start(gr),
        bp2 = end(gr)
    )
    return(peakinfo)
}

featureToGR <- function(feature){
    featureSplit <- stringr::str_split(paste0(feature), pattern = "_", n = 3, simplify = TRUE)
    gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
    return(gr)
}

makeCDS <- function(se,mdata, binarize = TRUE){
    peakinfo <- StringToGRanges(rownames(se), sep = c("_", "_"))
    peakinfo <- grToFeature(peakinfo)
    mat <- se
    if(binarize){
        mat@x[which(mat@x > 0)] <- 1
    }
    cellinfo <- mdata
    cellinfo$cells <- rownames(cellinfo)
    cds <-  suppressWarnings(newCellDataSet(mat,
                              phenoData = methods::new("AnnotatedDataFrame", data = cellinfo),
                              featureData = methods::new("AnnotatedDataFrame", data = peakinfo),
                              expressionFamily=negbinomial.size(),
                              lowerDetectionLimit=0))
    Biobase::fData(cds)$chr <- as.character(Biobase::fData(cds)$chr)
    Biobase::fData(cds)$bp1 <- as.numeric(as.character(Biobase::fData(cds)$bp1))
    Biobase::fData(cds)$bp2 <- as.numeric(as.character(Biobase::fData(cds)$bp2))
    cds <- cds[order(Biobase::fData(cds)$chr, Biobase::fData(cds)$bp1),]
    return(cds)
}


sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX) > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY) > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);
    } 

    return(cor);

  }'
)

#Cleaned up custom version of cicero_cds
custom_cicero_cds <- function(
                            cds,
                            reduced_coordinates,
                            k=50,
                            max_knn_iterations = 5000,
                            summary_stats = NULL,
                            size_factor_normalize = TRUE,
                            silent = FALSE) {
  
  start <- Sys.time()
  # Create a k-nearest neighbors map
  message("\nFNN k-nearest search...")
  nn_map <- FNN::knn.index(reduced_coordinates, k=(k-1)) # no data.frame wrapper
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- k * 2 # Compute once
  
  # function for sapply
  get_shared <- function(other, this_choice) {
    k2 - length(union(cell_sample[other,], this_choice))
  }
  
  while (length(good_choices) > 0 & it < max_knn_iterations) { # slow
    if(it %% 100 == 0) message(sprintf("%s of %s iterations", it, max_knn_iterations))
    it <- it + 1
    choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen,]

    others <- seq_len(nrow(cell_sample) - 1)
    this_choice <- cell_sample[nrow(cell_sample),]
    shared <- sapply(others, get_shared, this_choice = this_choice)
    
    if (max(shared) < .9 * k) {
      chosen <- new_chosen
    }
  }
  message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))

  cell_sample <- nn_map[chosen,]
  cell_sample_map <- lapply(seq_len(nrow(cell_sample)), function(x) rownames(reduced_coordinates)[cell_sample[x,]]) %>% Reduce("rbind",.) %>% data.frame
  rownames(cell_sample_map) <- rownames(cell_sample)

  if(!silent) {
    # Only need this slow step if !silent
    combs <- combn(nrow(cell_sample), 2)
    combs <- combs[,sample(seq_len(ncol(combs)),min(ncol(combs),10^6))] #sample 1 M because this really doesnt matter
    shared <- apply(combs, 2, function(x) {  #slow
      k2 - length(unique(as.vector(cell_sample[x,])))
    })
    
    message(paste0("\nOverlap QC metrics:\nCells per bin: ", k,
                   "\nMaximum shared cells bin-bin: ", max(shared),
                   "\nMean shared cells bin-bin: ", mean(shared),
                   "\nMedian shared cells bin-bin: ", median(shared)))
    
    if (mean(shared)/k > .1) warning("On average, more than 10% of cells are shared between paired bins.")
  }

  message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))
  message("\nMaking aggregated scATAC Matrix...")
  exprs_old <- Biobase::exprs(cds)
  
  new_exprs <- matrix(0, nrow = nrow(cell_sample), ncol = nrow(exprs_old))
  for(x in seq_len(nrow(cell_sample))){
    if(x %% 50 == 0){
        message(sprintf("%s of %s iterations : %s minutes since start", x, nrow(cell_sample), round(difftime(Sys.time(),start,units="mins"),1)))
    }
    new_exprs[x,] <- Matrix::rowSums(exprs_old[,cell_sample[x,]])
  }
  remove(exprs_old)

  message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))
  message("\nMaking aggregated CDS...")
  pdata <- Biobase::pData(cds)
  new_pcols <- "agg_cell"
  if(!is.null(summary_stats)) { 
    new_pcols <- c(new_pcols, paste0("mean_",summary_stats)) 
  }
  
  new_pdata <- plyr::adply(cell_sample,1, function(x) {
    sub <- pdata[x,]
    df_l <- list()
    df_l["temp"] <- 1
    for (att in summary_stats) {
      df_l[paste0("mean_", att)] <- mean(sub[,att])
    }
    data.frame(df_l)
  })
  
  new_pdata$agg_cell <- paste("agg", chosen, sep="")
  new_pdata <- new_pdata[,new_pcols, drop = FALSE] # fixes order, drops X1 and temp
  
  row.names(new_pdata) <- new_pdata$agg_cell
  row.names(new_exprs) <- new_pdata$agg_cell
  new_exprs <- as.matrix(t(new_exprs))
  
  fdf <- Biobase::fData(cds)
  new_pdata$temp <- NULL
  
  fd <- new("AnnotatedDataFrame", data = as.data.frame(fdf))
  pd <- new("AnnotatedDataFrame", data = as.data.frame(new_pdata))
  
  cicero_cds <-  new_cell_data_set(new_exprs,cell_metadata = as.data.frame(new_pdata),gene_metadata = as.data.frame(fdf))
  
  cicero_cds <- monocle3::detect_genes(cicero_cds,min_expr = .1)
  cicero_cds <- cicero_cds[Matrix::rowSums(exprs(cicero_cds)) != 0,] 
  cicero_cds <- estimate_size_factors(cicero_cds)
  #cicero_cds <- suppressWarnings(BiocGenerics::estimateDispersions(cicero_cds))

  if (any(!c("chr", "bp1", "bp2") %in% names(fData(cicero_cds)))) {
    fData(cicero_cds)$chr <- NULL
    fData(cicero_cds)$bp1 <- NULL
    fData(cicero_cds)$bp2 <- NULL
    fData(cicero_cds) <- cbind(fData(cicero_cds),
                               df_for_coords(row.names(fData(cicero_cds))))
  }
  message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))

  if (size_factor_normalize) {
    message("\nSize factor normalization...")
    new_cicero_expr <- t(t(exprs(cicero_cds))/pData(cicero_cds)$Size_Factor)
    cicero_cds <-  new_cell_data_set(new_cicero_expr,cell_metadata = as.data.frame(new_pdata[new_pdata$agg_cell%in%colnames(new_cicero_expr),,drop=FALSE]),gene_metadata = as.data.frame(fdf))
    cicero_cds <- monocle3::detect_genes(cicero_cds,min_expr = .1)
    cicero_cds <- cicero_cds[Matrix::rowSums(exprs(cicero_cds)) != 0,] 
    cicero_cds <- estimate_size_factors(cicero_cds)
    if (any(!c("chr", "bp1", "bp2") %in% names(fData(cicero_cds)))) {
      fData(cicero_cds)$chr <- NULL
      fData(cicero_cds)$bp1 <- NULL
      fData(cicero_cds)$bp2 <- NULL
      fData(cicero_cds) <- cbind(fData(cicero_cds),
                                 df_for_coords(row.names(fData(cicero_cds))))
    }
  }

  message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))
  
  return(list(ciceroCDS = cicero_cds, knnMap = cell_sample_map))

}

getGeneGTF <- function(file){
    #Import
    message("Reading in GTF...")
    importGTF <- rtracklayer::import(file)
    #Exon Info
    message("Computing Effective Exon Lengths...")
    exonGTF <- importGTF[importGTF$type=="exon",]
    exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
    exonReduced <- unlist(exonList, use.names=TRUE)
    mcols(exonReduced)$gene_id <- names(exonReduced)
    mcols(exonReduced)$widths <- width(exonReduced)
    exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
    exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>% 
        unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
    #Gene Info
    message("Constructing gene GTF...")
    geneGTF1 <- importGTF[importGTF$type=="gene",]
    geneGTF2 <- GRanges(
            seqnames=paste0("chr",seqnames(geneGTF1)),
            ranges=ranges(geneGTF1),
            strand=strand(geneGTF1),
            gene_name=geneGTF1$gene_name,
            gene_id=geneGTF1$gene_id
        ) %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
    mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
    return(geneGTF2)
}

####################################################
#Input Data
####################################################

#Specific Genome Libraries
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
bsgenome <- BSgenome.Mmusculus.UCSC.mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
orgdb <- org.Mm.eg.db

#Read input summarized experiment peaks x cells
cortex.atac <- readRDS('data/merged_scATAC_integrated.RDS')

mdata <- cortex.atac@meta.data
DefaultAssay(cortex.atac) <- 'MACS2peaks'
#Reduced Dimensions
dimred <- Reductions(cortex.atac,slot = 'harmony')
dimred <- as.data.frame(dimred@cell.embeddings)

se <- cortex.atac@assays$MACS2peaks@counts
row.names(se) <- gsub('-','_',row.names(se))
#Get ChromSizes
chromSizes <- seqlengths(bsgenome)[paste0("chr",c(1:19,"X"))]
genome <- data.frame(names(chromSizes),chromSizes)
rownames(genome) <- NULL



#Genes GTF File from 10x v3
gtfFile <- '/home/hpc/bonev/annotations/mm10/cellranger_rna/genes/genes.gtf'

#Window around TSS to be called promoter
tssWindow <- 2500

#Flanking distance from TSS in KB for Co-Accessibility
flank <- 500*10^3

#Correlation Cutoff for Co-Accessibility
corCutOff <- 0.35

#Run Cicero
obj <- makeCDS(se,mdata, binarize = TRUE)
obj <- detectGenes(obj)
obj <- estimateSizeFactors(obj)

### Default cicero  #####
ciceroObj <- make_cicero_cds(obj, k = 50, reduced_coordinates = dimred)
saveRDS(ciceroObj,file='results/default_ciceroObj.RDS')
######

#Cicero Object CDS
ciceroOut <- custom_cicero_cds(obj, k = 50, max_knn_iterations = 5000, reduced_coordinates = dimred[colnames(obj),])
ciceroObj <- ciceroOut[[1]]
saveRDS(ciceroObj, "results/scATAC/cicero_aggregated_accessibility_cds.RDS")

#Keep Cell Mappings!
cellMapKNN <- ciceroOut[[2]]
saveRDS(cellMapKNN, "results/scATAC/cicero_KNN_Groupings_cds.RDS")

#Compute Correlations
message("Computing grouped correlations...")
gr <- featureToGR(row.names(ciceroObj))
o <- suppressWarnings(as.matrix( findOverlaps(resize( resize(gr,1,"center"), 2*flank + 1, "center"), resize(gr,1,"center"), ignore.strand=TRUE) ))
o <- data.table::as.data.table(data.frame(i = matrixStats::rowMins(o), j = matrixStats::rowMaxs(o)))
o <- data.frame(o[!duplicated(o),])
o <- o[o[,1]!=o[,2],]

#Note
#Log2 Transform Prior to Computing Co-Accessibility
#This isnt a big deal since Cicero Groupings K=50
#And thus dynamic range is 0-50 for pearson correlation
#However we decided to proceed with a Log2 Transformed
#Aggregtate Matrix
logMat <- log2(edgeR::cpm(exprs(ciceroObj))+1)
o$cor <- rowCorCpp(o[,1], o[,2], logMat, logMat)
connections <- data.frame(
    Peak1 = row.names(ciceroObj)[o[,1]], 
    Peak2 = row.names(ciceroObj)[o[,2]], 
    coaccess = o[,3]
    )

#Annotate CDS
message("Annotating Cell Data Set...")

#Make TSS Window for genes
genes <- getGeneGTF(gtfFile) %>% resize(1,"start") %>% resize(tssWindow * 2 + 1, "center")
names(genes) <- genes$gene_name
geneDF <- data.frame(chromosome=seqnames(genes),start=start(genes),end=end(genes), gene=genes$gene_name)
obj <- annotate_cds_by_site(obj, geneDF)

#Prepare for Co-Accessibility
nSites <- Matrix::colSums(Biobase::exprs(obj))
names(nSites) <- row.names(Biobase::pData(obj))

#Cicero with Correlations
message("Calculating normalized gene activities...")
ciceroGA <- normalize_gene_activities(build_gene_activity_matrix(obj, connections, coaccess_cutoff = corCutOff), nSites)

seCicero <- SummarizedExperiment(
    assays = SimpleList(gA = ciceroGA),
    rowRanges = genes[rownames(ciceroGA),],
    colData = mdata
)

seCiceroLog <- SummarizedExperiment(
    assays = SimpleList(logGA = log2(10^6 * ciceroGA + 1)),
    rowRanges = genes[rownames(ciceroGA),],
    colData = mdata
)

#Save Output
saveRDS(connections, "results/scATAC/Peaks-Co-Accessibility.RDS")
saveRDS(seCicero, "results/scATAC/Cicero-Gene-Activity.RDS")
saveRDS(seCiceroLog, "results/scATAC/Cicero-Log2-Gene-Activity.RDS")

cortex.atac[['cicero_GA']] <- CreateAssayObject(counts = assay(seCicero))
saveRDS(cortex.atac,file = 'data/merged_scATAC_integrated_cicero.RDS')

