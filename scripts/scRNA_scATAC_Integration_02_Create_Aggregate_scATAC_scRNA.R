#Integration of scRNA and scATAC-seq data 
#With adaptations from Granja et al. (2019) 

#Multimodal profiling of the transcriptional regulatory landscape
#of developing mouse cortex identifies Neurog2 as a key epigenome remodeler
#Cite Noack et al. (2022) 

library(Seurat)
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
set.seed(1)

source('scripts/aux_functions.R')
####################################################
#Functions
####################################################


####################################################
#Input Data
####################################################

#Input Files
seATAC_file <- 'data/merged_scATAC_integrated_cicero.RDS'
ciceroObj <- "results/scATAC/cicero_aggregated_accessibility_cds.RDS"
ciceroObj <- readRDS(ciceroObj)
ciceroKNN_file <- "results/scATAC/cicero_KNN_Groupings_cds.RDS"
ciceroKNN_file <- readRDS(ciceroKNN_file)
ciceroATAC_file <- "results/scATAC/Cicero-Gene-Activity.RDS"
seRNA_file <- "data/merged_scRNA_unfiltered_IDs.RDS"
CCA_file <- "results/CCA_KNN_UMAP.RDS"
gtf_file <- '/home/hpc/bonev/annotations/mm10/cellranger_rna/genes/genes.gtf'

#Get Clusters information for each KNN Group top group/exp wins!
scATAC <- readRDS(seATAC_file)
KNN <- data.frame(ciceroKNN_file, stringsAsFactors=FALSE)
KNN <- apply(KNN,2,paste0)
KNNClusters <- apply(KNN, 2, function(x) as.data.frame(scATAC@active.ident)[x,"scATAC@active.ident"])
KNNGroups <- apply(KNN, 2, function(x) as.data.frame(scATAC$orig.ident)[x,"scATAC$orig.ident"])
KNN_Highest_Cluster <- lapply(seq_len(nrow(KNN)), function(x) names(sort(table(KNNClusters[x,]), decreasing=TRUE))[1]) %>% unlist
KNN_Highest_Experiment <- lapply(seq_len(nrow(KNN)), function(x) names(sort(table(KNNGroups[x,]), decreasing=TRUE))[1]) %>% unlist

####################################################
# scATAC-seq Merging from Cicero
####################################################

se <- SummarizedExperiment(
  assays = SimpleList(counts = assay(ciceroObj)),
  rowRanges = featureToGR(row.names(ciceroObj)),
  colData = DataFrame(
    row.names = colnames(assay(ciceroObj)), 
    clustATAC = KNN_Highest_Cluster, 
    groupATAC = KNN_Highest_Experiment
  ))
metadata(se)$knn <- KNN
metadata(se)$knnClust <- KNNClusters
metadata(se)$knnGroup <- KNNGroups
test <- grToFeature(rowRanges(seA),sep='_')
saveRDS(se, "results/scATAC-Merged-KNN-SVD.RDS")

####################################################
# scRNA-seq Merging from Cicero
####################################################
scRNA <- readRDS(seRNA_file)
CCA <- readRDS(CCA_file)[[2]]
CCA <- CCA[CCA$corCCA >= 0.5,] #Filter low correlation to CCA R > 0.5

#Get scRNA Matrix
scMat <- scRNA@assays$RNA@counts

#Create aggregated Matched RNA Matrix exclude non-mappings
matRNA <- matrix(NA, ncol = nrow(KNN), nrow = nrow(scMat))
for(i in seq_len(nrow(KNN))){
  if(i%%100==0) print(i)
  knnIdx <- paste0(t(KNN[i,]))
  rnaIdx <- CCA$y[match(knnIdx, CCA$x)]
  rnaIdx <- na.omit(rnaIdx)
  matRNA[, i] <- Matrix::rowSums(scMat[,rnaIdx])
}
colnames(matRNA) <- colnames(se)

#Create aggregated summarized experiment
seRNA <- SummarizedExperiment(
    assays = SimpleList(counts = matRNA)
  )
rownames(seRNA) <- rownames(scMat)
gtf <- getGeneGTF(gtf_file)

### Correct Genes for ucsc known gene coordinates
require(misha)
gsetroot('/home/hpc/bonev/trackdb/mm10')
genes <- gintervals.load("glp.intervals.ucscCanGenes")
genes$strand='*'
genes <- genes[match(gtf$gene_name,as.character(genes$geneName)),]
ranges(gtf)[!is.na(genes$start)] <- IRanges(start=genes$start[!is.na(genes$start)],end=genes$end[!is.na(genes$start)])  
#####
noMatch <- rownames(seRNA)[(is.na(match(rownames(seRNA),gtf$gene_name)))]
noMatch <- sub("\\..*", "", noMatch)
Indx <- match(rownames(seRNA),gtf$gene_name)
Indx[is.na(Indx)] <-  match(noMatch,gtf$gene_name)
gtfMatch <- gtf[Indx]
names(gtfMatch) <- rownames(seRNA)
rowRanges(seRNA) <- gtfMatch

#Use ATAC cluster info
colData(seRNA) <- DataFrame(row.names = colnames(seRNA), 
  clustATAC = colData(se)$clustATAC, 
  groupATAC = colData(se)$groupATAC
  )

saveRDS(seRNA, "results/scRNA-Merged-KNN-SVD.RDS")







