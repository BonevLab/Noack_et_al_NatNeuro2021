#Clustering and scATAC-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(Seurat)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(Rcpp)
set.seed(1)


source('scripts/config.R')
source('scripts/aux_functions.R')
####################################################
#Input Data
####################################################

cortex.rna <- readRDS("data/merged_scRNA_unfiltered_IDs.RDS")
cortex.atac <- readRDS('data/merged_scATAC_integrated_cicero.RDS')

#Parameters
nCCA <- 20
nVarGenes <- 5000
selectMethod <- "all"
atac_assay <- 'GeneBody'

DefaultAssay(cortex.atac) <- atac_assay
DefaultAssay(cortex.rna) <- 'RNA'
#####

cortex.rna <- FindVariableFeatures(cortex.rna,nfeatures = nVarGenes)
varGenes <- SelectIntegrationFeatures(c(cortex.rna,cortex.atac),nfeatures =  nVarGenes,fvf.nfeatures = nVarGenes)

cortex.rna$tech <- 'rna'
cortex.atac$tech <- 'atac'

CCA <- FindTransferAnchors(
  reference = cortex.rna,
  query = cortex.atac,features = varGenes,dims = 1:20,l2.norm=T,k.anchor = 20,k.filter=200,k.score = 30,max.features=500,
  reduction = 'cca',normalization.method = 'LogNormalize',query.assay = atac_assay,reference.assay = 'RNA'
)

predicted.labels <- TransferData(
  anchorset = CCA,
  refdata = cortex.rna@active.ident,
  weight.reduction = cortex.atac[['harmony']]
)

cortex.atac <- AddMetaData(object = cortex.atac, metadata = predicted.labels)
pdf(paste0('plots/scATAC/scATAC_integration_QC.pdf'),height=6,width=18)
hist(cortex.atac$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()
table(cortex.atac$prediction.score.max > 0.5)
write.table(cortex.atac$prediction.score.max,file = 'results/integration/prediction_score_max.tsv',quote=F,col.names=T,row.names=T,sep='\t')

refdata <- GetAssayData(cortex.rna, assay = "RNA", slot = "data")[varGenes, ]
imputation <- TransferData(anchorset = CCA, refdata = refdata, weight.reduction = cortex.atac[["harmony"]])
cortex.atac[["RNA"]] <- imputation
coembed <- merge(x = cortex.rna, y = cortex.atac)
# datasets
coembed <- ScaleData(coembed, features = varGenes, do.scale = FALSE)
coembed <- RunPCA(coembed, features = varGenes, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:20,n.epochs = 2000,spread=1.5,min.dist=0.5)
coembed$celltype <- coembed$predicted.id

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
pdf(paste0('plots/scATAC/scATAC_integration_coembedding.pdf'),height=6,width=18)
p1|p2 + theme_cowplot()
dev.off()

#Save objects
saveRDS(CCA, 'data/CCA_scRNA_scATAC.RDS')
saveRDS(coembed, 'data/coembed_scRNA_scATAC.RDS')

#Get CCA Matrix
#alignedCCA <- alignedCCA <- CCA@object.list[[1]]@reductions$cca     #or cca?
alignedCCA <- alignedCCA <- CCA@object.list[[1]]@reductions$cca@cell.embeddings

#KNN Search
#Alternatively for speed FNN::getknnx(query, reference, k = 1)
#We just used a simple function
matchedCells <- findNN(
  query = alignedCCA[grep('query',row.names(alignedCCA)),],
  reference = alignedCCA[grep('reference',row.names(alignedCCA)),], 
  method = "euclidean")

matchedCells$corCCA <- rowCorCpp(
  match(matchedCells$x, row.names(alignedCCA)), 
  match(matchedCells$y, row.names(alignedCCA)), 
  alignedCCA, alignedCCA)


matchedCells$x <- gsub('_query','',matchedCells$x)
matchedCells$y <- gsub('_reference','',matchedCells$y)

matchedCells$corVarRNA <- rowCorCpp(
  match(matchedCells$x, colnames(coembed@assays$RNA@data)), 
  match(matchedCells$y, colnames(coembed@assays$RNA@data)), 
  t(as.matrix(coembed@assays$RNA@data[coembed@assays$RNA@var.features,])), 
  t(as.matrix(coembed@assays$RNA@data[coembed@assays$RNA@var.features,])))


#-------------------------------------------------------
#UMAP
#-------------------------------------------------------


#Plot DF
plotDF <- as.data.frame(Embeddings(coembed,reduction = 'umap'))
rownames(plotDF) <- rownames(alignedCCA)
plotDF[grep('query',row.names(plotDF)),"protocol"] <- 'atac'
plotDF[grep('reference',row.names(plotDF)),"protocol"] <- 'rna'
#plotDF <- plotDF[sample(seq_len(nrow(plotDF)), nrow(plotDF)),, drop = FALSE]
library(ggplot2)
ggplot(plotDF, aes(UMAP_1,UMAP_2,color=protocol)) + geom_point() +
  theme_bw()
saveRDS(list(plotDF = plotDF, matchedCells = matchedCells), "results/CCA_KNN_UMAP.RDS")

#----------------------------
# ATAC PseudoBulk
#----------------------------

DefaultAssay(cortex.atac) <- 'MACS2peaks'

peaks <- featureToGR(row.names(cortex.atac),'-')
rna_IDs <- factor(x = gsub('_M','',cortex.atac$predicted.id),levels=c('NSC','IPC','PN1','PN2','PN3','CR','IN','MG','Mural'))

se <- SummarizedExperiment(
  assays = SimpleList(counts = cortex.atac@assays$MACS2peaks@counts), 
  rowRanges = peaks,
  colData=data.frame(Clusters= rna_IDs,Group=cortex.atac$orig.ident)
)

objPB <- createPseudoBulk(
  mat = assay(se), 
  groups=colData(se)$Clusters, 
  labels=colData(se)$Group, 
  ceiling = 1, 
  minCells = 100, 
  maxCells = 500, 
  minReps = 2, 
  prior.count = 3, 
  nSim = 250
)

cdPB <- DataFrame(
  row.names=colnames(objPB[[1]]), 
  Group = stringr::str_split(colnames(objPB[[1]]), "\\._\\.", simplify = TRUE)[,1],
  #  Group = levels(Idents(cortex.atac))[as.numeric(stringr::str_split(colnames(objPB[[1]]), "._.", simplify = TRUE)[,1])+1],
  Type  = stringr::str_split(stringr::str_split(colnames(objPB[[1]]), "\\._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,1],
  Cells = as.integer(stringr::str_split(stringr::str_split(colnames(objPB[[1]]), "\\._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,2]),
  Replace = stringr::str_split(stringr::str_split(colnames(objPB[[1]]), "\\._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,3]
)

sePB <- SummarizedExperiment(assays = SimpleList(counts = objPB[[1]]), rowRanges = rowRanges(se), colData = cdPB)
metadata(sePB)$matCS <- objPB[[2]] 

saveRDS(sePB, "results/integration/Cluster_PseudoBulk-Summarized-Experiment.RDS")

