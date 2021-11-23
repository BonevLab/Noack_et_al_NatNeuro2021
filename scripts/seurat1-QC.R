library(dplyr)
library(Seurat)
library(Hmisc)
library(SummarizedExperiment)
require(edgeR)
require(LSD)

source('scripts/plot_functions.R')
source('scripts/config.R')

batch_correction=FALSE
add_control=FALSE
min.cells <- 3
min.features <- 200
minGenes <- 1000
minRNA <- 2500
maxGenes <- 7000
#maxRNA <- 50000
maxMT <- 10

input_dirs <- paste0(rna_sample_dir,sample_names)
control_dirs <- c('scRNA/Telley2019_raw_counts.tsv')
control_names <- c('Telley2019')
control_reps <- c('rep1')
sample_list <- list()
control_list <- list()


#### Read in the data and create Seurat object #####
for (i in 1:length(input_dirs)){
  seurat.data <- Read10X(data.dir = input_dirs[i])
  seurat.data <- CreateSeuratObject(counts = seurat.data, project = sample_names[i], min.cells = min.cells, min.features = min.features)
  seurat.data[["percent.mt"]] <- PercentageFeatureSet(seurat.data, pattern = "^mt-")
  pdf(paste0('plots/scRNA/',sample_names[i],'_QC.pdf'),height=8,width=18)
  p1 <- VlnPlot(seurat.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01)
  p2 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(CombinePlots(plots = list(p1, p2,p3),ncol = 3))
  dev.off()
  sample_list[sample_names[i]] <- seurat.data
}

for (i in 1:length(control_dirs)){
  control_tsv <- read.table(control_dirs[i],sep='\t',header=T)
  seurat.data <- CreateSeuratObject(counts = control_tsv, project = control_names[i], min.cells = min.cells, min.features = min.features)
  seurat.data[["percent.mt"]] <- PercentageFeatureSet(seurat.data, pattern = "^mt-")
  pdf(paste0('plots/scRNA/',control_names[i],'_QC.pdf'),height=8,width=18)
  p1 <- VlnPlot(seurat.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01)
  p2 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(CombinePlots(plots = list(p1, p2,p3),ncol = 3))
  dev.off()
  control_list[control_names[i]] <- seurat.data
  }
################################################

#### Subset based on QC metrics #########
for (i in names(sample_list)){
  seurat.data <- sample_list[[i]]
  maxRNA <- quantile(seurat.data$nCount_RNA[seurat.data$nCount_RNA>minRNA],0.961)
  seurat.data <- subset(seurat.data, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxMT &nCount_RNA > minRNA&nCount_RNA < maxRNA&(nCount_RNA/nFeature_RNA)<10)
  pdf(paste0('plots/scRNA/',i,'_QCafterSubsetting.pdf'),height=8,width=18)
  p1 <- VlnPlot(seurat.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01)
  p2 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(CombinePlots(plots = list(p1, p2,p3),ncol = 3))
  dev.off()
  sample_list[[i]] <- seurat.data
}

for (i in names(control_list)){
  seurat.data <- control_list[[i]]
  seurat.data <- subset(seurat.data, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxMT &nCount_RNA > minRNA&nCount_RNA < maxRNA)
  pdf(paste0('plots/scRNA/',i,'_QCafterSubsetting.pdf'),height=8,width=18)
  p1 <- VlnPlot(seurat.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01)
  p2 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(CombinePlots(plots = list(p1, p2,p3),ncol = 3))
  dev.off()
  control_list[[i]] <- seurat.data
}

############################

##### Merging into one seurat object ##############

if (!batch_correction){
  #### Default merging and subsetting(no batch correction)
  if (!add_control){
    cortex <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)], add.cell.ids = names(sample_list), project = "RNA")
  } else {
    cortex <- merge(sample_list[[1]], y = c(sample_list[2:length(sample_list)],control_list), add.cell.ids = c(names(sample_list),names(control_list)), project = "RNA")
  }
  pdf('plots/scRNA/merged_QC.pdf',height=8,width=16)
  print(VlnPlot(cortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.001))
  print(FeatureScatter(cortex, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size = 0.001))
  print(FeatureScatter(cortex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.001))
  #print(CombinePlots(plots = list(p1, p2,p3),ncol = 3))
  dev.off()
  #############
} else {
  ###### Individual subsetting and batch correction #####
  if (!add_control){
    cortex.list <- sample_list
  } else {
    cortex.list <- c(sample_list,control_list)
  }
  for (i in 1:length(cortex.list)) {
    cortex.list[[i]] <- SCTransform(cortex.list[[i]], verbose = FALSE)
  }
  cortex.features <- SelectIntegrationFeatures(object.list = cortex.list, nfeatures = 2000)
  cortex.list <- PrepSCTIntegration(object.list = cortex.list, anchor.features = cortex.features, verbose = FALSE)
  cortex.anchors <- FindIntegrationAnchors(object.list = cortex.list, normalization.method = "SCT", anchor.features = cortex.features, verbose = FALSE)
  cortex.integrated <- IntegrateData(anchorset = cortex.anchors, normalization.method = "SCT", verbose = FALSE)
  cortex.integrated <- RunPCA(cortex.integrated, verbose = FALSE)
  cortex.integrated <- RunUMAP(cortex.integrated, dims = 1:50)
  plots <- DimPlot(cortex.integrated, group.by = c("orig.ident"), combine = FALSE)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
  CombinePlots(plots)
  cortex <- cortex.integrated
}

###### Correlation plots between replicates and samples #####
cortex_cpm <- as.data.frame(matrix(NA,nrow=nrow(cortex),ncol=length(unique(Idents(cortex)))))
colnames(cortex_cpm) <- unique(Idents(cortex))
for (i in unique(Idents(cortex))){
  cortex_sub <- subset(cortex,subset=orig.ident==i)
  cortex_cpm[,i] <- log2(edgeR::cpm(rowSums(cortex_sub@assays$RNA@counts), log = FALSE)+1)
}
png('plots/scRNA/correlation_QC.png',height=4200,width=4200,res = 300)
heatpairs(as.matrix(cortex_cpm[rowSums(cortex_cpm)>3,]),labels=colnames(cortex_cpm),colpal="bl2gr2rd")
dev.off()
############################################################


########### Cell Cycle Score identification but NO CORRECTION #####
s.genes <- capitalize(tolower(cc.genes$s.genes))
g2m.genes <- capitalize(tolower(cc.genes$g2m.genes))
cortex <- CellCycleScoring(cortex, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cortex$CC.Difference <- cortex$S.Score - cortex$G2M.Score

saveRDS(cortex,'data/merged_scRNA_filtered.RDS')







