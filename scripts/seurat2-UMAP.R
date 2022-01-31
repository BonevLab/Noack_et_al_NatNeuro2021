library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
set.seed(1234)

source('scripts/plot_functions.R')
source('scripts/config.R')

normalizations=c('VST','SCT')
nDims <- c(20,25,30,35)
resolutions <- c(0.8,1,1.2)
umap_methods <- 1

##### Initialize objects #####
cortex <- readRDS('data/merged_scRNA_filtered.RDS')


#### Some Bustools inspired QC plots
pdf('plots/scRNA/UMIs_vs_Genes.pdf',height=8,width=18)
print(ggplot(cortex@meta.data, aes(nCount_RNA, nFeature_RNA,color=orig.ident)) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(x = "Total UMI counts per cell", y = "Number of genes detected"))
dev.off()

### Removing replicate information #####
#cortex$orig.ident <- gsub('_rep1','',cortex$orig.ident)
#cortex$orig.ident <- gsub('_rep2','',cortex$orig.ident)

#### Normal transformation ####
DefaultAssay(cortex) <- 'RNA'
cortex <- NormalizeData(cortex)
cortex <- FindVariableFeatures(cortex, selection.method = "vst")
top50 <- head(VariableFeatures(cortex), 50)
pdf('plots/scRNA/VST_VarFeatures.pdf',height=8,width=18)
plot1 <- VariableFeaturePlot(cortex)
LabelPoints(plot = plot1, points = top50, repel = TRUE)
dev.off()
cortex <- ScaleData(cortex,features =rownames(cortex))
cortex <- RunPCA(cortex,npcs = 50)
pdf('plots/scRNA/VST_PCA.pdf',height=24,width=16)
DimHeatmap(cortex, dims = 1:20, cells = 500, balanced = TRUE)
DimHeatmap(cortex, dims = 20:40, cells = 500, balanced = TRUE)
ElbowPlot(cortex,ndims = 50)
dev.off()
write.table(cortex[["pca"]]@feature.loadings,'results/VST_PCA.txt',quote=F,row.names=T,col.names=T,sep='\t')
cortex <- FindNeighbors(cortex, dims = 1:20,force.recalc = T)
cortex <- FindClusters(cortex, resolution = 0.6,n.start = 20,n.iter = 30,algorithm = 1)
cortex <- RunUMAP(cortex, dims = 1:20,min.dist = 0.5,spread = 1.5,n.components = 2L)
pdf('plots/scRNA/VST_UMAP.pdf',height=8,width=16)
DimPlot(cortex, reduction = "umap",label = T,group.by='orig.ident',cols=sample_colors(length(unique(cortex$orig.ident))))
DimPlot(cortex, reduction = "umap",label = T,repel=T,cols=cluster_colors(length(unique(Idents(cortex)))))
dev.off()
for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
  pdf(paste0('plots/scRNA/VST_',set,'.pdf'),height=16,width=16)
  print(FeaturePlot(cortex, features = get(set),cols = gene_colors(3),pt.size = 0.1,label = T))
  dev.off()
}
pdf(paste0('plots/scRNA/VST_','Features','.pdf'),height=16,width=16)
print(FeaturePlot(cortex, features = c('nCount_RNA','nFeature_RNA','percent.mt'),cols = gene_colors(3),pt.size = 0.1,label = T))
print(VlnPlot(cortex,'nFeature_RNA',pt.size = 0.001))
print(VlnPlot(cortex,'nCount_RNA',pt.size = 0.001))
print(VlnPlot(cortex,'percent.mt',pt.size = 0.001))
dev.off()
cortexVST <- cortex

###### SCT Transform ######

cortex <- SCTransform(cortex, verbose = TRUE)
cortex <- RunPCA(cortex, npcs = 40, verbose = FALSE)
pdf('plots/scRNA/SCT_PCA.pdf',height=24,width=16)
DimHeatmap(cortex, dims = 1:20, cells = 500, balanced = TRUE)
DimHeatmap(cortex, dims = 20:40, cells = 500, balanced = TRUE)
ElbowPlot(cortex,ndims = 40)
dev.off()
write.table(cortex[["pca"]]@feature.loadings,'results/SCT_PCA.txt',quote=F,row.names=T,col.names=T,sep='\t')
cortex <- FindNeighbors(cortex, dims = 1:20,force.recalc = T)
cortex <- FindClusters(cortex, resolution = 0.6,n.start = 20,n.iter = 20,algorithm = 1)
cortex <- RunUMAP(cortex, dims = 1:20,min.dist = 0.5,spread = 1.5,n.components = 2L)
pdf('plots/scRNA/SCT_UMAP.pdf',height=8,width=16)
DimPlot(cortex, reduction = "umap",label = T,group.by='orig.ident',cols=sample_colors(length(unique(cortex$orig.ident))))
DimPlot(cortex, reduction = "umap",label = T,repel=T,cols=cluster_colors(length(unique(Idents(cortex)))))
dev.off()
for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
  pdf(paste0('plots/scRNA/SCT_',set,'.pdf'),height=16,width=16)
  print(FeaturePlot(cortex, features = get(set),cols = gene_colors(3),pt.size = 0.1,label = T))
  dev.off()
}
pdf(paste0('plots/scRNA/SCT_','Features','.pdf'),height=16,width=16)
print(FeaturePlot(cortex, features = c('nCount_RNA','nFeature_RNA','percent.mt'),cols = gene_colors(3),pt.size = 0.1,label = T))
dev.off()
cortexSCT <- cortex

######### DATA Exploration ####################

for (normalization in normalizations){
  if(normalization=='VST'){cortex <- cortexVST} else {cortex <- cortexSCT}
  DefaultAssay(cortex) <- ifelse(normalization=='VST','RNA','SCT')
  for (nDim in nDims){
    cortex <- RunUMAP(cortex, dims = 1:nDim,min.dist = 0.5,spread = 1.5,n.epochs=500)
    cortex <- FindNeighbors(cortex, dims = 1:nDim)
    for (resolution in resolutions){
      for (umap_method in umap_methods){
        cortex <- FindClusters(cortex, resolution = resolution,algorithm = umap_method)
        pdf(paste0('plots/scRNA/UMAP_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'.pdf'),height=8,width=16)
        print(DimPlot(cortex, reduction = "umap",label = T,group.by='orig.ident',cols=sample_colors(length(unique(cortex$orig.ident)))))
        print(DimPlot(cortex, reduction = "umap",label = T,repel=T,cols=cluster_colors(length(unique(Idents(cortex))))))
        dev.off()
        for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
          pdf(paste0('plots/scRNA/UMAP_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_',set,'.pdf'),height=16,width=16)
          print(FeaturePlot(cortex, features = get(set),cols = gene_colors(3),pt.size = 0.1,label = T))
          dev.off()
        }    
        pdf(paste0('plots/scRNA/UMAP_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_Features3.pdf'),height=16,width=16)
        print(VlnPlot(cortex,'nFeature_RNA',pt.size = 0.001))
        print(VlnPlot(cortex,'nCount_RNA',pt.size = 0.001))
        print(VlnPlot(cortex,'percent.mt',pt.size = 0.001))
        dev.off()
      }
    }
  }
}


########################################################
####### Choose values for subsequent analysis ##########
########################################################

normalization <- 'VST'
cluster_method <- 1
resolution <- 0.3
nDim <- 20
doHarmony=TRUE

if(normalization=='VST'){cortex <- cortexVST} else {cortex <- cortexSCT}
DefaultAssay(cortex) <- ifelse(normalization=='VST','RNA','SCT')
cortex <- RunHarmony(object = cortex,group.by.vars = 'orig.ident',reduction = 'pca',assay.use = DefaultAssay(cortex),project.dim = FALSE)
cortex <- FindNeighbors(cortex,reduction = ifelse(doHarmony,'harmony','pca'), dims = 1:nDim,force.recalc = T)
cortex <- FindClusters(cortex, resolution = resolution,n.start = 100,n.iter = 500,algorithm = cluster_method)
cortex <- RunUMAP(cortex, dims = 1:nDim,reduction = ifelse(doHarmony,'harmony','pca'),n.components = 2L,min.dist=0.5,spread = 1,n.epochs = 500)
pdf(paste0('plots/scRNA/chosen_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_UMAP.pdf'),height=16,width=16)
DimPlot(cortex, reduction = "umap",label = T,group.by='orig.ident',cols=c('darkgreen','darkred'))
DimPlot(cortex, reduction = "umap",label = T,repel=T,cols=cluster_colors(length(levels(cortex))))
dev.off()
for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
  pdf(paste0('plots/scRNA/chosen_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_',set,'Features.pdf'),height=16,width=16)
  print(FeaturePlot(cortex, features = get(set),cols = gene_colors(3),pt.size = 0.1,label = T))
  dev.off()
}    
pdf(paste0('plots/scRNA/chosen_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_Features2.pdf'),height=12,width=12)
print(DotPlot(cortex, features = c(diff_markers,layer_markers,temp_markers,other_markers,regional_markers),cols = gene_colors(3)) + RotatedAxis())
dev.off()
pdf(paste0('plots/scRNA/chosen_',normalization,'_dims',nDim,'_res',resolution,'_umap',umap_method,'_Features3.pdf'),height=16,width=16)
p1 <- VlnPlot(cortex,'nFeature_RNA',pt.size = 0.001)
p2 <- VlnPlot(cortex,'nCount_RNA',pt.size = 0.001)
p3 <- VlnPlot(cortex,'percent.mt',pt.size = 0.001)
print(CombinePlots(list(p1,p2,p3)))
dev.off()

saveRDS(cortex,file='data/merged_scRNA_filtered_umap.RDS')


