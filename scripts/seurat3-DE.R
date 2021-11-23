library(dplyr)
library(Seurat)
library(ggplot2)
library(MAST)

source('scripts/plot_functions.R')
source('scripts/config.R')

DE_Method='MAST'
top_DE <- 40
resolution <- 0.2

##### Initialize objects #####
cortex <- readRDS('data/merged_scRNA_filtered_umap.RDS')


###### Store Initial cluster Plot #####################
pdf('plots/scRNA/UMAPforDE.pdf',height=8,width=8)
DimPlot(cortex, reduction = "umap",label = T,repel=T,pt.size = 0.2)
print(p)
dev.off()


DefaultAssay(cortex) <- 'RNA'
cortex.markers <- FindAllMarkers(cortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)
saveRDS(cortex.markers,'data/scRNA_DE_markers.RDS')

top40 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
top20 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(top40,file='results/scRNA_top40_DEmarkers.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(cortex) <- 'RNA'
pdf('plots/scRNA/DE_top20.pdf',height=40,width=20)
DoHeatmap(cortex, features = top20$gene,raster = T) + NoLegend()
dev.off()
pdf('plots/scRNA/DE_top40.pdf',height=40,width=20)
DoHeatmap(cortex, features = top40$gene,raster = T) + NoLegend()
dev.off()

pdf('plots/scRNA/CellCycle_VlnPLot.pdf',height=10,width=20)
print(VlnPlot(cortex, c("S.Score","G2M.Score"),pt.size = 0.001))
dev.off()        

cluster_annot <- read.table('results/RNA_clusters_unfiltered.tsv',header=T,row.names = "ClusterID")
new.cluster.ids <- as.vector(cluster_annot$CellType)
names(new.cluster.ids) <- as.factor(levels(cortex))
cortex <- RenameIdents(cortex, new.cluster.ids)

levels(cortex) <- c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3','CR','IN','MG','Mural')

pdf('plots/scRNA/UMAPwithIDs.pdf',height=8,width=8)
DimPlot(cortex, reduction = "umap",label = T,repel=T,pt.size = 0.2) + NoLegend()
dev.off()
saveRDS(cortex,'data/merged_scRNA_unfiltered_IDs.RDS') 

DefaultAssay(cortex) <- 'RNA'

cortex.markers <- FindAllMarkers(cortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)
saveRDS(cortex.markers,'data/scRNA_DE_markers_IDs.RDS')

top40 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
top20 <- cortex.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(top40,file='results/scRNA_top40_DEmarkers_IDs.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(cortex) <- 'RNA'
pdf('plots/scRNA/DE_top20.pdf',height=40,width=20)
DoHeatmap(cortex, features = top20$gene,raster = T) + NoLegend()
dev.off()
pdf('plots/scRNA/DE_top40.pdf',height=40,width=20)
DoHeatmap(cortex, features = top40$gene,raster = T) + NoLegend()
dev.off()
#idents_to_keep <- cluster_annot
#idents_to_keep <- idents_to_keep[c(1,2,3,4,5,6,7,8,9,10,11,13,15,17,18,20,22,23),]
#idents_to_keep$CellType <- droplevels(idents_to_keep$CellType)
#cortex_filtered <- subset(cortex,idents = idents_to_keep$CellType)
#saveRDS(cortex_filtered,'data/merged_scRNA_filtered_IDs.RDS')  

#### Normal transformation ####

#cortex <- cortex_filtered
#cortex <- RunUMAP(cortex, dims = 1:20,min.dist = 0.5,spread = 1.5,n.components = 2L,n.epochs = 2000)
#pdf('plots/scRNA/filtered_UMAP.pdf',height=8,width=16)
#DimPlot(cortex, reduction = "umap",label = T,group.by='orig.ident',cols=sample_colors(length(unique(cortex$orig.ident))))
#DimPlot(cortex, reduction = "umap",label = T,repel=T,cols=cluster_colors(length(unique(Idents(cortex)))))
#dev.off()
#for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
#  pdf(paste0('plots/scRNA/unfiltered_',set,'.pdf'),height=16,width=16)
#  print(FeaturePlot(cortex, features = get(set),cols = gene_colors(3),pt.size = 0.1,label = T))
#  dev.off()
#}
#pdf(paste0('plots/scRNA/unfiltered_','Features','.pdf'),height=16,width=16)
#print(FeaturePlot(cortex, features = c('nCount_RNA','nFeature_RNA','percent.mt'),cols = gene_colors(3),pt.size = 0.1,label = T))
#print(VlnPlot(cortex,'nFeature_RNA',pt.size = 0.001))
#print(VlnPlot(cortex,'nCount_RNA',pt.size = 0.001))
#print(VlnPlot(cortex,'percent.mt',pt.size = 0.001))
#dev.off()



saveRDS(cortex_filtered,'data/merged_scRNA_filtered_IDs.RDS')  


