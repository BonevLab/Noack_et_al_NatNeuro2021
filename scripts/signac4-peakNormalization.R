library(Signac)
library(Seurat)
library(reticulate)
library(ggplot2)
library(LSD)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table)
library(cowplot)
#library(future)
#plan("multiprocess", workers = 4)

fragment.path <- "data/fragments_filtered.tsv.gz"
source('scripts/config.R')
source('scripts/plot_functions.R')
source('scripts/aux_functions.R')
doHarmony=TRUE
#-----------------------------------------
# NORMALIZATION AND DIMENSIONAL REDUCTION
#-----------------------------------------
if (!file.exists('data/merged_scATAC_MACSpeaks_200k.RDS')){
  cortex.atac <- readRDS('data/merged_scATAC_filtered.RDS')
  unionPeaks <- readRDS('data/unionPeaks.RDS')
  #### Reduce size of the object by setting Bins assays to NULL
  
  peaks <- FeatureMatrix(
    fragments = fragment.path,
    features = unionPeaks,
    cells = colnames(cortex.atac),
    chunk = 50
  )
  
  cortex.atac[['MACS2peaks']] <- CreateAssayObject(counts = peaks)
  
  saveRDS(cortex.atac,'data/merged_scATAC_MACSpeaks_200k.RDS')
} else {
  cortex.atac <- readRDS('data/merged_scATAC_MACSpeaks_200k.RDS')
}
DefaultAssay(cortex.atac) <- 'MACS2peaks'

cortex.atac <- BinarizeCounts(cortex.atac)

#### Perform initial clustering to estimate total counts
assay_type <- 'MACS2peaks'
TFIDF_method <- 3
cutoff <- 'q0'
max_dims <- 20
resolution <- 1

mat <- cortex.atac@assays$MACS2peaks@counts
clusterSums <- groupSums(mat = mat, groups = paste0("C",cortex.atac$TopBins_snn_res.1), sparse = TRUE)
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3)


#### Select the top 20000 MACS2 peaks with highest coverage 
nFeatures <- 25000
varPeaks <- mat[head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nFeatures),]
cortex.atac[['TopVarPeaks']] <- CreateAssayObject(counts = varPeaks)

cortex.atac@assays$Bins5000 <- NULL
cortex.atac@assays$peaks <- NULL

#cortex.atac@assays$MACS2peaks@var.features <- row.names(varPeaks)
assay_types <- c('TopVarPeaks')
TFIDF_methods <- c(1,3)
cutoffs <- c('q0')
max_dims <- c(20,25)
resolutions <- c(0.6,0.8,1,1.2)

assay_type <- 'TopVarPeaks'
DefaultAssay(cortex.atac) <- assay_type
TFIDF_method <- 3
cutoff <- 'q0'
max_dims <- 20
resolution <- 0.9

#cortex.atac <- BinarizeCounts(cortex.atac)

for (assay_type in assay_types){
  
  DefaultAssay(cortex.atac) <- assay_type
  
  for (TFIDF_method in TFIDF_methods){
    
    cortex.atac <- RunTFIDF(cortex.atac,method = TFIDF_method,scale.factor = 10000,verbose = T)
    
    for (cutoff in cutoffs){
      
      cortex.atac <- FindTopFeatures(cortex.atac,min.cutoff = cutoff)
      cortex.atac <- RunSVD(
        object = cortex.atac,
        assay = assay_type,
        reduction.key = 'LSI_',
        reduction.name = 'lsi'
      )
      ElbowPlot(cortex.atac,ndims = 50,reduction = 'lsi')
      cortex.atac <- RunHarmony(
        object = cortex.atac,
        group.by.vars = 'batch',
        reduction = 'lsi',
        assay.use = assay_type,
        project.dim = FALSE
      )
      for (max_dim in max_dims){
        cortex.atac <- RunUMAP(object = cortex.atac, reduction = ifelse(doHarmony,'harmony','lsi'), dims = 1:max_dim,min.dist = 0.5,spread = 1.5,n.epochs = 2000,n.components = 2)
        cortex.atac <- FindNeighbors(object = cortex.atac, reduction = ifelse(doHarmony,'harmony','lsi'), dims = 1:max_dim,force.recalc = T)
        
        for (resolution in resolutions){
          cortex.atac <- FindClusters(object = cortex.atac, verbose = FALSE,resolution = resolution,n.start = 100,n.iter = 200)
          p1 <- DimPlot(object = cortex.atac, label = TRUE,group.by = 'orig.ident',cols=c('green','blue')) + ggtitle(paste0('N Var Features: ',length(cortex.atac@assays[[assay_type]]@var.features)))
          p2 <- DimPlot(object = cortex.atac, label = TRUE,cols=cluster_colors(length(unique(Idents(cortex.atac))))) + NoLegend()
          
          name <- paste0('VarPeaks_TFDF_',TFIDF_method,'_assay_',assay_type,'_cutoff_',cutoff,'_dims_',max_dim,'_res',resolution,'_features',length(cortex.atac@assays[[assay_type]]@var.features))
          pdf(paste0('plots/scATAC/',name,'.pdf'),height=8,width=16)
          print(CombinePlots(plots = list(p1, p2))) 
          dev.off()
          DefaultAssay(cortex.atac) <- 'GeneBody'
          for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
            pdf(paste0('plots/scATAC/',name,'_',set,'_GB.pdf'),height=10,width=10)
            print(FeaturePlot(cortex.atac, features = get(set),cols = gene_colors(3),pt.size = 0.1,min.cutoff = 'q5',label = T))
            dev.off()
          }
          DefaultAssay(cortex.atac) <- assay_type
        }
      }
    }
  }  
}


DefaultAssay(cortex.atac) <- assay_type
TFIDF_method <- 1
cutoff <- 'q0'
max_dims <- 20
resolution <- 0.3
saveRDS(cortex.atac,'data/merged_scATAC_MACSpeaks_200k.RDS')

########## Integrating with scRNA data to infer cluster names ##########
cortex.atac_rna <- readRDS("data/merged_scRNA_unfiltered_IDs.RDS")
cortex.atac <- readRDS('data/merged_scATAC_MACSpeaks_200k.RDS')
p1 <- DimPlot(cortex.atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(cortex.atac_rna,, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
pdf(paste0('plots/beforeIntegration.pdf'),height=8,width=16)
CombinePlots(plots = list(p1, p2)) + theme_cowplot()
dev.off()

DefaultAssay(cortex.atac) <- 'GeneBody'
DefaultAssay(cortex.atac_rna) <- 'RNA'
#####
varGenes <- FindVariableFeatures(
  object = cortex.atac_rna,
  nfeatures = 5000
)
varGenes <- varGenes@assays$RNA@var.features

transfer.anchors <- FindTransferAnchors(
  reference = cortex.atac_rna,
  query = cortex.atac,features = varGenes,dims = 1:20,
  reduction = 'cca',normalization.method = 'LogNormalize',query.assay = 'GeneBody',reference.assay = 'RNA'
)

saveRDS(transfer.anchors, 'data/transferAnchors_scRNA_scATAC.RDS')

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = cortex.atac_rna@active.ident,
  weight.reduction = cortex.atac[['lsi']]
)

cortex.atac <- AddMetaData(object = cortex.atac, metadata = predicted.labels)
pdf(paste0('plots/scATAC/scATAC_integration_QC.pdf'),height=6,width=18)
hist(cortex.atac$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()
table(cortex.atac$prediction.score.max > 0.5)

plot1 <- DimPlot(cortex.atac_rna, label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(cortex.atac, group.by='seurat_clusters',label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq orig')
plot3 <- DimPlot(cortex.atac, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq predicted')
pdf(paste0('plots/scATAC/scATAC_integrated2.pdf'),height=6,width=18)
CombinePlots(list(plot1,plot2,plot3), ncol = 3) + theme_cowplot()
dev.off()

DE_Method='t'
top_DE <- 40
resolution <- 0.2
DefaultAssay(cortex.atac) <- 'Prom'
cortex.atac <- ScaleData(cortex.atac)
cortex.atac.markers <- FindAllMarkers(cortex.atac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)

top40 <- cortex.atac.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
top20 <- cortex.atac.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(top40,file='results/scATAC/top40_promDEmarkers.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(cortex.atac) <- 'Prom'
pdf('plots/scATAC/promDE_top20.pdf',height=40,width=20)
DoHeatmap(cortex.atac, features = top20$gene,raster = T) + NoLegend()
dev.off()
pdf('plots/scATAC/promDE_top40.pdf',height=40,width=20)
DoHeatmap(cortex.atac, features = top40$gene,raster = T) + NoLegend()
dev.off()

DefaultAssay(cortex.atac) <- 'GeneBody'
cortex.atac <- ScaleData(cortex.atac)
cortex.atac.markers <- FindAllMarkers(cortex.atac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = DE_Method)

top40 <- cortex.atac.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
top20 <- cortex.atac.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(top40,file='results/scATAC/top40_genebodyDEmarkers.tsv',quote=F,col.names=T,row.names=F,sep='\t')
DefaultAssay(cortex.atac) <- 'GeneBody'
pdf('plots/scATAC/genebodyDE_top20.pdf',height=40,width=20)
DoHeatmap(cortex.atac, features = top20$gene,raster = T) + NoLegend()
dev.off()
pdf('plots/scATAC/genebodyDE_top40.pdf',height=40,width=20)
DoHeatmap(cortex.atac, features = top40$gene,raster = T) + NoLegend()
dev.off()

cluster_annot <- read.table('results/ATAC_clusters_unfiltered.tsv',header=T,row.names = "ClusterID")
new.cluster.ids <- as.vector(cluster_annot$CellType)
names(new.cluster.ids) <- as.factor(levels(cortex.atac))
cortex.atac <- RenameIdents(cortex.atac, new.cluster.ids)
levels(cortex.atac) <- c('NSC','IPC','PN1','PN2','CR','IN','MG+Mural')

pdf('plots/scATAC/UMAPwithIDs.pdf',height=8,width=8)
DimPlot(cortex.atac, reduction = "umap",label = T,repel=T,pt.size = 0.2) + NoLegend() + theme_cowplot()
dev.off()


cortex.atac$celltype <- Idents(cortex.atac)
UMAP_plot <- data.frame(cortex.atac[['umap']]@cell.embeddings, cortex.atac@meta.data)
colnames(UMAP_plot) <- c("x","y",colnames(UMAP_plot)[3:ncol(UMAP_plot)])
clustCol <- colnames(UMAP_plot)[grep("celltype",colnames(UMAP_plot))]


pdf("plots/scATAC/scATAC_integrated_ggplot.pdf")
p1 <- ggplot(UMAP_plot, aes(x=x,y=y,color=celltype)) + geom_point(size = 0.1) + 
  DarkTheme() + xlab("UMAP1") + ylab("UMAP2") + theme_cowplot()
print(p1)
dev.off()

saveRDS(cortex.atac,'data/merged_scATAC_integrated.RDS')


##### Export bigwigs files per cluster ######

dirClusters <- "results/scATAC_clusters/"
dir.create(dirClusters)
dirPeaks <- dirClusters
method <- "q"
cutoff <- 0.05
shift <- -75
extsize <- 150
genome_size <- 'mm'
all_cells <- Idents(cortex.atac)
all_reads <- readRDS('data/all_fragments.RDS')
insertions <- all_reads
insertions2 <- all_reads
insertions2$start <- insertions2$end-1
insertions$end <- insertions$start+1 
#  insertions <- rbind(insertions,insertions2)
insertions <- insertions[with(data = insertions, expr = order(chr, start)), ]
insertions2 <- insertions2[with(data = insertions2, expr = order(chr, start)), ]
rm(all_reads)

cat('source activate macs2;\n',file='scripts/macs.sh',append=F)
for (i in levels(all_cells)){
  cluster_cells <- names(all_cells[all_cells==i])
  i <- gsub('\\/','-',i)
  i <- gsub('\\?','',i)
  ClusterFragments(reads=insertions,
                   reads2=insertions2,
                   cells=cluster_cells,
                   output.path=paste0(dirClusters,gsub(" ", "_", i, fixed = TRUE),'.bed'),
                   assume.sorted = FALSE, verbose = TRUE
  )
  clusterBedj <- paste0(gsub(" ", "_", i, fixed = TRUE),".bed")
  cmdPeaks <- sprintf(
    "macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all -B --SPMR", 
    genome_size, 
    paste0(gsub(" ", "_", i, fixed = TRUE)), 
    paste0(dirClusters,clusterBedj), 
    dirPeaks
  )
  if (!is.null(shift) & !is.null(extsize)) {
    cmdPeaks <- sprintf("%s --shift %s --extsize %s", cmdPeaks, shift, extsize)
  }
  if (tolower(method) == "p") {
    cmdPeaks <- sprintf("%s -p %s", cmdPeaks, cutoff)
  }else {
    cmdPeaks <- sprintf("%s -q %s", cmdPeaks, cutoff)
  }
  #  message("Running Macs2...")
  # message(cmdPeaks)
  cmdBW <- sprintf(
    "\nsort -k1,1 -k2,2n %s > %s; bedGraphToBigWig %s %s %s", 
    paste0(dirPeaks,'/',gsub(" ", "_", i, fixed = TRUE),'_treat_pileup.bdg'), 
    paste0(dirPeaks,'/',gsub(" ", "_", i, fixed = TRUE),'.bdg'),
    paste0(dirPeaks,'/',gsub(" ", "_", i, fixed = TRUE),'.bdg'),
    '/home/hpc/bonev/annotations/mm10/mm10.chrom.sizes',
    paste0(dirPeaks,'/',gsub(" ", "_", i, fixed = TRUE),'.bw \n')
  ) 
  cat(cmdPeaks,file='scripts/macs.sh',append=T)
  cat(cmdBW,file='scripts/macs.sh',append=T)
}

### Calculate Variance #####
intraVar <- FindVariableFeatures(object=cortex.atac@assays$MACS2peaks@data)
saveRDS(intraVar,'results/scATAC/MACS2peaks_intraVar.RDS')

### import tracks in the database ######
source('/home/hpc/bonev/projects/hic/sc/config.R')
source(paste0(main_f,'scripts/main_functions.R'))
files_f <- list.files(dirClusters,pattern = '\\.bw')
for (file_f in files_f){
  import_track(paste0(dirClusters,track),'scATAC.E14',gsub('\\.bw','',file_f),gsub('\\.bw','',file_f),10)
}
##################################

