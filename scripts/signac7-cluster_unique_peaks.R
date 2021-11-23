library(Signac)
library(Seurat)
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(data.table)
library(edgeR)
set.seed(1)

source('scripts/config.R')
source('scripts/aux_functions.R')

#----------------------------
# Get Inputs
#----------------------------
cortex.atac <- readRDS('data/merged_scATAC_integrated_cicero.RDS')

which_cluster <- c("NSC","IPC","PN1","PN2")
cortex.atac <- subset(cortex.atac, idents=levels(cortex.atac)[(levels(cortex.atac)%in%which_cluster)])

DefaultAssay(cortex.atac) <- 'MACS2peaks'

peaks <- featureToGR(row.names(cortex.atac),'-')

se <- SummarizedExperiment(
  assays = SimpleList(counts = cortex.atac@assays$MACS2peaks@counts), 
  rowRanges = peaks,
  colData=data.frame(Clusters= Idents(cortex.atac),Group=cortex.atac$batch)
)
#----------------------------
# PseudoBulk
#----------------------------

objPB <- createPseudoBulk(
  mat = assay(se), 
  groups=Idents(cortex.atac), 
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

saveRDS(sePB, "results/scATAC/Cluster_PseudoBulk-Summarized-Experiment.RDS")

#----------------------------
# Unique Features
#----------------------------
params <- list(
  padj = 0.01,
  minSdRatio = 0.001,
  minLFC = 0.25,
  zCutoff = 1,
  breakPt = "last",
  groupMin = 25,
  maxGroupSize = max(floor(length(unique(colData(se)$Group))/3), 2),
  date = Sys.Date()
)

#Unique Features
uf <- uniqueFeatures(
  edgeR::cpm(assay(sePB),log=TRUE,prior.count=3)[,1:8],
  groups = colData(sePB)$Group[1:8],
  padj = params$padj,
  minSdRatio = params$minSdRatio,
  minLFC = params$minLFC,
  zCutoff = params$zCutoff,
  breakPt = params$breakPt,
  groupMin = params$groupMin,
  maxGroupSize = params$maxGroupSize,
  clusterCols = F
)

library(pheatmap)
pdf("plots/scATAC/Unique_Peaks.pdf",width=8,height=12)
m <- uf$groupMat
#m <- m[row.names(m)%in%row.names(logMat[rowMaxs(logMat)>2,]),]
#dimnames(m)[[2]] <- levels(Idents(cortex.atac))[as.numeric(stringr::str_split(dimnames(m)[[2]], "\\._\\.", simplify = TRUE)[,1])+1]
z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
z[z > 2] <- 2
z[z < -2] <- -2
pheatmap(z, cluster_cols=F, cluster_rows=F, show_rownames = FALSE,clustering_method = 'ward.D2')
dev.off()

saveRDS(uf, "results/scATAC/Unique_Peaks.RDS")
saveRDS(params, "results/scATAC/Unique_Peaks_Params.RDS")


library(misha)
require(tidyr)
gsetroot('/home/hpc/bonev/trackdb/mm10/')

tss <- gintervals.load("glp.intervals.ucscCanTSS")
peaks <- row.names(z)
peaks <- separate(as.data.frame(peaks),col = 'peaks', into = c("chrom", "start","end"),sep = "-")
peaks <- gintervals(as.character(peaks[,1]),as.numeric(peaks[,2]),as.numeric(peaks[,3]))
anno_peaks <- gintervals.neighbors(peaks,tss)
anno_peaks <- anno_peaks[abs(anno_peaks$dist)<5e5,]
anno <- unite(data = anno_peaks[,1:3],col = 'peaks',sep='-')
anno$gene <- as.character(anno_peaks[abs(anno_peaks$dist)<5e5,c('geneName')])
anno$dist <- anno_peaks[abs(anno_peaks$dist)<5e5,c('dist')]

############ Plot using complex heatmap ###############

require(ComplexHeatmap)
library(circlize)
require(dplyr)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


motifs_to_mark <- c("Hes5","Pax6","Hes1","Eomes","Btg2","Neurod6","Neurod2","Tle4","Sox5","Bcl11b","Fezf2","Satb2","Cux2","Gfap","Olig1","Gad2","Gad1")
motifs_to_mark <- c("Hes1","Eomes","Bcl11b","Olig1","Gad1")
#motifs_to_mark <- motifs_to_mark[c(3,4,6)]
reduced_anno <- anno[anno$gene%in%motifs_to_mark, ]
dim(reduced_anno)
#reduced_anno <- (reduced_anno %>% group_by(gene) %>% slice(1))
require(scales)
col.list <- hue_pal()(ncol(z))
#col.list <- rep(col.list,each=2)
#col.list <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(Idents(cortex.atac))))
names(col.list) <- colnames(z)
test <- as.data.frame(z)
test$n <- seq(1:nrow(z))
test <- merge(test,reduced_anno,by.x='row.names',by.y='peaks')
test$labels <- paste0(test$gene,':',test$dist)
#idxs <- c('Fezf2:−99404','Sox5:24831','Bcl11b:−74950','Sox5:−75708','Gad1:0','Gad2:−1152','Olig1:2745','Hes1:−5635','Hes1:−4379','Eomes:4005','Btg2:1308','Neurod2:1062')
#idxs <- c(96348,106484,113996)
#test <- test[test$n%in%idxs,]
ra = rowAnnotation(foo = anno_mark(at = test$n, labels = test$labels,labels_gp = gpar(fontsize = 16)))
col_names <- stringr::str_split(colnames(z), "\\._\\.", simplify = TRUE)[,1]
dup_indices <- !(duplicated(col_names))
mark_pos <- round(seq(1:length(col_names))[dup_indices])
ha = HeatmapAnnotation(foo = anno_mark(at = mark_pos, labels = col_names[dup_indices],labels_gp = gpar(fontsize = 16),labels_rot = 45),cluster=colnames(z),col = list(cluster=col.list),show_legend = F)


HM <- Heatmap(z,cluster_columns=T,cluster_rows=F,show_row_names = F,show_column_names = F,col = col_fun,right_annotation = ra,top_annotation = ha)

pdf('plots/scATAC/unique_peaks.pdf',width=12,height=24)
HM
dev.off()

############ Plot using complex heatmap ########################### Plot using complex heatmap ########################### Plot using complex heatmap ###############

###### Select AP specific peaks ####################
new_z <- as.matrix(z[93700:108100,])
peaks <- row.names(new_z)
peaks <- separate(as.data.frame(peaks),col = 'peaks', into = c("chrom", "start","end"),sep = "-")
peaks <- gintervals(as.character(peaks[,1]),as.numeric(peaks[,2]),as.numeric(peaks[,3]))
anno_peaks <- gintervals.neighbors(peaks,tss)
new_anno <- unite(data = anno_peaks[,1:3],col = 'peaks',sep='-')
new_anno$gene <- as.character(anno_peaks$geneName)
new_anno$dist <- as.numeric(anno_peaks$dist)

late_AP <- new_anno[new_anno$peaks%in%row.names(new_z[9500:nrow(new_z),]),]
early_AP <- new_anno[new_anno$peaks%in%row.names(new_z[3300:9500,]),]
write.table(early_AP,'results/uniquePeaks_earlyAP.tsv',sep='\t',col.names=T,row.names=F,quote=F)
write.table(late_AP,'results/uniquePeaks_lateAP.tsv',sep='\t',col.names=T,row.names=F,quote=F)

