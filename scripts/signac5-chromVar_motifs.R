library(Signac)
library(Seurat)
library(JASPAR2020)
library(JASPAR2018)
library(Matrix)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(Hmisc)
library(RColorBrewer)
require(misha)
library(BiocParallel)
library(chromVAR)
register(MulticoreParam(12, progressbar = TRUE))
#mm10_txdb <- loadDb("/home/hpc/bonev/annotations/mm10/mm10_txdb.sqlite")
set.seed(1234)

source('scripts/config.R')
source('scripts/aux_functions.R')
fragment.path <- "data/fragments_filtered.tsv.gz"
gsetroot('/home/hpc/bonev/trackdb/mm10')
plot=FALSE
change_motifNames=TRUE
##### Initialize the seurat object ########

cortex.atac <- readRDS('data/merged_scATAC_integrated.RDS')
DefaultAssay(cortex.atac) <- 'MACS2peaks'


cortex.atac <- SetFragments(
  object = cortex.atac,
  file = fragment.path
)


#### Plotting subsets 
if (plot){
  gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
  seqlevelsStyle(gene.coords) <- 'UCSC'
  gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
  
  selected_gene <- 'Eomes'
  expansion_u <- 3e4
  expansion_d <- 5e3
  
  region <- subset(gene.coords, symbol==selected_gene)
  
  p <- scCoveragePlot(
    object=cortex.atac,
    region=region,
    annotation = mm10_txdb,
    #  cells = c(test1,test2),
    window = 100,
    downsample = 1,
    extend.upstream = expansion_u,
    extend.downstream = expansion_d,
    idents = c("NSC","IPC","PN1","PN2"),
    sep = c("-", "-")
  )
}




# Get a list of motif position weight matrices from the JASPAR database
pwm2020 <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

pwm2018 <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)

pwm <- c(pwm2020,pwm2018[!name(pwm2018)%in%name(pwm2020)])      ### Add entries which were removed in 2020 release (very important neuronal genes)
saveRDS(pwm,'results/scATAC/combined_pwm.RDS')


# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(cortex.atac), sep = c("-", "-")),
  pwm = pwm,
  genome = 'mm10',
  sep = c("-", "-")
)

# Integrate ChIP data if available 
chips_f <- list.files('/home/hpc/bonev/data/hic/data/peaks/final',pattern = 'narrowPeak',full.names = T)
chip_matrix <- matrix(data=0,nrow=length(row.names(cortex.atac)),ncol=length(chips_f))
row.names(chip_matrix) <- row.names(cortex.atac)
colnames(chip_matrix) <- gsub('.narrowPeak','',basename(chips_f))
obj_features <- StringToGRanges(rownames(cortex.atac), sep = c("-", "-"))
for (i in seq_along(chips_f)){
  chip <- read.table(chips_f[i])
  colnames(chip)[1:3] <- c('chrom','start','end')
  chip_features <- makeGRangesFromDataFrame(chip)
  overlaps <- findOverlaps(obj_features,chip_features)
  chip_matrix[overlaps@from,i] <- 1
}
saveRDS(chip_matrix,'results/scATAC/chip_matrix.RDS')

# Compute chromvar deviations for chips ######

cortex.atac <- RunChromVAR(
  object = cortex.atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  new.assay.name = 'chip_chromvar',
  motif.matrix=chip_matrix,
  verbose=T
)

cortex.atac@assays$chip_chromvar <- cortex.atac@assays$chromvar

# Compute chromvar deviations for motifs ######

motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pwm
)

cortex.atac[['MACS2peaks']] <- AddMotifObject(
  object = cortex.atac[['MACS2peaks']],
  motif.object = motif
)

cortex.atac <- RegionStats(
  object = cortex.atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  sep = c("-", "-")
)

cortex.atac <- RunChromVAR(
  object = cortex.atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

########################################

DefaultAssay(cortex.atac) <- 'chromvar'
chromvar_matrix <- cortex.atac@assays$chromvar@data
if (change_motifNames){
  motif.names <- name(pwm)
  row.names(chromvar_matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),row.names(chromvar_matrix))])))
  cortex.atac@assays$chromvar@data <- chromvar_matrix
}
saveRDS(cortex.atac,'data/merged_scATAC_integrated.RDS')

stop()

#### Exploratory analysis below ##############


varMotifs <- chromvar_matrix[order(matrixStats::rowVars(chromvar_matrix), decreasing = TRUE),]
varMotifs <- as.matrix(varMotifs)
varMotifs <- varMotifs[,order(as.numeric(Idents(cortex.atac)))]

############ Using pheatmap ###############
require(pheatmap)
varMotifs[varMotifs<(-5)] <- -5
varMotifs[varMotifs>5] <- 5

breaksList = seq(-5, 5, by = 0.1)
p_colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList)) # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
clust_anno <- data.frame(cluster=as.factor(Idents(cortex.atac)))

hm <- pheatmap(varMotifs,color=p_colors,breaks = breaksList,cluster_cols=F,show_colnames = F,annotation_col = clust_anno)

pheatmap(varMotifs[hm$tree_row$order,],cluster_rows = F,cluster_cols = F,color=p_colors,breaks = breaksList,show_colnames = F,show_rownames = F,annotation_col = clust_anno,filename = 'plots/chromVar_noNames.png',width = 6,height=6)
pheatmap(varMotifs[hm$tree_row$order,],cluster_rows = F,cluster_cols = F,color=p_colors,breaks = breaksList,show_colnames = F,show_rownames = T,annotation_col = clust_anno,filename = 'plots/chromVar.png',width = 12,height=18,fontsize_row = 6)

############ Plot using complex heatmap ###############

require(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))

nMotifs <- 250

chromvar_matrix <- cortex.atac@assays$chip_chromvar@data
varMotifs <- chromvar_matrix[order(matrixStats::rowVars(chromvar_matrix), decreasing = TRUE),]
varMotifs <- as.matrix(varMotifs)
varMotifs <- varMotifs[,order(as.numeric(Idents(cortex.atac)))]

motifs_to_mark <- c('Neurog2','Neurod2','Bhlhe22','Ascl1','Tgif2','Tcf3',"Id4",
                    "Tbr1","Eomes","Mef2c","Tcf15","Elf1","Sp8","Pou3f3","Rfx3","Tead2","Nr2f1","Fos","Rel","Rbpj",
                    "Rorb","Pax3","Sox9","Sox10","Sox13","Jun","Nfix","Emx1","Emx2","Lhx2","Lhx9")
motifs_to_mark <- c('Neurog2',"Id4","Eomes","Mef2c","Lhx2",'Fezf2','Pou3f2',"Neurod1","Tal1:tcf3","Neurod2","Olig1",'Bhlhe22','Nfic')
require(scales)
col.list <- cluster_colors(length(levels(cortex.atac)))
#col.list <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(Idents(cortex.atac))))
names(col.list) <- levels(cortex.atac)
test <- as.data.frame(varMotifs)
test$n <- seq(1:nrow(test))
test <- test[row.names(test)%in%motifs_to_mark,]
ra = rowAnnotation(foo = anno_mark(at = test$n, labels = row.names(test),labels_gp = gpar(fontsize = 16)))
col_names <- Idents(cortex.atac)[order(as.numeric(Idents(cortex.atac)))]
dup_indices <- !(duplicated(col_names))
mark_pos <- round(seq(1:length(col_names))[dup_indices] + table(col_names)/2)
ha = HeatmapAnnotation(foo = anno_mark(at = mark_pos, labels = col_names[dup_indices],labels_gp = gpar(fontsize = 16),labels_rot = 45),cluster = Idents(cortex.atac)[order(as.numeric(Idents(cortex.atac)))],col = list(cluster=col.list),show_legend = F)


HM <- Heatmap(varMotifs,cluster_columns=F,show_row_dend = F,cluster_rows = T,show_row_names = T,show_column_names = F,col = col_fun,top_annotation = ha,right_annotation = ra)

pdf('plots/scATAC/chromVar_chip.pdf',width=12,height=10)
HM
dev.off()

############ Plot using complex heatmap ########################### Plot using complex heatmap ########################### Plot using complex heatmap ###############


################# Plot motif activity in umap coordinates ###################
genes <- c('Neurog2','Neurod2','Ascl1','Eomes','Pou3f2','Lhx2')
genes <- head(enriched.motifs$motif.name)
genes <- capitalize(tolower(as.vector(genes)))
pdf(paste0('plots/scATAC/TF_motifs.pdf'),height=16,width=16)
p1 <- FeaturePlot(
  object = cortex.atac,
  cols = gene_colors(10),
  features = genes,
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.01,ncol = 3
)
p1 + theme_cowplot()
dev.off()

###################################################################

differential.activity <- FindMarkers(
  object = cortex.atac,
  ident.1 = 'Early AP',
  ident.2 = 'SCPN',
  only.pos = FALSE,
  test.use = 'LR',
  min.pct = 0.25,
  logfc.threshold = 0.25,
  latent.vars = 'nCount_MACS2peaks'
)




# Find differentially accesible peaks between two clusters 

da_peaks <- FindMarkers(               
  object = cortex.atac,
  ident.1 = 'Early AP',
  ident.2 = 'SCPN',
  only.pos = FALSE,
  test.use = 'LR',
  min.pct = 0.25,
  logfc.threshold = 0.25,
  latent.vars = 'nCount_peaks'
)

#### Or for all clusters ###################
DefaultAssay(cortex.atac) <- 'MACS2peaks'
all_da_peaks <- FindAllMarkers(
  object=cortex.atac,
  only.pos = FALSE,
  test.use = 'LR',
  min.pct = 0.25,
  logfc.threshold = 0.25,
  latent.vars = 'nCount_MACS2peaks'
)
saveRDS(all_da_peaks,'results/scATAC/all_da_peaks.tsv',quote=F,sep='\t')

DefaultAssay(cortex.atac) <- 'chromvar'
all_da_motifs <- FindAllMarkers(
  object=cortex.atac,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.25,
  logfc.threshold = 0.25,
  latent.vars = 'nCount_MACS2peaks'
)
write.table(all_da_motifs,'results/scATAC/all_da_motifs_pos.tsv',quote=F,sep='\t',col.names=T,row.names=F)

top10 <- all_da_motifs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf('plots/scATAC/motif_clusters_pos.pdf',height=18,width=10)
DoHeatmap(cortex.atac, features = top10$gene,slot = 'data') + NoLegend()
dev.off()

sel_peaks <- readRDS('results/old/P2G_binaryMat.RDS')
peak_set_NPC <- sel_peaks[sel_peaks$`BP late.y`==1|sel_peaks$`early BP.y`==1|sel_peaks$`neurogenic early AP/BP.y`==1|sel_peaks$`Early AP.y`==1|sel_peaks$`late AP.y`==1,]
peak_set_CN <- sel_peaks[sel_peaks$CthPN.y==1|sel_peaks$SCPN.y==1|sel_peaks$`CPN E16.y`==1|sel_peaks$`CPN E18.y`==1,]
peak_set_NPC <- peak_set_NPC[order(peak_set_NPC$NPCscore-peak_set_NPC$CNscore,decreasing=T),]
peak_set_CN <- peak_set_CN[order(peak_set_CN$CNscore-peak_set_CN$NPCscore,decreasing=T),]
peaks_NPC <- separate(peak_set_NPC,col='Row.names', c("chrom1", "start1","end1"), "-")
peaks_NPC <- gintervals.2d(chroms1 = as.character(peaks_NPC$chrom1),starts1 = (as.numeric(peaks_NPC$start1)+5000),ends1 = (as.numeric(peaks_NPC$start1)+5001),chroms2 = peaks_NPC$gene_chr,starts2 = as.numeric(peaks_NPC$gene_start),ends2 = as.numeric(peaks_NPC$gene_start+1))
peaks_CN <- separate(peak_set_CN,col='Row.names', c("chrom1", "start1","end1"), "-")
peaks_CN <- gintervals.2d(chroms1 = as.character(peaks_CN$chrom1),starts1 = (as.numeric(peaks_CN$start1)+5000),ends1 = (as.numeric(peaks_CN$start1)+5001),chroms2 = peaks_CN$gene_chr,starts2 = as.numeric(peaks_CN$gene_start),ends2 = as.numeric(peaks_CN$gene_start+1))



# Test the differentially accessible peaks for overrepresented motifs
peak_set1 <- sel_peaks[sel_peaks$`neurogenic early AP/BP.y`==1,]
enriched.motifs1 <- FindMotifs(
  object = cortex.atac,
  features = as.character(peak_set1$Row.names)
)

peak_set2 <- sel_peaks[sel_peaks$`late AP.y`==1,]
enriched.motifs2 <- FindMotifs(
  object = cortex.atac,
  features = as.character(peak_set_CN$Row.names[1:500])
)
# Plot overrepresented motif sequences

motif.names <- name(pwm)
motif.names <- capitalize(tolower(as.vector(motif.names)))
names(motif.names) <- names(name(pwm))
pdf('plots/neuroAP_motifs.pdf',height=12,width=12)
MotifPlot(
  object = cortex.atac,
  motifs = enriched.motifs1$motif[1:60]
)
dev.off()
cortex <- readRDS('data/merged_scRNA_filtered_IDs.RDS')
features <- c('Lhx9','Lhx2','Ldb2','Tead3','Lsd1')
FeaturePlot(cortex,features = 'Lrfn5' ,cols = c('grey','red','black'),pt.size = 0.1,label = F,ncol = 3)
FeaturePlot(cortex.atac,features = 'Ascl1' ,cols = c('grey','red','black'),pt.size = 0.1,label = F,ncol = 3)

############### Calculate 3D genome interactions between accessible motifs ################
require(tidyr)
require(plyr)
require(dplyr)
motif.matrix <- cortex.atac@assays$MACS2peaks@misc$motif@data
pwm <- cortex.atac@assays$MACS2peaks@misc$motif@pwm
motif.names <- name(pwm)
colnames(motif.matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])))


peaks <- read.table('/home/hpc/bonev/data/hic/data/peaks/final/NPC_Neurog2.bed')
peaks <- gintervals(as.character(peaks[,1]),as.numeric(peaks[,2]),as.numeric(peaks[,3]))
peak_index <- 'Neurog2'
nPeaks <- 2000
interval_window <- 2000
min_dist <- 5e3
max_dist <- 2e6
expand=c(-1500,1500)
domains <- gintervals.load("hic.ncx_Hes5.ins_250_domains_expanded")
tracks <- c('hic.ES.score_k200','hic.ncx_Hes5.score_k200','hic.ncx_Dcx.score_k200')
peaks <- motif.matrix[motif.matrix[,grep(peak_index,colnames(motif.matrix))]==1,grep(peak_index,colnames(motif.matrix))]
peak_sums <- as.data.frame(rowSums(cortex.atac@assays$MACS2peaks@counts[names(peaks),]))
colnames(peak_sums) <- 'coverage'
peak_sums$peaks <- row.names(peak_sums)
peak_sums <- peak_sums[head(order(peak_sums[,1],decreasing=T),n=nPeaks),]
peak_sums <- separate(peak_sums,col='peaks', c("chrom", "start","end"), "-")
peaks <- gintervals(as.character(peak_sums$chrom),as.numeric(peak_sums$start),as.numeric(peak_sums$end))
peaks <- peaks[peaks$chrom %in% gintervals.all()$chrom,]
options(gmax.data.size=1e7)
options(gmultitasking=FALSE)
for (i in 1:length(tracks)){
  gvtrack.create(paste0('v_',tracks[i]),tracks[i],'max')
  gvtrack.iterator.2d(paste0('v_',tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
}

grid <- construct.grid(peaks,peaks,min_dist,max_dist)
dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
intra_grid <- gintervals.intersect(grid,dom)
inter_2d = construct.grid2(domains,domains,min_dist,max_dist)			
inter_grid <- gintervals.intersect(grid,inter_2d)

intra_scores <- gextract(paste0('v_',tracks),intervals = intra_grid,iterator = intra_grid,band = -c(2e6,5e3))
inter_scores <- gextract(paste0('v_',tracks),intervals = inter_grid,iterator = inter_grid,band = -c(2e6,5e3))

###########################################################################################################
