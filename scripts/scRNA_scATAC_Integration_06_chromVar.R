library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(matrixStats)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(Hmisc)
library(RColorBrewer)
require(misha)
library(BiocParallel)
set.seed(1234)
register(MulticoreParam(8, progressbar = TRUE))

source('scripts/aux_functions.R')
source('scripts/config.R')

cortex.rna <- readRDS('data/merged_scRNA_unfiltered_IDs.RDS')
expr_matrix <- cortex@assays$RNA@scale.data

change_motifNames=TRUE
target_TF <- 'Pax6'
chip_f <- '/home/hpc/bonev/data/hic/data/peaks/final/NPC_Pax6.narrowPeak'

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

### Identify most likely targets per TF ######
df_all <- read.table('results/integrated_P2G.tsv',sep='\t',header=T)
#df_all <- as.data.frame(mcols(readRDS("results/P2G-Links.RDS")$linksAll))
df <- df_all[abs(df_all$Correlation>0.35),]

###  Create Motif Matrix ##########

motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(df$peakName, sep = c("_", "_")),
  pwm = pwm,
  genome = 'mm10',
  sep = c("_", "_")
)
motif.names <- name(pwm)
colnames(motif.matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])))

### Import ChIP if exists ######
if (import_ChIP){
  chip <- read.table(chip_f)
  colnames(chip)[1:6] <- c('chrom','start','end','name','width','strand')
  df_features <- StringToGRanges(df$peakName, sep = c("_", "_"))
  chip_features <- makeGRangesFromDataFrame(chip)
  overlaps <- findOverlaps(df_features,chip_features)
  sel_df <- df[overlaps@from,]
} else {
  sel_df <- df[motif.matrix[,target_TF]==1,]
  }

# Scan the DNA sequence of each peak for the presence of each motif


linkageScore <- split(sel_df$Correlation^2, as.character(as.vector(sel_df$gene_name))) %>% 
  lapply(., sum) %>% 
  unlist(use.names=TRUE)

linkageScore <- linkageScore[order(linkageScore,decreasing = F)]
plot_df <- data.frame(order=1:length(linkageScore),linkageScore=linkageScore)
plot_df$expr_cor <- 0

for (i in 1:nrow(plot_df)){
  plot_df$expr_cor[i] <- cor.test(expr_matrix[row.names(plot_df)[i],],expr_matrix[target_TF,],method = 'spearman')$estimate 
}
plot_df <- plot_df[order(plot_df$expr_cor,plot_df$linkageScore,decreasing = F),]
plot_df$Gene <- row.names(plot_df)
#### Plot the ranked plot labeling top 10 genes ########
p <- ggplot(plot_df,aes(x=expr_cor,y=linkageScore)) + geom_point() +ggtitle(label = target_TF)
p1 <- p + ggrepel::geom_label_repel(
  data = rbind(tail(plot_df[plot_df$linkageScore>0.7,],10),plot_df['Rnd2',]), size = 3,
  aes(x=expr_cor,y=linkageScore,color=NULL,label=Gene))
name_f <- paste0('plots/',target_TF,'linkageScore.pdf')
pdf(name_f,width=6,height=6)
print(p1)
dev.off()


###########
cortex.atac <- readRDS('data/merged_scATAC_integrated_cicero.RDS')
fragment.path <- "data/fragments_filtered.tsv.gz"
DefaultAssay(cortex.atac) <- 'MACS2peaks'
cortex.atac <- SetFragments(
  object = cortex.atac,
  file = fragment.path
)

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

delta_lim <- -100
hic_lim <- -100
class <- 'posCor'
cluster='DL'
idents_to_keep <- c('AP1','AP2','AP3','BP','DL1','DL_CPN1','DL2','DL_CPN2','DL3')
domain='intraTAD'
#df <- df_all[abs(df_all$Correlation)>0.35&df_all$cluster==cluster,]
df <- df_all[df_all$class==class&df_all$cluster==cluster&df_all$domain==domain,]
df <- df[(df$deltaNPC_CN_HiC>delta_lim)&(rowMaxs(as.matrix(df[,grep('score',colnames(df))]))>hic_lim),]
subset.atac <- subset(cortex.atac,features=gsub('_','-',df$peakName),idents=Idents(cortex.atac)[Idents(cortex.atac)%in%idents_to_keep])


# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(subset.atac), sep = c("-", "-")),
  pwm = pwm,
  genome = 'mm10',
  sep = c("-", "-")
)

# Create a new Motif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pwm
)

# Add the Motif object to the assay
subset.atac[['MACS2peaks']] <- AddMotifObject(
  object = subset.atac[['MACS2peaks']],
  motif.object = motif
)

subset.atac <- RegionStats(
  object = subset.atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  sep = c("-", "-")
)

#### Option 1 Calculate motif enrichment / cluster ####

enriched.motifs <- FindMotifs_fisher(
  object = cortex.atac,
  features = gsub('_','-',df$peakName),
  background=length(df$peakName)
)

diff.genes <- FindMarkers(
  object = cortex.rna,ident.1 = 'BP',ident.2=NULL,only.pos = F
)

enriched.motifs$motif.name <- capitalize(tolower(as.vector(enriched.motifs$motif.name)))
enriched.motifs$fold.enrichment <- log2(enriched.motifs$fold.enrichment)
enriched.motifs$pvalue <- -log10(enriched.motifs$pvalue)
enriched.motifs$col <- 'grey80'
enriched.motifs$col[enriched.motifs$fold.enrichment>0.25&enriched.motifs$pvalue>=2] <- 'red'
enriched.motifs$col[enriched.motifs$fold.enrichment<(-0.25)&enriched.motifs$pvalue>=2] <- 'green'
p <- ggplot(enriched.motifs,aes(x=fold.enrichment,y=pvalue,col=col)) + geom_point() 
p <- p + scale_colour_manual(values = c("green", "grey80",'red'),labels=NULL) + xlim(-1.8,1.8)
p1 <- p + ggrepel::geom_label_repel(
  data = rbind(head(enriched.motifs[order(enriched.motifs$fold.enrichment,enriched.motifs$pvalue,decreasing=T),],5),head(enriched.motifs[order(enriched.motifs$fold.enrichment*(-1),enriched.motifs$pvalue,decreasing=T),],5)), size = 3,
  aes(x=fold.enrichment,y=pvalue,color=NULL,label=motif.name))
name_f <- paste0('plots/MotifEnr_',cluster,'.pdf')
pdf(name_f,height=6,width=6)
print(p1 + theme(legend.position = "none"))
dev.off()

#### Running ChromVar to get per-cell motif accessibility score

subset.atac <- RunChromVAR(
  object = subset.atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

nMotifs <- 250
DefaultAssay(subset.atac) <- 'chromvar'
chromvar_matrix <- subset.atac@assays$chromvar@data
if (change_motifNames){
  motif.names <- name(pwm)
  row.names(chromvar_matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),row.names(chromvar_matrix))])))
  subset.atac@assays$chromvar@data <- chromvar_matrix
}

varMotifs <- chromvar_matrix[head(order(matrixStats::rowVars(chromvar_matrix), decreasing = TRUE), nMotifs),]
varMotifs <- as.matrix(varMotifs)
idents <- Idents(subset.atac)
idents <- factor(idents,levels(idents)[c(6,4,7,3,5,2,1)])
subset.atac@active.ident <- idents
varMotifs <- varMotifs[,order(as.numeric(Idents(subset.atac)))]

require(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

nMotifs <- 250
DefaultAssay(subset.atac) <- 'chromvar'
chromvar_matrix <- subset.atac@assays$chromvar@data
varMotifs <- chromvar_matrix[head(order(matrixStats::rowVars(chromvar_matrix), decreasing = TRUE), nMotifs),]
varMotifs <- as.matrix(varMotifs)
varMotifs <- varMotifs[,order(as.numeric(Idents(subset.atac)))]

motifs_to_mark <- c('Neurog2','Neurod2','Bhlhe22','Ascl1','Tgif2','Tcf3',"Id4",
                    "Tbr1","Eomes","Mef2c","Tcf15","Elf1","Sp8","Pou3f3","Rfx3","Tead2","Nr2f1","Fos","Rel","Rbpj",
                    "Rorb","Pax6","Sox9","Sox10","Sox13","Jun","Nfix","Emx1","Emx2","Lhx2","Lhx9")
#motifs_to_mark <- c('Jun','Tead2','Nfix','Rfx2','Gli3','Pou3f2','Cux2','Tcf12(var.2)','Meis2')
require(scales)
col.list <- hue_pal()(length(levels(Idents(subset.atac))))
#col.list <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(Idents(subset.atac))))
names(col.list) <- levels(Idents(subset.atac))
test <- as.data.frame(varMotifs)
test$n <- seq(1:nrow(test))
test <- test[row.names(test)%in%motifs_to_mark,]
ra = rowAnnotation(foo = anno_mark(at = test$n, labels = row.names(test),labels_gp = gpar(fontsize = 16)))
col_names <- Idents(subset.atac)[order(as.numeric(Idents(subset.atac)))]
dup_indices <- !(duplicated(col_names))
mark_pos <- round(seq(1:length(col_names))[dup_indices] + table(col_names)/2)
ha = HeatmapAnnotation(foo = anno_mark(at = mark_pos, labels = col_names[dup_indices],labels_gp = gpar(fontsize = 16),labels_rot = 45),cluster = Idents(subset.atac)[order(as.numeric(Idents(subset.atac)))],col = list(cluster=col.list),show_legend = F)


HM <- Heatmap(varMotifs,cluster_columns=F,show_row_dend = F,cluster_rows = T,show_row_names = F,show_column_names = F,col = col_fun,top_annotation = ha,right_annotation = ra)

pdf(paste0('plots/scATAC/',cluster,'_',domain,'_hicLim',hic_lim,'_chromvar.pdf'),width=12,height=12)
HM
dev.off()

varMotifs <- varMotifs[row_order(HM),]
write.table(row.names(varMotifs),file = paste0('results/scATAC/',cluster,'_',domain,'_chromvar.tsv'),quote = F,sep='\t')
