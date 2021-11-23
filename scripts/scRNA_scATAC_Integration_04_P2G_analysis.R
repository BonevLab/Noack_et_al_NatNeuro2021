############################################

########### Plot Complex Heatmap #####################

library(Seurat)
library(Signac)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
library(plyr)
library(dplyr)
require(ComplexHeatmap)
library(circlize)
library(scales)
library(misha)

set.seed(1)

source('scripts/config.R')
source('scripts/aux_functions.R')
gsetroot('/home/hpc/bonev/trackdb/mm10')

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

#----------------------------
# Get Inputs
#----------------------------
outMatch <- readRDS('results/P2G-Links.RDS')
names(outMatch)[3:6] <- c('posCor','all','negCor','noCor')

scRNA <- readRDS('data/merged_scRNA_unfiltered_IDs.RDS') 
sePB_all <- readRDS("results/scATAC/Cluster_PseudoBulk-Summarized-Experiment.RDS")
scRNA$batch <- 'rep1'
scRNA$batch[grep('rep2',colnames(scRNA))] <- 'rep2'
motifs_to_mark <- c("Hes1","Neurog2","Eomes","Bcl11b")
#----------------------------
# Set Parameters
#----------------------------
integrated_list <- list()
for (sel_peaks_name in c('posCor','negCor','noCor')){
  
  selected_peaks <- outMatch[[sel_peaks_name]]
  
  #----------------------------
  # PseudoBulk ATAC
  #----------------------------
  
  #sePB <- subsetByOverlaps(sePB_all,selected_peaks,type='any')
  sePB <- sePB_all[match(gsub('_','-',selected_peaks$peakName),row.names(sePB_all)),]
  ##############################
  ###### Pseudobulk RNA
  #############################
  if (!file.exists("results/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")){
    se <- SummarizedExperiment(
      assays = SimpleList(counts = scRNA@assays$RNA@counts), 
      colData=data.frame(Clusters= Idents(scRNA),Group=scRNA$batch)
    )
    
    objSE <- createPseudoBulk(
      mat = assay(se), 
      groups=Idents(scRNA), 
      labels=colData(se)$Group, 
      ceiling = 1, 
      minCells = 100, 
      maxCells = 500, 
      minReps = 2, 
      prior.count = 3, 
      nSim = 250)
    
    cdPB <- DataFrame(
      row.names=colnames(objSE[[1]]), 
      Group = stringr::str_split(colnames(objSE[[1]]), "\\._\\.", simplify = TRUE)[,1],
      #  Group = levels(Idents(pbmc))[as.numeric(stringr::str_split(colnames(objSE[[1]]), "._.", simplify = TRUE)[,1])+1],
      Type  = stringr::str_split(stringr::str_split(colnames(objSE[[1]]), "._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,1],
      Cells = as.integer(stringr::str_split(stringr::str_split(colnames(objSE[[1]]), "\\._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,2]),
      Replace = stringr::str_split(stringr::str_split(colnames(objSE[[1]]), "\\._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,3]
    )
    
    seRNA <- SummarizedExperiment(assays = SimpleList(counts = objSE[[1]]), colData = cdPB)
    metadata(seRNA )$matCS <- objSE[[2]] 
    
    saveRDS(seRNA, "results/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")
  } else {
    seRNA <- readRDS("results/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")
  }
  
  ##################################
  ###### Subset unwanted clusters
  ##################################
 
  seRNA <- seRNA[row.names(seRNA)%in%selected_peaks$gene_name,]
 # seRNA <- seRNA[,grep('IN',colnames(seRNA),invert=T)]
 # seRNA <- seRNA[,grep('CR',colnames(seRNA),invert=T)]
 # seRNA <- seRNA[,grep('Blood',colnames(seRNA),invert=T)]
 # seRNA <- seRNA[,grep('GE',colnames(seRNA),invert=T)]
 # seRNA <- seRNA[,grep('MSN',colnames(seRNA),invert=T)]
 # seRNA <- seRNA[,grep('Pericytes',colnames(seRNA),invert=T)]
 # seRNA <- seRNA[,grep('Microglia',colnames(seRNA),invert=T)]
  
 # sePB <- sePB[,grep('IN',colnames(sePB),invert=T)]
 # sePB <- sePB[,grep('CR',colnames(sePB),invert=T)]
 # sePB <- sePB[,grep('Blood',colnames(sePB),invert=T)]
 # sePB <- sePB[,grep('GE',colnames(sePB),invert=T)]
 # sePB <- sePB[,grep('MSN',colnames(sePB),invert=T)]
 # sePB <- sePB[,grep('Pericytes',colnames(sePB),invert=T)]
 # sePB <- sePB[,grep('Microglia',colnames(sePB),invert=T)]
 # sePB <- sePB[,grep('E12',colnames(sePB),invert=T)]
  #----------------------------
  # Unique Features
  #----------------------------
  params <- list(
    padj = 1,
    minSdRatio = 0,
    minLFC = 0.0,
    zCutoff = 0,
    breakPt = "last",
    groupMin = 20,
    maxGroupSize = max(floor(length(unique(colData(sePB)$Group))/3), 2),
    date = Sys.Date()
  )
  
  #Unique Features
  uf <- uniqueFeatures(
    edgeR::cpm(assay(sePB),log=TRUE,prior.count=1),
    groups = colData(sePB)$Group,
    padj = params$padj,
    minSdRatio = params$minSdRatio,
    minLFC = params$minLFC,
    zCutoff = params$zCutoff,
    breakPt = params$breakPt,
    groupMin = params$groupMin,
    maxGroupSize = params$maxGroupSize,
    clusterCols = F
  )
  
  uf_rna <- uniqueFeatures(
    edgeR::cpm(assay(seRNA),log=TRUE,prior.count=1),
    groups = colData(seRNA)$Group,
    padj = params$padj,
    minSdRatio = params$minSdRatio,
    minLFC = params$minLFC,
    zCutoff = params$zCutoff,
    breakPt = params$breakPt,
    groupMin = params$groupMin,
    maxGroupSize = params$maxGroupSize,
    clusterCols = F
  )
  
  m <- uf$groupMat
  #m <- m[rowMax(m)>3,]
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
  z[z > 2] <- 2
  z[z < -2] <- -2
  
  m_rna <- uf_rna$groupMat
  z_rna <- sweep(m_rna - rowMeans(m_rna), 1, matrixStats::rowSds(m_rna), `/`)
  z_rna[z_rna > 2] <- 2
  z_rna[z_rna < -2] <- -2
  
  linksSig_IDs <- as.data.frame(cbind(as.vector(seqnames(selected_peaks)),as.numeric(start(selected_peaks )),as.numeric(end(selected_peaks ))))
  selected_peaks$IDs <- unite(linksSig_IDs,col='peaks_ID',sep='-')[,1]
  
  anno <- as.data.frame(mcols(selected_peaks),stringsAsFactors=FALSE)
  anno <- cbind(z[match(anno$IDs,row.names(z)),],anno)
  row.names(anno) <- 1:nrow(anno)
 # anno <- anno[,-1]
 # anno <- anno[complete.cases(anno),]
  z <- anno[,1:ncol(z)]
  anno <- anno[,(ncol(z)+1):ncol(anno)]
  
  comb_z <- as.data.frame(matrix(NA,nrow(z),ncol(z_rna)))
  row.names(comb_z) <- row.names(z)
  colnames(comb_z) <- colnames(z_rna)
  comb_z$gene_name <- 'unknown'
  
  for (i in 1:nrow(comb_z)){
    if (sum(row.names(z_rna)%in%anno$gene_name[i])>0){
      comb_z[i,1:(ncol(comb_z)-1)] <- z_rna[row.names(z_rna)%in%anno$gene_name[i],]
      comb_z$gene_name[i] <- anno$gene_name[i]
    }
  }
  
  ############ Plot using complex heatmap ###############
  
  
  #motifs_to_mark <- c('Eomes')
  #motifs_to_mark <- motifs_to_mark[c(3,4,6)]
  reduced_anno <- anno[anno$gene_name%in%motifs_to_mark, ]
  reduced_anno <- reduced_anno[order(reduced_anno$Correlation,decreasing = T),]
  dim(reduced_anno)
  reduced_anno_att <- reduced_anno %>% group_by(gene_name) %>% slice(1:2)
  reduced_anno <- reduced_anno[reduced_anno$peakName%in%reduced_anno_att$peakName,]
  col.list <- hue_pal()(ncol(z))
  #col.list <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(Idents(scATAC))))
  names(col.list) <- colnames(z)
  test <- as.data.frame(z)
  test$n <- seq(from=1,to=nrow(z),by=1)
  test <- merge(test,reduced_anno,by='row.names',sort=F)
  test$labels <- paste0(test$gene_name,':',test$distance)
  #idxs <- c('Fezf2:−99404','Sox5:24831','Bcl11b:−74950','Sox5:−75708','Gad1:0','Gad2:−1152','Olig1:2745','Hes1:−5635','Hes1:−4379','Eomes:4005','Btg2:1308','Neurod2:1062')
  #idxs <- c(50915,58107,66287,93244,96456,112978,112522,116861)
  #test <- test[test$n%in%idxs,]
  ra = rowAnnotation(foo = anno_mark(at = test$n, labels = test$labels,labels_gp = gpar(fontsize = 16)))
  col_names <- stringr::str_split(colnames(z), "\\._\\.", simplify = TRUE)[,1]
  dup_indices <- !(duplicated(col_names))
  mark_pos <- round(seq(1:length(col_names))[dup_indices])
  ha = HeatmapAnnotation(foo = anno_mark(at = mark_pos, labels = col_names[dup_indices],labels_gp = gpar(fontsize = 16),labels_rot = 45),cluster=colnames(z),col = list(cluster=col.list),show_legend = F)
  
  
  HM <- Heatmap(as.matrix(z),cluster_columns=F,cluster_rows=F,show_row_names = F,show_column_names = F,col = col_fun,right_annotation = ra,top_annotation = ha)
  
  pdf(paste0('plots/',sel_peaks_name,'_scATAC.pdf'),width=10,height=16)
  print(HM)
  dev.off()
  
  ######### Plot RNA #######################
  z <- as.matrix(comb_z[,-ncol(comb_z)])
  z <- z[,c('AP','BP','DL1','DL_CPN1','DL2','DL_CPN2','DL3','DL4')]
  #motifs_to_mark <- motifs_to_mark[c(3,4,6)]
  reduced_anno <- anno[anno$gene_name%in%motifs_to_mark, ]
  dim(reduced_anno)
  reduced_anno <- reduced_anno[order(reduced_anno$Correlation,decreasing = T),]
  dim(reduced_anno)
  reduced_anno_att <- reduced_anno %>% group_by(gene_name) %>% slice(1:2)
  reduced_anno <- reduced_anno[reduced_anno$peakName%in%reduced_anno_att$peakName,]
    #reduced_anno <- (reduced_anno %>% group_by(gene) %>% slice(1))
  col.list <- hue_pal()(ncol(z))
  #col.list <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(Idents(scATAC))))
  names(col.list) <- colnames(z)
  test <- as.data.frame(z)
  test$n <- seq(from=1,to=nrow(z),by=1)
  test <- merge(test,reduced_anno,by='row.names',sort=F)
  test$labels <- paste0(test$gene_name,':',test$distance)
  #idxs <- c('Fezf2:−99404','Sox5:24831','Bcl11b:−74950','Sox5:−75708','Gad1:0','Gad2:−1152','Olig1:2745','Hes1:−5635','Hes1:−4379','Eomes:4005','Btg2:1308','Neurod2:1062')
  #idxs <- c(50915,58107,66287,93244,96456,112978,112522,116861)
  #test <- test[test$n%in%idxs,]
  ra = rowAnnotation(foo = anno_mark(at = test$n, labels = test$labels,labels_gp = gpar(fontsize = 16)))
  col_names <- stringr::str_split(colnames(z), "\\._\\.", simplify = TRUE)[,1]
  dup_indices <- !(duplicated(col_names))
  mark_pos <- round(seq(1:length(col_names))[dup_indices])
  ha = HeatmapAnnotation(foo = anno_mark(at = mark_pos, labels = col_names[dup_indices],labels_gp = gpar(fontsize = 16),labels_rot = 45),cluster=colnames(z),col = list(cluster=col.list),show_legend = F)
  
  HM <- Heatmap(z,cluster_columns=F,cluster_rows=F,show_row_names = F,show_column_names = F,col = col_fun,right_annotation = ra,top_annotation = ha)
  
  pdf(paste0('plots/',sel_peaks_name,'_scRNA.pdf'),width=10,height=16)
  print(HM)
  dev.off()
  
  ###### Output mat table #####
  test <- as.data.frame(z)
  test$n <- seq(from=1,to=nrow(z),by=1)
  test <- merge(test,anno,by='row.names',sort=F)
  test$labels <- paste0(test$gene_name,':',test$distance)
  test2 <- merge(test,uf$binaryMat,by.x='IDs',by.y='row.names',sort=F)
  write.table(test2,file = paste0('results/P2G_binaryMat_',sel_peaks_name,'.tsv'),quote = F,col.names=T,row.names=F,sep='\t')
  #colnames(test2)[which(colnames(test2)=='DL1.y')] <- 'BP2'
  
  
  for (select_group in  c('NSC','IPC','PN')){
    indx <- intersect(colnames(test2)[grep(select_group,colnames(test2))],colnames(test2)[(ncol(test2)+1-ncol(m)):ncol(test2)])
    subset_peaks <- test2[rowSums(test2[,indx,drop=F])>=1,]
    
    name <- paste0('results/beds/',select_group,'_',sel_peaks_name,'_RA.bed')
    bed_f <- DataFrame(subset_peaks[,c('gene_chr','gene_start','gene_start')])
    colnames(bed_f) <- c('chrom','start','end')
    bed_f[subset_peaks$gene_strand=='+','end'] <- bed_f$end[subset_peaks$gene_strand=='+']+1
    bed_f[subset_peaks$gene_strand=='-','start'] <- bed_f$start[subset_peaks$gene_strand=='-']-1
    bed_f1 <- makeGRangesFromDataFrame(bed_f)
    rtracklayer::export.bed(bed_f1,con = name)
    
    ### Export data for further analysis
    df <- subset_peaks[,c('peakName','gene_name','distance','Correlation','ESscore','NPCscore','CNscore','deltaNPC_CN_HiC','domain')]
    df$cluster <- select_group
    df$class <- sel_peaks_name
    integrated_list[[paste0(sel_peaks_name,select_group,'_')]] <- df 
    bed_f2 <- separate(subset_peaks[,'peakName',drop=F],col='peakName',into=c('chrom','start','end'),sep='_')
    
    bed_f <- as.data.frame(bed_f)
    g_beds <- gintervals.2d(chroms1 = bed_f[,1],starts1 = bed_f[,2],ends1 = bed_f[,3],chroms2 = bed_f2[,1],starts2 = as.numeric(bed_f2[,2]),ends2 = as.numeric(bed_f2[,3]))
    g_beds_LA <- g_beds[,4:6]
    g_beds_RA <- g_beds[,1:3]
    colnames(g_beds_LA) <-  c('chrom','start','end')
    colnames(g_beds_RA) <-  c('chrom','start','end')
    name_RA <- gsub('\\.bed','',name)
    save(g_beds_RA,file = name_RA)
    
    name <- paste0('results/beds/',select_group,'_',sel_peaks_name,'_LA.bed')
    bed_f2 <- makeGRangesFromDataFrame(bed_f2)
    rtracklayer::export.bed(bed_f2,con = name)
    name_LA <- gsub('\\.bed','',name)
    save(g_beds_LA,file = name_LA)  

    name <- paste0('plots/HiC/',select_group,'_',sel_peaks_name,'.pdf')
    pdf(name,width=8,height=8)
    plot(density(subset_peaks$ESscore,na.rm=T),col='blue')
    lines(density(subset_peaks$NPCscore,na.rm=T),col='green')
    lines(density(subset_peaks$CNscore,na.rm=T),col='red')  
    dev.off()
  }
  ##########################
}

df <- integrated_list %>% Reduce("rbind",.)
write.table(df,'results/integrated_P2G.tsv',quote=F,col.names=T,row.names = F,sep='\t')


hic_cutoff <- 0.5
delta_hic_cutoff <- 20
cor_cutoff <- 0.35

outMatch <- readRDS('results/P2G-Links.RDS')
df <- outMatch$linksSig
df <- df[!is.na(df$deltaHiC),]
df$gene_end <- df$gene_start + 1
df$gene_end[df$gene_strand=='-'] <- df$gene_end[df$gene_strand=='-']-2
hic_values <- as.matrix(values(df)[,c('NPCscore','CNscore')])
highCor_highHiC <- df[(df$Correlation>0.35)&(rowMaxs(hic_values,na.rm=T)>hic_cutoff)&(df$deltaHiC>delta_hic_cutoff),]
highCor <- df[df$Correlation>cor_cutoff,]
peaks <- data.frame(highCor_highHiC$IDs,highCor_highHiC$gene_chr,highCor_highHiC$gene_start,highCor_highHiC$gene_end,highCor_highHiC$Correlation)
colnames(peaks) <- c('Peak1','gene_chr','gene_start','gene_end','coaccess')
peaks$Peak2 <- unite(peaks[,c('gene_chr','gene_start','gene_end')],col='Peak2',sep='_')[,1]
peaks <- peaks[,c('Peak1','Peak2','coaccess')]
write.csv(x = peaks, file = "results/cellOracle/highCor_highHiC_connections.csv")
peaks <- data.frame(highCor$IDs,highCor$gene_chr,highCor$gene_start,highCor$gene_end,highCor$Correlation)
colnames(peaks) <- c('Peak1','gene_chr','gene_start','gene_end','coaccess')
peaks$Peak2 <- unite(peaks[,c('gene_chr','gene_start','gene_end')],col='Peak2',sep='_')[,1]
peaks <- peaks[,c('Peak1','Peak2','coaccess')]
write.csv(x = peaks, file = "results/cellOracle/highCor_connections.csv")
ciceroObj <- readRDS("results/scATAC/cicero_aggregated_accessibility_cds.RDS")
all_peaks <- rownames(exprs(ciceroObj))
write.csv(x = all_peaks, file = "results/cellOracle/all_peaks.csv")


library(cowplot)
library(ggrepel)
library(ggpointdensity)

df_all <- read.table('results/integrated_P2G.tsv',sep='\t',header=T)

for (class in c('posCor','negCor','noCor')){
  for (cluster in c('AP','BP','DL')){
    for (domain in c('intraTAD','interTAD')){
      
    }
  }
}

class <- 'negCor'
cluster='DL'
domain='intraTAD'

df <- df_all[df_all$class==class&df_all$cluster==cluster&df_all$domain==domain,]
df$col <- 'grey80'
df$col[(df$NPCscore-df$CNscore)>10] <- 'green'
df$col[(df$CNscore-df$NPCscore)>10] <- 'red'

p1 <- ggplot( df, aes( x = NPCscore, y = CNscore )) + geom_pointdensity(adjust=5, shape = 20,size=0.5,alpha=1) + scale_color_gradient(low = 'black',high = 'red')             # scale_color_gradientn(colours=tim.colors(12))
p2 =p1 + theme_cowplot() + geom_abline(slope = 1,intercept = 0,linetype = 2) + xlim(0,100) + ylim(0,100)
g_name <- paste('plots/HiC/',class,'_',cluster,'_',domain,'.pdf')
pdf(g_name,height=4,width=5)
print(p3)
dev.off()




#p1 <- ggplot( df, aes( x = NPCscore, y = CNscore )) +
df <- df[order(abs(df$NPCscore-df$CNscore),df$NPCscore,decreasing=T),]
unique(df$gene_name)[1:200]
#  geom_point( aes( color = as.vector(col)), shape = 20,size=0.5,alpha=0.5 ) 
#p1 =p1 + theme_cowplot() + scale_color_manual(values = c("grey80"="grey80",'green'='green','red'='red'),guide=F) + geom_density_2d() + stat_density_2d()
#p2 <- ggExtra::ggMarginal(p = p1, type = 'density', margins = 'both',size = 5,groupColour = TRUE,fill = 'white',bins=50)
plot_gene <- 'Sox5'
plot_what <- 'CN'

n_genes <- 2
delta_lim <- 10
data_f <- subset(df, gene_name==plot_gene&(abs(NPCscore-CNscore)>delta_lim))

if(plot_what=='CN'){
  data_f <- data_f[order(data_f$CNscore,decreasing=T),]
  data_f <- data_f[c(1:n_genes),]
  nudge_y <- data_f$CNscore - 10
  nudge_x <- data_f$NPCscore - 10
  xlims <- c(0,30)
  ylims <- c(60,90)
} else {
  data_f <- data_f[order(data_f$NPCscore,decreasing=T),]
  data_f <- data_f[c(1:n_genes),]
  nudge_y <- data_f$CNscore - 30
  nudge_x <- data_f$NPCscore + 20
  xlims <- c(70,100)
  ylims <- c(0,35)  
}


p1 <- ggplot( df, aes( x = NPCscore, y = CNscore )) + geom_pointdensity(adjust=5, shape = 20,size=0.5,alpha=1) + scale_color_gradient(low = 'black',high = 'red')             # scale_color_gradientn(colours=tim.colors(12))
p2 =p1 + theme_cowplot() + geom_abline(slope = 1,intercept = 0,linetype = 2) + xlim(0,100) + ylim(0,100)
p3 <- p2 +  geom_label_repel(data         = data_f,
                             nudge_y       = nudge_y,
                             nudge_x = nudge_x,
                             aes(label=paste0(gene_name,':',distance)),
                             size          = 4,
                             box.padding   = 0.35,
                             point.padding = 0.5,
                             force         = 10,
                             segment.size  = 0.2,
                             segment.color = "black",
                             direction     = "y",
                             xlim=xlims,
                             ylim=ylims)
p3

g_name <- paste('plots/HiC/',class,'_',cluster,'_',domain,'_',plot_gene,'.pdf')
pdf(g_name,height=4,width=5)
print(p3)
dev.off()







mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T)
conns <- readRDS('results/cicero_conns.RDS')

#Check only for NSC enhancers 
mat <- mat[mat$NSC==1,]

conns <- conns[conns$Peak1%in%mat$peakName|conns$Peak2%in%mat$peakName,]
conns_da <- StringToGRanges(conns$Peak1, sep = c("_", "_"))
conns_ga <- StringToGRanges(conns$Peak2, sep = c("_", "_"))
mat_ga <- GRanges(seqnames = mat$gene_chr,ranges = IRanges(start = mat$gene_start,end = mat$gene_start+1),strand = mat$gene_strand)

test <- findOverlaps(mat_ga,conns_ga,select='first')
mat$genePeak <- NA
mat$genePeak[!is.na(test)] <- as.character(as.vector(conns$Peak2))[test[!is.na(test)]]
mat$ciceroConn <- NA
mat$ciceroConn <- conns$coaccess[match(paste0(mat$peakName,mat$genePeak),paste0(conns$Peak1,conns$Peak2))]
mat$ciceroConn2<- NA
mat$ciceroConn2 <- conns$coaccess[match(paste0(mat$genePeak,mat$peakName),paste0(conns$Peak1,conns$Peak2))]

cor.test(mat$Correlation,mat$ciceroConn,method = 'pearson')
test <- mat[!is.na(mat$ciceroConn),]

p <-  ggplot(as.data.frame(mat[,c('ciceroConn','Correlation')]), aes( x = ciceroConn, y = Correlation )) + geom_pointdensity(adjust=1.5,alpha=1,size=1) + scale_color_gradientn(colours = colorpalette("heat",8),name='')             # scale_color_gradientn(colours=tim.colors(12))
p<- p + theme_cowplot() + geom_abline(slope = 1,intercept = 0,linetype = 2) + xlab('Cicero Co-Accessibility') + ylab('ArchR Correlation')  
p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.3, 0.97),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.4), vjust = 8))
my_grob = grid.text(paste0('r=',round(cor.test(mat$ciceroConn,mat$Correlation)$estimate,3)), x=0.85,  y=0.1, gp=gpar(col="black", fontsize=14, fontface="bold"))
p1 <- p + annotation_custom(my_grob)

pdf('plots/CiceroVsCorrelation.pdf',height=6,width=6)
print(p1)
dev.off()








