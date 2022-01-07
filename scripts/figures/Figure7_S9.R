library(Seurat)
library(Signac)
library(ggplot2)
library(ggpubr)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(dplyr)
library(JASPAR2020)
library(Matrix)
library(TFBSTools)
library(Hmisc)
library(SummarizedExperiment)
require(edgeR)
require(LSD)
library(ggrepel)
library(ggpointdensity)
library(grid)
library(ComplexHeatmap)
library(plyr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(misha)
library(seqplots)

source('scripts/figures/config.R')
source('scripts/aux_functions.R')
source('scripts/figures/plot_functions.R')

source('scripts/hic/config.R')
source('scripts/hic/scripts/main_functions.R')
source('scripts/hic/scripts/aux_functions.R')
source('scripts/hic/scripts/plot_functions.R')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))

rna.object <- readRDS('data/merged_scRNA_unfiltered_IDs.RDS')
p2glinks <- readRDS('results/P2G-Links.RDS')
seRNA_all <- readRDS("results/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")
uf_rna <- suppressMessages(uniqueFeatures(
  edgeR::cpm(assay(seRNA_all),log=TRUE,prior.count=1),
  groups = colData(seRNA_all)$Group,
  padj = 1,
  minSdRatio = 0,
  minLFC = 0,
  zCutoff = 0,
  breakPt = "last",
  groupMin = 0,
  maxGroupSize = 2,
  clusterCols = F
))

#Aux functions
ep_vs_enh <- function(df,use_chip,chip_f,tracks,expand=c(-5000,5000),anno_mat,clusters=NULL,vars,inv_size=500,intraTAD=FALSE,dist_cutoff=5000,max_dist=5e5,motif_name="Neurog2(var.2)"){
  options(gmax.data.size=5e7)
  options(gmultitasking=F)
  for (i in 1:length(tracks)){
    gvtrack.create(paste0('v_',tracks[i]),tracks[i],'max')
    gvtrack.iterator.2d(paste0('v_',tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
  }
  anno_mat$cluster <- colnames(anno_mat[,1:7])[max.col(anno_mat[,1:7])]
  anno_mat$labels <- paste0(anno_mat$gene_name,':',anno_mat$distance)
  df$labels <- paste0(df$gene_name,':',df$distance)
  df$cluster <- anno_mat$cluster[match(df$labels,anno_mat$labels)]
  if(!is.null(clusters)){
    df <- df[df$cluster%in%clusters,]
  }
  if(use_chip){
    chip <- read.table(chip_f)
    colnames(chip)[1:3] <- c('chrom','start','end')
    chip <- chip[chip$chrom%in%gintervals.all()$chrom,]
    chip <- intervals.normalize(chip,500)
    enh_ranges <- StringToGRanges(df$peakName, sep = c("_", "_"))
    prom_ranges <- GRanges(seqnames=df$gene_chr,IRanges(start=df$gene_start-5000,end=df$gene_start+5000))
    chip_features <- makeGRangesFromDataFrame(chip)
    prom_overlaps <- GenomicRanges::countOverlaps(prom_ranges,chip_features)
    df$prom_motif <- prom_overlaps
    overlaps <- GenomicRanges::countOverlaps(enh_ranges,chip_features)
    df$motif <- overlaps
  } else {
    motif_list <- list()
    pwm <- readRDS(chip_f)
    motif.matrix <- CreateMotifMatrix(features = StringToGRanges(df$peakName, sep = c("_", "_")),pwm = pwm,genome = 'mm10',sep = c("-", "-"))
    motif.names <- name(pwm)
    colnames(motif.matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])))
    df$motif <- as.vector(motif.matrix[,motif_name])
    prom_ranges <- GRanges(seqnames=df$gene_chr,IRanges(start=df$gene_start-2000,end=df$gene_start+500))
    start(prom_ranges)[df$gene_strand=='-'] <- df$gene_start[df$gene_strand=='-']-500
    end(prom_ranges)[df$gene_strand=='-'] <- df$gene_start[df$gene_strand=='-']+2000
    motif.matrix <- CreateMotifMatrix(features = prom_ranges,pwm = pwm,genome = 'mm10',sep = c("-", "-"))
    motif.names <- name(pwm)
    colnames(motif.matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])))
    df$prom_motif <- as.vector(motif.matrix[,motif_name])
  }
  df <- as.data.frame(mcols(df))
  df$binding <- 'None'
  df$binding[df$motif>=1&df$prom_motif>=1] <- 'E-P'
  df$binding <- factor(df$binding,levels=c('None','P-P','E-E','E-P'))
  colnames(df) <- gsub('score','',colnames(df))
  df_sub <- df[abs(df$distance)>dist_cutoff,c('binding',vars)]
  ###
  df_sub <- df_sub[complete.cases(df_sub),]
  for (i in 1:2){
    feature <- c('motif','prom_motif')[i]
    features <- df[,feature,drop=F]
    features_M <- unique(row.names(features[features[,1]>=1,,drop=F]))
    features_M <- StringToGRanges(features_M, sep = c("_", "_"))
    features_M <- intervals.centers(gintervals(seqnames(features_M),start(features_M),end(features_M)))
    grid_M <- unique(construct.grid(features_M,features_M,dist_cutoff,max_dist))
    grid_M <- gintervals.canonic(grid_M)
    res_M <- gextract(paste0('v_',tracks),intervals = grid_M,iterator = grid_M,band = -c(max_dist+max(expand),dist_cutoff-max(expand)))
    res_M <- res_M[,grep('v_hic',colnames(res_M))]
    colnames(res_M) <- vars
    res_M$binding <- c('E-E','P-P')[i]
    df_sub <- rbind(df_sub,res_M)
  }
  df_sub <- df_sub[complete.cases(df_sub),]
  df_m <- melt(df_sub,measure.vars = vars)
  df_m$variable <- factor(df_m$variable,levels=vars)
  p <- ggplot(df_m,aes(x=binding,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =T,width=0.8)
  return(list(plot=p,df=df,df_sub=df_sub))
}


#Figure Functions

Figure7A <- function(files_f,features,cor_f,width=6,height=6,point.size=3,anno.size=5,leg.text.size=10,cols,out_f,recalculate_cor=F){
  df_all <- sapply(files_f,function(x){
    df <- get(load(x))
    return(c(median(df$E14_NSC$intra$v_score,na.rm = T),median(df$E14_IPC$intra$v_score,na.rm = T),median(df$E14_PN$intra$v_score,na.rm = T)))
  })
  mat_names <- gsub("\\.bed.*$", "", basename(files_f))
  mat_names <- gsub("-var.2-", "(var.2)", mat_names)
  mat_names <- gsub("-var.3-", "(var.3)", mat_names)
  colnames(df_all) <- mat_names
  df <- t(df_all)
  colnames(df) <- c('NSC','IPC','PN')
  row.names(df) <- gsub('Ctcf_For','Ctcf_ForRev',row.names(df))
  df_shuff <- df[grep('shuffled',row.names(df)),]
  df[,1] <- df[,1]-mean(df[grep('shuffled',row.names(df)),1])
  df[,2] <- df[,2]-mean(df[grep('shuffled',row.names(df)),2])
  df[,3] <- df[,3]-mean(df[grep('shuffled',row.names(df)),3])
  df_shuff[,1] <- df_shuff[,1]-mean(df_shuff[,1])
  df_shuff[,2] <- df_shuff[,2]-mean(df_shuff[,2])
  df_shuff[,3] <- df_shuff[,3]-mean(df_shuff[,3])
  plot_df_shuff <- data.frame(Max_score=rowMaxs(df_shuff),Max_var=rowVars(df_shuff),row.names = row.names(df_shuff))
  low_score_q <- quantile(plot_df_shuff$Max_score,0.05)
  high_score_q <- quantile(plot_df_shuff$Max_score,0.95)
  low_var_q <- quantile(plot_df_shuff$Max_var,0.05)
  high_var_q <- quantile(plot_df_shuff$Max_var,0.95)  
  if(recalculate_cor){
    chrom_mat <- atac.object@assays$chromvar@data
    rna_mat <- coembed.object@assays$RNA[,colnames(chrom_mat)]
    motifFactors <- sapply(strsplit(row.names(chrom_mat), "\\(|\\:|\\."), function(x) x[[1]])
    indx <- match(motifFactors,row.names(rna_mat))
    rna_mat <- rna_mat[indx[!is.na(indx)],]
    rna_names <- row.names(rna_mat)[match(motifFactors,row.names(rna_mat))]
    o <- data.frame(feature=row.names(chrom_mat),rna=rna_names,row=1:nrow(chrom_mat),cor=NA)
    indx <- o$row[!is.na(o$rna)]
    o$cor[indx] <- rowCorCpp(indx,indx,as.matrix(chrom_mat),as.matrix(rna_mat))
    cor_f <- o
  }
  df <- df[grep('shuffled',row.names(df),invert = T),]
  plot_df <- data.frame(Max_score=rowMaxs(df),Max_var=rowVars(df),row.names = row.names(df))
  plot_df$Cor <- cor_f$cor[match(row.names(plot_df),as.character(cor_f$feature))]
  plot_df$feature <- row.names(plot_df)
  plot_df <- plot_df[complete.cases(plot_df)|plot_df$feature=='Ctcf_ForRev',]
  sign_plotDF <- subset(plot_df,((Max_score>high_score_q)|(Max_score<low_score_q))&(Max_var>high_var_q))
  p <- ggplot(plot_df,aes(x=Max_score,y=Max_var))  + geom_point(color='grey',size = point.size*1.1) + geom_point(aes(fill = Cor),alpha=0.8, pch = I(21),size = point.size,data = sign_plotDF) 
  p <- p + scale_fill_gradientn(colours = cols,name='',breaks=c(-1,1),labels=c('-1','+1'),limits=c(-1,1)) + xlab('Max normalized Hi-C Score') + ylab('Max normalized Hi-C Variance')
  p <- p + theme(legend.direction = "horizontal",legend.justification=c(0,0), legend.position=c(0.1, 0.9),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.2), vjust = 9,size=leg.text.size))+ guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
  if (is.null(features)){
    features <- head(plot_df[order(plot_df$Max_score,decreasing=T),'feature'],20)
  }
  p <- p + ggrepel::geom_text_repel(
    data = plot_df, size = anno.size,box.padding = 0.8, min.segment.length = 0,max.iter = 10000,
    aes(x=Max_score,y=Max_var,color=NULL,label=ifelse(feature%in%features, as.character(feature), "")),force=10)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure7E <- function(mat,features,n_peaks,cor_f,width=6,height=6,point.size=3,anno.size=5,leg.text.size=10,cols,out_f,recalculate_cor=F){
  df <- data.frame(NSC=rep(NA,(ncol(mat)-4)),IPC=rep(NA,(ncol(mat)-4)),PN=rep(NA,(ncol(mat)-4)),row.names = colnames(mat)[1:(ncol(mat)-4)])
  for (i in 1:(ncol(mat)-4)){
    x <- mat[mat[,i]==1,(ncol(mat)-3):ncol(mat)]
    x <- head(x[order(x$maxAccess,decreasing=T),-ncol(x)],n_peaks)
    df[i,] <- as.vector(colMeans(x,na.rm=T))
  }
  mat <- mat[rowMaxs(as.matrix(mat[,1:(ncol(mat)-4)]))>0,(ncol(mat)-3):ncol(mat)]
  mat <- head(mat[order(mat$maxAccess,decreasing=T),],n_peaks*20)
  df_shuff <- data.frame(NSC=rep(NA,1000),IPC=rep(NA,1000),PN=rep(NA,1000))
  for (i in 1:1000){
    x <- mat[sample(row.names(mat),n_peaks),1:3]
    df_shuff[i,] <- as.vector(colMeans(x,na.rm=T))
  }
  df[,1] <- df[,1]-mean(df_shuff[,1])
  df[,2] <- df[,2]-mean(df_shuff[,2])
  df[,3] <- df[,3]-mean(df_shuff[,3])
  df_shuff[,1] <- df_shuff[,1]-mean(df_shuff[,1])
  df_shuff[,2] <- df_shuff[,2]-mean(df_shuff[,2])
  df_shuff[,3] <- df_shuff[,3]-mean(df_shuff[,3])
  plot_df_shuff <- data.frame(Min_meth=rowMins(as.matrix(df_shuff)),Max_var=rowSds(as.matrix(df_shuff)),row.names = row.names(df_shuff))
  low_score_q <- quantile(plot_df_shuff$Min_meth,0.05)
  high_score_q <- quantile(plot_df_shuff$Min_meth,0.95)
  low_var_q <- quantile(plot_df_shuff$Max_var,0.05)
  high_var_q <- quantile(plot_df_shuff$Max_var,0.95)  
  if(recalculate_cor){
    chrom_mat <- atac.object@assays$chromvar@data
    rna_mat <- coembed.object@assays$RNA[,colnames(chrom_mat)]
    motifFactors <- sapply(strsplit(row.names(chrom_mat), "\\(|\\:|\\."), function(x) x[[1]])
    indx <- match(motifFactors,row.names(rna_mat))
    rna_mat <- rna_mat[indx[!is.na(indx)],]
    rna_names <- row.names(rna_mat)[match(motifFactors,row.names(rna_mat))]
    o <- data.frame(feature=row.names(chrom_mat),rna=rna_names,row=1:nrow(chrom_mat),cor=NA)
    indx <- o$row[!is.na(o$rna)]
    o$cor[indx] <- rowCorCpp(indx,indx,as.matrix(chrom_mat),as.matrix(rna_mat))
    cor_f <- o
  }
  plot_df <- data.frame(Min_meth=rowMins(as.matrix(df)),Max_var=rowSds(as.matrix(df)),row.names = row.names(df))
  plot_df$Cor <- cor_f$cor[match(row.names(plot_df),as.character(cor_f$feature))]
  plot_df$feature <- row.names(plot_df)
  plot_df <- plot_df[complete.cases(plot_df)|plot_df$feature=='Ctcf',]
  sign_plotDF <- subset(plot_df,((Min_meth>high_score_q)|(Min_meth<low_score_q))&(Max_var>high_var_q))
  p <- ggplot(plot_df,aes(x=Min_meth,y=Max_var))  + geom_point(color='grey',size = point.size*1.1) + geom_point(aes(fill = Cor),alpha=0.8, pch = I(21),size = point.size,data = sign_plotDF) 
  p <- p + scale_fill_gradientn(colours = cols,name='',breaks=c(-1,1),labels=c('-1','+1'),limits=c(-1,1)) + xlab('Min normalized % CpG Methylation') + ylab('Max normalized CpG Methylation SD')
  p <- p + theme(legend.direction = "horizontal",legend.justification=c(0,0), legend.position=c(0.1, 0.9),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.2), vjust = 9,size=leg.text.size))+ guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
  if (is.null(features)){
    features <- head(plot_df[order(plot_df$Min_meth,decreasing=F),'feature'],20)
  }
  p <- p + ggrepel::geom_text_repel(
    data = plot_df, size = anno.size,box.padding = 0.8, min.segment.length = 0,max.iter = 10000,
    aes(x=Min_meth,y=Max_var,color=NULL,label=ifelse(feature%in%features, as.character(feature), "")),force=10)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure7HI <- function(mat,gene,out_f,features=NULL,anno.size=12,height=5,width=5,cols1,cols2){
  mat <- mat[mat$gene_name==gene]
  mat_conn <- as.matrix(mcols(mat[,c('NSCscore','IPCscore','PNscore')]))
  colnames(mat_conn) <- c('NSC','IPC','PN')
  row.names(mat_conn) <- mat$distance
  mat_conn <- mat_conn[order(mat$distance),]
  
  mat_meth <- as.matrix(mcols(mat[,c('distal.E14_NSC_10x','distal.E14_IPC_10x','distal.E14_PN_10x')]))
  colnames(mat_meth) <- c('NSC','IPC','PN')
  row.names(mat_meth) <- mat$distance
  mat_meth <- mat_meth[order(mat$distance),]
  
  indx <- complete.cases(mat_conn)&complete.cases(mat_meth)
  mat_conn <- mat_conn[indx,]
  mat_meth <- mat_meth[indx,]
  hm1 <- Heatmap(mat_conn,row_title = 'Hi-C Score', name = "Hi-C",cluster_rows = F,column_names_centered = T,show_column_names = T,cluster_columns = F,col=cols1,
                show_row_names = T,row_names_gp = gpar(fontsize=6),column_names_rot = 0,heatmap_legend_param = list(direction='horizontal',title='Hi-C Score',legend_width = unit(width/4, "inch"),at=c(0,50,100)))
  hm2 <- Heatmap(mat_meth,row_title = '% CpG Methylation', name = "% CpG",cluster_rows = F,column_names_centered = T,show_column_names = T,cluster_columns = F,col=cols2,
                 show_row_names = T,row_names_gp = gpar(fontsize=6),column_names_rot = 0,heatmap_legend_param = list(direction='horizontal',title='% CpG Methylation',legend_width = unit(width/4, "inch"),at=c(0,50,100)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  draw(hm1%v%hm2,heatmap_legend_side='bottom',ht_gap=unit(0.8,'inch'))
  dev.off()
}

Figure7K <- function(res,tracks,chip_f,out_f,height,width){
  res_neurog2_mut <- res[res$enh_type=='Neurog2-var2-mut',]
  res_neurog2_wt <- res[res$coord%in%res_neurog2_mut$coord&res$enh_type=='WTposCor',]
  for (track in tracks){
    gvtrack.create(paste0('v_',track),track,'avg')
  }
  df <- merge(res_neurog2_wt[,c('MPRA_name','coord','mad.score_NSC','mad.score_IPC','mad.score_PN')],res_neurog2_mut[,c('MPRA_name','coord','mad.score_NSC','mad.score_IPC','mad.score_PN')],by='coord')
  colnames(df) <- c('coord','MPRA_name_WT','WT_NSC','WT_IPC','WT_PN','MPRA_name_mut','Mut_NSC','Mut_IPC','Mut_PN')
  df_peaks <- stringr::str_split(df$coord, pattern ='_' , n = 3, simplify = TRUE)
  df_peaks <- gintervals(chroms = df_peaks[,1],starts = as.numeric(df_peaks[,2]),ends=as.numeric(df_peaks[,3]))
  res <- gextract(paste0('v_',tracks),df_peaks,iterator = df_peaks)
  res$intervalID <- paste0(res$chrom,'_',res$start,'_',res$end)
  chip <- read.table(chip_f)[,1:3]
  chip <- gintervals(chroms = chip[,1],starts = as.numeric(chip[,2]),ends=as.numeric(chip[,3]))
  test <- gintervals.neighbors(intervals.centers(res),intervals.centers(chip))
  res$Neurog2_chip <- 'unbound'
  res$Neurog2_chip[abs(test$dist)<=100] <- 'bound'
  res <- merge(df,res,by.x='coord',by.y='intervalID')
  res_df <- melt(res,id.vars = c('Neurog2_chip'),measure.vars = c('WT_NSC','WT_IPC','WT_PN','Mut_NSC','Mut_IPC','Mut_PN'))
  res_df$Neurog2_chip <- factor(res_df$Neurog2_chip,levels=c('bound','unbound'))
  res_df$type <- factor(stringr::str_split(res_df$variable, pattern ='_' , n = 2, simplify = TRUE)[,1],levels=c('WT','Mut'))
  res_df$condition <- factor(stringr::str_split(res_df$variable, pattern ='_' , n = 2, simplify = TRUE)[,2],levels=c('NSC','IPC','PN'))
  p <- ggplot(res_df,aes(x=condition,y=value,fill=Neurog2_chip)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =T,width=0.8) + scale_fill_manual(name='',values=c("#F1F1F1","#6B6B6B"))
  p <- p + coord_cartesian(ylim=c(-1,15)) + ylab('MPRA signal') + theme(legend.position = c(0.75,0.92)) + facet_wrap(~type) 
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p+xlab(''))
  dev.off()

}

#Supplementary Figure functions

FigureS9A <- function(results_df,binding='both',legend_title='',out_f,clusters,rna_expr,cutoff,features=NULL,anno.size=12,height=5,width=5,plot_motifs=T){
  results_df <- read.table(results_df,header=T)
  rna_expr <- rna_expr[match(results_df$gene,row.names(rna_expr)),]
  results_df <- cbind(results_df,rna_expr)
  results_df <- results_df[complete.cases(results_df),]
  results_df <- unique(results_df)
  if(binding!='both'){
    res_single <- results_df[results_df$binding=='single',]
    res_double <- results_df[results_df$binding=='double',]
    results_df <- results_df[results_df$binding==binding,]
  }
  results_df$cluster <- gsub('PN2','PN',results_df$cluster)
  res_single$cluster <- gsub('PN2','PN',res_single$cluster)
  colnames(results_df) <- gsub('PN2','PN',colnames(results_df))
  results_df <- results_df[results_df$cluster%in%clusters,]
  indx <- c()
  for (cluster in clusters){
    temp_df <- results_df[results_df$cluster==cluster,intersect(grep(cluster,colnames(results_df)),grep('PValue',colnames(results_df)))]
    min_indx <- which(rowMins(as.matrix(temp_df))<=0.05)
    max_indx <- which(rowMaxs(as.matrix(temp_df))>=0.95)
    temp_df <- temp_df[unique(c(min_indx,max_indx)),]
    temp_indx <- which(row.names(results_df)%in%row.names(temp_df))
    indx <- unique(c(indx,temp_indx))
  }
  results_df <- results_df[indx,]
  results_df <- results_df[rowMaxs(as.matrix(results_df[,clusters]),na.rm=T)>=cutoff,]
  
  mat <- as.matrix(results_df[,c('NSCmean','IPCmean','PNmean')])
  colnames(mat) <- clusters
  row.names(mat) <- results_df$motif
  s_mat <- t(scale(t(mat)))
  bind_mat <- as.matrix(mat-res_single[match(paste0(results_df$motif,':',results_df$cluster),paste0(res_single$motif,':',res_single$cluster)),c('NSCmean','IPCmean','PNmean')])
  row.names(bind_mat) <- row.names(mat)
  ###Prepare complex heatmap
  results_df$cluster <- factor(results_df$cluster,levels=clusters)
  ra1 = rowAnnotation(foo = anno_mark(at = which(paste0(results_df$cluster,':',row.names(s_mat))%in%features),side='right', labels = as.character(results_df$motif[which(paste0(results_df$cluster,':',row.names(s_mat))%in%features)]),labels_gp = gpar(fontsize = anno.size)))
  la1 <- rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                        labels = c("NSC", "IPC", "PN"), 
                                        labels_gp = gpar(col = "white", fontsize = 10)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm <- Heatmap(s_mat, name = "mat",left_annotation = la1,column_names_centered = T,show_row_dend = F,
                row_dend_reorder = TRUE,clustering_method_rows = 'ward.D2',row_split = factor(results_df$cluster,levels=clusters),cluster_columns = F,cluster_row_slices=F,right_annotation = ra1,
                show_row_names = F,row_names_gp = gpar(fontsize=6),column_names_rot = 0,row_title = NULL,heatmap_legend_param = list(direction='horizontal',title=legend_title,legend_width = unit(width/4, "inch")))
  draw(hm,heatmap_legend_side='bottom')
  dev.off()
  
}

FigureS7B <- function(results_df,binding='both',out_f,cluster,xcluster,order_by,rna_expr,point.size=1,cutoff,atac.object=atac.object,features=NULL,anno.size=12,height=5,width=5,plot_motifs=T,ylim=c(-30,30)){
  results_df <- read.table(results_df,header=T)
  rna_expr <- rna_expr[match(results_df$gene,row.names(rna_expr)),]
  results_df <- cbind(results_df,rna_expr)
  results_df <- results_df[complete.cases(results_df),]
  if(binding!='both'){
    results_df <- results_df[results_df$binding==binding,]
  }
  rna_indx <- ifelse(cluster=='PN','PN2',cluster)
  results_clust <- results_df[results_df$cluster==rna_indx,]
  results_clust <- results_clust[results_clust[,rna_indx]>=cutoff,]
  results_clust <- unique(results_clust)
  if (paste0('FC_',cluster,'vs',xcluster)%in%colnames(results_clust)){
    results_clust$FC <- results_clust[,paste0('FC_',cluster,'vs',xcluster)]
  } else {
    results_clust$FC <- results_clust[,paste0('FC_',xcluster,'vs',cluster)]*(-1)
  }
  results_clust$HiCscore <- results_clust[,paste0(cluster,'mean')]
  p_value <- results_clust[,grepl('PValue',colnames(results_clust))&grepl(cluster,colnames(results_clust))&grepl(xcluster,colnames(results_clust))]
  results_clust$FCsig <- 'grey'
  results_clust$FCsig[(p_value<=0.05|p_value>=0.95)&results_clust$FC>mean(results_clust$FC[(p_value>0.05&p_value<0.95)])] <- 'red'
  results_clust$FCsig[(p_value<=0.05|p_value>=0.95)&results_clust$FC<mean(results_clust$FC[(p_value>0.05&p_value<0.95)])] <- 'blue'
  p_value <- results_clust[,grepl('mean',colnames(results_clust))&grepl(cluster,colnames(results_clust))&!grepl('PValue',colnames(results_clust))]
  results_clust$Meansig <- ifelse((p_value<=0.05|p_value>=0.95),'yes','no')
  results_clust <- results_clust[order(results_clust$HiCscore,results_clust$FC,decreasing=T),]
  p <- ggplot(results_clust,aes(x=HiCscore,y=FC,shape=binding)) + geom_point(size=point.size,col='grey',data=results_clust[results_clust$FCsig=='grey',])+ geom_point(size=point.size,col='red',data=results_clust[results_clust$FCsig=='red',]) + geom_point(size=point.size,col='blue',data=results_clust[results_clust$FCsig=='blue',]) + scale_shape_manual(values=c(19,4))  
  if (is.null(features)){
    test <- results_clust[results_clust$HiCscore>=30,]
    features <- unique(c(as.character(head(results_clust[order(results_clust$HiCscore,decreasing=T),'motif'],10)),as.character(head(test[order(test$FC,decreasing=T),'motif'],10))))
  }
  p_labels <- results_clust$motif
  p <- p + geom_text_repel(
    data = results_clust, size = anno.size,box.padding = 0.8,segment.alpha = 0.8, min.segment.length = 0,max.iter = 10000,
    aes(col=binding,label=ifelse(motif%in%features, as.character(p_labels), "")),force=10) + scale_color_manual(values=setNames(c("black","grey50"), levels(results_clust$binding)))
  p <- p + theme(legend.position = 'none') + xlab(cluster) + ylab(paste0(cluster,' - ',xcluster))
  if (plot_motifs){
    motif.names <- unlist(GetMotifObject(atac.object)@motif.names)
    motif.names <-  capitalize(tolower(motif.names))
    p1 <- MotifPlot(
      object = atac.object,
      motifs = names(motif.names)[match(features,motif.names)],
      use.names=T,
      ncol=1
    )
    p1 <- p1 +theme_cowplot(font_size = 12) + theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line=element_blank())
    p <- p+ p1 + plot_layout(ncol=2,widths=c(4,1))
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p+ylim(ylim))
  dev.off()
  return(results_clust)
}



FigureS7C <- function(df_all,sig.mat_all=NULL,plot_density=F,annot_TF=null,point.size=1,clusters=c('NSC','IPC','PN'),domain='intraTAD',features=NULL,anno.size=12,height,width,cols,out_f,boxplot_only=F){
  df_all$cluster <- colnames(df_all[,1:7])[max.col(df_all[,1:7])]
  df_all$labels <- paste0(df_all$gene_name,':',df_all$distance)
  plot_list <- list()
  features1 <- c()
  for (cluster in clusters){ 
    df <- df_all[grep(cluster,df_all$cluster),]
    df <- df[df$cluster!='PN1',]
    if(!is.null(domain)){
      df <- df[df$domain==domain,]
    }
    yscore <- paste0(cluster,'score')
    for (xcluster in clusters[grep(cluster,clusters,invert=T)]){
      xscore <- paste0(xcluster,'score')
      if (plot_density){
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_pointdensity(adjust=5, shape = 20,size=point.size,alpha=1) + scale_color_gradientn(colours = c('black','red','orange'),name='')
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='black') + xlim(0,100) + ylim(0,100)
      } else {
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_point(size=point.size,alpha=1,col='darkblue') 
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='red') + xlim(0,100) + ylim(0,100)
      }
      if(!is.null(sig.mat_all)){
        sig.mat <- sig.mat_all[(sig.mat_all$contrast==paste0(cluster,'-',xcluster))|(sig.mat_all$contrast==paste0(xcluster,'-',cluster)),]
        if (sig.mat$contrast[1]!=paste0(cluster,'-',xcluster)){
          sig.mat$logFC <- sig.mat$logFC*(-1)
        }
        sig.mat_idx <- paste0(sig.mat$GA_anchor,':',sig.mat$DA_anchor)
        df_idx <- paste0(df$peakName,':',df$gene_name)
        df_idx[df$peakCenter>df$gene_start] <- paste0(df$gene_name,':',df$peakName)[df$peakCenter>df$gene_start]
        df$diffHiC_PValue <- sig.mat$PValue[match(df_idx,sig.mat_idx)]
        df$logFC <- sig.mat$logFC[match(df_idx,sig.mat_idx)]
        df$type <- factor('notSig',levels=c('posSig','notSig','negSig'))
        df$type[df$logFC>0&df$diffHiC_PValue<0.05] <- 'posSig'
        df$type[df$logFC<0&df$diffHiC_PValue<0.05] <- 'negSig'
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_point(size=point.size,alpha=1,col='grey',data=df[df$type=='notSig',]) + geom_point(size=point.size,alpha=1,col='blue',data=df[df$type=='negSig',]) + geom_point(size=point.size,alpha=1,col='red',data=df[df$type=='posSig',])           # scale_color_gradientn(colours=tim.colors(12))
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='black') + xlim(0,100) + ylim(0,100)
      }
      if(!is.null(annot_TF)){
        sel_df <- readRDS('results/HiC/posCor_EnhancerMotifs.RDS')
        df$Enh_Motif <- sel_df[match(df$peakName,sel_df$peakName),annot_TF]
        prom_sel_df <- readRDS('results/HiC/posCor_PromMotifs.RDS')
        df$Prom_Motif <- prom_sel_df[match(df$gene_name,prom_sel_df$gene_name),annot_TF]
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_point(size=point.size,alpha=1,col='grey',data=df[df$Enh_Motif==0&df$Prom_Motif==0,]) + geom_point(size=point.size,alpha=1,col='green',data=df[df$Enh_Motif==1,]) + geom_point(size=point.size,alpha=1,col='red',data=df[df$Prom_Motif==1,])  + geom_point(size=point.size,alpha=1,col='yellow',data=df[df$Prom_Motif==1&df$Enh_Motif==1,])          # scale_color_gradientn(colours=tim.colors(12))
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='black') + xlim(0,100) + ylim(0,100)
      }
      df <- df[order(df[,yscore]-df[,xscore],df[,yscore],decreasing=T),]
      if (is.null(features)){
        df1 <- as.data.frame(df %>% group_by(gene_name) %>% dplyr::slice(2))
        df1 <- df1[order(df1[,yscore]-df1[,xscore],df1[,yscore],decreasing=T),]
        df1 <- df1[df1[,yscore]>0&df1[,xscore]>0,]
        features1 <- c(head(df1$labels,5),features1)
      } else {
        features1 <- features
      }
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(0,0), legend.position=c(0.85, 0.05),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_blank())+ guides(color = guide_colorbar(barwidth = 3, barheight = 1))
      p <- p + geom_text_repel(
        data = df, size = anno.size,box.padding = 0.5,segment.alpha = 0.5, min.segment.length = 0,max.iter = 20000,
        aes(color=NULL,label=ifelse(labels%in%features1, as.character(labels), "")),force=10)
      if(!boxplot_only){plot_list[[paste0(cluster,'vs',xcluster)]] <- p}
    }
    df <- melt(df[,c('NSCscore','IPCscore','PNscore')])
    df$variable <- gsub('score','',df$variable)
    df$variable <- factor(df$variable,levels=c('NSC','IPC','PN'))
    p1 <- ggplot(df,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    p1 <- p1 + scale_fill_manual(values=cols) + xlab('') + ylab('Hi-C Score') + theme(legend.position = "none")
    plot_list[[paste0(cluster,'_sign')]] <- p1 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','PN'),c('IPC','PN')),label = "p.format",method='wilcox')
  }
  p <- Reduce("+", plot_list) 
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height*(length(clusters)),width=width*(length(clusters)))
  if(!boxplot_only){
    print(p + theme_cowplot()+ plot_layout(nrow=length(clusters),ncol=length(clusters),widths = c(rep(3,length(clusters)-1),1)))
  } else {
    print(p + theme_cowplot()+ plot_layout(nrow=length(clusters),ncol=1,heights = rep(3,length(clusters))))
  }
  dev.off()
  return(features1)
}


### Plot Figures 

res <- Figure7A(files_f=list.files('/home/hpc/bonev/projects/hic/sc/data/cis_decay/',full.names = T,pattern='score'),
          features=c('Neurod1','Neurog2(var.2)','Pou3f2','Pou3f3','Emx2','Eomes','Sox2','Pax6','Lhx2','Sox4','Sox8','Prrx1','Ctcf_ForRev','Dmrta2','Tgif2','Rfx4','Tfap4','Tead2','Fosb::Jun'),
          cor_f=readRDS('results/integration/motif_rna_correlation.RDS'),
          width=6,height=6,point.size=3,anno.size=5,leg.text.size=10,
          cols=colorRampPalette(c('blue','white','red'))(101),out_f='Figure7A',recalculate_cor=F)
#Figure7B
pdf('plots/figures/Figure7B.pdf',width=6.5,height=3.5)
layout(matrix(c(1:6,7,7),nrow=2,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4,4),respect = T)
for (cell in cells[1:3]){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='Pou3f2.bed',intervals2='Pou3f2.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='Neurog2-var.2-.bed',intervals2='Neurog2-var.2-.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(3,2,3,2.5))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()
#Figure7C
p <- VlnPlot(rna.object, features = c('Pou3f2','Neurog2'),cols=RNA_cluster_colors[c(1,3,5:7)],idents = c('NSC','IPC','PN1','PN2','PN3'),pt.size = 0,combine = F)
pdf('plots/figures/Figure7C.pdf',height=4,width=6)
print(p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('')+ p[[2]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5)) + xlab(''))
dev.off()
#Figure7D
pwm_f <- 'results/scATAC/combined_pwm.RDS'
res <- ep_vs_enh(df=p2glinks$posCor,tracks=score_tracks[1:3],use_chip=T,chip_f="data/beds/NPC_Neurog2.bed",motif_name='Neurog2(var.2)',anno_mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T),vars=c('NSC','IPC','PN'),dist_cutoff=1e4,intraTAD=F)
p <- res$plot + scale_fill_manual(name='',values=ATAC_cluster_colors[c(1,2,4)]) + xlab('') + ylab('Hi-C score') + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.97, 0.97))
pdf('plots/figures/Figure7D.pdf',height=5,width=5.5)
print(p)
dev.off()
#Figure7E
Figure7E(mat=readRDS('results/methylation/motif_Meth.RDS'),n_peaks = 5000,
         features=c('Neurog2(var.2)','Pou3f2','Eomes','Sox2','Tfap2c(var.2)','Ctcf','Nrf1','Insm1','Neurod2','Dmrta2','Tead2','Fosb::jun'),
         cor_f=readRDS('results/integration/motif_rna_correlation.RDS'),
         width=6,height=6,point.size=3,anno.size=5,leg.text.size=10,
         cols=colorRampPalette(c('blue','white','red'))(101),out_f='Figure7E',recalculate_cor=F)
#Figure7F
tracks <- paste0('results/methylation/',c('NSC','IPC','PN'),'_methylation_CpG_10x.bw')
res1 <- getPlotSetArray(tracks=tracks,features="results/HiC/motif_beds_top5000/Nrf1.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
res2 <- getPlotSetArray(tracks=tracks,features="results/HiC/motif_beds_top5000/Neurog2-var.2-.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
res3 <- getPlotSetArray(tracks=tracks,features="results/HiC/motif_beds_top5000/Neurod2.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
pdf('plots/figures/Figure7F_1.pdf',height=5,width=5)
plotAverage(plotset=res1, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
pdf('plots/figures/Figure7F_2.pdf',height=5,width=5)
plotAverage(plotset=res2, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
              colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
pdf('plots/figures/Figure7F_3.pdf',height=5,width=5)
plotAverage(plotset=res3, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
#Figure7G
mat=p2glinks$posCor
plotMisha(targetGene='Eomes',out_f='Figure7G',upstream=4.5e5,downstream=5e4,plot.dgram=T,
          chipNames=c('NSC','IPC','PN','Neurog2'),window_scale=1.8,pointCEX=1,chipRes =10,conditions=score_tracks[1:3],binSize=5e3,radius=1.5e4,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN2','chipseq_RPM.NPC_Neurog2'),methTracksToExtract=c('methylation.E14_NSC_10x','methylation.E14_IPC_10x','methylation.E14_PN_10x'),methColors=ATAC_cluster_colors[c(1,2,4)],
          chipColors=c(ATAC_cluster_colors[c(1,2,4)],'black'),scoreTrackToExtract =score_tracks[1:3],arcIntervals=mat[mat$gene_name=='Eomes'],arcColors=colorRampPalette(archr_colors[[29]]),
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE,MPRA=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,MPRA=0.3,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))
#Figure7H
Figure7HI(mat=p2glinks$posCor,gene='Eomes',out_f='Figure7HI',features=NULL,anno.size=12,height=8,width=2.5,cols1=colorpalette('matlablike2',11),cols2=colorpalette('matlablike2',11))
#Figure7J
res <- vroom::vroom('results/immunoMPRA_res.tsv')
mpra_coords <- res$coord[res$enh_type=='WTposCor']
mpra_coords <- stringr::str_split(paste0(mpra_coords), pattern ='_' , n = 3, simplify = TRUE)
mpra_coords <- gintervals(as.character(mpra_coords[,1]),as.numeric(mpra_coords[,2]),as.numeric(mpra_coords[,3]))
mpra_coords <- intervals.normalize(mpra_coords,270)
mpra_tracks <- c("mpra.WT_NSC","mpra.Neurog2mut_NSC","mpra.WT_IPC","mpra.Neurog2mut_IPC","mpra.WT_PN","mpra.Neurog2mut_PN")
plotMisha(targetGene='chr9,118327840,118328340',outDir='plots/figures/',out_f='Figure7J',upstream=500,downstream=500,
          chipNames=c('','','',''),plot.dgram=F,mpra_coords=mpra_coords,mpra_tracks=mpra_tracks,mpra_colors=colorpalette('matlablike2',11),mpraNames=c('WT','MUT','WT','MUT','WT','MUT'),
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN1','chipseq_RPM.NPC_Neurog2'),
          chipColors=c(ATAC_cluster_colors[1:3],'black'),imgPlotScale=4,
          methTracksToExtract=c('methylation.E14_NSC_10x','methylation.E14_IPC_10x','methylation.E14_PN_10x'),methColors=ATAC_cluster_colors[c(1,2,4)],methNames=c('','',''),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=TRUE,MPRA=TRUE,axis=FALSE,scATAC=FALSE, genes=FALSE,ideogram=FALSE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.6,MPRA=0.2,meth=0.6, domains=0.15, genes=0.7, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))
pdf('plots/figures/Figure7J_scalebar.pdf',height=4,width=2)
par(mar=c(1,4,1,4))
image.scale(as.matrix(NA),zlim=c(-0.8,10,12), col=colorpalette('matlablike2',100),axis.pos=4,adj=1,cex.axis = 1.5)
dev.off()
#Figure7K
Figure7K(res=vroom::vroom('results/immunoMPRA_res.tsv'),tracks=c('chipseq_RPM.NPC_Neurog2'),chip_f='data/beds/NPC_Neurog2.bed',out_f='Figure7K',height=4,width=5)

####Plot Supplementary Figure #####
FigureS9A(results_df='results/HiC/posCor_motif_HiCscore.tsv',binding='double',
         out_f='FigureS9A',clusters=c('NSC','IPC','PN'),cutoff=6,legend_title='Scaled Hi-C score',
         features=c('NSC:Emx2','NSC:Pax6','NSC:Lhx2','NSC:Sox2','NSC:Prrx1','NPC:Pou3f2','IPC:Insm1','IPC:Neurog2(var.2)','IPC:Eomes','IPC:Neurod1','IPC:Pou3f2','PN:Tbr1','PN:Foxp2','PN:Neurod2','PN:Neurod1'),
         anno.size=10,height=8,width=6,plot_motifs=F,
         rna_expr=uf_rna$groupMat)
#FigureS9B
pdf('plots/figures/FigureS9B.pdf',width=6.5,height=5)
layout(matrix(c(1:9,10,10,10),nrow=3,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4,4,4),respect = T)
for (cell in cells[1:3]){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='NPC_Pax6.narrowPeak',intervals2='NPC_Pax6.narrowPeak',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='NPC_Neurog2.narrowPeak',intervals2='NPC_Neurog2.narrowPeak',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='NPC_Eomes_top3k.bed',intervals2='NPC_Eomes_top3k.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(7,2,7,2.5))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()
#FigureS9C
tracks <- paste0('results/scATAC_clusters/',c('NSC','IPC','PN2'),'.bw')
res1 <- getPlotSetArray(tracks=tracks,features="data/beds/NPC_Pax6.narrowPeak",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
res2 <- getPlotSetArray(tracks=tracks,features="data/beds/NPC_Neurog2.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
res3 <- getPlotSetArray(tracks=tracks,features="data/beds/NPC_Eomes_top3k.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
pdf('plots/figures/FigureS9C_1.pdf',height=4,width=4)
plotAverage(plotset=res1, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = 'Average Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
pdf('plots/figures/FigureS9C_2.pdf',height=4,width=4)
plotAverage(plotset=res2, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = 'Average Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
pdf('plots/figures/FigureS9C_3.pdf',height=4,width=4)
plotAverage(plotset=res3, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = 'Average Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
#FigureS9D
pdf('plots/figures/FigureS9D.pdf',width=6.5,height=3.5)
layout(matrix(c(1:6,7,7),nrow=2,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4,4),respect = T)
for (cell in cells[1:3]){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='Tead2.bed',intervals2='Tead2.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='Fosb::jun.bed',intervals2='Fosb::jun.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(3,2,3,2.5))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()
#FigureS9E
tracks <- paste0('results/methylation/',c('NSC','IPC','PN'),'_methylation_CpG_10x.bw')
res1 <- getPlotSetArray(tracks=tracks,features="results/HiC/motif_beds_top5000/Tead2.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
res2 <- getPlotSetArray(tracks=tracks,features="results/HiC/motif_beds_top5000/Fosb::jun.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
pdf('plots/figures/FigureS9E_1.pdf',height=5,width=5)
plotAverage(plotset=res1, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
pdf('plots/figures/FigureS9E_2.pdf',height=5,width=5)
plotAverage(plotset=res2, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
#FigureS9F
FigureS9A(results_df='results/methylation/posCor_motif_HiCscore.tsv',binding='double',
         out_f='FigureS9F',clusters=c('NSC','IPC','PN'),cutoff=6,legend_title='Scaled DNA Methylation',
         features=c('NSC:Emx2','NSC:Sox2','NSC:Pax6','IPC:Neurog2(var.2)','IPC:Eomes','IPC:Insm1','PN:Tbr1','PN:Foxp2','PN:Neurod2','PN:Neurod1'),
         anno.size=10,height=8,width=6,plot_motifs=F,
         rna_expr=uf_rna$groupMat)
#FigureS9G

tracks <- paste0('results/methylation/',c('NSC','IPC','PN'),'_methylation_CpG_10x.bw')
res1 <- getPlotSetArray(tracks=tracks,features="data/beds/NPC_Neurog2_top3k_noTss.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
res2 <- getPlotSetArray(tracks=tracks,features="data/beds/CN_NeuroD2_top3k_noTss.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20)
pdf('plots/figures/FigureS9G_1.pdf',height=5,width=5)
plotAverage(plotset=res1, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()
pdf('plots/figures/FigureS9G_2.pdf',height=5,width=5)
plotAverage(plotset=res2, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = ATAC_cluster_colors[c(1,2,4)], pointsize = 12)
axis(side = 1,labels=c('-2000',0,'+2000'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()

plotMisha(targetGene='Rnd2',out_f='FigureS9H',upstream=2e4,downstream=17e4,radius=5e3,
          chipNames=c('NSC','IPC','PN','Neurog2'),window_scale=1.2,pointCEX=0.5,chipRes =10,pointCEX=3,conditions=score_tracks[1:3],binSize=5e3,radius=2e3,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN2','chipseq_RPM.NPC_Neurog2'),methTracksToExtract=c('methylation.E14_NSC_10x','methylation.E14_IPC_10x','methylation.E14_PN_10x'),methColors=ATAC_cluster_colors[c(1,2,4)],
          chipColors=c(ATAC_cluster_colors[c(1,2,4)],'black'),scoreTrackToExtract =score_tracks[1:3],arcIntervals=mat[mat$gene_name=='Rnd2'],arcColors=colorRampPalette(archr_colors[[29]]),
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE,meth=FALSE,MPRA=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))
