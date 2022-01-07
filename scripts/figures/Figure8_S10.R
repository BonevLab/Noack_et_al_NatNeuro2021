library(Signac)
library(Seurat)
library(edgeR)
library(JASPAR2020)
library(Matrix)
library(TFBSTools)
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(dplyr)
library(Hmisc)
library(Matrix)
library(motifmatchr)
library(BSgenome)
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
library(doParallel)
registerDoParallel(8)

source('scripts/figures/config.R')
source('scripts/aux_functions.R')
source('scripts/figures/plot_functions.R')

source('scripts/hic/config.R')
source('scripts/hic/scripts/main_functions.R')
source('scripts/hic/scripts/aux_functions.R')
source('scripts/hic/scripts/plot_functions.R')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))

iue_cols <- c("#008B45FF","#631879FF")
iue_labels <- c('GFP','Neurog2')
iue_rep_labels <- paste0(rep(iue_labels,each=2),'_',c('rep1','rep2','rep1','rep2')) 

access_tracks <- c('data/IUE/GFP_IUE24h_GpC_cov10x.bw','data/IUE/NGN2_IUE24h_GpC_cov10x.bw')
access_misha_tracks <- c("methylation.GFP_IUE24h_GpC_cov10x","methylation.NGN2_IUE24h_GpC_cov10x")

meth_tracks <- c('data/IUE/GFP_IUE24h_CpG_cov10x.bw','data/IUE/NGN2_IUE24h_CpG_cov10x.bw')
meth_misha_tracks <- c("methylation.GFP_IUE24h_CpG_cov10x","methylation.NGN2_IUE24h_CpG_cov10x")

object <- readRDS('data/merged_scRNA_unfiltered_IDs.RDS')
p2glinks <- readRDS('results/IUE24h_P2G-Links.RDS')

plot_pair <- function(borders,tss,extra_tracks,ins_tracks,features=NULL,width,height,plot_what='ins',plot_cells=c('NSC','PN'),point.size=1,anno.size=4,out_f,mode='scatter',cols,ylab_n){
  for (track in c(ins_tracks,extra_tracks)){
    gvtrack.create(paste0('v_',track),track,'avg')
  }
  if(!is.null(extra_tracks)){
    df_borders <- gextract(c(paste0('v_',ins_tracks,'*(-1)'),paste0('v_',extra_tracks)),intervals = borders,iterator=borders)
  } else {
    df_borders <- gextract(c(paste0('v_',ins_tracks,'*(-1)')),intervals = borders,iterator=borders)
  }
  
  df_borders$labels <- paste0(df_borders$chrom,':',df_borders$start,'-',df_borders$end)
  df_borders$genes <- gintervals.neighbors(df_borders,tss)$geneName
  if(mode=='scatter'){
    xscore <- colnames(df_borders)[intersect(grep(plot_cells[1],colnames(df_borders)),grep(plot_what,colnames(df_borders)))]
    yscore <- colnames(df_borders)[intersect(grep(plot_cells[2],colnames(df_borders)),grep(plot_what,colnames(df_borders)))]
    p <- ggplot(df_borders, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_point(size=point.size,alpha=1,col='darkblue') 
    p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='red') + xlab(plot_cells[1]) + ylab(plot_cells[2])
    p <- p + geom_text_repel(
      data = df_borders, size = anno.size,box.padding = 0.8, min.segment.length = 0,max.iter = 10000,
      aes(label=ifelse(genes%in%features, as.character(genes), "")),force=10)
    pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
    print(p)
    dev.off()
  } else{
    colnames(df_borders)[4:(3+length(plot_cells))] <- plot_cells
    df_borders <- df_borders[,c('labels',plot_cells)]
    df <- melt(df_borders,id.vars = 'labels')
    p <- ggplot(df,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
    p <- p + scale_fill_manual(values=cols) + xlab('') + ylab(ylab_n) + theme(legend.position = "none")
    return(p)
  }
}

Figure8C <- function(file_f,conditions,reps,which_genes,cols){
  mat <- read.csv(file_f,header=T,row.names = 'Gene')
  df <- mat[which_genes,]
  df$Gene <- factor(which_genes,levels=which_genes)
  df <- melt(df,varnames = c('Gene'),value.name = 'FC')
  df$rep <- as.factor(gsub(".*_","",df$variable))
  df$Condition <- factor(rep(conditions,each=4*length(which_genes)),levels=conditions)
  df_summ <- data_summary(df, varname="FC",groupnames=c('Condition','Gene'))
  p <- ggplot(df_summ, aes(x=Gene, y=FC, fill=Condition)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black", width=.8,position=position_dodge(.85)) + xlab('') + ylab('Relative Fold Change') +
    geom_errorbar(aes(ymin=FC-sd, ymax=FC+sd), width=.2,position=position_dodge(.85))
  p <- p + theme(legend.position='none') + geom_point(data = df,aes(x=Gene, y=FC, fill=Condition),position = position_jitterdodge(dodge.width =.85,jitter.width = 0.25,seed = 42),size=2) 
  p <- p + stat_compare_means(data=df,aes(group=Condition),label = "p.format",method = 't.test',paired = T) 
  return(p)
}

Figure7A <- function(object,clusters,count_matrix,cols=colorRampPalette(rev(colorpalette('ylorrd')))(100)){
  object_cpm <- as.data.frame(matrix(NA,nrow=nrow(object),ncol=length(clusters)*2))
  colnames(object_cpm) <- paste0(rep(clusters,each=2),c('_rep1','_rep2'))
  row.names(object_cpm) <- row.names(object)
  object_cts <- object_cpm
  i=1
  for (s in clusters){
    for (k in unique(object$orig.ident)){
      object_sub <- object[,intersect(grep(s,Idents(object)),grep(k,object$orig.ident))]
      object_cpm[,i] <- cpm(rowSums(object_sub@assays$RNA@counts),log=T,prior.count = 1)
      object_cts[,i] <- rowSums(object_sub@assays$RNA@counts,na.rm=T)
      
      i <- i+1
    }
  }
  cts <- fread(count_matrix)
  colnames(cts)[1] <- 'gene_name'
  cts <- as.data.frame(cts)
  for (i in 2:ncol(cts)){
    cts[,i] <- cpm(cts[,i],log=T,prior.count = 1)
  }
  cts <- cts[,c(1,2,4,3,5)]
  df <- merge(object_cpm,cts,by.x='row.names',by.y='gene_name')
  res <- cor(df[,-1],method = 'spearman')
  cor_dist <- dist(1-res)
  res <- round(res, 2)
  annotation_df <- data.frame(condition=factor(c(rep(clusters,each=2),rep(iue_labels,each=2)),levels=c(clusters,iue_labels)))
  colnames(annotation_df) <- c('condition')
  row.names(annotation_df) <- colnames(res)
  
  print(pheatmap(res,color=cols,breaks=seq(0.5,1,length=101),clustering_distance_rows=cor_dist,clustering_distance_cols=cor_dist,show_colnames = F,
                 angle_col = 0,display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE,legend_breaks=c(0.5,0.75,1)))
  
}

#Figure functions#

Figure8C <- function(file_f,conditions,reps,which_genes,cols){
  mat <- read.csv(file_f,header=T,row.names = 'Gene')
  df <- mat[which_genes,]
  df$Gene <- factor(which_genes,levels=which_genes)
  df <- melt(df,varnames = c('Gene'),value.name = 'FC')
  df$rep <- as.factor(gsub(".*_","",df$variable))
  df$Condition <- factor(rep(conditions,each=4*length(which_genes)),levels=conditions)
  df_summ <- data_summary(df, varname="FC",groupnames=c('Condition','Gene'))
  p <- ggplot(df_summ, aes(x=Gene, y=FC, fill=Condition)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black", width=.8,position=position_dodge(.85)) + xlab('') + ylab('Relative Fold Change') +
    geom_errorbar(aes(ymin=FC-sd, ymax=FC+sd), width=.2,position=position_dodge(.85))
  p <- p + theme(legend.position='none') + geom_point(data = df,aes(x=Gene, y=FC, fill=Condition),position = position_jitterdodge(dodge.width =.85,jitter.width = 0.25,seed = 42),size=2) 
  p <- p + stat_compare_means(data=df,aes(group=Condition),label = "p.format",method = 't.test',paired = T) 
  return(p)
}

Figure8J <- function(score_f,obs_f=NULL,min_dist,max_dist,TAD,labels,cols,p_stats,add_pvalue=T,ylim=NULL){
  if(is.null(obs_f)){
    df <- get(load(score_f))
    df_e <- cbind(df[[4]]$intra$v_score,df[[5]]$intra$v_score,df[[1]]$intra$v_score,df[[2]]$intra$v_score,df[[3]]$intra$v_score)
    colnames(df_e) <- labels
    #colnames(df_e) <- names(df)
    df_o <- melt(df_e)
    df_o$Var2 <- factor(df_o$Var2,levels=labels)
    p <- ggplot(df_o,aes(x=Var2,y=value,fill=Var2)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
    p <- p + scale_fill_manual(values=cols) + xlab('') + ylab('Hi-C score') + theme(legend.position = "none") + ylim(c(0,100)) + ggtitle('')
    if(add_pvalue){
      p <- p + stat_compare_means(comparisons = list(p_stats),label = "p.format",method = 'wilcox',paired = T,label.y = 100,tip.length = c(0.03,0.01)) 
    }
    return(p)
  } else{
    df <- get(load(obs_f))
    dist = 10^seq(4,8,by=0.1)
    min_ind <- findInterval(min_dist,dist)+1
    max_ind <- findInterval(max_dist,dist)
    if(TAD=='intra'){
      df_obs <- cbind(df[[5]][[1]]$obs_intra,df[[5]][[2]]$obs_intra,df[[5]][[3]]$obs_intra,df[[6]][[1]]$obs_intra,df[[6]][[2]]$obs_intra,df[[6]][[3]]$obs_intra,df[[1]][[1]]$obs_intra,df[[1]][[2]]$obs_intra,df[[1]][[3]]$obs_intra,df[[3]][[1]]$obs_intra,df[[3]][[2]]$obs_intra,df[[3]][[3]]$obs_intra,df[[4]][[1]]$obs_intra,df[[4]][[2]]$obs_intra,df[[4]][[3]]$obs_intra)
      df_exp <- cbind(df[[5]][[1]]$exp_intra,df[[5]][[2]]$exp_intra,df[[5]][[3]]$exp_intra,df[[6]][[1]]$exp_intra,df[[6]][[2]]$exp_intra,df[[6]][[3]]$exp_intra,df[[1]][[1]]$exp_intra,df[[1]][[2]]$exp_intra,df[[1]][[3]]$exp_intra,df[[3]][[1]]$exp_intra,df[[3]][[2]]$exp_intra,df[[3]][[3]]$exp_intra,df[[4]][[1]]$exp_intra,df[[4]][[2]]$exp_intra,df[[4]][[3]]$exp_intra)
    } else if(TAD=='inter') {
      df_obs <- cbind(df[[5]][[1]]$obs_inter,df[[5]][[2]]$obs_inter,df[[5]][[3]]$obs_inter,df[[6]][[1]]$obs_inter,df[[6]][[2]]$obs_inter,df[[6]][[3]]$obs_inter,df[[1]][[1]]$obs_inter,df[[1]][[2]]$obs_inter,df[[1]][[3]]$obs_inter,df[[3]][[1]]$obs_inter,df[[3]][[2]]$obs_inter,df[[3]][[3]]$obs_inter,df[[4]][[1]]$obs_inter,df[[4]][[2]]$obs_inter,df[[4]][[3]]$obs_inter)
      df_exp <- cbind(df[[5]][[1]]$exp_inter,df[[5]][[2]]$exp_inter,df[[5]][[3]]$exp_inter,df[[6]][[1]]$exp_inter,df[[6]][[2]]$exp_inter,df[[6]][[3]]$exp_inter,df[[1]][[1]]$exp_inter,df[[1]][[2]]$exp_inter,df[[1]][[3]]$exp_inter,df[[3]][[1]]$exp_inter,df[[3]][[2]]$exp_inter,df[[3]][[3]]$exp_inter,df[[4]][[1]]$exp_inter,df[[4]][[2]]$exp_inter,df[[4]][[3]]$exp_inter)
    } else{
      df_obs1 <- cbind(df[[5]][[1]]$obs_intra,df[[5]][[2]]$obs_intra,df[[5]][[3]]$obs_intra,df[[6]][[1]]$obs_intra,df[[6]][[2]]$obs_intra,df[[6]][[3]]$obs_intra,df[[1]][[1]]$obs_intra,df[[1]][[2]]$obs_intra,df[[1]][[3]]$obs_intra,df[[3]][[1]]$obs_intra,df[[3]][[2]]$obs_intra,df[[3]][[3]]$obs_intra,df[[4]][[1]]$obs_intra,df[[4]][[2]]$obs_intra,df[[4]][[3]]$obs_intra)
      df_exp1 <- cbind(df[[5]][[1]]$exp_intra,df[[5]][[2]]$exp_intra,df[[5]][[3]]$exp_intra,df[[6]][[1]]$exp_intra,df[[6]][[2]]$exp_intra,df[[6]][[3]]$exp_intra,df[[1]][[1]]$exp_intra,df[[1]][[2]]$exp_intra,df[[1]][[3]]$exp_intra,df[[3]][[1]]$exp_intra,df[[3]][[2]]$exp_intra,df[[3]][[3]]$exp_intra,df[[4]][[1]]$exp_intra,df[[4]][[2]]$exp_intra,df[[4]][[3]]$exp_intra)
      df_obs2 <- cbind(df[[5]][[1]]$obs_inter,df[[5]][[2]]$obs_inter,df[[5]][[3]]$obs_inter,df[[6]][[1]]$obs_inter,df[[6]][[2]]$obs_inter,df[[6]][[3]]$obs_inter,df[[1]][[1]]$obs_inter,df[[1]][[2]]$obs_inter,df[[1]][[3]]$obs_inter,df[[3]][[1]]$obs_inter,df[[3]][[2]]$obs_inter,df[[3]][[3]]$obs_inter,df[[4]][[1]]$obs_inter,df[[4]][[2]]$obs_inter,df[[4]][[3]]$obs_inter)
      df_exp2 <- cbind(df[[5]][[1]]$exp_inter,df[[5]][[2]]$exp_inter,df[[5]][[3]]$exp_inter,df[[6]][[1]]$exp_inter,df[[6]][[2]]$exp_inter,df[[6]][[3]]$exp_inter,df[[1]][[1]]$exp_inter,df[[1]][[2]]$exp_inter,df[[1]][[3]]$exp_inter,df[[3]][[1]]$exp_inter,df[[3]][[2]]$exp_inter,df[[3]][[3]]$exp_inter,df[[4]][[1]]$exp_inter,df[[4]][[2]]$exp_inter,df[[4]][[3]]$exp_inter)
      df_obs <- df_obs1+df_obs2
      df_exp <- df_exp1+df_exp2
    }
    res <- data.frame(condition=rep(labels,each=3),rep=paste0(rep('rep_',15),1:3),obs_exp=log2(colSums(df_obs[min_ind:max_ind,],na.rm=T)/colSums(df_exp[min_ind:max_ind,],na.rm=T)))
    res$condition <- factor(res$condition,levels=labels)
    df_summ <- data_summary(res, varname="obs_exp",groupnames=c('condition'))
    p <- ggplot(df_summ, aes(x=condition, y=obs_exp, fill=condition)) + scale_fill_manual(values = cols) +
      geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('log2(obs/exp)') +
      geom_errorbar(aes(ymin=obs_exp-sd, ymax=obs_exp+sd), width=.2,position=position_dodge(.9))
    p <- p + geom_point(data = res,aes(x=condition, y=obs_exp, fill=condition,shape=rep),position = position_jitterdodge(dodge.width =.2,jitter.width = 0.1,seed = 42),size=2) 
    if(add_pvalue){
      p <- p + stat_compare_means(data = res,comparisons = list(p_stats),label = "p.format",method = 't.test',paired = T,label.y = df_summ$obs_exp[2]*1.05,tip.length = c(0.05,0.01)) 
    }
  }
}

Figure8K <- function(df,use_chip,chip_f,anno_mat,clusters,vars,inv_size=500,intraTAD=FALSE,dist_cutoff=5000){
  if(use_chip){
    chip <- read.table(chip_f)
    colnames(chip)[1:3] <- c('chrom','start','end')
    chip <- chip[chip$chrom%in%gintervals.all()$chrom,]
    chip <- intervals.normalize(chip,500)
    enh_ranges <- StringToGRanges(df$peakName, sep = c("_", "_"))
    prom_ranges <- GRanges(seqnames=df$gene_chr,IRanges(start=df$gene_start-5000,end=df$gene_start+5000))
    chip_features <- makeGRangesFromDataFrame(chip)
    prom_overlaps <- GenomicRanges::countOverlaps(prom_ranges,chip_features)
    df$prom_Neurog2_ChIP <- prom_overlaps
    overlaps <- GenomicRanges::countOverlaps(enh_ranges,chip_features)
    df$Neurog2_ChIP <- overlaps
  } else {
    motif_list <- list()
    pwm <- readRDS(chip_f)
    motif.matrix <- CreateMotifMatrix(features = StringToGRanges(df$peakName, sep = c("_", "_")),pwm = pwm,genome = 'mm10',sep = c("-", "-"))
    motif.names <- name(pwm)
    colnames(motif.matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])))
    df$Neurog2_motif <- as.vector(motif.matrix[,"Neurog2(var.2)"])
    motif_list[['Enhancer']] <- motif.matrix
    
    prom_ranges <- GRanges(seqnames=df$gene_chr,IRanges(start=df$gene_start-2000,end=df$gene_start+500))
    start(prom_ranges)[df$gene_strand=='-'] <- df$gene_start[df$gene_strand=='-']-500
    end(prom_ranges)[df$gene_strand=='-'] <- df$gene_start[df$gene_strand=='-']+2000
    motif.matrix <- CreateMotifMatrix(features = prom_ranges,pwm = pwm,genome = 'mm10',sep = c("-", "-"))
    motif.names <- name(pwm)
    colnames(motif.matrix) <- capitalize(tolower(as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])))
    df$prom_Neurog2_motif <- as.vector(motif.matrix[,"Neurog2(var.2)"])
    motif_list[['Prom']] <- motif.matrix
  }
  anno_mat$cluster <- colnames(anno_mat[,1:7])[max.col(anno_mat[,1:7])]
  anno_mat$labels <- paste0(anno_mat$gene_name,':',anno_mat$distance)
  df$labels <- paste0(df$gene_name,':',df$distance)
  df$cluster <- anno_mat$cluster[match(df$labels,anno_mat$labels)]
  if(use_chip){
    df <- as.data.frame(mcols(df))
    df$Neurog2_binding <- 'none'
    df$Neurog2_binding[df$Neurog2_ChIP>=1&df$prom_Neurog2_ChIP>=1] <- 'both'
    df$Neurog2_binding[(df$Neurog2_ChIP>=1&df$prom_Neurog2_ChIP==0)] <- 'distal'
    df$Neurog2_binding[(df$Neurog2_ChIP==0&df$prom_Neurog2_ChIP>=1)] <- 'promoter'
    df$Neurog2_binding <- factor(df$Neurog2_binding,levels=c('none','distal','promoter','both'))
  } else {
    df <- as.data.frame(mcols(df))
    df$Neurog2_binding <- 'no'
    df$Neurog2_binding[df$Neurog2_motif==1&df$prom_Neurog2_motif==1] <- 'both'
    df$Neurog2_binding[(df$Neurog2_motif==1&df$prom_Neurog2_motif==0)] <- 'enh'
    df$Neurog2_binding[(df$Neurog2_motif==0&df$prom_Neurog2_motif==1)] <- 'prom'
    df$Neurog2_binding <- factor(df$Neurog2_binding,levels=c('no','enh','prom','both'))
  }
  df_sub <- df[abs(df$distance)>dist_cutoff,]
  colnames(df_sub) <- gsub('GFP_score','GFP',colnames(df_sub))
  colnames(df_sub) <- gsub('NGN2_score','Neurog2',colnames(df_sub))
  df_sub <- df_sub[complete.cases(df_sub[,vars]),]
  df_m <- melt(df_sub,measure.vars = vars)
  if(intraTAD){
    df_m <- df_m[df_m$domain=='intraTAD',]
  }
  df_m$variable <- factor(df_m$variable,levels=vars)
  p <- ggplot(df_m,aes(x=Neurog2_binding,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =T,width=0.8)
  p1 <- p + stat_compare_means(aes(group=variable),label = "p.format",method = 'wilcox',paired = T)
  return(list(plot=p,p1=p1,df=df,df_sub=df_sub))
}

Figure8L <- function(fc,df,cols,vars='RNA_fc',id.vars='Neurog2_binding'){
  df$RNA_fc <- 2^fc[match(df$gene_name,row.names(fc)),'log2FoldChange']
  df$CpG_fc <- df$distal.NGN2_IUE24h_CpG_cov10x-df$distal.GFP_IUE24h_CpG_cov10x
  df$GpC_fc <- df$distal.NGN2_IUE24h_GpC_cov10x-df$distal.GFP_IUE24h_GpC_cov10x
  df$HiC_fc <- df$Neurog2-df$GFP
  df_m <- melt(df,measure.vars = vars,id.vars = id.vars)
  df_m$variable <- factor(df_m$variable,levels=vars)
  p <- ggplot(df_m,aes(x=Neurog2_binding,y=value,fill=Neurog2_binding)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
  return(p)
}

#Supplementary Figure Functions #

FigureS10B <- function(bin_res,out_f,heigh,width){
  df<-read.table(bin_res, header=T, sep='\t')
  stat.test <- df %>%
    group_by(Zone) %>%
    t_test(Value ~ Condition) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = "Zone", dodge = 0.8)
  p<-ggbarplot(df, x = "Zone", y = "Value",
               fill= "Condition", add.params = list(shape = "Condition"), palette = c("#008B45FF","#631879FF"), xlab=FALSE, ylab='GFP+ cells in bins (in %)',
               position = position_dodge(0.8), add = c("mean_sd",'jitter')) + stat_pvalue_manual(stat.test,  label = "p", tip.length = 0.01)+ labs(fill = "", shape="")
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS10C <- function(expr_df,cluster,out_f,width,height){
  marker_genes <- FindMarkers(object,ident.1 = cluster,ident.2 = NULL,only.pos = T)
  expr_df <- log10(expr_df[row.names(expr_df)%in%row.names(marker_genes),]+1)  #log10 FPKM +1 transformation
  df_summ <- ddply(res$df[,c('gene_name','Neurog2_binding')],.(gene_name),function(x){return(table(x$Neurog2_binding))})
  df_summ$bound <- 'unbound'
  df_summ$bound[rowSums(df_summ[,3:5])>=1] <- 'bound'
  expr_df$binding <- df_summ$bound[match(row.names(expr_df),df_summ$gene_name)]
  expr_df <- expr_df[complete.cases(expr_df),]
  expr_df$gene_name <- row.names(expr_df)
  res_df <- melt(expr_df,id.vars = c('gene_name','binding'))
  res_df$binding <- factor(res_df$binding,levels=c('bound','unbound'))
  res_df$variable <- factor(res_df$variable,levels=c('GFP','Neurog2','NSC','IPC','PN1'))
  res_df$xc <- 1
  res_df$xc[res_df$variable=='PN1'] <- 5.5
  res_df$xc[res_df$variable=='IPC'] <- 4.5
  res_df$xc[res_df$variable=='NSC'] <- 3.5
  res_df$xc[res_df$variable=='Neurog2'] <- 2
  p <- ggplot(res_df,aes(x=xc,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =T,width=0.8) + scale_fill_manual(name='',values=c(iue_cols,rep('#aaa9ad',3)))
  p <- p + ylab('log10(FPKM+1)') + xlab('') + scale_x_discrete(labels='') + theme(legend.position = 'none') + facet_wrap(~binding) 
  message('Bound p-value:',wilcox.test(expr_df$NSC[expr_df$binding=='bound'],expr_df$IPC[expr_df$binding=='bound'])$p.value)
  message('Unbound p-value:',wilcox.test(expr_df$NSC[expr_df$binding=='unbound'],expr_df$IPC[expr_df$binding=='unbound'])$p.value)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

#####Plot Figures####
#Figure8C
p1 <- Figure8C(file_f="results/IUE_qPCR.tsv",conditions=iue_labels,reps=4,which_genes="Neurog2",cols=iue_cols)
p2 <- Figure8C(file_f="results/IUE_qPCR.tsv",conditions=iue_labels,reps=4,which_genes=c("Pax6","Eomes"),cols=iue_cols)
pdf('plots/figures/Figure8C.pdf',height=5,width=2,useDingbats=FALSE)
print(p1+p2)
dev.off()
#Figure8D
res <- getPlotSetArray(tracks=access_tracks,features="data/beds/NPC_Neurog2.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20)
pdf('plots/figures/Figure8D.pdf',height=6,width=6)
plotAverage(plotset=res, labels = iue_labels, xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = 'Average Accessibility (% GpC Methylation)',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = iue_cols, pointsize = 12)
axis(side = 1,labels=c('-1000','-500',0,'+500','+1000'),at=c(-1000,-500,0,500,1000),pos=7.5,tick = T)
dev.off()
#Figure8E
res <- extract_lin(regions="data/beds/NPC_Neurog2.bed",window=200,tracks=access_misha_tracks)
p <- misha_lin_boxplot(df=res,condition_names=iue_labels,cols=iue_cols,ylab_n='Accessibility (%GpC Methylation)')
p <- p + coord_cartesian(ylim=c(0,80)) + scale_y_continuous(breaks = seq(0, 80, by = 20)) + stat_compare_means(comparisons = list(iue_labels),label = "p.format",method='wilcox',paired = T,label.y = 77,tip.length = c(0.03,0.01))
pdf('plots/figures/Figure8E.pdf',height=5,width=4)
print(p)
dev.off()
#Figure8F
res <- getPlotSetArray(tracks=meth_tracks,features='data/beds/NPC_Neurog2.bed',refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 50)
pdf('plots/figures/Figure8F.pdf',height=6,width=6)
plotAverage(plotset=res, labels = iue_labels, xlim = NULL,
            ylim = c(5,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = iue_cols, pointsize = 12)
axis(side = 1,labels=c('-1000','-500',0,'+500','+1000'),at=c(-1000,-500,0,500,1000),pos=5,tick = T)
dev.off()
#Figure8G
res <- extract_lin(regions="data/beds/NPC_Neurog2.bed",window=200,tracks=meth_misha_tracks)
p <- misha_lin_boxplot(df=res[complete.cases(res),],condition_names=iue_labels,cols=iue_cols,ylab_n='%CpG Methylation')
p <- p + coord_cartesian(ylim=c(0,55)) + scale_y_continuous(breaks = seq(0, 50, by = 10)) + stat_compare_means(comparisons = list(iue_labels),label = "p.format",method='wilcox',paired = T,label.y = 52,tip.length = c(0.01,0.03))
pdf('plots/figures/Figure8G.pdf',height=5,width=3)
print(p)
dev.off()
#Figure8H
res <- extract_lin(regions='results/beds/NPC_Neurog2_GFP_meth20.bed',window=200,tracks=meth_misha_tracks)
p <- misha_lin_boxplot(df=res[complete.cases(res),],condition_names=iue_labels,cols=iue_cols,ylab_n='%CpG Methylation')
p <- p  + stat_compare_means(comparisons = list(iue_labels),label = "p.format",method='wilcox',paired = T,label.y = 100,tip.length = c(0.01,0.03))
pdf('plots/figures/Figure8H.pdf',height=5,width=3)
print(p)
dev.off()
#Figure8I
pdf('plots/figures/Figure8I.pdf',width=4.5,height=4)
layout(matrix(c(1:3),nrow=1,ncol=3,byrow=F),widths = c(4,4,2),heights=c(4),respect = T)
for (cell in c('GFP_IUE24h','NGN2_IUE24h')){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='NPC_Neurog2.bed',intervals2='NPC_Neurog2.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(2.5,1.5,2.5,4))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()
#Figure8J
p <- Figure8J(obs_f='/home/hpc/bonev/projects/hic/sc/data/cis_decay/NPC_Neurog2.bed_NPC_Neurog2.bed.1D.10000',min_dist = 1e5,max_dist = 4e6,TAD='intra',labels=c('GFP','Neurog2','NSC','IPC','PN'),cols=c(iue_cols,rep('#aaa9ad',3)),p_stats=c('GFP','Neurog2'))
p <- p + coord_cartesian(ylim=c(0.5,1)) + theme(legend.position = 'none')
pdf('plots/figures/Figure8J.pdf',height=4,width=5.5,useDingbats=FALSE)
print(p)
dev.off()
#Figure8K
res <- Figure8K(df=p2glinks$posCor,use_chip=T,chip_f="data/beds/NPC_Neurog2.bed",anno_mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T),vars=c('GFP','Neurog2'),dist_cutoff=1e4,intraTAD=F)
p <- res$plot + scale_fill_manual(name='',values=iue_cols) + xlab('') + ylab('Hi-C score') + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.97, 0.97)) + ylim(c(0,100))
pdf('plots/figures/Figure8K.pdf',height=4,width=5)
print(p)
dev.off()
#Figure8L
deseq2_fc <- read.table('data/IUE/Neurog2vsGFP_RNA_deseq2.tsv')
p <- Figure8L(fc=deseq2_fc,df=res$df_sub,vars='RNA_fc',id.vars='Neurog2_binding')
p <- p + coord_cartesian(ylim=c(0.25,2.5)) + scale_fill_grey(start = 1, end = .01) + xlab('') + ylab('Expression fold change (Neurog2/GFP)') + theme(legend.position = "none")
pdf('plots/figures/Figure8L.pdf',height=4,width=4)
print(p)
dev.off()  
#Figure8M
plotMisha(targetGene='Eomes',out_f='Figure8M',upstream=4.5e5,downstream=5e4,chipYlim <- matrix(c(0,2,0,100,0,100,0,100,0,100),nrow = 3,ncol = 2,byrow = T),
          chipNames=c('Neurog2','GFP CpG','Neurog2 CpG','GFP GpC','Neurog2 GpC'),window_scale=1.8,chipRes =20,pointCEX=1,conditions=c("hic.GFP_IUE24h.score_k100","hic.NGN2_IUE24h.score_k100"),binSize=5e3,radius=1.5e4,plot.dgram=FALSE,
          chipTracksToExtract=c('chipseq_RPM.NPC_Neurog2',"methylation.GFP_IUE24h_CpG_cov10x","methylation.NGN2_IUE24h_CpG_cov10x","methylation.GFP_IUE24h_GpC_cov10x","methylation.NGN2_IUE24h_GpC_cov10x"),methTracksToExtract=c("methylation.GFP_IUE24h_CpG_cov10x","methylation.NGN2_IUE24h_CpG_cov10x","methylation.GFP_IUE24h_GpC_cov10x","methylation.NGN2_IUE24h_GpC_cov10x"),methColors=rep(iue_cols,each=2),
          chipColors=c('black',rep(iue_cols,2)),scoreTrackToExtract =c("hic.GFP_IUE24h.score_k100","hic.NGN2_IUE24h.score_k100"),arcIntervals=p2glinks$posCor[p2glinks$posCor$gene_name=='Eomes'],arcColors=colorRampPalette(archr_colors[[29]]),
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

#Supplementary Figure S10
FigureS10B(bin_res='data/IUE/GFP_bin_counting.txt',out_f='FigureS10B',heigh=4,width=5)
FigureS10C(expr_df=read.table('data/IUE/rna_fpkm_average.tsv'),cluster='IPC',out_f='FigureS10C',width=5,height=4)
#FigureS10D
res <- getPlotSetArray(tracks=access_tracks,features="data/beds/NPC_CTCF.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20)
pdf('plots/figures/FigureS10D.pdf',height=6,width=6)
plotAverage(plotset=res, labels = iue_labels, xlim = NULL,
            ylim = c(4,25), main = NULL, xlab = "", ylab = 'Average Accessibility (% GpC Methylation)',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = iue_cols, pointsize = 12)
axis(side = 1,labels=c('-1000','-500','Ctcf','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=4,tick = T)
dev.off()
#FigureS10E 
res <- extract_MotifPos(pwm=readRDS('results/scATAC/combined_pwm.RDS'),chip_f="data/beds/NPC_Tbr2.bed",motif.name='Eomes') #Generate and save Eomes/Neurog2 motif positions
write.table(res,'results/beds/Eomes_ChIP_motifs.bed',quote=F,sep='\t',col.names=F,row.names=F)
res1 <- getPlotSetArray(tracks=access_tracks,features="results/beds/Neurog2_ChIP_motifs.bed",refgenome='mm10',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 10)
res2 <- getPlotSetArray(tracks=access_tracks,features="results/beds/Eomes_ChIP_motifs.bed",refgenome='mm10',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 10)
pdf('plots/figures/FigureS10E.pdf',height=6,width=6)
plotAverage(plotset=res1, labels = iue_labels, xlim = NULL,
            ylim = c(5,40), main = NULL, xlab = "", ylab = 'Average Accessibility (% GpC Methylation)',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = iue_cols, pointsize = 12)
axis(side = 1,labels=c('-1000','-500','Neurog2(var.2)','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=5,tick = T)
dev.off()
#FigureS10F
pdf('plots/figures/FigureS10F.pdf',height=6,width=6)
plotAverage(plotset=res2, labels = iue_labels, xlim = NULL,
            ylim = c(5,40), main = NULL, xlab = "", ylab = 'Average Accessibility (% GpC Methylation)',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = iue_cols, pointsize = 12)
axis(side = 1,labels=c('-1000','-500','Neurog2(var.2)','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=5,tick = T)
dev.off()
#FigureS10G
res <- getPlotSetArray(tracks=meth_tracks,features="data/beds/NPC_CTCF.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20)
pdf('plots/figures/FigureS10G.pdf',height=6,width=6)
plotAverage(plotset=res, labels = iue_labels, xlim = NULL,
            ylim = c(15,80), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = iue_cols, pointsize = 12)
axis(side = 1,labels=c('-1000','-500','Ctcf','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=15,tick = T)
dev.off()
#FigureS10H
plot_rankMatrix(file_f=paste0(main_f,'analysis/compartments/','rankMatrix_250kb_ranks',100),out_f='plots/figures/FigureS10H_1.pdf',zlim=c(-1,1),cells=cells[4:5],col=wide_red_blue_pal(1000),plot_chip=FALSE)
#FigureS10I
res_comp <- comp_boundaries(type='eigen',cells=cells[4:5], write_borders=F)
res_comp$condition <- iue_labels
df <- melt(res_comp)
df$condition <- factor(df$condition,levels=iue_labels)
df$variable <- factor(df$variable,levels=c('A-B','B-A'))
p <- ggplot(df, aes(x=condition,y=value,fill=variable)) +  geom_bar(position="stack", stat="identity",colour="black") + scale_fill_manual(name='',values = c("#F0F0F0","#525252"))
p <- p + ylab('Compartment transitions') + xlab('') + ylim(c(0,1300)) + theme(legend.position=c(0.2,0.95),legend.direction = 'horizontal')
pdf('plots/figures/FigureS10I.pdf',height=4,width=3)
print(p)
dev.off()
#FigureS10J
res1 <- averageTrack(tracks=c(meth_misha_tracks[1]),regions="hic.GFP_IUE24h.ins_250_domains_expanded",bins=100,anchor_middle=F)
res2 <- averageTrack(tracks=c(meth_misha_tracks[2]),regions="hic.NGN2_IUE24h.ins_250_domains_expanded",bins=100,anchor_middle=F)
zlim=c(-1.5,1)
blue_white_pal = colorRampPalette(c("purple", "navy", "blue", "#87FFFF", "white"))
white_red_pal = colorRampPalette(c("white","#FF413D", "black", "orange", "yellow"))
averageTAD_colors <- c(blue_white_pal(length(seq(zlim[1],0.01,by=0.01))),'white',white_red_pal(length(seq(0.01,zlim[2],by=0.01))))
plot_averageTAD(tracks=all_tracks,fig_name='plots/figures/FigureS10J.pdf',add_matrix=list(GFP_IUE24h=res1,NGN2_IUE24h=res2),mat_lim=c(60,90),cells=cells[4:5],path='/home/hpc/bonev/projects/hic/sc/analysis/averageTAD/',file_f='IUE_averageTAD',stats_f=F,plot_what='',zlim=zlim,flip=T,z_colors=averageTAD_colors,height=2.5,width=4.5)
#FigureS10K
p <- plot_pair(borders=readRDS('results/HiC/unity_borders.RDS'),tss=gintervals.load(tss_f),
               extra_tracks=NULL,ins_tracks=ins_tracks[c(4,5)],
               plot_what='ins',plot_cells=c('GFP','NGN2'),mode = 'boxplot',cols = c(iue_cols,rep('#aaa9ad',2)),ylab_n = 'Insulation score')
p <- p + coord_cartesian(ylim=c(1.2,3.7)) + scale_y_continuous(n.breaks=6)
pdf('plots/figures/FigureS10K.pdf',height=4,width=2)
print(p)
dev.off()
#FigureS10L
run_cisDecay(tracks=all_tracks[10:13],out_f='plots/figures/FigureS10L.pdf',cells=cells[4:5],path=paste0(main_f,'analysis/cis_decay/'),log_base=10,file_f='IUE_cisDecay',maxDist=1e8,colors=iue_cols,alpha=0.3)
#FigureS10M
pdf('plots/figures/FigureS10M.pdf',width=4.5,height=4)
layout(matrix(c(1:3),nrow=1,ncol=3,byrow=F),widths = c(4,4,2),heights=c(4),respect = T)
for (cell in c('GFP_IUE24h','NGN2_IUE24h')){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='NPC_Pax6',intervals2='NPC_Pax6',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(3,1,3,3))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()
#FigureS10N
p <- Figure8J(obs_f='/home/hpc/bonev/projects/hic/sc/data/cis_decay/NPC_Pax6_NPC_Pax6.1D.10000',min_dist = 1e5,max_dist = 4e6,TAD='intra',labels=c('GFP','Neurog2','NSC','IPC','PN'),cols=c(iue_cols,rep('#aaa9ad',3)),p_stats=c('GFP','Neurog2'))
p <- p + coord_cartesian(ylim=c(0.5,1))+xlim(c('GFP','Neurog2'))
pdf('plots/figures/FigureS10N.pdf',height=4,width=4,useDingbats=FALSE)
print(p)
dev.off()



