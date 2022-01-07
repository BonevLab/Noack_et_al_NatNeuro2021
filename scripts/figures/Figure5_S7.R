library(Seurat)
library(LSD)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(readr)

set.seed(123)
source('scripts/figures/config.R')
source('scripts/aux_functions.R')
source('scripts/figures/plot_functions.R')

source('scripts/hic/config.R')
source('scripts/hic/scripts/main_functions.R')
source('scripts/hic/scripts/aux_functions.R')
source('scripts/hic/scripts/plot_functions.R')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))

meth_tracks <- c("methylation.E14_NSC_10x","methylation.E14_IPC_10x","methylation.E14_PN_10x")
rna.object <- readRDS('data/merged_scRNA_unfiltered_IDs.RDS')

Figure5E <- function(input,anno.size=12,features=NULL,width,height,out_f){
  set.seed(123)
  input_con <- input$constant
  input <- input$differential
  input[,4:6] <- input[,4:6]*(-1)
  input <- gintervals.neighbors(input,gintervals.load(tss_f)[,1:5])
  input$rowVars <- rowVars(as.matrix(input[,4:6]))
  input <- input[order(input$rowVars,decreasing=T),]
  row.names(input) <- 1:nrow(input)
  mat <- as.matrix(input[,4:6])
  colnames(mat) <- c('NSC','IPC','PN')
  clust <- kmeans(mat,centers = 4,nstart = 1000,iter.max = 100000)
  input$cluster <- clust$cluster
  la1 <- rowAnnotation(foo = anno_block(gp = gpar(fill = c('darkblue','blue','red','darkgreen')),labels = c("C1", "C2", "C3","C4"),labels_gp = gpar(col = "white", fontsize = 14)))
  ra1 = rowAnnotation(foo = anno_mark(at = which(input$geneName%in%features),side='right', labels = paste0(input$geneName,':',input$dist)[input$geneName%in%features],labels_gp = gpar(fontsize = anno.size)))
  hm1 <- Heatmap(mat,name='Insulation',column_title = '',cluster_rows =  F,show_row_dend = F,row_split = clust$cluster,row_title = NULL,column_names_rot = 0,column_names_centered = T,row_dend_reorder = FALSE,row_order=clust$order,cluster_columns = F,show_column_names = T,show_row_names = F,right_annotation = ra1,left_annotation = la1,
                 heatmap_legend_param=list(direction = "horizontal",title='Insulation Score',legend_width = unit(width/2, "inch")))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  draw(hm1, merge_legends=F, ht_gap = unit(0.2, "inch"),heatmap_legend_side='bottom',annotation_legend_side = "bottom") 
  dev.off()
}


Figure5H <- function(object,features,cols,point.size,out_f,height,width,plot_filled=F,anno.size=10,theme=NULL,direction='vertical'){
  plot_list <- list()
  for (feature in features){
    plot.data <- extractFeatures(object,features=feature,cells_toInclude = c('all'),cells_toExclude = 'none')
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(color=cols[1],size = point.size,data = plot.data[plot.data$feature==0,]) + geom_point(aes(color=feature),size = point.size,data = plot.data[plot.data$feature>0,]) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_color_gradientn(colours=cols,name='',breaks=c(0,max(plot.data$feature)),labels=c(0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1, 1),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
      p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
    } else {
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=0.1) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colours=cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(min(plot.data$feature,na.rm=T),round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1.05, 0.98),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 7.5))
      p <- p + guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
    }
    if(!is.null(theme)){
      p <- p+theme
      p <- p + annotate("text",x = -Inf, y = Inf, label = feature,size=anno.size,hjust = -0.25, vjust = 1.5)
    } else {
      p <- p + ggtitle(feature)
    }
    plot_list[[feature]] <- p
  }
  p <- Reduce("+", plot_list) 
  if (direction=='vertical'){
    height=length(features)*height
    p_layout <- plot_layout(nrow=length(features),ncol=1)
  } else {
    width=length(features)*width
    p_layout <- plot_layout(ncol=length(features),nrow=1)
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p+p_layout)
  dev.off()
}


FigureS7B <- function(file_f,conditions,reps,which_genes,out_f,width,height,cols,ylim=NULL){
  mat <- read.csv(file_f,header=T,row.names = 'Gene')
  s_names <- paste0(rep(conditions,each=3),reps)
  df <- mat[which_genes,]
  df$Gene <- factor(which_genes,levels=which_genes)
  df <- melt(df,varnames = c('Gene'),value.name = 'FC')
  df$Condition <- factor(rep(conditions,each=3*length(which_genes)),levels=conditions)
  df_summ <- data_summary(df, varname="FC",groupnames=c('Condition','Gene'))
  p <- ggplot(df_summ, aes(x=Gene, y=FC, fill=Condition)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Relative Fold Change') +
    geom_errorbar(aes(ymin=FC-sd, ymax=FC+sd), width=.2,position=position_dodge(.9))
  p <- p + theme(legend.position='none') + geom_point(data = df,aes(x=Gene, y=FC, fill=Condition),position = position_jitterdodge(dodge.width =.9,jitter.width = 0.25,seed = 42),size=2) 
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
  if(!is.null(ylim)){
    print(p+ylim(ylim))
  } else {
    print(p)
  }
  dev.off()
}

FigureS7CD <- function(path,conditions,reps,which_control,ylim,out_f,width,height,cols){
  s_names <- paste0(rep(conditions,each=3),reps)
  all_files <- list.files(path = path,pattern = 'bt2_PE_report.txt',full.names = T,recursive = T,include.dirs = T)
  all_files <- all_files[grep(which_control,all_files)]
  df <- matrix(NA,nrow=3,ncol=3)
  colnames(df) <- conditions
  row.names(df) <- reps
  for (cond in seq_along(conditions)){
    for (rep in reps){
      file_l <- all_files[grep(paste0(conditions[cond],rep),all_files)]
      df_l <- suppressWarnings(readLines(file_l))
      df[rep,cond] <-  round(parse_number(df_l)[26]/(parse_number(df_l)[26]+parse_number(df_l)[31])*100,3)
    }
  }
  df <- melt(df,varnames = c('Rep','Condition'),value.name = 'CpG_Methylation')
  df_summ <- data_summary(df, varname="CpG_Methylation",groupnames=c("Condition"))
  p <- ggplot(df_summ, aes(x=Condition, y=CpG_Methylation, fill=Condition)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('% CpG Methylation') +
    geom_errorbar(aes(ymin=CpG_Methylation-sd, ymax=CpG_Methylation+sd), width=.2,position=position_dodge(.9))
  p <- p + theme(legend.position='none') + geom_point(data = df,position=position_jitter(w = 0.1, h = 0),size=2) + annotate("text", x = c(1.25,2.25,3.25), y = df_summ$CpG_Methylation*1.05, label = round(df_summ$CpG_Methylation,2)) + scale_y_continuous(expand=c(0,0),limits=ylim)
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height,useDingbats = F)
  print(p)
  dev.off()
}


FigureS7F <- function(cor_file,cells=c('NSC','IPC','PN'),cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),cluster_cols=c("#F51111","#16A810","#4861F0"),out_f,width,height){
  res <- read.table(cor_file)
  row.names(res) <- paste0(rep(cells,each=3),'_',1:3)
  colnames(res) <- row.names(res)
  distance <- dist(1-res)
  res <- round(res, 2)
  annotation_df <- data.frame(cell_Type=factor(rep(cells,each=3),levels=cells))
  colnames(annotation_df) <- c('Condition')
  row.names(annotation_df) <- colnames(res)
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
  pheatmap(res,color=cols,breaks=seq(0.85,1,length=101),clustering_distance_rows=distance,clustering_distance_cols=distance,clustering_method='ward.D',show_colnames = F,annotation_colors=list(Condition=c(NSC=cluster_cols[1],IPC=cluster_cols[2],PN=cluster_cols[3])),
           angle_col = 0,display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE,legend_breaks=c(0.85,0.9,0.95,1))
  dev.off()
}

FigureS7G<-function(meth_misha_rep_tracks,height,width,out_f,cluster_cols,peaks,window_size) {
  res <- extract_lin(regions=peaks,window=window_size,tracks=meth_misha_rep_tracks)
  colnames(res)[4:12] <- c(paste0('NSC rep',1:3),paste0('IPC rep',1:3),paste0('PN rep',1:3))
  res<-na.omit(res)
  df<-as.matrix(res[,4:12])
  res<-cor(df,method='pearson')
  distance <- dist(1-res)
  res <- round(res, 2)
  cells<-c('NSC','IPC','PN')
  annotation_df <- data.frame(cell_Type=factor(rep(cells,each=3),levels=cells))
  colnames(annotation_df) <- c('Condition')
  row.names(annotation_df) <- colnames(res)
  cols=colorRampPalette(rev(colorpalette('ylorrd')))(100)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(pheatmap(res,color=cols,breaks=seq(0.85,1,length=101),clustering_distance_rows=distance,clustering_distance_cols=distance,clustering_method='ward.D',show_colnames = F,annotation_colors=list(Condition=c(NSC=cluster_cols[1],IPC=cluster_cols[2],PN=cluster_cols[3])),
                 angle_col = 0,display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE,legend_breaks=c(0.85,0.9,0.95,1))
  )
  dev.off()
}

FigureS7H <- function(cells,levels_names=c('NSC','IPC','PN'),out_f,height=4,width=3){
  res_comp <- comp_boundaries(type='eigen',cells=cells, write_borders=F)
  res_comp$condition <- gsub('E14_','',row.names(res_comp))
  df <- melt(res_comp)
  df$condition <- factor(df$condition,levels=levels_names)
  df$variable <- factor(df$variable,levels=c('A-B','B-A'))
  p <- ggplot(df, aes(x=condition,y=value,fill=variable)) +  geom_bar(position="stack", stat="identity",colour="black") + scale_fill_manual(name='',values = c("#F0F0F0","#525252"))
  p <- p + ylab('Compartment transitions') + xlab('') + theme(legend.position=c(0.75,0.95))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS7IJ <- function(borders,tss,extra_tracks,ins_tracks,features=NULL,width,height,plot_what='ins',plot_cells=c('NSC','PN'),point.size=1,anno.size=4,out_f,mode='scatter',cols,ylab_n){
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


igureS5B <- function(input=readRDS('/home/hpc/bonev/projects/hic/sc/analysis/insulation/res_0.8.RDS'),all_b,tss=gintervals.load(tss_f)){
  diff <- input$differential
  const <- input$constant
  diff$dist <- gintervals.neighbors(diff,tss)$dist
  const$dist <- gintervals.neighbors(const,tss)$dist
  all_b$dist <- gintervals.neighbors(all_b,tss)$dist
}

figureS5K <- function(input,cols=ATAC_cluster_colors[1:3],width,height){
  input$labels <- paste0(input$chrom,':',input$start,'-',input$end)
  df_ins <- melt(input[,c(grep('labels',colnames(input)),grep('ins',colnames(input)),grep('cluster',colnames(input)))],id.vars = c('labels','cluster'))
  df_ins$variable <- gsub('v_hic.E14_','',df_ins$variable)
  df_ins$variable <- gsub('.ins_250','',df_ins$variable)
  df_ins$variable <- factor(df_ins$variable,levels=c('NSC','IPC','PN'))
  df_meth <- melt(input[,c(grep('labels',colnames(input)),grep('meth',colnames(input)),grep('cluster',colnames(input)))],id.vars = c('labels','cluster'))
  df_meth$variable <- gsub('v_methylation.E14_','',df_meth$variable)
  df_meth$variable <- gsub('_10x','',df_meth$variable)
  df_meth$variable <- factor(df_meth$variable,levels=c('NSC','IPC','PN'))
  df_expr <- melt(input[abs(input$dist)<=5e4,c(grep('labels',colnames(input)),which(colnames(input)%in%c('NSC','IPC','PN')),grep('cluster',colnames(input)))],id.vars = c('labels','cluster'))
  df_expr$variable <- factor(df_expr$variable,levels=c('NSC','IPC','PN'))
  plot_list <- list()
  for (i in 1:max(input$cluster)){
    p1 <- ggplot(subset(df_ins,cluster==i),aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    p1 <- p1 + scale_fill_manual(values=cols) + xlab('')  + scale_y_continuous(n.breaks =4,limits = c(1.75,3.5)) + ylab(ifelse(i==1,'Insulation Score','')) + theme(legend.position = "none")
    plot_list[[paste0('ins_',i)]] <- p1 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','PN'),c('IPC','PN')),label = "p.signif",method='wilcox',method.args = list(paired=TRUE)) + ggtitle(paste0('Cluster',i))
    p2 <- ggplot(subset(df_meth,cluster==i),aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    p2 <- p2 + scale_fill_manual(values=cols) + xlab('') + scale_y_continuous(breaks =c(0,50,100),limits = c(0,120)) + ylab(ifelse(i==1,'% CpG Methylation','')) + theme(legend.position = "none")
    plot_list[[paste0('meth_',i)]] <- p2 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','PN'),c('IPC','PN')),label = "p.signif",method='wilcox',method.args = list(paired=TRUE)) + ggtitle(' ')
  #  p3 <- ggplot(subset(df_expr,cluster==i),aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
 #   p3 <- p3 + scale_fill_manual(values=cols) + xlab('') + ylab('log2(FPKM+1)') + theme(legend.position = "none")
  #  plot_list[[paste0('expr_',i)]] <- p3 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','PN'),c('IPC','PN')),label = "p.signif",method='wilcox',method.args = list(paired=TRUE)) + ggtitle(' ')
  }
  p <- Reduce("+", plot_list) 
  pdf('plots/figures/figureS5K.pdf',width=width*(length(unique(input$cluster))),height=height*2)
  print(p + theme_cowplot()+ plot_layout(ncol=length(unique(input$cluster)),nrow=2,byrow = F))
  dev.off()
  
}


#### Plot Figures ######

#Figure 5B
plotBinned(extra_tracks=meth_tracks,fig_name='plots/figures/Figure5B_1.pdf',cells='E14_NSC',init_cell="E14_NSC",chr='chr3',binSize=2e5,path=paste0(main_f,'analysis/compartments/'),file_f=paste0(main_f,'analysis/compartments/chr3_200kb'),plot_what='obs',balance=T,width=4,height=4.5)
plotBinned(extra_tracks=meth_tracks,fig_name='plots/figures/Figure5B_2.pdf',cells='E14_IPC',init_cell="E14_NSC",chr='chr3',binSize=2e5,path=paste0(main_f,'analysis/compartments/'),file_f=paste0(main_f,'analysis/compartments/chr3_200kb'),plot_what='obs',balance=T,width=4,height=4.5)
plotBinned(extra_tracks=meth_tracks,fig_name='plots/figures/Figure5B_3.pdf',cells='E14_PN',init_cell="E14_NSC",chr='chr3',binSize=2e5,path=paste0(main_f,'analysis/compartments/'),file_f=paste0(main_f,'analysis/compartments/chr3_200kb'),plot_what='obs',balance=T,width=4,height=4.5)
pdf('plots/figures/Figure5B_scaleBar.pdf',width=0.75,height=4.5)
par(mar=c(0.5,0.5,0.5,1.5))
image.scale(as.matrix(1),zlim=c(0,0.001628), col=colorRampPalette(c("white","orange","red","darkRed"))(1000),axis.pos=4,adj=1,cex.axis=1,axis_lab=c('','',''))
dev.off()

#Figure 5C
zlim=c(-1.5,1)
blue_white_pal = colorRampPalette(c("purple", "navy", "blue", "#87FFFF", "white"))
white_red_pal = colorRampPalette(c("white","#FF413D", "black", "orange", "yellow"))
averageTAD_colors <- c(blue_white_pal(length(seq(zlim[1],0.01,by=0.01))),'white',white_red_pal(length(seq(0.01,zlim[2],by=0.01))))
res1 <- averageTrack(tracks="methylation.E14_NSC_10x",regions="hic.E14_NSC.ins_250_domains_expanded",bins=100,anchor_middle=F)
res2 <- averageTrack(tracks="methylation.E14_IPC_10x",regions="hic.E14_IPC.ins_250_domains_expanded",bins=100,anchor_middle=F)
res3 <- averageTrack(tracks="methylation.E14_PN_10x",regions="hic.E14_PN.ins_250_domains_expanded",bins=100,anchor_middle=F)
plot_averageTAD(tracks=all_tracks,fig_name='plots/figures/Figure5C.pdf',add_matrix=list(E14_NSC=res1,E14_IPC=res2,E14_PN=res3),mat_lim=c(60,90),cells=cells[1:3],path='/home/hpc/bonev/projects/hic/sc/analysis/averageTAD/',file_f='averageTAD',stats_f=F,plot_what='',zlim=zlim,flip=T,z_colors=averageTAD_colors,height=2.5,width=4.5)

#Figure 5D
plot_rankMatrix(file_f=paste0(main_f,'analysis/compartments/','E14_rankMatrix_250kb_ranks',100),out_f='plots/figures/Figure5D.pdf',zlim=c(-1,1),cells=cells[1],col=wide_red_blue_pal(1000),plot_chip=FALSE)

#Figure 5E
Figure5E(input=readRDS('/home/hpc/bonev/projects/hic/sc/analysis/insulation/res_0.8.RDS'),anno.size=12,features=c('Flrt2','Gli3','Gas1','Myt1l','Cxcl12'),width=6,height=8,out_f='Figure5E')

#Figure 5F
plotMisha(targetGene='Gas1',out_f='Figure5F',upstream=4e5,downstream=3.5e5,
          chipNames=c('NSC','IPC','PN'),window_scale=2,chipRes =50,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN2'),pointCEX=1.5,
          chipColors=ATAC_cluster_colors[c(1,2,4)],scoreTrackToExtract =score_tracks[1:3],
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.6, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

#Figure 5G
plotMisha(object=atac.object,targetGene='Flrt2',out_f='Figure5G',upstream=4e6,downstream=2.5e6,
          chipNames=c('NSC','IPC','PN'),window_scale=2,pointCEX=0.5,chipRes =200,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN2'),
          chipColors=ATAC_cluster_colors[c(1,2,4)],scoreTrackToExtract =score_tracks[1:3],
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))
pdf('plots/figures/Figure5G_scaleBar.pdf',height=6,width=1.2)
par(mar=c(0.5,0.5,0.5,4))
image.scale(as.matrix(1),zlim=c(-100,100), col=col.scores[1:201],axis.pos=4,adj=1,cex.axis = 1.2)
dev.off()

#Figure 5H
Figure5H(object=rna.object,c('Gas1','Flrt2'),cols=gene_colors(50),point.size=2,out_f='Figure5H',width=5,height=5,plot_filled=T,direction='vertical',anno.size=8,theme=theme_border)

### Supplementary Figure 7

FigureS7B(file_f="data/E14_methylHiC_FACS_qPCR.csv",which_genes=c('Hes1','Pax6','Sox2'),conditions=c('NSC','IPC','PN'),reps=1:3,out_f='FigureS7B_1',width=6,height=5,cols=c("#F51111","#16A810","#4861F0"))
FigureS7B(file_f="data/E14_methylHiC_FACS_qPCR.csv",which_genes=c('Eomes'),conditions=c('NSC','IPC','PN'),reps=1:3,out_f='FigureS7B_2',width=2.5,height=5,cols=c("#F51111","#16A810","#4861F0"))
FigureS7B(file_f="data/E14_methylHiC_FACS_qPCR.csv",which_genes=c('Btg2'),conditions=c('NSC','IPC','PN'),reps=1:3,out_f='FigureS7B_3',width=2.5,height=5,cols=c("#F51111","#16A810","#4861F0"))
FigureS7B(file_f="data/E14_methylHiC_FACS_qPCR.csv",which_genes=c('Sox5','Tubb3'),conditions=c('NSC','IPC','PN'),reps=1:3,out_f='FigureS7B_4',width=4,height=5,cols=c("#F51111","#16A810","#4861F0"))
FigureS7CD(path='/home/hpc/bonev/projects/hic/e14_fnoack/juicer/',conditions=c('NSC','IPC','PN'),reps=1:3,which_control='lambda',ylim=c(0,105),out_f='FigureS7C',width=3.5,height=6,cols=c("#F51111","#16A810","#4861F0"))
FigureS7CD(path='/home/hpc/bonev/projects/hic/e14_fnoack/juicer/',conditions=c('NSC','IPC','PN'),reps=1:3,which_control='puc19',ylim=c(0,105),out_f='FigureS7D',width=3.5,height=6,cols=c("#F51111","#16A810","#4861F0"))
#Figure S7E
run_cisDecay(tracks=all_tracks[1:9],out_f='plots/figures/FigureS7E.pdf',cells=cells[1:3],path=paste0(main_f,'analysis/cis_decay/'),log_base=10,file_f='cisDecay',maxDist=1e8,colors=c("#F71911B1","#13A810BC","#485CF2BC"),alpha=0.3,disp_colors=c('#FAD3D03A','#B0FFB156','#66E8FF48'))
#Figure S7F
FigureS7F(cor_file='results/HiC/hic_correlation.tsv',cells=c('NSC','IPC','PN'),cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),out_f='FigureS7F',width=8,height=8)
#FigureS7G
meth_misha_rep_tracks <- c("methylation.NSC_rep1_methylation_CpG_5x","methylation.NSC_rep2_methylation_CpG_5x","methylation.NSC_rep3_methylation_CpG_5x",
                           "methylation.IPC_rep1_methylation_CpG_5x","methylation.IPC_rep2_methylation_CpG_5x","methylation.IPC_rep3_methylation_CpG_5x",
                           "methylation.PN_rep1_methylation_CpG_5x","methylation.PN_rep2_methylation_CpG_5x","methylation.PN_rep3_methylation_CpG_5x")
FigureS7G(meth_misha_rep_tracks=meth_misha_rep_tracks,height=8, width=8,out_f='FigureS7G',cluster_cols=c("#F51111","#16A810","#4861F0"),peaks='data/unionPeaks.bed',window_size=500)
#Figure S7H
FigureS7H(cells=cells[1:3],levels_names=c('NSC','IPC','PN'),out_f='FigureS7H',height=4,width=3)
#Figure S7I
p <- FigureS7IJ(borders=readRDS('results/HiC/unity_borders.RDS'),tss=gintervals.load(tss_f),
               extra_tracks=NULL,ins_tracks=ins_tracks[1:3],
               plot_what='ins',plot_cells=c('NSC','IPC','PN'),mode = 'boxplot',cols = c("#F51111","#16A810","#4861F0"),ylab_n = 'Insulation score')
p <- p + coord_cartesian(ylim=c(1.2,3.7)) + stat_compare_means(comparisons = list(c('NSC','IPC'),c('IPC','PN'),c('NSC','PN')),label = "p.format",method='wilcox',label.y = c(3.25,3.4,3.6),tip.length = c(0)) + scale_y_continuous(n.breaks=6)
pdf('plots/figures/FigureS7I_1.pdf',height=4,width=2)
print(p)
dev.off()
p <- FigureS7IJ(borders=intervals.normalize(gintervals.load(tss_f),1000),tss=gintervals.load(tss_f),
               extra_tracks=NULL,ins_tracks=ins_tracks[1:3],
               plot_what='ins',plot_cells=c('NSC','IPC','PN'),mode = 'boxplot',cols = c("#F51111","#16A810","#4861F0"),ylab_n = 'Insulation score')
p <- p + coord_cartesian(ylim=c(1.2,3.7)) + stat_compare_means(comparisons = list(c('NSC','IPC'),c('IPC','PN'),c('NSC','PN')),label = "p.format",method='wilcox',label.y = c(3.25,3.4,3.6),tip.length = c(0)) + scale_y_continuous(n.breaks=6)
pdf('plots/figures/FigureS7I_2.pdf',height=4,width=2)
print(p)
dev.off()

#Figure S7J
FigureS7IJ(borders=readRDS('results/HiC/unity_borders.RDS'),tss=gintervals.load(tss_f),
          out_f='FigureS7J_1',extra_tracks=meth_tracks[c(2,1,3)],ins_tracks=ins_tracks[1:3],features=c('Flrt2','Gas1'),width=4,height=4,
          plot_what='ins',plot_cells=c('NSC','PN'),point.size=0.5,anno.size=4)
FigureS7IJ(borders=readRDS('results/HiC/unity_borders.RDS'),tss=gintervals.load(tss_f),
          out_f='FigureS7J_2',extra_tracks=meth_tracks[c(2,1,3)],ins_tracks=ins_tracks[1:3],features=NULL,width=4,height=4,
          plot_what='meth',plot_cells=c('NSC','PN'),point.size=0.5,anno.size=4)
FigureS7IJ(borders=intervals.normalize(gintervals.load(tss_f),1000),tss=gintervals.load(tss_f),
          out_f='FigureS7J_3',extra_tracks=meth_tracks[c(2,1,3)],ins_tracks=ins_tracks[1:3],features=c('Flrt2','Gas1'),width=4,height=4,
          plot_what='ins',plot_cells=c('NSC','PN'),point.size=0.5,anno.size=4)
FigureS7IJ(borders=intervals.normalize(gintervals.load(tss_f),1000),tss=gintervals.load(tss_f),
          out_f='FigureS7J_4',extra_tracks=meth_tracks[c(2,1,3)],ins_tracks=ins_tracks[1:3],width=4,height=4,
          plot_what='meth',plot_cells=c('NSC','PN'),point.size=0.5,anno.size=4)

library(seqplots)
tracks <- paste0('results/methylation/',c('NSC','IPC','PN'),'_methylation_CpG_10x.bw')
tracks2 <- paste0('results/scATAC_clusters/',c('NSC','IPC','PN2'),'.bw')
res_meth <- getPlotSetArray(tracks=tracks,features="data/beds/NPC_CTCF.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10)
res_atac <- getPlotSetArray(tracks=tracks2,features="data/beds/NPC_CTCF.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10)

pdf('plots/figures/FigureS7K_1.pdf',height=5,width=5)
plotAverage(plotset=res_meth, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(10,70), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = c("#F51111","#16A810","#4861F0"), pointsize = 12)
axis(side = 1,labels=c('-1000','CTCF','+1000'),at=c(-1000,0,1000),pos=10,tick = T)
dev.off()

pdf('plots/figures/FigureS7K_2.pdf',height=5,width=5)
plotAverage(plotset=res_atac, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,1.2), main = NULL, xlab = "", ylab = 'Average Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = c("#F51111","#16A810","#4861F0"), pointsize = 12)
axis(side = 1,labels=c('-1000','CTCF','+1000'),at=c(-1000,0,1000),pos=0,tick = T)
dev.off()

plotMisha(object=atac.object,targetGene='Cxcl12',out_f='FigureS7L',upstream=7e5,downstream=8e5,
          chipNames=c('NSC','IPC','PN'),window_scale=2,chipRes =100,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN2'),pointCEX=1.5,
          chipColors=ATAC_cluster_colors[c(1,2,4)],scoreTrackToExtract =score_tracks[1:3],
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.6, domains=0.15, genes=1.2, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

figureS5K(input=diff_res,cols=ATAC_cluster_colors[1:3],width=4,height=4.2)





