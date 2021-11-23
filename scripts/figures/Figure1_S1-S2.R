library(Seurat)
library(ggplot2)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(dplyr)
library(Hmisc)
library(SummarizedExperiment)
require(edgeR)
require(LSD)
library(ggrepel)
library(ggpointdensity)
library(ggrastr)
library(grid)
library(ComplexHeatmap)
library(org.Mm.eg.db)
library(clusterProfiler)

## Set main folder as current dir ####

source('scripts/figures/config.R')
source('scripts/figures/plot_functions.R')

object <- readRDS('data/merged_scRNA_unfiltered_IDs.RDS')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))

Figure1B <- function(object,cols,out_f='Figure1B',point.size,anno.size,key.size,height,width,plot_filled=F,theme=NULL,stroke=0.1,rows=1){
  plot.data <- extractFeatures(object,features='seurat_clusters',cells_toInclude = c('all'),cells_toExclude = 'none')
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=labels)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = labels), pch = I(21),size = point.size,stroke=stroke,alpha=1) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0,(1.03+0.025*rows)),plot.margin = margin(height/10,0.1,0.1,0.1,unit='inches'), legend.box = "horizontal",legend.spacing.x = unit(0.1, 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = rows)) 
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  if(!is.null(theme)){
    print(p + theme)}
  else {
    print(p)
  }
  dev.off()
  }

Figure1C <- function(object,features,cols,point.size,height,width,plot_filled=F,anno.size=10,theme=NULL,direction='vertical',stroke=0.1,out_f){
  plot_list <- list()
  for (feature in features){
    plot.data <- extractFeatures(object,features=feature,cells_toInclude = c('all'),cells_toExclude = 'none')
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(color=cols[1],size = point.size,data = plot.data[plot.data$feature==0,]) + geom_point(aes(color=feature),size = point.size,data = plot.data[plot.data$feature>0,]) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_color_gradientn(colours=cols,name='',breaks=c(0,max(plot.data$feature)),labels=c(0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1, 1),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
      p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
    } else {
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colours=cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(min(plot.data$feature,na.rm=T),round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1.05, 0.98),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 7.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))
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

Figure1F <- function(object,pgraph=NULL,feature,cols,point.size,line.size,out_f,width,height,plot_filled=F,theme=NULL,stroke=0.1){
  plot.data <- extractFeatures(object,features=feature,cells_toInclude = c('all'),cells_toExclude = 'none')
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") 
    p <- p + scale_color_gradientn(colours = cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(round(min(plot.data$feature,na.rm=T)),round(max(plot.data$feature,na.rm=T))))
    p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 1),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.4), vjust = 8))
    p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
    p <- p + scale_fill_gradientn(colours = cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(round(min(plot.data$feature,na.rm=T)),round(max(plot.data$feature,na.rm=T))))
    p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1.02, 0.98),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 7.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))
    p <- p + guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
  }
  if (!is.null(pgraph)){
    edge_df <- readRDS(pgraph)
    p <- p + geom_segment(aes_string(x="source_prin_graph_dim_1",y="source_prin_graph_dim_2",xend="target_prin_graph_dim_1",yend="target_prin_graph_dim_2"),
                   size=line.size,color=I('black'),linetype="solid",na.rm=TRUE,data=edge_df)
  }
  if(!is.null(theme)){
    p <- p+theme
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
  }

Figure1G <- function(mat,object,features,cols,cluster_cols,anno.size,out_f,width,height){
  ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%(features)), labels = row.names(mat)[which(row.names(mat)%in%unique(c(features)))],labels_gp = gpar(fontsize = anno.size)))
  col.list <- cluster_cols
  names(col.list) <- levels(object)
  idents_indx <- Idents(object)
  idents_indx <- idents_indx[match(colnames(mat),names(idents_indx))]
  
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm <- Heatmap(mat, name = "Pseudotime",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, right_annotation = ra,col = cols,top_annotation = ha,
          heatmap_legend_param=list(labels = c(0,100),at=c(0,1), direction = "vertical",title='% Max',legend_height = unit(height/6, "inch")))
  draw(hm, merge_legends=T, ht_gap = unit(height/2, "inch")) 
  dev.off()
}

Figure1H <- function(object,features,out_f,cols,point.size,alpha,anno.size,key.size,height,width,plot_filled=F,direction='horizontal'){
  plot_list <- list()
  for (feature in features){
    plot.data <- extractFeatures(object,features=c(feature,'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    plot.data <- plot.data[!is.na(plot.data$pseudotime),]
    plot.data$labels <- droplevels(plot.data$labels)
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=pseudotime,y=feature,color=labels)) + geom_point(size = point.size,alpha=alpha) + ylab("Expression") + xlab("Pseudotime") + scale_colour_manual(name='',values=as.character(cols)) + geom_smooth(aes(x=pseudotime,y=feature),inherit.aes = F,colour='black')  
      p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) + ggtitle(feature)
    } else {
      p <- ggplot(plot.data, aes(x=pseudotime,y=feature)) + geom_point(aes(fill = labels),colour='black', pch = I(21),size = point.size,alpha=alpha) + ylab("Expression") + xlab("Pseudotime") + scale_fill_manual(name='',values=as.character(cols)) + geom_smooth(aes(x=pseudotime,y=feature),inherit.aes = F,colour='black')
      p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = 1)) + ggtitle(feature)
    }
    plot_list[[feature]] <- p
  }
  if(direction=='vertical'){
    p <- Reduce("+", plot_list)
    p <- p + guide_area() + plot_layout(guides="collect",heights = c(rep(9,length(features)),1),ncol=1,nrow=length(features)+1)
    height=length(features)*height
  } else {
    p <- Reduce("|", plot_list) 
    p <- guide_area()/p + plot_layout(guides="collect",heights = c(1, 9))
    width=length(features)*width
  }
  
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + theme_cowplot())
  dev.off()
}


#### Figure S1 functions ########

FigureS1A <- function(object,out_f,height,width){
  object_cpm <- as.data.frame(matrix(NA,nrow=nrow(object),ncol=length(unique(object$orig.ident))))
  colnames(object_cpm) <- unique(object$orig.ident)
  for (s in unique(object$orig.ident)){
    object_sub <- object[,object$orig.ident==s]
    object_cpm[,s] <- cpm(rowSums(object_sub@assays$RNA@counts),log=T,prior.count = 1)
  }
  p <-  ggplot(object_cpm, aes( x = E14_rep1, y = E14_rep2 )) + geom_pointdensity(adjust=1.5,alpha=1,size=1) + scale_color_gradientn(colours = colorpalette("heat",8),name='',breaks=c(1,3000),labels=c('min','max'))             # scale_color_gradientn(colours=tim.colors(12))
  p<- p + theme_cowplot() + geom_abline(slope = 1,intercept = 0,linetype = 2) + xlab('Rep1 log2(CPM+1)') + ylab('Rep2 log2(CPM+1)')  
  p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.3, 0.97),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.4), vjust = 8))
  my_grob = grid.text(paste0('r=',round(cor.test(object_cpm[,1],object_cpm[,2])$estimate,3)), x=0.85,  y=0.1, gp=gpar(col="black", fontsize=14, fontface="bold"))
  p1 <- p + annotation_custom(my_grob)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p1)
  dev.off()
}

FigureS1B <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  p <- ggplot(plot.data, aes(x=reps,fill=labels)) + geom_bar(stat="count",fill=cols,color='black') + scale_fill_manual(values = as.character(cols))
  p <- p + ylab('Number of Cells Passing Filter') + xlab('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + ggtitle(paste0('Total = ',nrow(plot.data))))
  dev.off()
}

FigureS1C <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  p <- ggplot(plot.data, aes(x=reps,y=feature,fill=reps)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = cols) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of UMIs / cell') + ggtitle(paste0('Median = ',median(plot.data$feature)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS1D <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=labels,y=feature,fill=labels)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = as.character(cols)) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of UMIs / cell') + ggtitle('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS1E <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  p <- ggplot(plot.data, aes(x=reps,y=feature,fill=reps)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = cols) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of Genes / cell') + ggtitle(paste0('Median = ',median(plot.data$feature)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS1F <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=labels,y=feature,fill=labels)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = as.character(cols)) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of Genes / cell') + ggtitle('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS1G <- function(object,features,cols,out_f,height,width){
  p <- DotPlot(object=object, features = features) + RotatedAxis()
  p <-p + scale_color_gradientn(colours = cols) + xlab('') + ylab('')
  p <- p + scale_y_discrete(limits=rev(levels(object))) + scale_x_discrete(limits=features)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS1H <- function(object,cols,point.size,alpha=1,anno.size,key.size,out_f,height,width,plot_filled=F,theme,stroke=0.1){
  plot.data <- extractFeatures(object,features='orig.ident',cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  plot.data <- plot.data[sample(row.names(plot.data),nrow(plot.data)),]
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=reps)) + geom_point(size = point.size,alpha=alpha) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = reps), pch = I(21),size = point.size,alpha=alpha,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0,0.95), legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  }
  if(!is.null(theme)){
    p <- p+theme
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS1I <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  p <- ggplot(plot.data, aes(x=reps,y=1,fill=labels)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='cluster',values = as.character(cols))
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS1J <- function(object,out_f){
  plot_list <- list()
  for (cluster in levels(object)){
    marker_genes <- FindMarkers(object,ident.1 = cluster,ident.2 = NULL,only.pos = T)
    res <- enrichGO_wrapper(row.names(marker_genes))
    p <- dotplot(res, showCategory=10,orderBy = "GeneRatio", x='GeneRatio',title=cluster,font.size=10)
    plot_list[[cluster]] <- p
  }
  p1 <- Reduce("+",plot_list[1:2])
  pdf(paste0('plots/figures/',out_f,'_1.pdf'),height=5,width=12)
  print(p1+plot_layout(ncol=2,nrow = 1,byrow = F))
  dev.off()
  p1 <- Reduce("+",plot_list[3:length(plot_list)])
  pdf(paste0('plots/figures/',out_f,'_2.pdf'),height=12,width=24)
  print(p1+plot_layout(ncol=3,nrow = 3,byrow = T))
  dev.off()
}


#### Figure S2 functions ########


FigureS2C <- function(object,features,out_f,point.size,cols,idents,height,width){
  p <- VlnPlot(object,features = features,pt.size = point.size,idents=idents,cols=cols)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + NoLegend() + ggtitle('') + xlab('') + ylab('Pseudotime') + theme(axis.text.x = element_text(angle=0,hjust = 0.5),title = element_text('')))
  dev.off()
}

###########################
#### Plot Figures #########
###########################

### Figure 1 

Figure1B(object,cols=RNA_cluster_colors,out_f='Figure1B',point.size=2,anno.size=14,key.size=4,height=6,width=5.5,plot_filled=T,theme = theme_border,rows=2,stroke=0.2)
Figure1C(object=object,c('Pax6','Eomes','Tubb3'),out_f='Figure1C',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
## Figure 1E is generated externally using scVelo/scanpy ####
Figure1F(object,pgraph=NULL,feature='pseudotime',out_f='Figure1F',cols=brewer.pal(9,'Purples'),point.size=2,width=6,height=6,plot_filled=T,theme=theme_border,stroke=0.2)
Figure1G(mat=as.matrix(readRDS('results/Monocle3/pd_genes_GAMmat.RDS')),object=object,
         features=c('Fhl1','Fos','Id2','Id4','Tfap2c','Prom1','Rbbp7','Hmga2','Chd7','Lrfn5','Sox2','Hes1','Hes5','Fezf2','Sox5','Eomes','Neurog2','Tle4','Bcl11b','Rnd2','Mapt','Pax6','Dcx','Fabp7','Nes','Neurod1','Neurod2','Satb2','Neurod6','Vim','Glast','Camk2b'),
         anno.size=12,out_f='Figure1G',
         cols=rev(brewer.pal(n = 9, name = "RdYlBu")),
         cluster_cols=RNA_cluster_colors,width=8,height=10)
Figure1H(object,out_f='Figure1H',features=c('Fhl1','Chd7','Lrfn5'),cols=RNA_cluster_colors,point.size=1,anno.size=12,key.size=4,alpha=0.5,width=6,height=4,plot_filled=F,direction='vertical')

### Supplementary Figure 1 ###

FigureS1A(object,out_f='FigureS1A',height=6,width=6)
FigureS1B(object,features='seurat_clusters',cols=rep_colors,out_f='FigureS1B',height=6,width=3)
FigureS1C(object,features='nCount_RNA',cols=rep_colors,out_f='FigureS1C',height=6,width=3)
FigureS1D(object,features='nCount_RNA',cols=RNA_cluster_colors,out_f='FigureS1D',height=6,width=8)
FigureS1E(object,features='nFeature_RNA',cols=rep_colors,out_f='FigureS1E',height=6,width=3)
FigureS1F(object,features='nFeature_RNA',cols=RNA_cluster_colors,out_f='FigureS1F',height=6,width=8)
FigureS1G(object,
          features=c('Hes1','Eomes','Mki67','Neurod1','Cntn2','Satb2','Sox5','Bcl11b','Tle4','Mapt','Reln','Gad2','Ptprc','Pdgfrb'),
          cols=c('grey','red','black'),out_f='FigureS1G',height=5,width=7)
FigureS1H(object,cols=rep_colors,point.size=2,alpha=1,anno.size=14,key.size=6,out_f='FigureS1H',height=6,width=6,plot_filled=T,stroke=0.2,theme=theme_border)
FigureS1I(object,features='seurat_clusters',cols=RNA_cluster_colors,out_f='FigureS1I',height=4,width=4)
FigureS1J(object,out_f='FigureS1J')


### Supplementary Figure 2  ###
  
# Figure S2A
Figure1C(object=object,c('Fbxo32','Ctgf','Cyr61'),out_f='FigureS2A_1',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
Figure1C(object=object,c('Moxd1','Hopx','Fam107a','Mt3'),out_f='FigureS2A_2',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
Figure1C(object=object,c('Cryab','Nr4a1','Foxj1'),out_f='FigureS2A_3',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
Figure1C(object=object,c('Ppp1r17','Penk','Neurog1'),out_f='FigureS2A_4',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
## Figure S2B is generated externally using scVelo ####
FigureS2C(object,features='pseudotime',out_f='FigureS2C',point.size=0,idents = c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3'),cols = RNA_cluster_colors,width=6,height=6)
Figure1C(object=object,c('Hes1','Id4','Hes5','Neurog2','Eomes','Neurod2','Rnd2','Mapt'),out_f='FigureS2D',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
Figure1H(object,out_f='FigureS2E',features=c('Hes1','Id4','Hes5','Neurog2','Eomes','Neurod2','Rnd2','Mapt'),cols=RNA_cluster_colors,point.size=1,anno.size=12,key.size=4,alpha=0.5,width=4,height=4,plot_filled=F)

