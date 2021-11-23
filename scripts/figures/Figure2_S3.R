library(Seurat)
library(Signac)
library(ggplot2)
library(ggpubr)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(dplyr)
library(Hmisc)
library(Matrix)
library(SummarizedExperiment)
require(edgeR)
require(LSD)
library(ggrepel)
library(ggsci)
library(ggpointdensity)
library(grid)
library(ComplexHeatmap)
library(GenomicAlignments)

source('scripts/figures/config.R')
source('scripts/aux_functions.R')
source('scripts/figures/plot_functions.R')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))

atac.object <- readRDS('data/merged_scATAC_integrated_cicero.RDS')
DefaultAssay(atac.object) <- 'MACS2peaks'
which_cluster <- c('NSC','IPC','PN1','PN2')
object.sub <- subset(atac.object, idents=levels(Idents(atac.object))[(levels(Idents(atac.object))%in%which_cluster)])
fragment.path <- "data/fragmentsNorm_filtered.tsv.bgz"
### Load ArchR object ####
require(ArchR)
addArchRThreads(threads = 12)
archr_obj <- loadArchRProject("archr/seurat_atac3/",showLogo = F)

##### Read promoter coordinates from gtf file ####
gtfFile <- 'data/genes.gtf'
tssWindow <- 5000
gtf <- getGeneGTF(gtfFile)
prom.coords <- gtf %>% resize(1,"start") %>% resize(tssWindow * 2 + 1, "center")
names(prom.coords) <- prom.coords$gene_name
tss.coords <- resize(prom.coords,width = 1, "center")
seRNA_all <- readRDS("results/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")
gtf <- unique(gtf[gtf$gene_name%in%row.names(seRNA_all)])
gene.length=gtf$exonLength[match(row.names(seRNA_all),gtf$gene_name)]
uf_rpkm <- edgeR::rpkmByGroup(assay(seRNA_all)[,c(1,2,5,6,9:14)],group=factor(c('NSC','NSC','IPC','IPC','PN1','PN1','PN2','PN2','PN3','PN3'),levels=c('NSC','IPC','PN1','PN2','PN3')),
                              gene.length=gene.length,log=F,prior.count=0)
#####################

Figure2A <- function(object,cols,point.size,anno.size,key.size,out_f,height,width,plot_filled=F,theme=NULL,stroke=0.1,rows=1){
  plot.data <- extractFeatures(object,features='seurat_clusters',cells_toInclude = c('all'),cells_toExclude = 'none')
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=labels)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = labels), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0,(1.03+0.025*rows)),plot.margin = margin(height/10,0.1,0.1,0.1,unit='inches'), legend.box = "horizontal",legend.spacing.x = unit(0.001, 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  if(!is.null(theme)){
    print(p + theme)}
  else {
    print(p)
  }
  dev.off()
}

Figure2C <- function(mat,nVar,object,features,clust_method='ward.D',cols,cluster_cols,anno.size,out_f,width,height,expr=NULL){
  if (!is.null(expr)){
    tf_name <- sapply(strsplit(as.character(row.names(mat)), "\\(|\\:|\\."), function(x) x[[1]])
    expr <- expr[match(tf_name,row.names(expr)),]
    mat <- mat[complete.cases(expr),]
  }
  sub_mat <- mat[row.names(mat)%in%features,]
  mat <- head(mat[order(matrixStats::rowVars(mat), decreasing = TRUE),],nVar)
  mat <- rbind(mat,sub_mat[!row.names(sub_mat)%in%row.names(mat),])
  idents_indx <- Idents(object)
  idents_indx <- idents_indx[match(colnames(mat),names(idents_indx))]
  
  intraMean <- groupMeans(mat, groups = idents_indx, sparse = F, na.rm = TRUE)
  rowClust <- hclust(dist(intraMean,method='euclidean'),method = clust_method)
  mat <- mat[rowClust$order,]
  intraMean <- intraMean[rowClust$order,]
  ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%(features)), labels = row.names(mat)[which(row.names(mat)%in%unique(c(features)))],labels_gp = gpar(fontsize = anno.size)))
  col.list <- cluster_cols
  names(col.list) <- levels(object)
  
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm <- Heatmap(mat,cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 3, right_annotation = ra,col = cols,top_annotation = ha,
                heatmap_legend_param=list(direction = "vertical",title='',legend_height = unit(height/6, "inch")))
  draw(hm, merge_legends=T, ht_gap = unit(height/2, "inch")) 
  dev.off()
  return(intraMean)
}

Figure2D <- function(mat,object,features,clust_method='ward.D',cols,cluster_cols,anno.size,out_f,width,height){
  row.names(mat) <- gsub('NeuroD2','Neurod2',row.names(mat))
  row.names(mat) <- gsub('NPC-','',row.names(mat))
  row.names(mat) <- gsub('CN-','',row.names(mat))
  row.names(mat) <- gsub('-peaks','',row.names(mat))
  row.names(mat) <- gsub('Tbr2','Eomes',row.names(mat))
  row.names(mat) <- gsub('Brn2','Pou3f2',row.names(mat))
  mat <- mat[row.names(mat)%in%features,]
  mat <- mat[match(features,row.names(mat)),]
  idents_indx <- Idents(object)
  idents_indx <- idents_indx[match(colnames(mat),names(idents_indx))]
  
  intraMean <- groupMeans(mat, groups = idents_indx, sparse = F, na.rm = TRUE)
  ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%(features)), labels = row.names(mat)[which(row.names(mat)%in%unique(c(features)))],labels_gp = gpar(fontsize = anno.size)))
  col.list <- cluster_cols
  names(col.list) <- levels(object)
  
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm <- Heatmap(mat,cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10, right_annotation = ra,col = cols,top_annotation = ha,
                heatmap_legend_param=list(direction = "vertical",title='',legend_height = unit(height/6, "inch")))
  draw(hm, merge_legends=T, ht_gap = unit(height/2, "inch")) 
  dev.off()
  return(intraMean)
}

Figure2E <- function(object,features,assay='Prom',out_f,min.cutoff=NA,max.cutoff=NA,cols1,cols2,point.size,height,width,plot_filled=F,theme=NULL,stroke=0.1){
  plot_list1 <- list()
  plot_list2 <- list()
  for (feature_f in features){
    DefaultAssay(object) <- assay
    plot.data <- extractFeatures(object,features=feature_f,min.cutoff=min.cutoff,max.cutoff=max.cutoff,cells_toInclude = c('all'),cells_toExclude = 'none')
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colors=cols1,name='',breaks=c(min(plot.data$feature),max(plot.data$feature)),labels=c(round(min(plot.data$feature)),round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.975),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
      p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
    } else {
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colors=cols1,name='',breaks=c(min(plot.data$feature),max(plot.data$feature)),labels=c(round(min(plot.data$feature)),round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.975),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
      p <- p + guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
    }
    if(!is.null(theme)){
      p <- p + theme_border}
    plot_list1[[feature_f]] <- p + ggtitle(paste0(feature_f,' ',assay))
    
    DefaultAssay(object) <- 'chromvar'
    plot.data <- extractFeatures(object,features=feature_f,cells_toInclude = c('all'),cells_toExclude = 'none')
    cols_lower <- colorpalette(cols2,round(abs(min(plot.data[, 3],na.rm=T)))*2)[1:round(abs(min(plot.data[, 3],na.rm=T)))]
    cols_upper <- colorpalette(cols2,round(abs(max(plot.data[, 3],na.rm=T)))*2)[round(abs(max(plot.data[, 3],na.rm=T))+1):(round(abs(max(plot.data[, 3],na.rm=T)))*2)]
    cols <- c(cols_lower,cols_upper)
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_color_gradientn(colors=cols,name='',breaks=c(min(plot.data$feature),0,max(plot.data$feature)),labels=c(round(min(plot.data$feature)),0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.975),legend.background =element_blank(),legend.key = element_blank())
      p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
    } else {
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colors=cols,name='',breaks=c(min(plot.data$feature),0,max(plot.data$feature)),labels=c(round(min(plot.data$feature)),0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.975),legend.background =element_blank(),legend.key = element_blank())
      p <- p + guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
    }
    if(!is.null(theme)){
      p <- p + theme_border}
    plot_list2[[feature_f]] <- p + ggtitle(paste0(feature_f,' Motif'))
  }
  p1 <- Reduce("/", plot_list1) 
  p2 <- Reduce("/", plot_list2)
  height=length(features)*height
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width*2)
  print(p1|p2)
  dev.off()
}

Figure2F <- function(mat_f="results/scATAC/Unique_Peaks.RDS",which_clusters,tss.coords,distance=5000,scale=T,scaleMax=NULL,object,features,cols,cluster_cols,anno.size,out_f,width,height){
  mat <- readRDS(mat_f)$groupMat
  if (length(which_clusters)!=length(colnames(mat))){
    mat <- mat[,grep(paste0(which_clusters,collapse='|'),colnames(mat))]
  }
  if (scale) {
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
  }
  if (is.numeric(scaleMax)){
    mat[mat > scaleMax] <- scaleMax
    mat[mat < -scaleMax] <- -scaleMax
  }
  n_clusters <- ncol(mat)
  if(is.data.frame(tss.coords)){
    tss.coords <- makeGRangesFromDataFrame(tss.coords[,1:3])
  }
  obj_features <- StringToGRanges(rownames(mat), sep = c("-", "-"))
  overlaps <- distanceToNearest(obj_features,tss.coords,ignore.strand=T)
  mat <- as.data.frame(mat)
  mat$nearestGene <- tss.coords$gene_name[overlaps@to]
  mat$dist <- start(resize(obj_features,width = 1,'center'))[overlaps@from] - start(tss.coords)[overlaps@to]
  ### Invert sign of distance for - strand genes ####
  mat$dist[as.vector(strand(tss.coords)=='-')[overlaps@to]] <- mat$dist[as.vector(strand(tss.coords)=='-')[overlaps@to]]*(-1)

  mat_d <- mat[abs(mat$dist)>distance,]
  mat_p <- mat[abs(mat$dist)<=500,]
  
  mat_d$labels <- paste0(mat_d$nearestGene,':',mat_d$dist)
  mat_p$labels <- mat_p$nearestGene
  
  col.list <- cluster_cols
  names(col.list) <- colnames(mat)[1:length(which_clusters)]

  
  ha = HeatmapAnnotation(cluster = colnames(mat_d[1:n_clusters]) ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ra_d = rowAnnotation(foo = anno_mark(at = which(mat_d$labels%in%features), labels = mat_d$labels[which(mat_d$labels%in%features)],labels_gp = gpar(fontsize = anno.size)))
  ra_p = rowAnnotation(foo = anno_mark(at = which(mat_p$labels%in%features), labels = mat_p$labels[which(mat_p$labels%in%features)],labels_gp = gpar(fontsize = anno.size)))
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm_d <- Heatmap(as.matrix(mat_d[,1:n_clusters]),name='Distal',row_title = paste0('Distal: ',nrow(mat_d)),cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10, right_annotation = ra_d,col = cols,top_annotation = ha,
                height=unit(height*3/4-0.5, "inch"),heatmap_legend_param=list(direction = "vertical",title='',legend_height = unit(height/6, "inch")))
  hm_p <- Heatmap(as.matrix(mat_p[,1:n_clusters]),name='Promoter',row_title = paste0('Promoter: ',nrow(mat_p)),cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10,
                  height=unit(height/4, "inch"),right_annotation = ra_p,col = cols,show_heatmap_legend = F)
  hm <- hm_d %v% hm_p
  draw(hm,main_heatmap='Distal',merge_legends=T, ht_gap = unit(0.1, "inch")) 
  dev.off()
  return(list(Distal=mat_d,Promoter=mat_p))
}

Figure2G <- function(mat='results/scATAC/MACS2peaks_intraVar.RDS',prom.coords,cols,out_f,height,width){
  mat <- readRDS(mat)
  obj_features <- StringToGRanges(rownames(mat), sep = c("-", "-"))
  if(is.data.frame(prom.coords)){
    prom.coords <- makeGRangesFromDataFrame(prom.coords[,1:3])
  }
  overlaps <- findOverlaps(obj_features,prom.coords)
  mat$type <- 'Distal'
  mat$type[overlaps@from] <- 'Promoter'
  mat$type <- factor(mat$type, levels=c('Promoter','Distal'))
  p <- ggplot(mat,aes(x=type,y=vst.variance.standardized,fill=type)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8) + ylim(c(0.7,1.4))
  p <- p + scale_fill_manual(values=cols) + xlab('') + ylab('Standardized Variance') + theme(legend.position = "none")
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + stat_compare_means(comparisons = list(c('Promoter','Distal')),label = "p.format",label.y = 1.3,tip.length = c(0.01,0.05)) )
  dev.off()
}

FigureS3A <- function(object,out_f,height,width){
  object_cpm <- as.data.frame(matrix(NA,nrow=nrow(object),ncol=length(unique(object$orig.ident))))
  colnames(object_cpm) <- unique(object$orig.ident)
  for (s in unique(object$orig.ident)){
    object_sub <- object[,object$orig.ident==s]
    object_cpm[,s] <- cpm(rowSums(object_sub@assays$MACS2peaks@counts),log=T,prior.count = 1)
  }
  p <- ggplot(object_cpm, aes( x = E14_rep1, y = E14_rep2 )) + geom_pointdensity(adjust=1,alpha=1,size=1) + scale_color_gradientn(colours = colorpalette("heat",8),name='',breaks=c(0.001,0.27),labels=c('min','max'))             # scale_color_gradientn(colours=tim.colors(12))
  p <- p + theme_cowplot() + geom_abline(slope = 1,intercept = 0,linetype = 2) + xlab('Rep1 log2(CPM+1)') + ylab('Rep2 log2(CPM+1)') + xlim(c(-5,5)) + ylim(c(-5,5))  
  p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.3, 0.97),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.4), vjust = 8))
  my_grob = grid.text(paste0('r=',round(cor.test(object_cpm[,1],object_cpm[,2])$estimate,3)), x=0.85,  y=0.1, gp=gpar(col="black", fontsize=14, fontface="bold"))
  p1 <- p + annotation_custom(my_grob)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p1)
  dev.off()
}

FigureS3B <- function(files_f,filterFrags,max_filterFrags,filterTSS,max_filterTSS,out_f,height,width){
  plot_list <- list()
  for (file_f in files_f){
    mat <- read.table(file_f)
    mat <- mat[mat$uniqueFrags > 1000&mat$enrichment > 1,]
    p <-  ggplot(mat, aes( x = uniqueFrags, y = enrichment )) + xlim(c(0,150000)) + geom_pointdensity(adjust=1,alpha=1,size=1) + scale_color_gradientn(colours = colorpalette("heat",8),name='',breaks=c(1,1000),labels=c('min','max'))             # scale_color_gradientn(colours=tim.colors(12))
    p<- p + xlab('Unique Fragments') + ylab('TSS Enrichment') + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) 
    p <- p + geom_vline(xintercept=c(filterFrags,max_filterFrags),linetype = 2) + geom_hline(yintercept=c(filterTSS,max_filterTSS),linetype = 2) + ylim(c(1,27))
    p <-  p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 1),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1,-0.4), vjust = 8))
    plot_list[[basename(file_f)]] <- p + guides(color = guide_colourbar(barwidth = 2, barheight = 1))
  }
  p <- Reduce("|", plot_list) 
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS3C <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  p <- ggplot(plot.data, aes(x=reps,fill=labels)) + geom_bar(stat="count",fill=cols,color='black') + scale_fill_manual(values = as.character(cols))
  p <- p + ylab('Number of Cells Passing Filter') + xlab('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + ggtitle(paste0('Total = ',nrow(plot.data))))
  dev.off()
}

FigureS3D <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  p <- ggplot(plot.data, aes(x=reps,y=feature,fill=reps)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = cols) + geom_boxplot(width=0.1, fill="white",outlier.size = 0) + theme(legend.position=c('none')) + xlab('') + ylab('Unique Fragments / cell') + ggtitle(paste0('Median = ',median(plot.data$feature)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS3E <- function(archr_obj,logFile,cols,out_f,height,width){
  df <- plotFragmentSizes(ArchRProj = archr_obj,maxSize = 1000,returnDF = TRUE,logFile = logFile )
  plotDF <- data.frame(df)
  p <- ggplot(plotDF, aes(fragmentSize, fragmentPercent,color=sampleName)) + 
    geom_line(size = 0.75) +
    xlab("Fragment Size (bp)") +
    ylab("Percentage of Fragments") + theme(legend.position=c(0.8, 0.95)) + 
    scale_color_manual(values = cols,name='') +
    scale_y_continuous(limits = c(0, max(plotDF$fragmentPercent)*1.05), expand = c(0,0)) +
    scale_x_continuous(limits = c(min(plotDF$fragmentSize), 500),expand = expansion(mult = c(0,0.05)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS3F <- function(object,cols,point.size,alpha=1,anno.size,key.size,height,width,out_f,plot_filled=F,theme,stroke=0.1){
  plot.data <- extractFeatures(object,features='orig.ident',cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  plot.data <- plot.data[sample(row.names(plot.data),nrow(plot.data)),]
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=reps)) + geom_point(size = point.size,alpha=alpha) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = reps), pch = I(21),size = point.size,alpha=alpha,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0.7,0.95), legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  }
  if(!is.null(theme)){
    p <- p+theme
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS3G <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$reps <- gsub('E14_r','R',plot.data$reps)
  p <- ggplot(plot.data, aes(x=reps,y=1,fill=labels)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='cluster',values = as.character(cols))
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS3H <- function(object,assay,features,cols,height,width,out_f){
  p <- DotPlot(object=object,assay = assay, features = features) + RotatedAxis()
  p <-p + scale_color_gradientn(colours = cols) + xlab('') + ylab('')
  p <- p + scale_y_discrete(limits=rev(levels(object))) + scale_x_discrete(limits=features)
  p$guides$colour$title <- 'Average Accessibility'
  p$guides$size$title <- 'Percent Accessible'
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS3I <- function(object,features,assays=c('Prom','GeneBody','cicero_GA'),out_f,min.cutoff=NA,max.cutoff=NA,cols1,cols2,point.size,height,width,plot_filled=F,theme=NULL,stroke=0.1,direction='horizontal'){
  plot_list <- list()
  for (assay in assays){
    for (feature_f in features){
      DefaultAssay(object) <- assay
      plot.data <- extractFeatures(object,features=feature_f,min.cutoff=min.cutoff,max.cutoff=max.cutoff,cells_toInclude = c('all'),cells_toExclude = 'none')
      if(!plot_filled){
        p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") 
        p <- p + scale_fill_gradientn(colors=cols1,name='',breaks=c(min(plot.data$feature),max(plot.data$feature)),labels=c(round(min(plot.data$feature)),round(max(plot.data$feature))))
        p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.975),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
        p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
      } else {
        p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
        p <- p + scale_fill_gradientn(colors=cols1,name='',breaks=c(min(plot.data$feature),max(plot.data$feature)),labels=c(round(min(plot.data$feature)),round(max(plot.data$feature))))
        p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.975),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
        p <- p + guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
      }
      if(!is.null(theme)){
        p <- p + theme_border}
      plot_list[[paste0(assay,'_',feature_f)]] <- p + ggtitle(paste0(feature_f,' ',assay))
    }
  }
  p <- Reduce("+", plot_list)
  if(direction=='horizontal'){
  p1 <- p + plot_layout(ncol=length(assays),nrow=length(features))
  height=length(features)*height
  width=length(assays)*width
  } else {
    p1 <- p + plot_layout(nrow=length(assays),ncol=length(features))
    width=length(features)*width
    height=length(assays)*height
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width)
  print(p1)
  dev.off()
}

#### Figure 2 ############

Figure2A(object=atac.object,cols=ATAC_cluster_colors,point.size=2,anno.size=14,key.size=4,out_f='Figure2A',height=6,width=5.5,plot_filled=T,theme = theme_border,stroke=0.2)
plotMisha(object=atac.object,targetGene='Dll1',out_f='Figure2B',upstream=10e4,downstream=8e4,plotClusters=c('NSC','IPC','PN1','PN2','CR','IN','MG+Mural'),
          chipTracksToExtract=c('scATAC.E14.NSC','scATAC.E14.IPC','scATAC.E14.PN1','scATAC.E14.PN2','scATAC.E14.CR','scATAC.E14.IN','scATAC.E14.MG_Mural'),
          cluster_cols=ATAC_cluster_colors,chipColors=ATAC_cluster_colors,sample_cells=1000,chipNames=c('NSC','IPC','PN1','PN2','CR','IN','MG+Mural'),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.4, domains=0.15, genes=0.7, arcs=0.7, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,axis=FALSE,scATAC=TRUE, genes=TRUE,ideogram=TRUE))
mat <- Figure2C(mat=as.matrix(object.sub@assays$chromvar@data[,order(as.numeric(Idents(object.sub)))]),
         nVar=100,object=object.sub,expr = uf_rpkm[rowMaxs(uf_rpkm,na.rm=T)>=1,],
         features=c('Sox2','Fos::jun(var.2)','Lhx2','Id4','Tead2','Tbr2','Pax6','Eomes','Neurog2(var.2)','Neurod2','Neurod6','Mef2c','Tbr1','Fezf2','Satb2'),
         anno.size=12,clust_method='ward.D',out_f='Figure2C',
         cols=colorRamp2(seq(-10,10,length.out=length(heatmap_colors)), heatmap_colors),
         cluster_cols=ATAC_cluster_colors[1:4],width=8,height=10)
mat <- Figure2D(mat=as.matrix(object.sub@assays$chip_chromvar@data[,order(as.numeric(Idents(object.sub)))]),
                object=object.sub,
                features=c('Sox2','Pax6','Neurog2','Eomes','Neurod2','Tbr1'),
                anno.size=12,clust_method='ward.D',out_f='Figure2D',
                cols=colorRamp2(seq(-10,10,length.out=length(heatmap_colors)), heatmap_colors),
                cluster_cols=ATAC_cluster_colors[1:4],width=6,height=3)
Figure2E(object=atac.object,assay = 'GeneBody',out_f='Figure2E',features=c('Sox2','Eomes','Neurod2'),cols1=rev(colorpalette('ylgnbu')),cols2='matlablike',min.cutoff='q5',max.cutoff='q95',point.size=1,width=4,height=4,plot_filled=T,stroke=0.1,theme = theme_border)
res <- Figure2F(mat_f="results/scATAC/Unique_Peaks.RDS",which_clusters = levels(atac.object),
         tss.coords=tss.coords,scale=T,scaleMax=2,
         features=c('Gad2','Gad2:21900','Hes5','Eomes','Hes5:-15937','Neurog2','Neurog2:13176','Neurog2:72610','Neurog2:65117','Neurod2','Neurod2:59306','Bcl11b','Bcl11b:105239','Reln','Reln:-7935','Pdgfrb:13405','Pdgfrb'),
         cols=colorRamp2(seq(-2,2,length.out=length(heatmap_colors)), heatmap_colors),
         cluster_cols=ATAC_cluster_colors,anno.size=12,out_f='Figure2F',width=8,height=12)
Figure2G(mat='results/scATAC/MACS2peaks_intraVar.RDS',prom.coords=prom.coords,cols=pal_npg('nrc')(2),out_f='Figure2G',height=3.5,width=3.5)
plotMisha(object=atac.object,targetGene='Fhl1',out_f='Figure2H',upstream=1e5,downstream=1e5,
          plotClusters=c('NSC'),chipNames=c('NSC','IPC','PN1','PN2','CR','IN'),
          chipTracksToExtract=c('scATAC.E14.NSC','scATAC.E14.IPC','scATAC.E14.PN1','scATAC.E14.PN2','scATAC.E14.CR','scATAC.E14.IN'),
          cluster_cols=ATAC_cluster_colors[1],sample_cells=1000,chipColors=ATAC_cluster_colors,chipNames=c('NSC','IPC','PN1','PN2','CR','IN','MG+Mural'),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=TRUE, arcs=FALSE, rna=FALSE, chip=TRUE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.6, domains=0.15, genes=0.7, arcs=0.7, axis=0.4,ideogram=0.2,scATAC=2,anno=0.15))

#### Figure S3 ############

FigureS3A(object=atac.object,out_f='FigureS3A',height=6,width=6)
FigureS3B(files_f=list.files(path = "results/scATAC/",pattern = "_FilterCells.txt",full.names = T),filterFrags=10000,max_filterFrags=120000,filterTSS=9,max_filterTSS=25,out_f='FigureS3B',height=4,width=8)
FigureS3C(object=atac.object,features='seurat_clusters',cols=rep_colors,out_f='FigureS3C',height=6,width=3)
FigureS3D(object=atac.object,features='passed_filters',cols=rep_colors,out_f='FigureS3D',height=6,width=3)
FigureS3E(archr_obj=archr_obj,logFile='./archr/ArchRLogs/FragmentSize.log',cols=c('darkred','darkblue'),out_f='FigureS3E',height=4,width=6) 
FigureS3F(object=atac.object,cols=rep_colors,point.size=2,alpha=1,anno.size=14,key.size=6,height=6,width=6,plot_filled=T,stroke=0.2,out_f='FigureS3F',theme=theme_border)
FigureS3G(object=atac.object,features='seurat_clusters',cols=ATAC_cluster_colors,out_f='FigureS3G',height=4,width=4)
FigureS3H(object=atac.object,assay='GeneBody',out_f='FigureS3H',
          features=c('Hes1','Sox2','Eomes','Neurod2','Satb2','Sox5','Lrfn5','Mapt','Reln','Gad2','Pdgfrb'),
          cols=c('grey','red','black'),height=5,width=7)
FigureS3I(object=atac.object,out_f='FigureS3I',assays='GeneBody',features=c('Pax6','Eomes','Tubb3'),cols1=rev(colorpalette('ylgnbu')),cols2='matlablike',min.cutoff=NA,max.cutoff=NA,point.size=1,width=4,height=4,plot_filled=T,stroke=0.2,theme = theme_border)
FigureS3H(object=atac.object,assay='Prom',out_f='FigureS3J',
          features=c('Hes1','Sox2','Eomes','Neurod2','Satb2','Sox5','Lrfn5','Mapt','Reln','Gad2','Pdgfrb'),
          cols=c('grey','red','black'),height=5,width=7)

