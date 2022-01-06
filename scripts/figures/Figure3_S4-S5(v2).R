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
library(ggpointdensity)
library(grid)
library(Rcpp)
library(ComplexHeatmap)
library(Rcpp)
library(htmlwidgets)
require(webshot)
library(monaLisa)
library(TFBSTools)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(GenomicRanges)
## Set main folder as current dir ####

source('scripts/figures/config.R')
source('scripts/aux_functions.R')
source('scripts/figures/plot_functions.R')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))

p2glinks <- readRDS('results/P2G-Links.RDS')
names(p2glinks)[3:6] <- c('posCor','all','negCor','noCor')
rna.object <- readRDS('data/merged_scRNA_unfiltered_IDs.RDS')
atac.object <- readRDS('data/merged_scATAC_integrated_cicero.RDS')
DefaultAssay(atac.object) <- 'MACS2peaks'
coembed.object <- readRDS('data/coembed_scRNA_scATAC.RDS')
Idents(coembed.object) <- coembed.object$celltype
levels(coembed.object) <- c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3','CR','IN','MG','Mural')

seRNA_all <- readRDS("results/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")
sePB_all <- readRDS("results/scATAC/Cluster_PseudoBulk-Summarized-Experiment.RDS")

uf_rna <- suppressMessages(uniqueFeatures(
  edgeR::cpm(assay(seRNA_all),log=TRUE,prior.count=1),
  groups = colData(seRNA_all)$Group,
  padj = 1,
  minSdRatio = 0,
  minLFC = 0,
  zCutoff = 0,
  breakPt = "last",
  groupMin = 20,
  maxGroupSize = 2,
  clusterCols = F
))
expr <- uf_rna$groupMat

Figure3A <- function(selected_peaks,scaleMax=NULL,seRNA_all,sePB_all,out_f,features,padj=1,minLFC=0,minSdRatio=0,zCutoff=0,which_clusters=unique(sePB_all$Group),cluster_cols1,cluster_cols2,which_clusters_rna=unique(seRNA_all$Group),scale=T,cols1,cols2,anno.size=12,height,width){
  seRNA <- seRNA_all[row.names(seRNA_all)%in%selected_peaks$gene_name,]
  sePB <- sePB_all[match(gsub('_','-',selected_peaks$peakName),row.names(sePB_all)),]
  uf <- suppressMessages(uniqueFeatures(
    edgeR::cpm(assay(sePB),log=TRUE,prior.count=1),
    groups = colData(sePB)$Group,
    padj = padj,
    minSdRatio = minSdRatio,
    minLFC = minLFC,
    zCutoff = zCutoff,
    breakPt = "last",
    groupMin = 20,
    maxGroupSize = 2,
    clusterCols = F
  ))
  uf_rna <- suppressMessages(uniqueFeatures(
    edgeR::cpm(assay(seRNA),log=TRUE,prior.count=1),
    groups = colData(seRNA)$Group,
    padj = padj,
    minSdRatio = minSdRatio,
    minLFC = minLFC,
    zCutoff = zCutoff,
    breakPt = "last",
    groupMin = 0,
    maxGroupSize = 2,
    clusterCols = F
  ))
  mat <- uf$groupMat
  indx <- rep(TRUE,nrow(mat))
  if (length(unique(which_clusters))!=length(unique(colnames(mat)))){
    mat <- mat[,colnames(mat)%in%which_clusters]
    bin_mat <- uf$binaryMat
    indx <- rowMaxs(bin_mat[,colnames(bin_mat)%in%which_clusters])>0
    mat <- mat[indx,]
  }
  if (scale) {
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
  }
  if (is.numeric(scaleMax)){
    mat[mat > scaleMax] <- scaleMax
    mat[mat < -scaleMax] <- -scaleMax
  }
  n_clusters <- ncol(mat)
  mat_r <- uf_rna$groupMat
  if (length(unique(which_clusters_rna))!=length(unique(colnames(mat_r)))){
    mat_r <- mat_r[,colnames(mat_r)%in%which_clusters_rna]
  }
  if (scale) {
    mat_r <- sweep(mat_r - rowMeans(mat_r), 1, matrixStats::rowSds(mat_r), `/`)
  }
  if (is.numeric(scaleMax)){
    mat_r[mat_r > scaleMax] <- scaleMax
    mat_r[mat_r < -scaleMax] <- -scaleMax
  }
  n_clusters_r <- ncol(mat_r)
  linksSig_IDs <- as.data.frame(cbind(as.vector(seqnames(selected_peaks)),as.numeric(start(selected_peaks)),as.numeric(end(selected_peaks ))))
  selected_peaks$IDs <- unite(linksSig_IDs,col='peaks_ID',sep='-')[,1]
  anno <- as.data.frame(mcols(selected_peaks),stringsAsFactors=FALSE)
  row.names(anno) <- paste0(anno$gene_name,':',anno$distance)
  anno <- anno[uf$rowOrder,]
  anno_orig <- anno
  anno <- anno[indx,]
  row.names(mat) <- paste0(anno$gene_name,':',anno$distance)
  comb_mat <- mat_r[match(anno$gene_name,row.names(mat_r)),]
  col.list1 <- cluster_cols1[1:n_clusters]
  names(col.list1) <- colnames(mat) 
  col.list2 <- cluster_cols2[1:n_clusters_r]
  names(col.list2) <- colnames(comb_mat)
  ha1 = HeatmapAnnotation(ATAC_clusters = colnames(mat[,1:n_clusters]) ,col = list(ATAC_clusters=col.list1),show_legend = T,show_annotation_name=F,annotation_legend_param=list(ATAC_clusters=list(direction='vertical')))
  ha2 = HeatmapAnnotation(RNA_clusters = colnames(comb_mat[,1:n_clusters_r]) ,col = list(RNA_clusters=col.list2),show_legend = T,show_annotation_name=F,annotation_legend_param=list(RNA_clusters=list(direction='vertical')))
  ha1 = HeatmapAnnotation(foo = anno_mark(at = 1:n_clusters, labels = colnames(mat[,1:n_clusters]),labels_gp = gpar(fontsize = 12),labels_rot = 0),ATAC_clusters=colnames(mat[,1:n_clusters]),col = list(ATAC_clusters=col.list1),show_legend = F,show_annotation_name=F,annotation_legend_param=list(ATAC_clusters=list(direction='vertical')))
  ha2 = HeatmapAnnotation(foo = anno_mark(at = 1:n_clusters_r, labels = colnames(comb_mat[,1:n_clusters_r]),labels_gp = gpar(fontsize = 12),labels_rot = 0),RNA_clusters=colnames(comb_mat[,1:n_clusters_r]),col = list(RNA_clusters=col.list2),show_legend = F,show_annotation_name=F,annotation_legend_param=list(RNA_clusters=list(direction='vertical')))
  la1 <- columnAnnotation(foo = anno_text(colnames(mat), location = 0.5,rot=0, just = "center",gp = gpar(fill = col.list1, col = "black", border = "black", fontsize = 12,height = unit(height/50,'inch'))))
  la2 <- columnAnnotation(foo = anno_text(colnames(comb_mat), location = 0.5,rot=0, just = "center",gp = gpar(fill = col.list2, col = "black", border = "black", fontsize = 12),height = unit(height/50,'inch')))
  if(!is.null(features)){
    ra1 = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%features),side='left', labels = row.names(mat)[which(row.names(mat)%in%features)],labels_gp = gpar(fontsize = anno.size)))
  } else {
    ra1 = NULL
  }
  hm1 <- Heatmap(mat,name='Accessibility',column_title = 'Accessibility',cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10, left_annotation = ra1,col = cols1,top_annotation = la1,
                 heatmap_legend_param=list(direction = "horizontal",title='Accessibility Z-score',legend_width = unit(width/6, "inch")))
  hm2 <- Heatmap(comb_mat,name='Expression',column_title = 'Expression',cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10,col = cols2,top_annotation = la2,
                 heatmap_legend_param=list(direction = "horizontal",title='Expression Z-score',legend_width = unit(width/6, "inch")))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  draw(hm1+hm2, merge_legends=F, ht_gap = unit(0.2, "inch"),heatmap_legend_side='bottom',annotation_legend_side = "bottom") 
  dev.off()
  anno_orig$peakChr <- as.character(seqnames(selected_peaks[uf$rowOrder]))
  anno_orig$peakCenter <- as.numeric(start(selected_peaks[uf$rowOrder])+end(selected_peaks[uf$rowOrder]))/2
  return(cbind(uf$binaryMat,anno_orig,uf$groupMat))
}

Figure3C_1 <- function(mat_f,out_f,which_clusters=NULL,expr=NULL,background_others=F,object,features=NULL,cols=c('green','grey80','red'),logFC=log2(1.5),logP=5,height=4,width=4,point.size=2,anno.size=8,plot_motifs=F,n_cols=NULL){
  plot_list <- list()
  motif_list <- list()
  for (which_cluster in which_clusters){
    if(!is.null(which_cluster)){
      mat <- mat_f[mat_f[,which_cluster]==1,]
    }
    if(background_others){
      background <- gsub('_','-',mat_f[mat_f[,which_cluster]!=1,'IDs'])
    } else {
      background <- length(mat$IDs)
    }
    enriched.motifs <- FindMotifs_fisher(
      object = object,
      features = gsub('_','-',mat$IDs),
      background=10000
    )
    enriched.motifs$motif.name <- capitalize(tolower(as.vector(enriched.motifs$motif.name)))
    enriched.motifs$fold.enrichment <- log2(enriched.motifs$fold.enrichment)
    enriched.motifs$pvalue <- -log10(enriched.motifs$pvalue)
    enriched.motifs$col <- 'grey80'
    enriched.motifs$col[enriched.motifs$fold.enrichment>logFC&enriched.motifs$pvalue>=logP] <- 'red'
    enriched.motifs$col[enriched.motifs$fold.enrichment<(-logFC)&enriched.motifs$pvalue>=logP] <- 'green'
    enriched.motifs <- enriched.motifs[order(abs(enriched.motifs$fold.enrichment),enriched.motifs$pvalue,decreasing = T),]
    if (!is.null(expr)){
      tf_name <- sapply(strsplit(as.character(enriched.motifs$motif.name), "\\(|\\:|\\."), function(x) x[[1]])
      expr <- expr[match(tf_name,row.names(expr)),]
      enriched.motifs <- enriched.motifs[complete.cases(expr),]
    }
    p <- ggplot(enriched.motifs,aes(x=fold.enrichment,y=pvalue)) + geom_point(aes(fill = col), pch = I(21),size = point.size) 
    p <- p + scale_fill_manual(values = cols,labels=NULL) + xlab(expression(Log[2]~Fold~Change)) + ylab(expression(-Log[10]~(P)))
    if (is.null(features)){
      sel_df <- rbind(head(enriched.motifs[order(enriched.motifs$fold.enrichment,enriched.motifs$pvalue,decreasing=T),],10),head(enriched.motifs[order(enriched.motifs$fold.enrichment*(-1),enriched.motifs$pvalue,decreasing=T),],10))
      features <- sel_df$motif.name
    }
    p <- p + ggrepel::geom_text_repel(
      data = enriched.motifs, size = anno.size,seed = 42,
      box.padding =0.8, min.segment.length = 0,max.iter = 10000,
      aes(x=fold.enrichment,y=pvalue,color=NULL,label=ifelse(motif.name%in%features, as.character(motif.name), "")),force=10)
    enriched.motifs <- enriched.motifs[!is.infinite(enriched.motifs$fold.enrichment),]
    p <- p + theme(legend.position = "none") + ggtitle(which_cluster) + xlim(c(-max(abs(enriched.motifs$fold.enrichment),na.rm=T),max(abs(enriched.motifs$fold.enrichment),na.rm=T)))
    plot_list[[which_cluster]] <- p
    motif_list[[which_cluster]] <- enriched.motifs
  }
  p <- Reduce('+',plot_list)
  if (plot_motifs){
    plot_motifs <- enriched.motifs$motif[match(features,enriched.motifs$motif.name)]
    plot_motifs <- plot_motifs[!is.na(plot_motifs)]
    p1 <- MotifPlot(
      object = object,
      motifs = plot_motifs,
      use.names=T,
      ncol=round(length(features)/6)
    )
    p1 <- p1 +theme_cowplot(font_size = 12) + theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.margin = margin(0,0,0,0),
      axis.line=element_blank())
    widths=c(rep(4,length(which_clusters)),2)
    heights=c(rep(4,length(which_clusters)),5)
    p <- p+ p1 + plot_layout(ncol=length(which_clusters)+1,widths=widths,heights=4)
    width=width*length(which_clusters)+width/2
  } else {
    if(is.null(n_cols)){
      n_cols <- length(which_clusters)
    }
    width=width*n_cols
    height=height*(round(length(which_clusters)/n_cols))
    p <- p+ plot_layout(ncol=n_cols,byrow = T,nrow = round(length(which_clusters)/n_cols))
  }
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width,useDingbats=FALSE)
  print(p)
  dev.off()
  return(motif_list)
}

Figure3C_2 <- function(object,out_f,features,n_cols,height,width){
  motifs <- as.vector(GetMotifObject(object)@motif.names)
  motifs <- data.frame(names=capitalize(tolower(motifs)),motifs=names(motifs))
  plot_motifs <- motifs[match(features,motifs$names),]
  plot_motifs <- plot_motifs[!is.na(plot_motifs$motifs),]
  p1 <- MotifPlot(
    object = object,
    motifs = as.character(plot_motifs$motifs),
    use.names=T,
    ncol=n_cols)
  p1 <- p1 +theme_cowplot(font_size = 12) + theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line=element_blank())
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p1)
  dev.off()
}

Figure3D <- function(mat,mat_f,mat_rna,anno_mat,which_clusters=NULL,object=coembed.object,features,cols1,cols2,cluster_cols,anno.size,width,height){
  mat$labels <- paste0(mat$gene_name,':',mat$distance)
  mat_f <- mat_f[row.names(mat_f)%in%(gsub('_','-',mat$peakName)),]
  pd <- object$pseudotime
  row.names(mat_f) <- mat[match(row.names(mat_f),gsub('_','-',mat$peakName))]$labels
  mat_rna <- mat_rna[match(mat$gene_name[match(row.names(mat_f),mat$labels)],row.names(mat_rna)),]
  mat_rna <- mat_rna[,colnames(mat_f)]
  o <- data.frame(feature=row.names(mat_f),rna=row.names(mat_rna))
  o$cor <- rowCorCpp(1:nrow(o),1:nrow(o), as.matrix(mat_f), as.matrix(mat_rna))
  o$maxPD_f <- pd[match(colnames(mat_f)[max.col(mat_f)],names(pd))]
  o$maxPD_rna <- pd[match(colnames(mat_rna)[max.col(mat_rna)],names(pd))]
  o$PDdiff <- o$maxPD_f-o$maxPD_rna               #Pseudotime difference peak(Enhancer Accessibility - Gene expression) 
  col_levels <- levels(rna.object)[1:7]
  col.list <- cluster_cols[1:length(col_levels)]
  names(col.list) <- col_levels
  idents_indx <- factor(object$celltype,levels=col_levels)
  idents_indx <- idents_indx[match(colnames(mat_f),names(idents_indx))]
  anno_mat$labels <- paste0(anno_mat$gene_name,':',anno_mat$distance)
  anno_mat <- anno_mat[match(row.names(mat_f),anno_mat$labels),]
  if (is.null(which_clusters)){
    which_clusters <- colnames(anno_mat)[1:7]
  }
  anno_mat <- anno_mat[rowMaxs(as.matrix(anno_mat[,which_clusters]),na.rm=T)>0,]
  mat_f <- mat_f[row.names(mat_f)%in%anno_mat$labels,]
  o <- o[o$feature%in%anno_mat$labels,]
  colnames(mat_f) <- as.numeric(idents_indx)
  mat2 <- bin.matrix.rows(mat_f,bin.size = 25)
  mat2_idents <- levels(idents_indx)[as.numeric(colnames(mat2))]
  o$subgroup <- mat2_idents[max.col(mat2)]
  mat_f <- mat_f[!is.na(o$subgroup),]
  o <- o[!is.na(o$subgroup),]
  o$subgroup <- gsub('_M','',o$subgroup)
  o$subgroup[1:4500] <- 'NSC'
  o$subgroup[4501:5200] <- 'IPC'  
  o$subgroup[5201:6500] <- 'PN1' 
  o$subgroup[6501:nrow(o)] <- 'PN2' 
  o$subgroup <- factor(o$subgroup, levels=c('NSC','IPC','PN1','PN2'))
  panel_fun = function(index, nm) {
    pushViewport(viewport(xscale = c(-30,30), yscale = c(0, 2)))
    grid.rect()
    grid.xaxis(gp = gpar(fontsize = 8))
    grid.boxplot(o$PDdiff[index], pos = 1, direction = "horizontal",outline=F,gp = gpar(fill = "grey"))
    p_diff <- wilcox.test(o$PDdiff[index],mu = 0, alternative = "less")$p.value
    if(p_diff<0.001){p_diff <- scientific(p_diff, digits = 2)}
    if(as.numeric(p_diff)>0.05){p_diff <- ''} else {p_diff <- paste0('; p=',p_diff)}
    grid.text(paste0('M: ',round(median(o$PDdiff[index],na.rm=T),2),p_diff), 0, y = 1.9,
              just = "top", default.units = "native", gp = gpar(fontsize = 12))
    popViewport()
  }
  anno_box = anno_zoom(align_to = o$subgroup, which = "row", panel_fun = panel_fun, 
                       size = unit(height/12, "inch"), gap = unit(0.4, "inch"), width = unit(width/6, "inch"))
  
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat_f)%in%features),side='left', labels = row.names(mat_f)[row.names(mat_f)%in%features],labels_gp = gpar(fontsize = anno.size)))
  pdf('plots/figures/Figure3D.pdf',height=height,width=width)
  hm1 <- Heatmap(mat_f, name = "Pseudotime",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, left_annotation = ra,col = cols1,top_annotation = ha,raster_quality = 6,
                 heatmap_legend_param=list(labels = c(0,100),at=c(0,1), direction = "vertical",title='% Max',legend_height = unit(height/6, "inch")),width=unit(5.5*width/10, "inch"))
  hm2 <- Heatmap(o$PDdiff, name = "PDdiff",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE,col=cols2,right_annotation = rowAnnotation(foo = anno_box),
                 heatmap_legend_param=list(labels = c(-5,0,5),at=c(-5,0,5), direction = "vertical",title='dPD',legend_height = unit(height/6, "inch")),width=unit(0.6*width/10, "inch"))
  
  draw(hm1+hm2, merge_legends=T, ht_gap = unit(0.2, "inch")) 
  dev.off()
  return(o)
}

Figure3E <- function(PWM,scRNA_expression,posCOR_pairs,pseudotime,FDR,out_f,expression_cutoff,height,width)
{
  pwms <-readRDS(PWM)
  pwms<- toPWM(pwms, type=c("log2probratio", "prob"), pseudocounts=0.8,bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  pwms<-pwms[!duplicated(names(pwms))]
  expression<-scRNA_expression[,c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3')] 
  expression=expression[rowMaxs(expression)>=expression_cutoff,] 
  rownames(expression)<-toupper(rownames(expression))
  pos_cor_pairs<-read.table(posCOR_pairs, sep='\t', header=T)
  pos_cor_pairs<-subset(pos_cor_pairs, pos_cor_pairs$IPC ==1 | pos_cor_pairs$PN1 ==1 ) #select IPC and PN enhancer
  pos_cor_pairs<-pos_cor_pairs[,c("gene_chr","gene_start","gene_strand","gene_name","gene_id","peakName","labels")]
  pseudo_t<-read.table(pseudotime, sep='\t',header=T)
  pos_cor_pairs<-na.omit(merge(pos_cor_pairs,pseudo_t[,c(1,6)], by.x='labels', by.y='feature', all.x=T))
  bins <- bin(x =pos_cor_pairs$PDdiff, binmode = "equalWidth", minAbsX=2, nBins=3)
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, StringToGRanges(pos_cor_pairs$peakName,sep = c("_","_")))
  se <- calcBinnedMotifEnr(seqs = seqs, bins = bins, motifs = pwms,verbose = T,method = "R")
  sel <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR
  seSel <- se[sel, ]
  #subsetting for expressed TFs 
  expression<-scRNA_expression[,c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3')] 
  expression=expression[rowMaxs(expression)>=expression_cutoff,] 
  rownames(expression)<-toupper(rownames(expression))
  seSel@elementMetadata$gene<-sapply(strsplit(as.character( seSel@elementMetadata$motif.name), "\\(|\\:|\\."), function(x) x[[1]])
  seSel<-seSel[rowData(seSel)$gene %in% row.names(expression),]
  plotMotifHeatmaps(x = seSel, which.plots = c("log2enr","negLog10Padj","negLog10P"), width = 2.0,cluster = TRUE, maxEnr = 2, maxSig = 10, show_motif_GC = TRUE)
  SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm)
  range(SimMatSel)
  hcl <- hclust(dist(assay(seSel,"log2enr")))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  plotMotifHeatmaps(x = seSel, which.plots = c("log2enr"), width = 2.5,
                    cluster = hcl, maxEnr = 2, maxSig = 10,
                    show_dendrogram = TRUE, show_seqlogo = TRUE,
                    width.seqlogo = 1.5)
  dev.off()
}

Figure3F <- function(out_f,PWM,expression_cutoff,pseudotime,motif_name,chip_f,posCOR_pairs,FDR_cutoff,width,height)
{
  ####filter PWMs
  expr<-expr[,c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3')]
  expr=expr[rowMaxs(expr)>=expression_cutoff,]
  rownames(expr)<-toupper(rownames(expr))
  pwms <-readRDS(PWM)
  pwms@metadata$gene_name<-sapply(strsplit(as.character(name(pwms)), "\\(|\\:|\\."), function(x) x[[1]])
  pwms<-pwms[pwms@metadata$gene_name %in% rownames(expr)]
  pwms<- toPWM(pwms, type=c("log2probratio", "prob"), pseudocounts=0.8,bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  pwms<-pwms[!duplicated(names(pwms))]
  pos_cor_pairs<-read.table(posCOR_pairs, sep='\t', header=T)
  pos_cor_pairs<-subset(pos_cor_pairs, pos_cor_pairs$IPC ==1 | pos_cor_pairs$PN1 ==1 ) #select IPC and PN! enhancer
  pos_cor_pairs<-pos_cor_pairs[,c("gene_chr","gene_start","gene_strand","gene_name","gene_id","peakName","labels")]
  res <- extract_MotifPos(pwm=readRDS(PWM),chip_f=chip_f,motif.name=motif_name)
  colnames(res)[1:3] <- c('chrom','start','end')
  chip <- makeGRangesFromDataFrame(res[,1:3])
  chip <- GenomicRanges::resize(chip,501,fix='center')
  chip_features <- makeGRangesFromDataFrame(chip)
  enh_ranges <- StringToGRanges(pos_cor_pairs$peakName, sep = c("_", "_"))
  prom_ranges <- GRanges(seqnames=pos_cor_pairs$gene_chr,IRanges(start=pos_cor_pairs$gene_start-5000,end=pos_cor_pairs$gene_start+5000))
  prom_overlaps <- GenomicRanges::countOverlaps(prom_ranges,chip_features)
  pos_cor_pairs$prom_ChIP <- prom_overlaps
  overlaps <- GenomicRanges::countOverlaps(enh_ranges,chip_features)
  pos_cor_pairs$enh_ChIP <- overlaps
  pseudo_t<-read.table(pseudotime, sep='\t',header=T)
  pos_cor_pairs<-merge(pos_cor_pairs,pseudo_t[,c(1,6)], by.x='labels', by.y='feature', all.x=T)
  subgroup<-subset(pos_cor_pairs, pos_cor_pairs$prom_ChIP==1 | pos_cor_pairs$enh_ChIP==1)
  bins <- bin(x =subgroup$PDdiff, binmode = "equalWidth", minAbsX=2, nBins=3)
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, StringToGRanges(subgroup$peakName,sep = c("_","_")))
  se <- calcBinnedMotifEnr(seqs = seqs, bins = bins, motifs = pwms,verbose = T,method = "R")
  sel <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR_cutoff
  seSel <- se[sel, ]
  SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm)
  #range(SimMatSel)
  hcl <- hclust(dist(assay(seSel,"log2enr")))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  plotMotifHeatmaps(x = seSel, which.plots = c("log2enr"), width = 2.5,
                    cluster = hcl, maxEnr = 2, maxSig = 10,
                    show_dendrogram = TRUE, show_seqlogo = TRUE,
                    width.seqlogo = 1.5)
  dev.off()
}


Figure3G <- function(object,mat,out_f,enhancers=NULL,features,motifs,cols,height,width,cluster_cols=cluster_colors,anno.size=12,direction='vertical'){
  mat$labels <- paste0(mat$gene_name,':',mat$distance)
  res_list <- list()
  for (i in seq_along(features)){
    DefaultAssay(object) <- 'MACS2peaks'
    distal <- as.character(mat[mat$gene_name==features[i],'IDs'])
    plot.data_d <- extractFeatures(object,features=c(distal,'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    plot.data_d <- plot.data_d[order(plot.data_d$pseudotime),-c(1,2,ncol(plot.data_d)-1,ncol(plot.data_d))]
    plot.data_d$agg <- rowSums(plot.data_d[,1:length(distal)],na.rm=T)
    colnames(plot.data_d)[1:length(distal)] <- paste0('enh',1:length(distal))
    plot.data_d <- plot.data_d[complete.cases(plot.data_d),]
    for (s in 1:length(distal)){
      plot.data_d[,s] <- scales::rescale(predict(mgcv::gam(formula = as.formula(paste0(colnames(plot.data_d)[s],' ~ s(pseudotime, bs = "cs")')),method = "REML", data=plot.data_d)),to=c(0,100))
    }
    plot.data_d$agg <- scales::rescale(predict(mgcv::gam(formula = agg ~ s(pseudotime, bs = "cs"),method = "REML", data=plot.data_d)),to=c(0,100))
    colnames(plot.data_d)[1:length(distal)] <- as.character(mat[mat$gene_name==features[i],'labels'])
    plot.data_d <- t(plot.data_d)
    
    plot.data_agg <- t(as.matrix(plot.data_d['agg',]))
    colnames(plot.data_agg) <- colnames(plot.data_d)
    row.names(plot.data_agg) <- 'aggregated'
    plot.data_d <- plot.data_d[-c((nrow(plot.data_d)-1):nrow(plot.data_d)),]
    
    DefaultAssay(object) <- 'RNA'
    plot.data_e <- extractFeatures(object,features=c(features[i],'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    plot.data_e <- plot.data_e[match(colnames(plot.data_d),row.names(plot.data_e)),]
    plot.data_e$rescaled <- scales::rescale(predict(mgcv::gam(formula = feature ~ s(pseudotime, bs = "cs"),method = "REML", data=plot.data_e)),to=c(0,100))
    plot.data_e <- t(as.matrix(plot.data_e$rescaled,nrow=1))
    colnames(plot.data_e) <- colnames(plot.data_d)
    row.names(plot.data_e) <- 'expression'
    DefaultAssay(object) <- 'chromvar'
    plot.data_m <- extractFeatures(object,features=c(motifs[i],'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    plot.data_m <- plot.data_m[match(colnames(plot.data_d),row.names(plot.data_m)),]
    plot.data_m$rescaled <- scales::rescale(predict(mgcv::gam(formula = feature ~ s(pseudotime, bs = "cs"),method = "REML", data=plot.data_m)),to=c(0,100))
    plot.data_m <- t(as.matrix(plot.data_m$rescaled,nrow=1))
    colnames(plot.data_m) <- colnames(plot.data_d)
    row.names(plot.data_m) <- 'motif'
    ## Create heatmap
    if(is.null(enhancers)){
      enhancers_anno <- row.names(plot.data_d)
    } else {
      enhancers_anno <- enhancers
    }
    col_levels <- levels(object)[1:7]
    col.list <- cluster_cols[1:length(col_levels)]
    names(col.list) <- col_levels
    idents_indx <- factor(c(as.character(Idents(rna.object)),as.character(object$celltype[!is.na(object$celltype)])),levels=col_levels)
    names(idents_indx) <- c(names(Idents(rna.object)),names(object$celltype[!is.na(object$celltype)]))
    idents_indx <- idents_indx[match(colnames(plot.data_d),names(idents_indx))]
    ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = ifelse(i==length(features),TRUE,FALSE),show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
    ra = rowAnnotation(foo = anno_mark(at = which(row.names(plot.data_d)%in%enhancers_anno),side='right', labels = row.names(plot.data_d)[row.names(plot.data_d)%in%enhancers_anno],labels_gp = gpar(fontsize = anno.size)))
    if (i==length(features)){
      hm1 <- Heatmap(as.matrix(plot.data_d),cluster_rows = T,row_title='Enhancers',show_row_dend = F,name = features[i],cluster_columns = F,show_column_names = F,show_row_names = F,col = cols,right_annotation = ra,top_annotation = ha,height=unit(height*0.7, "inch"),use_raster = T,raster_quality = 10,raster_device = 'CairoJPEG',
                     heatmap_legend_param=list(labels = c(0,100),at=c(0,100), direction = "vertical",title='% Max',legend_height = unit(height/3, "inch")),width=unit(6*width/10, "inch"))
    } else {
      hm1 <- Heatmap(as.matrix(plot.data_d),cluster_rows = T,row_title='Enhancers',show_row_dend = F,name = features[i],cluster_columns = F,show_column_names = F,show_row_names = F,col = cols,right_annotation = ra,top_annotation = ha,height=unit(height*0.7, "inch"),use_raster = T,raster_quality = 10,raster_device = 'CairoJPEG',
                     show_heatmap_legend = F)
    }
    hm2 <- Heatmap(as.matrix(plot.data_agg),cluster_rows = F,row_title='',cluster_columns = F,show_column_names = F,show_row_names = F,col = cols,show_heatmap_legend = F,height=unit(height*0.05, "inch"),use_raster = T,raster_quality = 10,raster_device = 'CairoJPEG',right_annotation = rowAnnotation(foo = anno_text('Aggregate', gp = gpar(fontsize = anno.size))))  
    hm3 <- Heatmap(as.matrix(plot.data_e),cluster_rows = F,row_title='',cluster_columns = F,show_column_names = F,show_row_names = F,col = cols,show_heatmap_legend = F,height=unit(height*0.05, "inch"),use_raster = T,raster_quality = 10,raster_device = 'CairoJPEG',right_annotation = rowAnnotation(foo = anno_text('Expression', gp = gpar(fontsize = anno.size))))
    hm4 <- Heatmap(as.matrix(plot.data_m),cluster_rows = F,row_title='',cluster_columns = F,show_column_names = F,show_row_names = F,col = cols,show_heatmap_legend = F,height=unit(height*0.05, "inch"),use_raster = T,raster_quality = 10,raster_device = 'CairoJPEG',right_annotation = rowAnnotation(foo = anno_text('Motif', gp = gpar(fontsize = anno.size))))
    
    ht_list = hm1 %v% hm2%v% hm3%v% hm4
    res_list[[i]] <- ht_list
  }
  grid_layout <- grid.layout(nrow=ifelse(direction=='horizontal',1,length(features)),ncol=ifelse(direction=='horizontal',length(features),1))
  if(direction=='horizontal'){
    width=width*length(features)
  } else {
    height=height*length(features)
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  grid.newpage()
  pushViewport(viewport(layout=grid_layout))
  for (k in seq_along(features)){
    if (direction=='horizontal'){
      pushViewport(viewport(layout.pos.row=1, layout.pos.col=k))
    } else {
      pushViewport(viewport(layout.pos.row=k, layout.pos.col=1))
    }
    draw(res_list[[k]],merge_legends=T,newpage=F,column_title =features[k])
    popViewport()
  }
  dev.off()
}

FigureS4A <- function(rna.object,atac.object, prediction.score, cols, point.size, anno.size, key.size, height, width, plot_filled=F,out_f, theme = NULL,rows=1,stroke=0.2) { 
  atac.object$predicted.id[atac.object$prediction.score.max<0.5] <- NA
  plot.data <- extractFeatures(atac.object,features='predicted.id',cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$feature <- factor(plot.data$feature, levels = c(levels(rna.object),'N/A'))
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols),na.value='grey')  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0,(1.03+0.025*rows)),plot.margin = margin(height/10,0.1,0.1,0.1,unit='inches'), legend.box = "horizontal",legend.spacing.x = unit(0.001, 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = rows)) 
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  if(!is.null(theme)){
    print(p + theme)}
  else {
    print(p)
  }
  dev.off()
}

FigureS4B <- function(atac.object,width,out_f,height) {
  plot.data <- as.data.frame(atac.object$prediction.score.max)
  colnames(plot.data)[1] <- 'Prediction_Score'
  p <- ggplot(plot.data, aes(x=Prediction_Score)) + geom_histogram(binwidth = 0.02,aes(y=2*..density..), colour="black", fill="white")+geom_density(aes(y=2*..density..),alpha=.2, fill="#FF6666") 
  p <- p + xlim(c(0,1)) + xlab('Prediction Score') + ylab('% Total') + scale_y_continuous(limits = c(0,5.2),expand=c(0,0)) + xlim(c(0,1))
  p <- p + geom_vline(xintercept = 0.5,linetype='dashed',color='black',size=1.5)
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
  print(p) 
  dev.off()
}

FigureS4C <- function(atac.object,coembed.object, linkCol,out_f,width=600,height=1200) { 
  celltypes <- as.data.frame(atac.object@active.ident, row.names = rownames(atac.object@active.ident))
  cluster.idents <- data.frame(atac.object@meta.data$seurat_clusters)
  celltype <- data.frame(celltypes, rownames(celltypes), cluster.idents)
  colnames(celltype) <- c("celltype_N", "Barcode", "celltype")
  
  predicted <- as.data.frame(coembed.object$predicted.id)
  predicted <- data.frame(predicted, rownames(predicted))
  colnames(predicted) <- c("Prediction", "Barcode")
  
  order <- rownames(celltypes)
  predicted <- predicted[match(order, predicted$Barcode),]
  
  cluster_predicted <- predicted %>% 
    mutate(prediction = ifelse(Prediction == 'PN1', "0",
                               ifelse(Prediction %in% c("IPC", "IPC_M"), "1",
                                      ifelse(Prediction %in% c("PN2", "PN3"), "3",
                                             ifelse(Prediction %in% c("NSC", "NSC_M"), "2",
                                                    ifelse(Prediction == "IN", "4", 
                                                           ifelse(Prediction =="CR", "5",
                                                                  ifelse(Prediction %in% c("MG", "Mural"), "6", "unassigned"))))))))
  
  
  cluster_predicted$Barcode <- NULL
  celltypes$Barcode <- NULL
  
  DF <- cbind(cluster_predicted, celltype)
  
  DF <- DF %>% 
    mutate(links = ifelse(prediction == celltype, "Correct",
                          ifelse(prediction != celltype, "Wrong", 
                                 ifelse(is.na(prediction), "unassigned", "unassigned"))))
  
  DF$prediction <- as.factor(DF$prediction)
  DF$links <- as.factor(DF$links)
  
  link_color <- data.frame(link = c("Correct", "Wrong", "unassigned"),
                           color = linkCol)
  node_color <- data.frame(node = c('NSC', 'NSC_M', 'IPC', 'IPC_M', 'PN1', 'PN2', 'PN3', 'CR', 'IN', 'MG', "Mural", "'MG+Mural"),
                           color = c(unique(c(RNA_cluster_colors, ATAC_cluster_colors), "#F78F4ADF" )))
  
  p <- plot_sankey(true_labels = DF$Prediction, prediction = DF$celltype_N, links = DF$links,
                   custom_link_color = link_color,  custom_node_color = node_color)
  saveWidget(p, file = paste0('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/plots/figures/',out_f,".html"))
  webshot(paste0('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/plots/figures/',out_f,".html"), file = sprintf("/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/plots/figures/%s.png", out_f), vwidth = width, vheight = height)
}

FigureS4D <- function(object,cols,point.size,alpha=1,anno.size,key.size,out_f,height,width,plot_filled=F,theme,stroke=0.1){
  plot.data <- extractFeatures(object,features='tech',cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data <- plot.data[sample(row.names(plot.data),nrow(plot.data)),]
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size,alpha=alpha) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,alpha=alpha,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0.05,0.95), legend.box = "horizontal",legend.spacing.x = unit(0.1, 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  }
  if(!is.null(theme)){
    p <- p+theme
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS4E <- function(hitobject, cat_names=c("Singleome", "Multiome"), col = c("#440154ff", '#21908dff'), out_f, height, width) {
  
  pdf(paste0(out_f,'.pdf'),height=height,width=width)
  draw.pairwise.venn(area1=as.numeric(queryLength(hitobject)), area2 = as.numeric(subjectLength(hitobject)), cross.area = length(hitobject) , category=cat_names, 
                     lwd = 1, col= col, fill=c(alpha(col[1],0.3), alpha(col[2],0.3)),
                     cex=0.8, fontfamily = 'sans', cat.pos = c(-135,135), cat.dist = c(0.05, 0.05), cat.cex = 0.5)
  dev.off()
  
  
}


FigureS4F <- function(p2g_multi, p2g_single, hitobject, out_f, title, xlab='Multiome Correlation', ylab='Singleome Correlation', height, width) {
  
  hits_m <- as.data.frame(hitobject)
  corr_single <- data.frame(single_corr=p2g_single$posCor$Correlation, m_gene = p2g_single$posCor$gene_name)
  corr_multi <- data.frame(multi_corr=p2g_multi$posCor$Correlation, s_gene=p2g_multi$posCor$gene_name)
  vec_single <- hits_m$queryHits
  vec_single <- hits_m$queryHits
  corr_single <- corr_single[vec_single,]
  vec_multi <- hits_m$subjectHits
  corr_multi <- corr_multi[vec_multi,]
  df <- cbind(corr_multi, corr_single)
  p <- ggplot(df, aes(x=multi_corr, y=single_corr)) + geom_pointdensity(adjust=1,alpha=1,size=1) + scale_color_gradientn(colours = colorpalette("heat",8),name='', breaks=c(1000,15000),labels=c('min','max'))  
  
  p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.34, 0.97),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.4), vjust = 8))
  #p <- p+stat_cor(method = "pearson", label.x = 0.35, label.y = -0.70)  
  my_grob = grid.text(paste0('r=',round(cor.test(df[,1],df[,3])$estimate,3)), x=0.85,  y=0.1, gp=gpar(col="black", fontsize=14, fontface="bold"))
  p1 <- p + annotation_custom(my_grob)
  jpeg(paste0(out_f,'.jpeg'), type = 'cairo', width = width, height=height, pointsize = 1)
  print(p1)
  dev.off()
}



FigureS4G <- function(mat,cols,height,width){
  res <- list()
  for (s in c('posCor','negCor','noCor')){
    df <- as.data.frame(mcols(mat[[s]])) %>% dplyr::count(gene_name)
    df$type <- s
    res[[s]] <- as.data.frame(df)
  }
  df <- Reduce('rbind',res)
  df$type <- gsub('noCor','control',df$type)
  df$type <- factor(df$type,levels=c('posCor','negCor','control'))
  df_med <- df %>% dplyr::group_by(type) %>% summarise(median=median(n))
  p <- ggplot(df,aes(x=type,y=n,fill=type)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
  p <- p + scale_fill_manual(values=cols) + xlab('') + ylab('N distal peaks per gene') + theme(legend.position = "none")
  p <- p + stat_compare_means(comparisons = list(c('posCor','negCor'),c('posCor','control'),c('negCor','control')),label = "p.format",method='wilcox')
  p <- p +  annotate("text", x = c(1.45,2.45,3.45), y = df_med$median, label = df_med$median,size=5) 
  pdf('plots/figures/FigureS3J.pdf',height=height,width=width)
  print(p)
  dev.off()
}



FigureS4H <- function(p2glinks, xlabel, ylabel, posCol, negCol, controlCol,out_f, width, height) {
  df_dist_pos <- data.frame(Distance= p2glinks$posCor$distance, Significance = rep('PosCor', length(p2glinks$posCor)))
  df_dist_neg <- data.frame(Distance = p2glinks$negCor$distance, Significance = rep('NegCor', length(p2glinks$negCor)))
  df_dist_control <- data.frame(Distance = p2glinks$noCor$distance,  Significance = rep('Control', length(p2glinks$noCor))) 
  dist_Df <- rbind(df_dist_pos, df_dist_control, df_dist_neg)
  dist_Df$Significance <- factor(dist_Df$Significance,levels=c('PosCor','NegCor','Control'))
  dist_Df_abs <- dist_Df %>% mutate(., absDistance = abs(Distance)
                                    , kb = ifelse(absDistance >500, round((absDistance/1000), digits = 0), 10))
  dist_Df_abs$cat <- cut(dist_Df_abs$kb, c(0,50,100, 150, 200, 250, 300, 350, 400, 450, 500))
  p <- ggplot(dist_Df_abs, aes(x=cat)) + 
    geom_bar(aes(fill=Significance), position = 'dodge', width = 0.7)+ 
    scale_x_discrete(labels = c('50', '100', '150', '200', '250',
                                '300', '350', '400', '450', '500')) + 
    scale_fill_manual("", values = c("PosCor" = posCol, "NegCor" = negCol, "Control" = controlCol)) +
    xlab(xlabel) +
    ylab(ylabel) + theme(legend.position=c(0.8,0.95)) + scale_y_continuous(expand=c(0,0),n.breaks = 6)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS4I <- function(controlTable,posCorTable,negCorTable,xlabel, ylabel, col1, col2, col3, width, height){
  Df_positive <- posCorTable %>% mutate(Correlation = rep("PosCor", nrow(.)))
  Df_negative <- negCorTable %>% mutate(Correlation = rep("NegCor", nrow(.)))
  Df_control <- controlTable %>% mutate(Correlation = rep("Control", nrow(.)))
  Final <- rbind(Df_positive, Df_negative, Df_control)
  Final$Correlation <- factor(Final$Correlation,levels=c('PosCor','NegCor','Control'))
  
  p <- ggplot(Final, aes(col = Correlation)) +
    stat_smooth(aes(x = position, y = norm.value)) + 
    scale_color_manual('',values = c('PosCor' = col1, 'NegCor' = col2, 'Control' = col3)) + 
    xlab(xlabel) +
    ylab(ylabel) + theme(legend.position = c(0.75,0.9))
  pdf('plots/figures/FigureS3I.pdf', height = height, width = width)
  print(p)
  dev.off()
  
}


###############
#major changes here
###############

#Figure S5A + plotting
p <- VlnPlot(rna.object, features = c('Neurog2'),cols=RNA_cluster_colors[c(1,3,5:7)],idents = c('NSC','IPC','PN1','PN2','PN3'),pt.size = 0,combine = F,same.y.lims = F)
p1 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('')+ ggtitle('Neurog2') #+ ylab('Normalized expression')
p <- VlnPlot(atac.object,assay = 'chromvar', features = c('Neurog2(var.2)'),cols=ATAC_cluster_colors[c(1:4)],idents = c('NSC','IPC','PN1','PN2'),pt.size = 0,combine = F)
p2 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('') + ggtitle('Neurog2 motif') + ylab('Chromvar Deviation Score')
pdf('plots/figures/FigureS5Apdf',width=5,height=4)
print(p1+p2)
dev.off()

#Figure S5B + plotting
p <- VlnPlot(rna.object, features = c('Tead2','Yap1'),cols=RNA_cluster_colors[c(1,3,5:7)],idents = c('NSC','IPC','PN1','PN2','PN3'),pt.size = 0,combine = F)
p1 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('')+ p[[2]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5)) + xlab('')
p <- VlnPlot(atac.object,assay = 'chromvar', features = c('Tead2'),cols=ATAC_cluster_colors[c(1:4)],idents = c('NSC','IPC','PN1','PN2'),pt.size = 0,combine = F)
p2 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('') + ggtitle('Tead2 motif') + ylab('Chromvar Deviation Score')
pdf('plots/figures/FigureS5B.pdf',width=9,height=6)
print(p1+p2)
dev.off()

#Figure S5C + plotting
p <- VlnPlot(rna.object, features = c('Fos','Jun'),cols=RNA_cluster_colors[c(1,3,5:7)],idents = c('NSC','IPC','PN1','PN2','PN3'),pt.size = 0,combine = F,same.y.lims = F)
p1 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('')+ p[[2]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5)) + xlab('')
p <- VlnPlot(atac.object,assay = 'chromvar', features = c('Fosb::jun'),cols=ATAC_cluster_colors[c(1:4)],idents = c('NSC','IPC','PN1','PN2'),pt.size = 0,combine = F)
p2 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('') + ggtitle('Fosb::jun motif') + ylab('Chromvar Deviation Score')
pdf('plots/figures/FigureS5C.pdf',width=9,height=6)
print(p1+p2)
dev.off()

FigureS5D <- function(archr_obj=archr_obj,flank=500,motif=NULL,regions=NULL,normalize_by='Subtract',file_f,which_clusters=levels(atac.object),cols=cluster_colors(length(levels(atac.object))),plot_bias=T,anno.size=12,key.size=4,width,height){
  if (is.character(regions)){
    chip <- read.table(regions)
    colnames(chip)[1:3] <- c('chrom','start','end')
    regions <- makeGRangesFromDataFrame(chip)
  }
  if(!is.null(motif)){
    motifPositions <- getPositions(archr_obj)
    names(motifPositions) <- capitalize(tolower(names(motifPositions)))
    name_arch <- names(motifPositions)[grep(gsub('\\(|\\:|\\)','.',motif),names(motifPositions))]
    positions <- motifPositions[name_arch,]
    if(!is.null(regions)){
      subset_positions <- GRangesList(subsetByOverlaps(positions[[1]],regions))
      names(subset_positions) <- names(positions)[1]
      positions <- subset_positions
    }
  } else {
    positions <- GRangesList(regions)
    name_arch <- NULL
  }
  seFoot <- suppressMessages(getFootprints(
    ArchRProj = archr_obj, flank = flank,
    positions = positions, 
    groupBy = "Clusters",useGroups = which_clusters,threads = 1,verbose = F,logFile = '/dev/null'
  ))
  p <- ggFootprint(
    seFoot = seFoot,
    name = name_arch[1],
    pal = cols,
    smoothWindow = 20,
    flank = flank,
    flankNorm = 50,
    baseSize = 10,
    normMethod = 'Subtract',
    which_clusters=which_clusters
  )
  
  p1 <- p$ggFoot + theme(legend.text = element_text(size=anno.size),legend.justification=c(0,0), legend.position=c(0.89, 0.75), legend.box = "vertical",legend.spacing.x = unit(0.1, 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),ncol = 1)) + ggtitle(ifelse(is.null(motif),'',motif))
  p2 <- p$ggBias + theme(legend.position="none") + ylab('') + scale_y_continuous(n.breaks=3) + xlab('')
  if (plot_bias){
    p <- p1+ p2 + plot_layout(nrow=2,heights=c(4,1))
  } else {
    p <- p1
  }
  pdf(paste0('plots/figures/',file_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}


FigureS5E <- function(mat,mat_expr,mat_prom_acc,anno_mat,which_clusters=NULL,object=coembed.object,features,cols1,cols2,cluster_cols,anno.size,width,height){
  pd <- object$pseudotime
  mat_expr <- mat_expr[row.names(mat_expr)%in%row.names(mat_prom_acc)&row.names(mat_expr)%in%row.names(mat_gb_acc),]
  mat_prom_acc <- mat_prom_acc[row.names(mat_prom_acc)%in%row.names(mat_expr),]
  mat_prom_acc <- mat_prom_acc[match(row.names(mat_expr),row.names(mat_prom_acc)),]
  mat_gb_acc <- mat_gb_acc[row.names(mat_gb_acc)%in%row.names(mat_expr),]
  mat_gb_acc <- mat_gb_acc[match(row.names(mat_expr),row.names(mat_gb_acc)),]
  mat_expr <- mat_expr[,colnames(mat_prom_acc)]
  o <- data.frame(feature=row.names(mat_expr),rna=row.names(mat_prom_acc))
  o$exprvsprom_cor <- rowCorCpp(1:nrow(o),1:nrow(o), as.matrix(mat_expr), as.matrix(mat_prom_acc))
  o$exprvsgb_cor <- rowCorCpp(1:nrow(o),1:nrow(o), as.matrix(mat_expr), as.matrix(mat_gb_acc))
  o$promvsgb_cor <- rowCorCpp(1:nrow(o),1:nrow(o), as.matrix(mat_prom_acc), as.matrix(mat_gb_acc))
  
  o$maxPD_expr <- pd[match(colnames(mat_expr)[max.col(mat_expr)],names(pd))]
  o$maxPD_promAcc <- pd[match(colnames(mat_prom_acc)[max.col(mat_prom_acc)],names(pd))]
  o$maxPD_GeneBodyAcc <- pd[match(colnames(mat_gb_acc)[max.col(mat_gb_acc)],names(pd))]
  
  o$PD_PromExpr <- o$maxPD_promAcc - o$maxPD_expr               #Pseudotime difference Prom accessibility - expression) 
  o$PD_GeneBodyExpr <- o$maxPD_GeneBodyAcc - o$maxPD_expr
  o$PD_PromGeneBody <- o$maxPD_promAcc - o$maxPD_GeneBodyAcc
  o_orig <- o
  col_levels <- levels(rna.object)[1:7]
  col.list <- cluster_cols[1:length(col_levels)]
  names(col.list) <- col_levels
  idents_indx <- factor(object$celltype,levels=col_levels)
  idents_indx <- idents_indx[match(colnames(mat_expr),names(idents_indx))]
  anno_mat$labels <- paste0(anno_mat$gene_name,':',anno_mat$distance)
  #anno_mat <- anno_mat[match(row.names(mat_expr),anno_mat$gene_name),]
  if (is.null(which_clusters)){
    which_clusters <- colnames(anno_mat)[1:7]
  }
  anno_mat <- anno_mat[rowMaxs(as.matrix(anno_mat[,which_clusters]),na.rm=T)>0,]
  mat_expr <- mat_expr[row.names(mat_expr)%in%anno_mat$gene_name,]
  o <- o[o$feature%in%anno_mat$gene_name,]
  colnames(mat_expr) <- as.numeric(idents_indx)
  mat2 <- bin.matrix.rows(mat_expr,bin.size = 25)
  mat2_idents <- levels(idents_indx)[as.numeric(colnames(mat2))]
  o$subgroup <- mat2_idents[max.col(mat2)]
  mat_expr <- mat_expr[!is.na(o$subgroup),]
  mat_expr <- mat_expr[row.names(mat_expr)%in%o$feature,]
  o <- o[o$feature%in%row.names(mat_expr),]
  mat_prom_acc <- mat_prom_acc[row.names(mat_prom_acc)%in%row.names(mat_expr),]
  mat_gb_acc <- mat_gb_acc[row.names(mat_gb_acc)%in%row.names(mat_expr),]
  o <- o[!is.na(o$subgroup),]
  o$subgroup <- gsub('_M','',o$subgroup)
  
  o$subgroup[1:1368] <- 'NSC'
  o$subgroup[1369:1550] <- 'IPC'  
  o$subgroup[1551:1870] <- 'PN1' 
  o$subgroup[1871:nrow(o)] <- 'PN2' 
  # o$subgroup[which(o$subgroup=='IPC')[1]:(which(o$subgroup=='IPC')[1]+250)] <- o$subgroup[1]        ### Adjust the starting row of IPC to account for cells in transition
  # o$subgroup[grep('PN3',o$subgroup)] <- 'PN2'
  o$subgroup <- factor(o$subgroup, levels=c('NSC','IPC','PN1','PN2'))
  panel_fun = function(index, nm) {
    pushViewport(viewport(xscale = c(-70,40), yscale = c(0, 2)))
    grid.rect()
    grid.xaxis(gp = gpar(fontsize = 8))
    grid.boxplot(o$PD_PromExpr[index], pos = 1, direction = "horizontal",outline=F,gp = gpar(fill = "grey"))
    p_diff <- wilcox.test(o$PD_PromExpr[index],mu = 0, alternative = "less")$p.value
    if(p_diff<0.001){p_diff <- scientific(p_diff, digits = 2)}
    if(as.numeric(p_diff)>0.05){p_diff <- ''} else {p_diff <- paste0('; p=',p_diff)}
    grid.text(paste0('M: ',round(median(o$PD_PromExpr[index],na.rm=T),2),p_diff), 0, y = 1.9,
              just = "top", default.units = "native", gp = gpar(fontsize = 10))
    popViewport()
  }
  anno_box = anno_zoom(align_to = o$subgroup, which = "row", panel_fun = panel_fun, 
                       size = unit(height/12, "inch"), gap = unit(0.4, "inch"), width = unit(width/6, "inch"))
  
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat_expr)%in%features),side='left', labels = row.names(mat_expr)[row.names(mat_expr)%in%features],labels_gp = gpar(fontsize = anno.size)))
  pdf('plots/figures/Figure3D_All.pdf',height=height,width=width)
  hm0 <- Heatmap(mat_gb_acc, name = "GeneBody Accessibility",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, left_annotation = ra,col = cols1,top_annotation = ha,raster_quality = 6,
                 heatmap_legend_param=list(labels = c(0,100),at=c(0,1), direction = "vertical",title='% Max',legend_height = unit(height/6, "inch")),width=unit(2*width/10, "inch"))
  hm1 <- Heatmap(mat_prom_acc, name = "Prom Accessibility",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE,col = cols1,top_annotation = ha,raster_quality = 6,
                 heatmap_legend_param=list(labels = c(0,100),at=c(0,1), direction = "vertical",title='% Max',legend_height = unit(height/6, "inch")),width=unit(2*width/10, "inch"))
  hm2 <- Heatmap(mat_expr, name = "Gene expression",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, col = cols1,top_annotation = ha,raster_quality = 6,
                 heatmap_legend_param=list(labels = c(0,100),at=c(0,1), direction = "vertical",title='% Max',legend_height = unit(height/6, "inch")),width=unit(2*width/10, "inch"))
  
  # hm2 <- Heatmap(o$cor, name = "Correlation",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE,col=cols2,
  #                 heatmap_legend_param=list(labels = c(0.35,1),at=c(0.35,1), direction = "vertical",title='Cor',legend_height = unit(height/6, "inch")),width=0.5*width/12)
  hm3 <- Heatmap(o$PD_PromExpr, name = "PD_Prom-Expr",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE,col=cols2,right_annotation = rowAnnotation(foo = anno_box),
                 heatmap_legend_param=list(labels = c(-5,0,5),at=c(-5,0,5), direction = "vertical",title='dPD',legend_height = unit(height/6, "inch")),width=unit(0.35*width/10, "inch"))
  
  draw(hm0+hm1+hm2+hm3, merge_legends=T, ht_gap = unit(0.2, "inch")) 
  dev.off()
  return(o_orig)
}


FigureS5F <- function(mat_f,mat_rna,anno_mat,object=coembed.object,subset_mat=T,features,cols1,cols2,cols3,cluster_cols,anno.size,width,height){
  mat_rna <- mat_rna[,colnames(mat_f)]
  mat_f$labels <- sapply(strsplit(as.character(row.names(mat_f)), "\\(|\\.|\\:|\\)"), function(x) x[[1]])
  pd <- object$pseudotime
  mat_rna <- mat_rna[match(mat_f$labels,row.names(mat_rna)),]
  o <- data.frame(feature=row.names(mat_f),rna=row.names(mat_rna),row=1:nrow(mat_f),cor=NA,maxPD_f=NA,maxPD_rna=NA)
  indx <- o$row[!is.na(o$rna)]
  o$cor[indx] <- rowCorCpp(indx,indx,as.matrix(mat_f[,-ncol(mat_f)]),as.matrix(mat_rna))
  o$maxPD_f[indx] <- pd[match(colnames(mat_f[,-ncol(mat_f)])[max.col(mat_f[,-ncol(mat_f)])],names(pd))][indx]
  o$maxPD_rna[indx] <- pd[match(colnames(mat_rna)[max.col(mat_rna)],names(pd))][indx]
  o$PDdiff <- o$maxPD_f-o$maxPD_rna
  col_levels <- levels(object)[1:7]
  col.list <- cluster_cols[1:length(col_levels)]
  names(col.list) <- col_levels
  idents_indx <- factor(object$celltype,levels=col_levels)
  idents_indx <- idents_indx[match(colnames(mat_f[,-ncol(mat_f)]),names(idents_indx))]
  anno_mat <- anno_mat[match(mat_f$labels,anno_mat$gene_name),]
  subgroup <- factor(colnames(anno_mat)[max.col(anno_mat[,1:4])],levels=colnames(anno_mat)[1:4])
  o$subgroup <- subgroup
  if (subset_mat){
    mat_f <- mat_f[complete.cases(o),]
    o <- o[complete.cases(o),]
  }
  ha = HeatmapAnnotation(cluster = idents_indx ,col = list(cluster=col.list),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat_f)%in%features),side='left', labels = row.names(mat_f)[row.names(mat_f)%in%features],labels_gp = gpar(fontsize = anno.size)))
  pdf('plots/figures/Figure3E.pdf',height=height,width=width)
  hm1 <- Heatmap(as.matrix(mat_f[,-ncol(mat_f)]), name = "Pseudotime",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 6, left_annotation = ra,col = cols1,top_annotation = ha,
                 heatmap_legend_param=list(labels = c(0,100),at=c(0,1), direction = "vertical",title='% Max',legend_height = unit(height/6, "inch")),width=unit(6*width/10, "inch"))
  hm2 <- Heatmap(o$cor, name = "Correlation",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE,col=cols2, raster_quality = 10,
                 heatmap_legend_param=list(labels = c(-1,0,1),at=c(-1,0,1), direction = "vertical",title='Cor',legend_height = unit(height/6, "inch")),width=0.6*width/10)
  hm3 <- Heatmap(o$PDdiff, name = "PDdiff",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10,col=cols3,
                 heatmap_legend_param=list(labels = c(-5,0,5),at=c(-5,0,5), direction = "vertical",title='dPD',legend_height = unit(height/6, "inch")),width=unit(0.6*width/10, "inch"))
  
  draw(hm1+hm2+hm3, merge_legends=T, ht_gap = unit(0.2, "inch")) 
  dev.off()
  return(o)
}

FigureS5H <- function(object,out_f,mat,enhancers,aggregate_enh=T,features,motifs,rescale=F,cols,key.size,height,width,cluster_cols=cluster_colors,anno.size=12,direction='vertical'){
  mat$labels <- paste0(mat$gene_name,':',mat$distance)
  plot_list <- list()
  for (i in seq_along(features)){
    DefaultAssay(object) <- 'MACS2peaks'
    if (aggregate_enh){
      distal <- as.character(mat[mat$gene_name==enhancers[i],'IDs'])
      plot.data_d <- extractFeatures(object,features=c(distal,'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
      plot.data_d$feature <- rowSums(plot.data_d[,3:(2+length(distal))],na.rm=T)
      plot.data_d <- plot.data_d[,c('UMAP_1','UMAP_2','feature','pseudotime','labels','reps')]
    } else {
      distal <- as.character(mat[mat$labels==enhancers[i],'IDs'])
      plot.data_d <- extractFeatures(object,features=c(distal,'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    }
    DefaultAssay(object) <- 'RNA'
    plot.data_e <- extractFeatures(object,features=c(features[i],'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    
    DefaultAssay(object) <- 'chromvar'
    plot.data_m <- extractFeatures(object,features=c(motifs[i],'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    plot.data_d$labels <- 'Enhancer Accessibility'
    plot.data_e$labels <- 'Gene Expression'
    plot.data_m$labels <- 'Motif Accessibility'
    plot.data_d <- plot.data_d[!is.na(plot.data_d$pseudotime),]
    plot.data_d$rescaled <- scales::rescale(predict(mgcv::gam(formula = feature ~ s(pseudotime, bs = "cs"),method = "REML", data=plot.data_d)),to=c(0,1))
    plot.data_e <- plot.data_e[!is.na(plot.data_e$pseudotime),]
    plot.data_e$rescaled <- scales::rescale(predict(mgcv::gam(formula = feature ~ s(pseudotime, bs = "cs"),method = "REML", data=plot.data_e)),to=c(0,1))
    plot.data_m <- plot.data_m[!is.na(plot.data_m$pseudotime),]
    plot.data_m$rescaled <- scales::rescale(predict(mgcv::gam(formula = feature ~ s(pseudotime, bs = "cs"),method = "REML", data=plot.data_m)),to=c(0,1))
    plot.data <- rbind(plot.data_d,plot.data_e,plot.data_m)
    if(rescale){
      plot.data$feature <- plot.data$rescaled*100
      ylab_f <- '% Max'
    } else {
      ylab_f <- ''
    }
    p <- ggplot(plot.data, aes(x=pseudotime,y=feature,color=labels)) + ylab("") + xlab("Pseudotime") + scale_colour_manual(name='',values=as.character(cols)) + geom_smooth(aes(x=pseudotime,y=feature,color=labels),inherit.aes = F)         #+ geom_point(plot.data[sample(row.names(plot.data[plot.data$feature!=0,]),1000),],mapping=aes(x=pseudotime,y=feature,color=labels),size = 0.5,alpha=0.5,inherit.aes = T)
    p <- p + geom_vline(mapping=aes(xintercept=pseudotime,color=labels),data=plot.data[plot.data$rescaled==1,],linetype = "longdash",show.legend=F,alpha=0.5,lwd=1)
    p <- p + theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal") + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) + ggtitle(features[i])
    p <- p + ylab(ylab_f) + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),plot.margin = margin(0,0,0.5,0))
    plot_list[[features[i]]] <- p + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(limits = c(-5,105),expand = c(0,0)) + geom_segment(aes(x=0,y=-5,xend=max(pseudotime),yend=-5),col='black',arrow=arrow(length=unit(0.2,"cm"),type = 'closed'))
  }
  plot.data <- FetchData(object,vars=c('celltype','pseudotime'))
  colnames(plot.data)[1] <- 'labels'
  plot.data <- plot.data[complete.cases(plot.data),]
  plot.data <- plot.data[!plot.data$labels%in%c('MG','Mural','IN',"MG+Mural","CR"),]
  plot.data$labels <- factor(plot.data$labels,levels=levels(object))
  plot.data$labels <- droplevels(plot.data$labels)
  p1 <- ggplot(plot.data) + geom_rect(mapping=aes(xmin=pseudotime-0.1,xmax=pseudotime+0.1,ymin=0,ymax=1,fill=labels)) + xlab('Pseudotime')
  p1 <- p1 + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.line.y=element_blank(),plot.margin = unit(c(0,0,0,0), "cm")) + scale_fill_manual(name='',values=as.character(cluster_cols[1:length(levels(plot.data$labels))])) 
  p1 <- p1 + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(plot.data$labels))*10), 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  if(direction=='vertical'){
    p <- Reduce("/", plot_list) 
    p <- guide_area()/p + plot_layout(guides="collect",heights = c(1,3*length(features)))
    p <- guide_area()/(p)/p1 + plot_layout(guides = 'collect',ncol=1,heights=c(0.2,length(features)*height,0.2))
    pdf(paste0('plots/figures/',out_f,'.pdf'),height=length(features)*height,width=width)
  } else {
    p <- Reduce("|", plot_list) 
    p <- guide_area()/p + plot_layout(guides="collect",heights = c(1, 9))
    #(p)/guide_area() + plot_layout(guides = 'collect',nrow=2,heights=c(8,0.5),widths=c(8))
    pdf(paste0('plots/figures/',out_f,'.pdf'),width=length(features)*width,height=height)
  }
  print(p)
  dev.off()
}



### Plot Figures #######

mat <- Figure3A(selected_peaks=p2glinks$posCor,scaleMax=F,seRNA=seRNA_all,sePB=sePB_all,out_f='Figure3A',
                which_clusters_rna = c('NSC','IPC','PN1','PN2','PN3'),which_clusters = c('NSC','IPC','PN1','PN2'),
                features=c("Sox2:434656","Emx2:-36039","Tfap2c:177606","Dmrta2:231841","Eomes:-159192","Eomes:-249213","Neurog2:13167","Rnd2:148829",'Rnd2:9546',"Gli3-241669","Gli3:81089",'Pax6:405689',"Neurod6:34596","Sox5:297066","Sox5:-98202"),
                cols1=colorRamp2(seq(-2,2,length.out=length(heatmap_colors)), heatmap_colors),cluster_cols1=RNA_cluster_colors[-c(2,4)],cluster_cols2=RNA_cluster_colors[-c(2,4)],
                cols2=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))),anno.size=12,height=12,width=10)
#write.table(mat,file = paste0('results/P2G_binaryMat_posCor.tsv'),quote = F,col.names=T,row.names=F,sep='\t')
plotMisha(object=atac.object,targetGene='Rnd2',out_f='Figure3B',upstream=2e4,downstream=17e4,
          chipNames=c('NSC','IPC','PN1','PN2','CR','IN','MG+Mural'),
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN1','scATAC.E14_PN2','scATAC.E14_CR','scATAC.E14_IN','scATAC.E14_MG_Mural'),
          chipColors=ATAC_cluster_colors,arcIntervals=p2glinks$posCor,arcColors=colorRampPalette(archr_colors[[29]]),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.6, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

motifs <- Figure3C_1(mat_f=read.table('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/P2G_binaryMat_posCor.tsv',header=T),expr=expr[rowMaxs(expr)>=6,],
                     out_f='Figure3C_1',which_clusters=c('NSC','IPC','PN1','PN2'),background_others = F,object=atac.object,
                     features=c('Neurog2(var.2)','Neurod2','Eomes','Tfap2c(var.2)','Tead2','Sox2','Lhx2','Fosb::jun','Insm1','Mef2c','Pou3f2','Ctcf'),
                     cols=c("blue",'grey80','red'),plot_motifs=F,n_cols =2,
                     logFC=0.25,logP=2,height=4,width=5,point.size=3,anno.size=5)

Figure3C_2(object=atac.object,out_f='Figure3C_2',features=c('Neurog2(var.2)','Neurod2','Eomes','Tfap2c(var.2)','Tead2','Sox2','Lhx2','Fosb::jun','Insm1','Mef2c','Ctcf'),
           n_cols=5,height=2.2,width=10)

o <- Figure3D(mat=p2glinks$posCor,anno_mat=read.table('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/P2G_binaryMat_posCor.tsv',header=T),which_clusters=c("NSC","IPC","PN1","PN2"),
              mat_f=as.data.frame(readRDS('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/Monocle3_int/pd_P2G_GAMmat.RDS')),
              mat_rna=as.data.frame(readRDS('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/Monocle3_int/pd_genes_GAMmat.RDS')),
              object=coembed.object,
              features=c("Sox2:434656","Emx2:-36039","Tfap2c:177606","Dmrta2:231841","Eomes:-159192","Eomes:-249213","Neurog2:13167","Rnd2:148829",'Rnd2:9546',"Gli3-241669","Gli3:81089","Pax6:50535",'Pax6:405689',"Neurod6:34596","Sox5:297066","Sox5:-98202"),
              cols1=colorRamp2(seq(0,1,length.out=length(heatmap_colors)), heatmap_colors),
              cols2=colorRamp2(seq(-5,5,length.out=11), rev(colorpalette('rdbu',11))),
              cluster_cols=RNA_cluster_colors,
              anno.size=10,width=12,height=10)
#write.table(o,file = paste0('results/posCor_P2G_PD.tsv'),quote = F,col.names=T,row.names=T,sep='\t')

Figure3E(PWM='/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/scATAC/combined_pwm.RDS', 
         scRNA_expression=expr,
         posCOR_pairs='/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/P2G_binaryMat_posCor.tsv',
         pseudotime='/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/PD_Enh_RNA.tsv',
         FDR=1,out_f='Figure3E',
         expression_cutoff=3,height=12,width=8)

Figure3F(out_f='Figure3F',PWM='/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/scATAC/combined_pwm.RDS',expression_cutoff=3,
         pseudotime='/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/PD_Enh_RNA.tsv',
         motif_name="Neurog2(var.2)",chip_f="/home/hpc/bonev/data/hic/data/peaks/final/NPC_Neurog2.bed",
         posCOR_pairs='/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/P2G_binaryMat_posCor.tsv',FDR_cutoff=0.2, width=6,height=8)

Figure3G(object=coembed.object,mat=read.table('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/results/P2G_binaryMat_posCor.tsv',header=T),
         cluster_cols=RNA_cluster_colors,out_f='Figure3G',
         enhancers=c('Sox2:434656','Sox2:-206704','Insm1:88500','Insm1:-37679','Eomes:-159192','Eomes:-249213','Eomes:-150098','Id4:-312757','Id4:-219895','Neurod2:8330','Neurod2:41915'),
         features=c('Sox2','Id4','Eomes','Insm1','Neurod2'),direction='horizontal',
         motifs=c('Sox2','Id4','Eomes','Insm1','Neurod2'),
         cols=colorRamp2(seq(0,100,length.out=length(heatmap_colors)), heatmap_colors),height=7,width=6,anno.size=12)



## Supplementary Figure 4

FigureS4A(rna.object=rna.object,atac.object=atac.object, prediction.score = 0.5, cols=RNA_cluster_colors, point.size=2,anno.size=14,out_f='FigureS4A',key.size=4,height=6,width=5.5,rows=2,plot_filled=T,theme = theme_border,stroke=0.2)
FigureS4B(atac.object=atac.object,out_f='FigureS4B',width=4, height=4)
FigureS4C(atac.object,coembed.object,linkCol=rep('grey',3),out_f='FigureS4C',width=500,height=800,zoom=2)
FigureS4D(object=coembed.object,rep_colors,point.size=1.5,alpha=0.8,anno.size=14,key.size=6,out_f='FigureS4D',height=6,width=6,plot_filled=T,stroke=0.1,theme=theme_border)
FigureS4E(hitobject = hits, out_f='FigureS4E',height=3, width=3.45 )
FigureS4F(p2g_multi = outMatch_m, p2g_single = outMatch_s, hitobject = hits, out_f='FigureS4F', height=520, width=480 )
FigureS4G(mat=p2glinks,cols=c('darkred','darkblue','grey'),height=6,width=3.5)
FigureS4H(p2glinks, xlabel = "Distance (KB)", ylabel = "N connections", posCol = "darkred", negCol = "darkblue", controlCol = "grey",out_f='FigureS4G', width = 6, height=6)
#FigureS4H(p2glinks, ylabels = "% of all", col1 = "darkred", col2 = "darkblue", xlabels = "Correlation Type", width = 6, height = 6)
FigureS4I(controlTable=readRDS('/home/faye/E14_Analysis/Figures/controlTable.RDS'),
           posCorTable=readRDS('/home/faye/E14_Analysis/Figures/posCorTable.RDS'),
           negCorTable=readRDS('/home/faye/E14_Analysis/Figures/negCorTable.RDS'),
           xlabel="Genomic Distance (KB)", ylabel='Mean Accessibility score', col1 = 'darkred', col2 = 'darkblue', col3 = 'grey', width = 5, height = 5.5)
#Figure S4J
mat <- Figure3A(selected_peaks=p2glinks$noCor,scaleMax=F,seRNA=seRNA_all,sePB=sePB_all,out_f='FigureS4J',
                 which_clusters_rna = c('NSC','IPC','PN1','PN2','PN3'),which_clusters = c('NSC','IPC','PN1','PN2'),
                 features=NULL,
                 cols1=colorRamp2(seq(-2,2,length.out=length(heatmap_colors)), heatmap_colors),cluster_cols1=RNA_cluster_colors[-c(2,4)],cluster_cols2=RNA_cluster_colors[-c(2,4)],
                 cols2=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))),anno.size=12,height=12,width=10)
 write.table(mat,file = paste0('results/P2G_binaryMat_noCor.tsv'),quote = F,col.names=T,row.names=F,sep='\t')
 #Figure S4K
 mat <- Figure3A(selected_peaks=p2glinks$negCor,scaleMax=F,seRNA=seRNA_all,sePB=sePB_all,out_f='FigureS4K',
                 which_clusters_rna = c('NSC','IPC','PN1','PN2','PN3'),which_clusters = c('NSC','IPC','PN1','PN2'),
                 features=NULL,
                 cols1=colorRamp2(seq(-2,2,length.out=length(heatmap_colors)), heatmap_colors),cluster_cols1=RNA_cluster_colors[-c(2,4)],cluster_cols2=RNA_cluster_colors[-c(2,4)],
                 cols2=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))),anno.size=12,height=12,width=10)
 write.table(mat,file = paste0('results/P2G_binaryMat_negCor.tsv'),quote = F,col.names=T,row.names=F,sep='\t')


#plotMisha(object=atac.object,targetGene='Neurod2',out_f='Figure3H',upstream=5.5e4,downstream=5.5e4,
#          chipNames=c('NSC','IPC','PN1','PN2','Neurog2'),
#          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN1','scATAC.E14_PN2','chipseq_RPM.NPC_Neurog2'),
#         chipColors=c(ATAC_cluster_colors[1:4],'black'),arcIntervals=mat[mat$gene_name=='Neurod2'&(mat$distance==41915|mat$distance==8330)],arcColors=colorRampPalette(archr_colors[[29]]),
#         plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
#          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.6, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

#Figure S5D
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))
require(ArchR)
archr_obj <- loadArchRProject("/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/archr/seurat_atac3/",showLogo = F)

FigureS5D(archr_obj=archr_obj,flank=500,motif='Neurog2(var.2)',regions=NULL,
                   normalize_by='Subtract',file_f='FigureS3L',cols=ATAC_cluster_colors[1:4],
                   which_clusters=c('NSC','IPC','PN1','PN2'),plot_bias=T,anno.size=12,key.size=4,width=6,height=6)

FigureS5D(archr_obj=archr_obj,flank=500,motif='Tead2',regions=NULL,
                   normalize_by='Subtract',file_f='FigureS3L_2',cols=ATAC_cluster_colors[1:4],
                   which_clusters=c('NSC','IPC','PN1','PN2'),plot_bias=T,anno.size=12,key.size=4,width=6,height=6)

FigureS5D(archr_obj=archr_obj,flank=500,motif='Fos..jun.var.2_327',regions=NULL,
                   normalize_by='Subtract',file_f='FigureS3L_3',cols=ATAC_cluster_colors[1:4],
                   which_clusters=c('NSC','IPC','PN1','PN2'),plot_bias=T,anno.size=12,key.size=4,width=6,height=6)

FigureS5D(archr_obj=archr_obj,flank=500,motif='Ctcf',regions=NULL,
                   normalize_by='Subtract',file_f='FigureS3L_4',cols=ATAC_cluster_colors[1:4],
                   which_clusters=c('NSC','IPC','PN1','PN2'),plot_bias=T,anno.size=12,key.size=4,width=6,height=6)


o <- FigureS5E(mat=p2glinks$posCor,anno_mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T),which_clusters=c("NSC","IPC","PN1","PN2"),
               mat_expr=as.data.frame(readRDS('results/Monocle3_int/pd_genes_GAMmat.RDS')),
               mat_prom_acc=as.data.frame(readRDS('results/Monocle3_int/pd_PromAccess_GAMmat.RDS')),
               mat_gb_acc=as.data.frame(readRDS('results/Monocle3_int/pd_GeneBodyAccess_GAMmat.RDS')),
               object=coembed.object,
               features=c("Sox2","Emx2","Tfap2c","Dmrta2","Eomes","Neurog2","Rnd2","Gli3","Pax6","Neurod6","Sox5"),
               cols1=colorRamp2(seq(0,1,length.out=length(heatmap_colors)), heatmap_colors),
               cols2=colorRamp2(seq(-5,5,length.out=11), rev(colorpalette('rdbu',11))),
               cluster_cols=RNA_cluster_colors,
               anno.size=10,width=12,height=10)
write.table(o,file = paste0('results/promAcc_expr_PD.tsv'),quote = F,col.names=T,row.names=T,sep='\t')

o_m <- FigureS5F(anno_mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T),
                 mat_f=as.data.frame(readRDS('results/Monocle3_int/pd_motifs_mat.RDS')),
                 mat_rna=as.matrix(readRDS('results/Monocle3_int/pd_genes_GAMmat.RDS')),
                 object=coembed.object,subset_mat=T,
                 features=c('Pax6','Lhx2','Fosb::jun','Hes1','Hes5','Zeb1','Id4','Sox2','Neurog2(var.2)','Neurod2','Dmrta2','Tead2','Eomes','Tfap2c(var.2)','Xbp1','Insm1','Mef2c'),
                 cols1=colorRamp2(seq(0,1,length.out=length(heatmap_colors)), heatmap_colors),
                 cols2=colorRamp2(seq(-1,1,length.out=11),rev(colorpalette('puor',11))),
                 cols3=colorRamp2(seq(-5,5,length.out=11), rev(colorpalette('rdbu',11))),
                 cluster_cols=RNA_cluster_colors,
                 anno.size=10,width=10,height=5)
plotMisha(object=atac.object,targetGene='Neurod2',out_f='FigureS5G',upstream=5.5e4,downstream=5.5e4,
          chipNames=c('NSC','IPC','PN1','PN2','Neurog2'),
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN1','scATAC.E14_PN2','chipseq_RPM.NPC_Neurog2'),
          chipColors=c(ATAC_cluster_colors[1:4],'black'),arcIntervals=mat[mat$gene_name=='Neurod2'&(mat$distance==41915|mat$distance==8330)],arcColors=colorRampPalette(archr_colors[[29]]),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.6, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

FigureS5H(object=coembed.object,mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T),cluster_cols=RNA_cluster_colors,
          enhancers=c('Sox2','Id4','Eomes','Insm1','Neurod2'),aggregate_enh=T,
          features=c('Sox2','Id4','Eomes','Insm1','Neurod2'),direction='horizontal',
          motifs=c('Sox2','Id4','Eomes','Insm1','Neurod2'),out_f='FigureS5H',
          rescale=T,cols=c('darkblue','darkgreen','darkred'),key.size=2,height=4,width=6,anno.size=12)


