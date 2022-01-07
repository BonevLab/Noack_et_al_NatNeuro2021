library(plyr)
library(dplyr)
library(stringr)
library(MPRAnalyze)
library(ggplot2)
library(ggpubr)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(Hmisc)
library(Matrix)
library(matrixStats)
require(LSD)
library(ggrepel)
library(ggpointdensity)
library(grid)
library(Rcpp)
library(ComplexHeatmap)
library(monaLisa)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SummarizedExperiment)

theme_set(theme_cowplot())

source('scripts/figures/config.R')
source('scripts/aux_functions.R')
source('scripts/figures/plot_functions.R')

source('scripts/hic/config.R')
source('scripts/hic/scripts/main_functions.R')
source('scripts/hic/scripts/aux_functions.R')
source('scripts/hic/scripts/plot_functions.R')

featureToGR <- function(feature,pattern="_"){
  featureSplit <- stringr::str_split(paste0(feature), pattern =pattern , n = 3, simplify = TRUE)
  gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
  return(gr)
}

getColsByBin <- function(b) {
  res <- c(rev(colorpalette('brbg',nlevels(b)-1)),'grey')
  names(res) <- levels(b)
  return(res)
}

#Temporarily overwrite monaLisa functions to get the colors right for Figure4G
plotMotifHeatmaps <- function(x,
                              which.plots = c("negLog10P", "pearsonResid", "negLog10Padj", "log2enr"),
                              width = 4,
                              col.enr = c("#053061","#2166AC","#4393C3","#92C5DE",
                                          "#D1E5F0","#F7F7F7","#FDDBC7","#F4A582",
                                          "#D6604D","#B2182B","#67001F"),
                              col.sig = c("#F0F0F0","#D9D9D9","#BDBDBD","#969696",
                                          "#737373","#525252","#252525","#000000"),
                              col.gc = c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B",
                                         "#74C476","#41AB5D","#238B45","#006D2C",
                                         "#00441B"),
                              maxEnr = NULL,
                              maxSig = NULL,
                              highlight = NULL,
                              cluster = FALSE,
                              show_dendrogram = FALSE,
                              show_motif_GC = FALSE,
                              show_seqlogo = FALSE,
                              width.seqlogo = 1.5,
                              use_raster = FALSE,
                              na_col = "white", 
                              ...) {
  
  b <- metadata(x)$bins
  bincols <- getColsByBin(b)
  if (identical(cluster, TRUE)) {
    clAssayName <- "pearsonResid"
    clAssay <- assay(x, clAssayName)
    allNA <- rowSums(is.na(clAssay)) == ncol(clAssay)
    if (any(allNA)) {
      warning("removing motifs without finite values in '",
              clAssayName, "': ",
              paste(rownames(clAssay)[allNA], collapse = ", "))
      x <- x[!allNA, ]
      clAssay <- clAssay[!allNA, ]
    }
    clres <- hclust(dist(clAssay))
  } else if (identical(cluster, FALSE)) {
    clres <- FALSE
  } else if (is(cluster, "hclust")) {
    clres <- cluster
  } else {
    stop("'cluster' must be either TRUE, FALSE or an hclust-object.")
  }
  hmBin <- HeatmapAnnotation(df = data.frame(bin = colnames(x)), name = "bin",
                             col = list(bin = bincols),
                             show_annotation_name = FALSE,
                             which = "column", width = unit(width,"inch"),
                             annotation_height = unit(width / 16, "inch"),
                             show_legend = FALSE)
  tmp <- matrix(if (!is.null(highlight)) as.character(highlight) else rep(NA, nrow(x)),
                ncol = 1, dimnames = list(unname(rowData(x)$motif.name), NULL))
  hmSeqlogo <- NULL
  if (show_seqlogo) {
    pfms <- rowData(x)$motif.pfm
    maxwidth <- max(sapply(TFBSTools::Matrix(pfms), ncol))
    grobL <- lapply(pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")
    hmSeqlogo <- HeatmapAnnotation(
      logo = anno_seqlogo(grobL = grobL, which = "row",
                          space = unit(0.5, "mm"),
                          width = unit(width.seqlogo, "inch")),
      show_legend = FALSE, show_annotation_name = FALSE, which = "row")
  }
  hmMotifs <- Heatmap(matrix = tmp, name = "names",
                      width = unit(if (!is.null(highlight)) .2 else 0, "inch"),
                      na_col = NA, col = c("TRUE" = "green3", "FALSE" = "white"),
                      cluster_rows = clres, show_row_dend = show_dendrogram,
                      cluster_columns = FALSE, show_row_names = TRUE,
                      row_names_side = "left", show_column_names = FALSE,
                      show_heatmap_legend = FALSE, left_annotation = hmSeqlogo)
  
  assayNameMap1 <- c(negLog10P = "P value",
                     negLog10Padj = "adj. P value",
                     pearsonResid = "Pearson residual",
                     log2enr = "log2 enrichment")
  assayNameMap2 <- c(negLog10P = "P value (-log10)",
                     negLog10Padj = "adj. P value (-log10)",
                     pearsonResid = "Pearson residual (o-e)/sqrt(e)",
                     log2enr = "enrichment (log2)")
  L <- list(labels = hmMotifs)
  if (show_motif_GC) {
    tmp <- as.matrix(rowData(x)[, "motif.percentGC", drop = FALSE])
    hmPercentGC <- Heatmap(matrix = tmp, name = "Percent G+C",
                           width = unit(0.2, "inch"), na_col = NA,
                           col = colorRamp2(breaks = c(0, seq(20, 80, length.out = 254), 100),
                                            colors = colorRampPalette(col.gc)(256)),
                           cluster_rows = FALSE, cluster_columns = FALSE,
                           show_row_names = FALSE, show_column_names = FALSE,
                           show_heatmap_legend = TRUE,
                           heatmap_legend_param = list(color_bar = "continuous"),
                           use_raster = use_raster)
    L <- c(L, list("percentGC" = hmPercentGC))
  }
  ret <- c(L, lapply(which.plots, function(w) {
    dat <- assay(x, w)
    if ((w == "pearsonResid") | (w == "log2enr")) {
      rng <- c(-1, 1) * if (is.null(maxEnr)) quantile(abs(dat), .995, na.rm = TRUE) else maxEnr
      cols <- col.enr
    } else {
      rng <- c(0, if (is.null(maxSig)) quantile(dat, .995, na.rm = TRUE) else maxSig)
      cols <- col.sig
    }
    Heatmap(matrix = dat,
            name = assayNameMap1[w],
            width = unit(width,"inch"),
            column_title = assayNameMap2[w],
            col = colorRamp2(breaks = seq(rng[1], rng[2], length.out = 256),
                             colors = colorRampPalette(cols)(256)),
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            ##column_names_side = "bottom", column_names_max_height = unit(1.5,"inch"),
            top_annotation = hmBin, show_heatmap_legend = TRUE,
            heatmap_legend_param = list(color_bar = "continuous"),
            use_raster = use_raster,
            na_col = na_col, 
            ...)
  }))
  names(ret)[seq(length(ret) - length(which.plots) + 1L, length(ret))] <- which.plots
  show(Reduce(ComplexHeatmap::add_heatmap, ret))
  invisible(ret)
}


getColsByBin <- function(b) {
  res <- c(rev(colorpalette('brbg',nlevels(b)-1)),'grey')
  names(res) <- levels(b)
  return(res)
}



# Import data
p2glinks <- readRDS('results/P2G-Links.RDS')
mpra_annot <- vroom::vroom('results/MPRA_anno.tsv')
res_bulk <- vroom::vroom('results/MPRA18k_Bulk_res.tsv')
res <- vroom::vroom('results/immunoMPRA_res.tsv')
res_bulk$enh_type[grep('scrContr',res_bulk$enh_type)] <-'scrContr'
#####

Figure4B <- function(df,CRE_groups,CRE_groupNames,cols,ylab_n){
  df <- df[,c('enh_type','statistic')]
  df <- df[df$enh_type%in%CRE_groups,]
  df$enh_type <- factor(df$enh_type,levels=CRE_groups)
  levels(df$enh_type) <- CRE_groupNames
  colnames(df) <- c('CRE','MPRA')
  p <- ggplot(df,aes(x=CRE,y=MPRA,fill=CRE)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
  p <- p + scale_fill_manual(values=cols) + xlab('') + ylab(ylab_n) + theme(legend.position = "none")
  p <- p + coord_cartesian(ylim=c(0,1.5)) + stat_compare_means(comparisons = list(c('scrambled','noCor'),c('scrambled','posCor'),c('noCor','posCor')),label = "p.format",method='wilcox',label.y = c(1.1,1.3,1.5),tip.length = c(0))
  pdf(out_f,height=height,width = width)
  print(p)
  dev.off()
}

Figure4C <- function(df,anno,CRE_groups,CRE_groupNames,cols,ylab_n,max_sig=2,features=NULL,anno.size=2,point.size=2,out_f,width,height){
  df <- df[,c('enh_type','statistic','MPRA_name','pval.mad')]
  df <- df[df$enh_type%in%CRE_groups,]
  df$enh_type <- factor(df$enh_type,levels=CRE_groups)
  levels(df$enh_type) <- CRE_groupNames
  df$labels <- anno$labels[match(df$MPRA_name,anno$MPRA_name)]         ### IS not 100% accurate because of double enh usage
  colnames(df) <- c('CRE','MPRA',"name",'P.Value')
  df <- df[!is.na(df$MPRA),]
  df <- df[order(df$MPRA,decreasing=T),]
  df$order <- seq(1:nrow(df))
  df$order <- factor(df$order,levels=as.numeric(df$order))
  df$sig <- -log10(df$P.Value)
  df$sig[is.infinite(df$sig)] <- max(df$sig[!is.infinite(df$sig)],na.rm=T)
  df$sig[df$sig>max_sig] <- max_sig
  p <- ggplot(df,aes(x=order,y=MPRA,color=sig)) + geom_point(size = point.size) + ylab(ylab_n) + xlab('')
  p <- p + scale_color_gradientn(colours = cols,name='',breaks=c(min(df$sig,na.rm=T),max(df$sig,na.rm=T)),labels=c(round(min(df$sig,na.rm=T)),round(max(df$sig,na.rm=T)))) + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
  p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.95))
  p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1)) 
  if(!is.null(features)){
    sub_anno <- anno[anno$labels%in%features&anno$type=='WTposCor',]
    df$labels <- sub_anno$labels[match(df$name,sub_anno$MPRA_name)]
    p <- p + geom_text_repel(
      data = df[df$labels%in%features,], size = anno.size,seed = 42,
      box.padding =0.8, min.segment.length = 0,max.iter = 10000,inherit.aes = FALSE,
      aes(x=order,y=MPRA,color=NULL,label=labels))
  }                                                                                                                                                                          
  p <- rasterize(p, layers='Point', dpi=300)
  pdf(out_f,height=height,width = width)
  print(p)
  dev.off()
}

Figure4E <- function(df,anno,clusters,CRE_groups,cols,which_anno_cluster,out_f,width,height){
  df <- df[df$enh_type%in%CRE_groups,]
  df$cluster <- anno[match(df$MPRA_name,anno$MPRA_name),which_anno_cluster]
  colnames(df) <- c('CRE','MPRA_name','NSC','IPC','PN','cluster')
  df <- df[df$cluster%in%clusters|df$CRE=='scrContr'|df$CRE=='WTnoCor',]
  df$cluster[df$CRE=='scrContr'] <- 'scrambled'
  df$cluster[df$CRE=='WTnoCor'] <- 'noCor'
  df$cluster <- factor(df$cluster,levels=c('scrambled','noCor',clusters))
  simple_df <- df[,-c(1:2)]
  mat <- ddply(simple_df,.(cluster),function(x){
    return(colMedians(as.matrix(x[,1:3])))
  })
  colnames(mat) <- c('cluster','NSC','IPC','PN')
  row.names(mat) <- mat$cluster
  mat <- round(as.matrix(mat[,-1]),2)
  hm <- Heatmap(mat,name = "MPRA signal",cluster_rows = F,cluster_columns = F,show_row_names = T,border = 'black',column_names_rot = 0,column_names_centered = T,col = cols,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
  )
  pdf(out_f,width=width,height=height)
  print(hm)
  dev.off()
}

Figure4F <- function(df,anno,clusters,CRE_groups,which_anno_cluster,sig_df,cols,out_f,width,height){
  df <- df[df$enh_type%in%CRE_groups,]
  df$cluster <- anno[match(df$MPRA_name,anno$MPRA_name),which_anno_cluster]
  colnames(df) <- c('CRE','MPRA_name','NSC','IPC','PN','cluster')
  df <- df[df$cluster%in%clusters|df$CRE=='scrContr'|df$CRE=='WTnoCor',]
  df$cluster[df$CRE=='scrContr'] <- 'scrambled'
  df$cluster[df$CRE=='WTnoCor'] <- 'noCor'
  df$cluster <- factor(df$cluster,levels=c('scrambled','noCor',clusters))
  df <- df[df$CRE=='WTposCor',]
  sig_df <- sig_df[sig_df$MPRA_name%in%df$MPRA_name,]
  sig_df$sig <- 'no'
  sig_df$sig[rowMins(as.matrix(sig_df[,-c(1,ncol(sig_df))]),na.rm=T)<=0.05] <- 'yes'
  df_sig <- df[sig_df$sig=='yes',]
  mat <- t(scale(t(df_sig[,3:5])))
  mat_k <- cluster::pam(mat,6)
  
  la <- rowAnnotation(foo = anno_block(gp = gpar(fill = rev(colorpalette('brbg',6)))))
  hm <- Heatmap(df_sig[,3:5],name='MPRA signal',row_split = factor(mat_k$clustering,levels=c(5,3,1,2,4,6)),cluster_columns = F,show_row_names = F,cluster_row_slices = F,cluster_rows = F,left_annotation = la,
                col = cols,column_names_rot = 0,column_names_centered = T,border=T,row_title = NULL)
  pdf(out_f,width=width,height=height)
  print(hm)
  dev.off()
  # Save the clustering output to ensure reproducibility, renaming the clusters 1:6 for sig and 7 for non-sig
  df$Vars <- rowVars(as.matrix(df[,3:5]))
  df$sig <- 'no'
  df$sig[sig_df$sig=='yes'] <- 'yes'
  df$cluster <- 7
  df$cluster[sig_df$sig=='yes'] <- c(3,4,2,5,1,6)[mat_k$clustering]
  write.table(df,'results/immunoMPRA_sigDF_clustered.tsv',col.names=T,row.names=F,sep='\t',quote=F)
}


Figure4G <- function(df,pwms,cre_size,genome,mcparams,FDR.cutoff,enr.cutoff,out_f,height,width){
  pwms <- toPWM(pwms)
  pwms <- pwms[!duplicated(names(pwms))]
  peaks <- gsub('WTposCor_','',df$MPRA_name)
  peaks <- featureToGR(peaks)
  names(peaks) <- df$MPRA_name
  peaks <- resize(peaks,cre_size,fix='center')          #Change coords to actual MPRA coords
  peakseqs <- getSeq(genome, peaks)
  se <- calcBinnedMotifEnr(seqs = peakseqs, bins = factor(df$cluster), motifs = pwms,BPPARAM = mcparams)
  sel1 <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  sel2 <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  seSel <- se[sel1&sel2, ]
  #SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm,BPPARAM = mcparams)
  hcl <- hclust(dist(assay(seSel,"log2enr")), method = "single")
  pdf(out_f,height=height,width=width)
  plotMotifHeatmaps(x = seSel, which.plots = c("log2enr"), width = 1.8,cluster=hcl,
                    show_dendrogram = TRUE, show_seqlogo = TRUE,
                    width.seqlogo = 1.2)
  dev.off()
}

#Supplementary Functions

FigureS6B<-function(MPRA_data,MPRA_design,cutoff,groups,newnames,out_f,colours){
  file <- read.table(MPRA_data)
  colnames(file) <- c("seq", "number")
  file <- file[file$number > cutoff,]
  file$label <- str_split_fixed(file$seq, "_", 2)[,1]
  design <- read.table(MPRA_design, sep="\t")
  stats <- file %>% 
    group_by(label) %>% 
    dplyr::summarize(mean = mean(number), median = median(number))
  counter_design <- table(design$V2)
  counter_file <- table(file$label)
  counters <- data.frame(counter_file,counter_design[names(counter_design) %in% names(counter_file)],names = names(counter_file))
  counters <- counters[,c(2,4)]
  colnames(counters) <- c("obs", "exp")
  counters$fraction <- round(counters$obs / counters$exp, 3)
  stats <- cbind.data.frame(counters, stats)
  stats_subset <- stats[stats$label %in% groups,]
  mutated<-data.frame(obs = sum(stats[grep(stats$label,pattern='mut'),'obs']), exp= sum(stats[grep(stats$label,pattern='mut'),'exp']), 
                      fraction = (sum(stats[grep(stats$label,pattern='mut'),'obs'])/sum(stats[grep(stats$label,pattern='mut'),'exp'])),
                      label = 'Mutated motifs' , mean= NA, median= NA)
  stats_subset<-rbind(stats_subset,mutated)
  stats_subset$label <- factor(stats_subset$label, levels = groups, ordered=T)
  p <- ggplot(data = stats_subset, aes(x = label, y = fraction, fill = label)) + geom_bar(stat="identity",colour="black")+
    ylab("Fraction of recovered CREs") +  
    scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0, 0)) + scale_x_discrete(labels= newnames) +
    theme(axis.title.x=element_blank(),legend.title = element_blank(),axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), legend.position = 'none')+
    scale_fill_manual(values=colours)+geom_hline(yintercept=1.0, linetype="dashed",color = "blue", size=1)
  pdf(out_f, height= 5, width= 4)
  print(p)
  dev.off()
}

FigureS6C<-function(MPRA_data,MPRA_design,cutoff,groups,newnames,out_f,colours){
  file <- read.table(MPRA_data)
  design <- read.table(MPRA_design, sep="\t")
  colnames(file) <- c("seq", "number")
  file <- file[file$number > cutoff,]
  file$label <- str_split_fixed(file$seq, "_", 2)[,1]     
  file_subset <- file[file$label %in% groups,]
  mutated<-file[grep(file$label,pattern='mut'),]
  mutated$label<-groups[5]
  file_subset<-rbind(file_subset,mutated)
  file_subset$label <- factor(file_subset$label, levels =groups, ordered=T)
  p <- ggplot(file_subset, aes(x = label, y = log2(number), fill=label)) +scale_fill_manual(values=colours)+
    geom_violin() +
    scale_x_discrete(labels= newnames)+
    geom_boxplot(width=0.1, outlier.shape = NA) +
    theme(axis.title.x=element_blank(),legend.title = element_blank(),axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), legend.position = 'none')+ylab("log2 (Number of Barcodes per CRE)")+
    scale_y_continuous(breaks = seq(0,10,2)) 
  pdf(out_f,height= 5, width= 4)
  plot(p)
  dev.off()
}


FigureS6F<-function(MPRA_data,MPRA_design,groups,colours,newnames,out_f){
  lables<- read.table(MPRA_design, sep="\t")
  out <- strsplit(as.character(lables$V1),'_') 
  lables<-data.frame(lables, do.call(rbind, out))
  dnaCounts <- read.delim(MPRA_data,header=T,row.names = 'seq_id')
  #split in replicates
  rep1<-dnaCounts[1:(1+1)!=(1+1)]
  rep1<-subset(rep1,rowSums(rep1,na.rm=T)>0)
  out<-strsplit(as.character(row.names(rep1)),'_')
  rep1<-data.frame(rep1, do.call(rbind, out))[,c('X1','X2','X3')]
  rep2<-dnaCounts[1:(1+1)==(1+1)]
  rep2<-subset(rep2,rowSums(rep2,na.rm=T)>0)
  out<-strsplit(as.character(row.names(rep2)),'_')
  rep2<-data.frame(rep2, do.call(rbind, out))[,c('X1','X2','X3')]
  name_Vec<-unique(grep(lables$V2, pattern="",value = TRUE))
  rep_vec<-c("rep1","rep2")
  to_plot<-data.frame()
  for (k in 1:length(rep_vec)) {
    df<-data.frame()
    filtered_barcodes<-get(rep_vec[k])
    for (i in 1:length(name_Vec)) {
      term<-name_Vec[i]
      number<-as.data.frame(nrow(subset(filtered_barcodes,filtered_barcodes$X1 == term)) / nrow(subset(lables,lables$V2 == term)))*100
      colnames(number)<-"Percentage"
      number$type<-term
      df<-rbind(df,number)
    }
    df$replicate<-rep_vec[k]
    to_plot<-rbind(to_plot,df)
  }
  to_plot <- to_plot[to_plot$type %in% groups,]
  to_plot$type <- factor(to_plot$type,levels =groups)
  to_plot$Percentage<-to_plot$Percentage/100
  p<-ggplot(to_plot, aes(y=Percentage, x=type,fill=type))+geom_bar(stat = "summary", fun = "mean", colour='black')+geom_point(aes(y=Percentage,x=type),position = position_jitter(w = 0.1, h = 0), size=2)+
    ylab(label = "Fraction of recovered CREs")+scale_x_discrete(labels= newnames)+scale_fill_manual(values=colours)+
    theme(axis.title.x=element_blank(),legend.title = element_blank(),legend.position = 'none')+
    scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0, 0))+geom_hline(yintercept=1.0, linetype="dashed", color = "blue", size=1)
  pdf(file=out_f,height= 5, width= 3.5)
  print(p)
  dev.off()
}

FigureS6G<-function(MPRA_data,groups,out_f){
  dnaCounts <- read.delim(MPRA_data,header=T,row.names = 'seq_id')
  rep1<-dnaCounts[1:(1+1)!=(1+1)]
  rep1$barcodes<-rowSums(rep1>=1, na.rm = T)
  rep1$replicate<-'Replicate 1'
  rep2<-dnaCounts[1:(1+1)==(1+1)]
  rep2$barcodes<-rowSums(rep2>=1, na.rm = T)
  rep2$replicate<-'Replicate 2'
  to_plot<-rbind(rep1[,c('barcodes','replicate')],rep2[,c('barcodes','replicate')])
  out<-strsplit(as.character(row.names(to_plot)),'_')
  to_plot<-data.frame(to_plot, do.call(rbind, out))
  to_plot <- to_plot[to_plot$X1 %in% groups,]
  p <- ggplot(to_plot, aes(x = replicate, y = log2(barcodes), fill=replicate))+ geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    ylab("log2 (Number of Barcodes per CRE)") +
    theme(axis.title.x=element_blank(),legend.position = "none") +scale_y_continuous(breaks = seq(0,10,2)) 
  pdf(file=out_f,height= 5, width= 3)
  print(p)
  dev.off()
}

FigureS6HIJ<-function(dnaCounts,rnaCounts,removel,whattoplot,out_f){
  dnaCounts <- read.delim(dnaCounts,header=T,row.names = 'seq_id')
  rnaCounts <- read.delim(rnaCounts,header=T,row.names = 'seq_id')
  #split in replicates + normalize by total number of reads
  rep1_DNA<-dnaCounts[1:(1+1)!=(1+1)]
  rep1_DNA[is.na(rep1_DNA)]<-0
  rep1_RNA<-rnaCounts[1:(1+1)!=(1+1)]
  rep1_RNA[is.na(rep1_RNA)]<-0
  #removal of all barcodes that dont have at leats 1 read both in RNA and DNA
  if (removel==T) {
    rep1_RNA[rep1_DNA==0]<-0
    rep1_DNA[rep1_RNA==0]<-0
  }
  rep1_DNA<-rep1_DNA[]/sum(rep1_DNA)
  rep1_DNA$DNA_read_sum_rep1<-rowSums(rep1_DNA,na.rm=T)
  rep1_DNA$DNA_barcodes_rep1<-rowSums(rep1_DNA>0, na.rm = T)-1
  rep1_RNA<-rep1_RNA[]/sum(rep1_RNA)
  rep1_RNA$RNA_read_sum_rep1<-rowSums(rep1_RNA,na.rm=T)
  rep1_RNA$RNA_barcodes_rep1<-rowSums(rep1_RNA>0, na.rm = T)-1
  ###for replicate 2
  rep2_DNA<-dnaCounts[1:(1+1)==(1+1)]
  rep2_DNA[is.na(rep2_DNA)]<-0
  rep2_RNA<-rnaCounts[1:(1+1)==(1+1)]
  rep2_RNA[is.na(rep2_RNA)]<-0
  #removal of all barcodes that dont have at leats 1 read both in RNA and DNA
  if (removel==T) {
    rep2_RNA[rep2_DNA==0]<-0
    rep2_DNA[rep2_RNA==0]<-0
  }
  rep2_DNA<-rep2_DNA[]/sum(rep2_DNA)
  rep2_DNA$DNA_read_sum_rep2<-rowSums(rep2_DNA,na.rm=T)
  rep2_DNA$DNA_barcodes_rep2<-rowSums(rep2_DNA>0, na.rm = T)-1
  rep2_RNA<-rep2_RNA[]/sum(rep2_RNA)
  rep2_RNA$RNA_read_sum_rep2<-rowSums(rep2_RNA,na.rm=T)
  rep2_RNA$RNA_barcodes_rep2<-rowSums(rep2_RNA>0, na.rm = T)-1
  df<-as.data.frame(cbind(rep1_DNA[,c('DNA_read_sum_rep1','DNA_barcodes_rep1')],rep1_RNA[,c('RNA_read_sum_rep1','RNA_barcodes_rep1')],rep2_DNA[,c('DNA_read_sum_rep2','DNA_barcodes_rep2')],rep2_RNA[,c('RNA_read_sum_rep2','RNA_barcodes_rep2')]))
  row.names(df)<-row.names(dnaCounts)
  #####RNA and DNA correlations
  if (whattoplot=='dna') {
    p1<-ggplot(df,aes(x=DNA_read_sum_rep1,y=DNA_read_sum_rep2))+geom_point()+stat_cor(method='pearson')+xlab(label = 'Normalized DNA reads replicate 1')+ylab(label = 'Normalized DNA reads replicate 2')
    ggsave(file=out_f,plot=p1,units="in",dpi=320, device='png',height= 5, width= 5)
  }
  if (whattoplot=='rna') {
    p2<-ggplot(df,aes(x=RNA_read_sum_rep1,y=RNA_read_sum_rep2))+geom_point()+stat_cor(method='pearson')+xlab(label = 'Normalized RNA reads replicate 1')+ylab(label = 'Normalized RNA reads replicate 2')
    ggsave(file=out_f,plot=p2,units="in",dpi=320, device='png',height= 5, width= 5)
  }
  #####calculate the ration of sums
  if (whattoplot=='ratio') {
    df$ratio_sums_rep1<-(df$RNA_read_sum_rep1+1)/(df$DNA_read_sum_rep1+1)
    df$ratio_sums_rep2<-(df$RNA_read_sum_rep2+1)/(df$DNA_read_sum_rep2+1)
    p3<-ggplot(df,aes(x=log2(ratio_sums_rep1),y=log2(ratio_sums_rep2)))+geom_point()+stat_cor(method='pearson')+xlab(label = 'Log2(ratio of sums) replicate 1')+ylab(label = 'Log2(ratio of sums) replicate 2')
    ggsave(file=out_f,plot=p3,units="in",dpi=320, device='png',height= 5, width= 5)
  }
}

FigureS6K<-function(v_enhancer,MPRA_data,group,out_f){
  v_enhancer<-read.table(v_enhancer, sep='\t', header=F)
  gr_vista_enhancer<-makeGRangesFromDataFrame(v_enhancer, keep.extra.columns=F,ignore.strand=T,seqnames.field="V1",start.field="V2",end.field="V3")
  df<-read.table(MPRA_data, sep='\t', header=T)
  temp<-df[!is.na(df$chr) & !is.na(df$start),]
  gr_MPRA<-makeGRangesFromDataFrame(temp, keep.extra.columns=T,ignore.strand=T,seqnames.field="chr",start.field="start",end.field="end")
  overlap<-as.data.frame(gr_MPRA[gr_MPRA %over% gr_vista_enhancer])
  pos_cor<-subset(overlap, overlap$CDS_type == group)
  pos_cor$type<-'Vista enhancer'
  mut_cor<-subset(df, df$CDS_type == 'scrContr') 
  mut_cor$type<-'scrambled'
  for_ploting<-rbind(pos_cor[,c("CDS_type","alpha",'type')],mut_cor[,c("CDS_type","alpha",'type')])
  p<-ggplot(for_ploting,aes(x=type,y=alpha, fill=type))+
    stat_compare_means(ref.group="scrambled", method = "wilcox.test", paired = F,label= "p",size = 3.5,label.y = quantile(for_ploting$alpha,0.98,na.rm=T))+xlab(label = "")+
    geom_boxplot(width=0.8,outlier.shape = NA)+coord_cartesian(ylim = quantile(for_ploting$alpha, c(0.05, 0.98),na.rm=T))+ ylab(label = "MPRA signal")+
    scale_fill_manual(values=c("white",'grey80'))+theme(legend.position="none")
  pdf(out_f, width = 3, height = 5)
  print(p)
  dev.off()    
}

FigureS6L<-function(MPRA_data,annotation,group,out_f){
  anno<-read.table(annotation, sep='\t',header=T)
  df<-read.table(MPRA_data, sep='\t', header=T)
  df_subsetted<-subset(df, df$CDS_type == group[1] | df$CDS_type == group[2])
  df_subsetted$Enhancer_celltype<-anno$cluster[match(df_subsetted$name,anno$MPRA_name)]
  df_subsetted$Enhancer_celltype<-as.character(df_subsetted$Enhancer_celltype)
  df_subsetted$Enhancer_celltype[is.na(df_subsetted$Enhancer_celltype)] <-'scrambled'
  df_subsetted$Enhancer_celltype <- factor(df_subsetted$Enhancer_celltype,levels = c("scrambled", "NSC", "IPC", "PN1","PN2",'CR'))
  p<-ggplot(df_subsetted,aes(x=Enhancer_celltype,y=alpha, fill=Enhancer_celltype))+geom_boxplot(width=0.8,outlier.shape = NA)+
    ylab(label = "MPRA signal")+theme(legend.position="none")+scale_fill_manual(values=c("white",ATAC_cluster_colors[1:5]))+xlab(label = "")+coord_cartesian(ylim = quantile(df_subsetted$alpha, c(0.05, 0.96),na.rm=T))
  pdf(out_f)
  print(p)
  dev.off()
}

FigureS6O<-function(data,out_f){
  df<-read.table(data,sep='\t', header=T)
  df$Celltype <- factor(df$Celltype,levels = c("NSC","IPC","PN"))
  for (i in df$Gene){
    df_sub<-subset(df, df$Gene == i)
    p <- ggplot(df_sub, aes(x=Gene, y=Mean, fill=Celltype)) + scale_fill_manual(values = RNA_cluster_colors[c(1,3,6)]) +
      geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Relative Fold Change') +
      geom_errorbar(aes(ymin=Mean-STDEV, ymax=Mean+STDEV), width=.2,position=position_dodge(.9))
    p <- p + theme(legend.position='none') + geom_point(data = df_sub,aes(x=Gene, y=Value, fill=Celltype),position = position_jitterdodge(dodge.width =.9,jitter.width = 0.25,seed = 42),size=2)
    pdf(file=paste0(out_f,'_',i,'.pdf'), width=2,height=5)
    print(p)
    dev.off()
  }
}

FigureS6P<-function(MPRA_data_f,MPRA_design,cutoff,out_f){
  names<-c('scrambled',"WTnoCor",'WTposCor',"WTnegCor",'Mutated motifs')
  lables<- read.table(MPRA_design, sep="\t")
  out <- strsplit(as.character(lables$V1),'_') 
  lables<-data.frame(lables, do.call(rbind, out))
  celltypes<-c('PAX6','EOMES','TUBB3')
  to_plot<-data.frame()
  for (c in 1:length(celltypes)) {
    celltype<-celltypes[c]
    dnaCounts <- read.delim(paste0(MPRA_data_f,'immunoMPRA_',celltype,'_mpranalyze/dna_counts.tsv'),header=T,row.names = 'seq_id')
    rep1<-dnaCounts[1:(1+1)!=(1+1)]
    rep1$barcodes<-rowSums(rep1>=1, na.rm = T)
    rep1<-subset(rep1,rep1$barcodes >= cutoff)
    out<-strsplit(as.character(row.names(rep1)),'_')
    rep1<-data.frame(rep1, do.call(rbind, out))[,c('X1','X2','X3')]
    rep2<-dnaCounts[1:(1+1)==(1+1)]
    rep2$barcodes<-rowSums(rep2>=1, na.rm = T)
    rep2<-subset(rep2,rep2$barcodes >= cutoff)
    out<-strsplit(as.character(row.names(rep2)),'_')
    rep2<-data.frame(rep2, do.call(rbind, out))[,c('X1','X2','X3')]
    name_Vec<-unique(grep(lables$V2, pattern="",value = TRUE))
    rep_vec<-c("rep1","rep2")
    for (k in 1:length(rep_vec)) {
      df<-data.frame()
      filtered_barcodes<-get(rep_vec[k])
      for (i in 1:length(name_Vec)) {
        term<-name_Vec[i]
        number<-as.data.frame(nrow(subset(filtered_barcodes,filtered_barcodes$X1 == term)) / nrow(subset(lables,lables$V2 == term)))*100
        colnames(number)<-"Percentage"
        number$type<-term
        df<-rbind(df,number)
      }
      df$replicate<-rep_vec[k]
      df$celltype<-celltype
      to_plot<-rbind(to_plot,df)
    }
  }
  to_plot[grepl("scrContr",to_plot$type),'class']<-names[1]
  to_plot[grepl("WTposCor",to_plot$type),'class']<-names[3]
  to_plot[grepl("WTnoCor",to_plot$type),'class']<-names[2]
  to_plot[grepl("WTnegCor",to_plot$type),'class']<-names[4]
  to_plot[grepl("mut",to_plot$type),'class']<-names[5]
  to_plot<-na.omit(to_plot)
  to_plot<-to_plot %>% dplyr::group_by(class,celltype) %>% mutate(Mean=mean(Percentage))
  to_plot<-to_plot %>% dplyr::group_by(class,celltype,replicate) %>% mutate(Mean_replicate=mean(Percentage))
  to_plot$class <- factor(to_plot$class,levels = names)
  to_plot$celltype <- factor(to_plot$celltype,levels = c("PAX6", "EOMES", "TUBB3"))
  to_plot$Mean<-to_plot$Mean/100
  to_plot$Percentage<-to_plot$Percentage/100
  to_plot$Mean_replicate<-to_plot$Mean_replicate/100
  to_plot$type<-NULL
  to_plot$Percentage<-NULL
  to_plot<-unique(to_plot)
  p<-ggplot(to_plot, aes(y=Mean,x=class,fill=celltype))+geom_bar(stat="identity",position="dodge")+
    xlab(label = 'CRS Type')+ylab(label = paste0('Percentage recovery')) +
    ylab("Fraction of recovered CREs") +  
    scale_y_continuous(breaks = seq(0,1,0.1))+
    theme(axis.title.x=element_blank(),legend.title = element_blank(),axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))+
    scale_fill_manual(values=RNA_cluster_colors[c(1,3,6)])+geom_hline(yintercept=1.0, linetype="dashed",color = "blue", size=1)
  p<-p+ geom_point(data = to_plot,aes(x=class, y=Mean_replicate, fill=celltype),position = position_jitterdodge(dodge.width =.9,jitter.width = 0.25,seed = 42),size=2)+theme(legend.position='none') 
  pdf(file=out_f, width=4,height=5)
  print(p)
  dev.off()  
}

FigureS6Q<-function(MPRA_data_f,out_f){
  celltypes<-c('PAX6','EOMES','TUBB3')
  names<-c('NSC','IPC','PN')
  to_plot<-data.frame()
  for (c in 1:length(celltypes)) {
    celltype<-celltypes[c]
    dnaCounts <- read.delim(paste0(MPRA_data_f,'immunoMPRA_',celltype,'_mpranalyze/dna_counts.tsv'),header=T,row.names = 'seq_id')
    rep1<-dnaCounts[1:(1+1)!=(1+1)]
    rep1$barcodes<-rowSums(rep1>=1, na.rm = T)
    rep1$replicate<-'Rep1'
    rep2<-dnaCounts[1:(1+1)==(1+1)]
    rep2$barcodes<-rowSums(rep2>=1, na.rm = T)
    rep2$replicate<-'Rep2'
    df<-rbind(rep1[,c('barcodes','replicate')],rep2[,c('barcodes','replicate')])
    df$celltype<-names[c]
    to_plot<-rbind(to_plot,df)
  }
  to_plot$celltype <- factor(to_plot$celltype,levels = names)
  p <- ggplot(to_plot, aes(x = celltype, y = log2(barcodes), fill=celltype))+
    geom_violin() + scale_fill_manual(values =RNA_cluster_colors[c(1,3,6)])+geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_classic() +ylab("log2 (Number of Barcodes per CRE)") +
    theme(axis.title.x=element_blank(),legend.position = "none") +scale_y_continuous(breaks = seq(0,10,2)) 
  pdf(file=out_f, height = 5, width=3)
  print(p)
  dev.off() 
}


FigureS6R<-function(MPRA_data_f,removel,out_f){
  celltypes<-c('PAX6','EOMES','TUBB3')
  merged_df<-data.frame
  for (i in 1:length(celltypes)) {
    celltype<-celltypes[i]
    dnaCounts <- read.delim(paste0(MPRA_data_f,'immunoMPRA_',celltype,'_mpranalyze/dna_counts.tsv'),header=T,row.names = 'seq_id')
    rnaCounts <- read.delim(paste0(MPRA_data_f,'immunoMPRA_',celltype,'_mpranalyze/rna_counts.tsv'),header=T,row.names = 'seq_id')
    #split in replicates + normalize by total number of reads
    rep1_DNA<-dnaCounts[1:(1+1)!=(1+1)]
    rep1_DNA[is.na(rep1_DNA)]<-0
    rep1_RNA<-rnaCounts[1:(1+1)!=(1+1)]
    rep1_RNA[is.na(rep1_RNA)]<-0
    if (removel==T) {
      rep1_RNA[rep1_DNA==0]<-0
      rep1_DNA[rep1_RNA==0]<-0
    }
    rep1_DNA<-rep1_DNA[]/sum(rep1_DNA)
    rep1_DNA$DNA_read_sum_rep1<-rowSums(rep1_DNA,na.rm=T)
    rep1_DNA$DNA_barcodes_rep1<-rowSums(rep1_DNA>0, na.rm = T)-1
    rep1_RNA<-rep1_RNA[]/sum(rep1_RNA)
    rep1_RNA$RNA_read_sum_rep1<-rowSums(rep1_RNA,na.rm=T)
    rep1_RNA$RNA_barcodes_rep1<-rowSums(rep1_RNA>0, na.rm = T)-1
    ###for replicate 2
    rep2_DNA<-dnaCounts[1:(1+1)==(1+1)]
    rep2_DNA[is.na(rep2_DNA)]<-0
    rep2_RNA<-rnaCounts[1:(1+1)==(1+1)]
    rep2_RNA[is.na(rep2_RNA)]<-0
    if (removel==T) {
      rep2_RNA[rep2_DNA==0]<-0
      rep2_DNA[rep2_RNA==0]<-0
    }
    rep2_DNA<-rep2_DNA[]/sum(rep2_DNA)
    rep2_DNA$DNA_read_sum_rep2<-rowSums(rep2_DNA,na.rm=T)
    rep2_DNA$DNA_barcodes_rep2<-rowSums(rep2_DNA>0, na.rm = T)-1
    rep2_RNA<-rep2_RNA[]/sum(rep2_RNA)
    rep2_RNA$RNA_read_sum_rep2<-rowSums(rep2_RNA,na.rm=T)
    rep2_RNA$RNA_barcodes_rep2<-rowSums(rep2_RNA>0, na.rm = T)-1
    df<-as.data.frame(cbind(rep1_DNA[,c('DNA_read_sum_rep1','DNA_barcodes_rep1')],rep1_RNA[,c('RNA_read_sum_rep1','RNA_barcodes_rep1')],rep2_DNA[,c('DNA_read_sum_rep2','DNA_barcodes_rep2')],rep2_RNA[,c('RNA_read_sum_rep2','RNA_barcodes_rep2')]))
    row.names(df)<-row.names(dnaCounts)
    #####calculate the sum ratio 
    df$ratio_sums_rep1<-(df$RNA_read_sum_rep1+1)/(df$DNA_read_sum_rep1+1)
    df$ratio_sums_rep2<-(df$RNA_read_sum_rep2+1)/(df$DNA_read_sum_rep2+1)
    df$mean_ration_sums<-rowMeans(df[,c("ratio_sums_rep1", "ratio_sums_rep2")], na.rm=TRUE)
    df$name<-row.names(df)
    colnames(df)<-c(paste0(celltype,'_DNA_read_sum_rep1'),paste0(celltype,'_DNA_barcodes_rep1'),
                    paste0(celltype,'_RNA_read_sum_rep1'),paste0(celltype,'_RNA_barcodes_rep1'),
                    paste0(celltype,'_DNA_read_sum_rep2'),paste0(celltype,'_DNA_barcodes_rep2'),
                    paste0(celltype,'_RNA_read_sum_rep2'),paste0(celltype,'_RNA_barcodes_rep2'),
                    paste0(celltype,'_ratio_sums_rep1'),paste0(celltype,'_ratio_sums_rep2'),paste0(celltype,'_mean_ration_sums'),paste0(celltype,'_name')
    )
    ifelse(i==1,merged_df<-df,merged_df<-merge(merged_df,df,by.x=paste0(celltypes[1],'_name'), by.y=paste0(celltype,'_name'), all.x=T)) ###keep only CRS which are in all celltypes !!! otherwise all=T 
  }
  merged_df$CDS_type <- str_split_fixed(merged_df$PAX6_name,'_chr',2)[,1]
  merged_df[grepl("scrContr",merged_df$CDS_type),'CDS_type']<-'scrambled'
  out <- strsplit(as.character(merged_df$PAX6_name),'_') 
  out<-data.frame(do.call(rbind, out))
  colnames(out)<-c('name','chr','start','end')
  merged_df<-cbind(merged_df,out[,c(2,3,4)])
  merged_df[grepl("scrambled",merged_df$CDS_type),c('chr','start','end')]<-NA
  #####Generate PCA plot
  df<- merged_df[,c('PAX6_ratio_sums_rep1','PAX6_ratio_sums_rep2','EOMES_ratio_sums_rep1','EOMES_ratio_sums_rep2','TUBB3_ratio_sums_rep1','TUBB3_ratio_sums_rep2')]
  row.names(df)<-row.names(merged_df)
  df<-na.omit(df)
  df<-as.data.frame(as.matrix(t(df)))
  df_pca <- prcomp(df)
  df_pca_out<-as.data.frame(df_pca$x)
  df_pca_out$Celltype<-c('NSC','NSC','IPC','IPC','PN','PN')
  df_pca_out$replicate<-c('rep1','rep2','rep1','rep2','rep1','rep2')
  percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
  percentage <- paste(colnames(df_pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  df_pca_out$Celltype<- factor(df_pca_out$Celltype, levels = c("NSC","IPC","PN"))
  p<-ggplot(df_pca_out,aes(x=PC1,y=PC2,color=Celltype))+geom_point(size = 5)+xlab(percentage[1]) + ylab(percentage[2])+ 
    scale_color_manual(values = RNA_cluster_colors[c(1,3,6)])+theme(legend.position = c(0.75,0.9),legend.title = element_blank(),legend.text = element_text( size = 15))
  pdf(file=out_f, width=5,height=5)
  print(p)
  dev.off()
}

FigureS6ST<- function(MPRA_data,annotation,clusters=c('NSC','IPC','PN1','PN2'),groups,out_f){
  df=as.data.frame(vroom::vroom(MPRA_data))
  anno=as.data.frame(vroom::vroom(annotation))
  which_anno_cluster='rna_cluster_full' #rna_cluster_simple or 'cluster' or 'rna_cluster_full
  df <- df[,c('enh_type','MPRA_name','mad.score_NSC','mad.score_IPC','mad.score_PN')]
  df <- df[df$enh_type%in%groups,]
  df$cluster <- anno[match(df$MPRA_name,anno$MPRA_name),which_anno_cluster]
  colnames(df) <- c('CRE','MPRA_name','NSC','IPC','PN','cluster')
  df <- df[df$cluster%in%clusters|df$CRE=='scrContr'|df$CRE=='WTnoCor',]
  df$cluster[df$CRE=='scrContr'] <- 'scrambled'
  df$cluster[df$CRE=='WTnoCor'] <- 'noCor'
  df$cluster <- factor(df$cluster,levels=c('scrambled','noCor',clusters))
  simple_df <- df[,-c(1:2)]
  mat <- ddply(simple_df,.(cluster),function(x){
    return(colMedians(as.matrix(x[,1:3])))
  })
  colnames(mat) <- c('cluster','NSC','IPC','PN')
  row.names(mat) <- mat$cluster
  mat <- round(as.matrix(mat[,-1]),2)
  hm <- Heatmap(mat,name = "MPRA signal",cluster_rows = F,cluster_columns = F,show_row_names = T,border = 'black',column_names_rot = 0,column_names_centered = T,
                col = colorRamp2(c(0,0.13,0.15,seq(0.25,1.12,length.out=9)), rev(colorpalette('reds',12))),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))}
  )
  pdf(file=out_f, width= 4, height=5)
  print(hm)
  dev.off()
}



### Plot Figures #######

Figure4B(df=res_bulk,CRE_groups=c('scrContr','WTnoCor','WTposCor'),CRE_groupNames=c('scrambled','noCor','posCor'),
              cols=colorRampPalette(c('white','gray30'))(3),ylab_n ='MPRA signal',
              out_f='plots/figures/Figure4B.pdf',height=4,width=4)

Figure4C(df=res_bulk,anno=mpra_annot,CRE_groups=c('scrContr','WTnoCor','WTposCor'),
              CRE_groupNames=c('scrambled','noCor','posCor'),cols=rev(colorpalette('rdbu',11)),
              ylab_n='MPRA signal',max_sig=2,
              features=c("Sema3c:-74670","Sox11:-495899","Sstr2:-179704","Sox5:152614","Sema6d:-248213","Tfap2c:450163","Dmrta2:231841","Neurog2:13167","Rnd2:148829",'Rnd2:9546',"Neurod6:34596","Sox2:434656","Emx2:-36039","Gli3-241669","Gli3:81089",'Pax6:405689'),
              anno.size=3,point.size=1,out_f='plots/figures/Figure4C.pdf',height=4,width=6)

Figure4E(df=as.data.frame(res[,c('enh_type','MPRA_name','mad.score_NSC','mad.score_IPC','mad.score_PN')]),out_f='plots/figures/Figure4E.pdf',height=3.5,width=4,
         anno=as.data.frame(mpra_annot),clusters=c('NSC','IPC','PN1','PN2'),CRE_groups=c('scrContr','WTposCor'),
         which_anno_cluster='cluster',cols=colorRamp2(c(0,0.13,0.15,seq(0.25,1.12,length.out=9)), rev(colorpalette('reds',12))))

Figure4F(df=as.data.frame(res[,c('enh_type','MPRA_name','mad.score_NSC','mad.score_IPC','mad.score_PN')]),
                      sig_df=res[,c('MPRA_name','pval.mad_NSC','pval.mad_IPC','pval.mad_PN')],out_f='plots/figures/Figure4F.pdf',height=5,width=5,
                      anno=as.data.frame(mpra_annot),clusters=c('NSC','IPC','PN1','PN2'),CRE_groups=c('scrContr','WTposCor'),which_anno_cluster='cluster',
                      cols=colorRamp2(c(0,1,seq(1.5,15,length.out=10)), rev(colorpalette('reds',12))))

Figure4G(df=vroom::vroom('results/immunoMPRA_sigDF_clustered.tsv'),pwms=readRDS('data/mpra_combined_pwm.RDS'),
         cre_size=270,genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L),FDR.cutoff=0.9,enr.cutoff=0.5,
         out_f='plots/figures/Figure4G.pdf',height=5,width=7)

mpra_coords <- res$coord[res$enh_type=='WTposCor']
mpra_coords <- stringr::str_split(paste0(mpra_coords), pattern ='_' , n = 3, simplify = TRUE)
mpra_coords <- gintervals(as.character(mpra_coords[,1]),as.numeric(mpra_coords[,2]),as.numeric(mpra_coords[,3]))
mpra_coords <- intervals.normalize(mpra_coords,1000)      #For Visualisation in big regions
mpra_tracks <- c("mpra.WT_NSC","mpra.WT_IPC","mpra.WT_PN")

plotMisha(targetGene='Dll1',outDir='plots/figures/',out_f='Figure4I',upstream=2e5,downstream=8e4,imgPlotScale=8,
          chipNames=c('','','','',''),plot.dgram=F,mpra_coords=mpra_coords,mpra_tracks=mpra_tracks,mpra_colors=colorpalette('matlablike2',11),mpraNames=c('','',''),
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN1','scATAC.E14_PN2','chipseq_RPM.NPC_Neurog2'),mpra_ylim=c(-0.66,6.24),
          chipColors=c(ATAC_cluster_colors[1:4],'black'),arcIntervals=p2glinks$posCor,arcColors=colorRampPalette(archr_colors[[29]]),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE,MPRA=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.4,MPRA=0.3,meth=0.6, domains=0.15, genes=0.7, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))


mpra_coords <- intervals.normalize(mpra_coords,270) #Make coords exact
mpra_tracks <- c("mpra.WT_NSC","mpra.Neurog2mut_NSC","mpra.WT_IPC","mpra.Neurog2mut_IPC","mpra.WT_PN","mpra.Neurog2mut_PN")

plotMisha(targetGene='chr17,15247517,15248617',outDir='plots/figures/',out_f='Figure4J_1',upstream=500,downstream=500,
          chipNames=c('','','','',''),plot.dgram=F,mpra_coords=mpra_coords,mpra_tracks=mpra_tracks,mpra_colors=colorpalette('matlablike2',11),mpraNames=c('WT','MUT','WT','MUT','WT','MUT'),
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN1','scATAC.E14_PN2','chipseq_RPM.CN_NeuroD2'),
          chipColors=c(ATAC_cluster_colors[1:4],'black'),imgPlotScale=2.5,mpra_ylim=c(-0.66,6.24),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,MPRA=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=FALSE,ideogram=FALSE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.6,MPRA=0.2,meth=0.6, domains=0.15, genes=0.7, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

plotMisha(targetGene='chr17,15440466,15440966',outDir='plots/figures/',out_f='Figure4J_2',upstream=500,downstream=500,
          chipNames=c('','','','',''),plot.dgram=F,mpra_coords=mpra_coords,mpra_tracks=mpra_tracks,mpra_colors=colorpalette('matlablike2',11),mpraNames=c('WT','MUT','WT','MUT','WT','MUT'),
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN1','scATAC.E14_PN2','chipseq_RPM.CN_NeuroD2'),
          chipColors=c(ATAC_cluster_colors[1:4],'black'),imgPlotScale=2.5,mpra_ylim=c(-0.66,6.24),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,MPRA=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=FALSE,ideogram=FALSE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.6,MPRA=0.2,meth=0.6, domains=0.15, genes=0.7, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

pdf('plots/figures/Figure4I_scalebar.pdf',height=4,width=2)
par(mar=c(1,4,1,4))
image.scale(as.matrix(NA),zlim=c(-0.66,6.24), col=colorpalette('matlablike2',100),axis.pos=4,adj=1,cex.axis = 1.5)
dev.off()

#Plot Supplementary 

FigureS6B(MPRA_data='data/MPRA/association_combined_both_filtered.tsv',
          MPRA_design='data/MPRA/MPRA_labels.tsv',
          cutoff=1,colours=c("white",'grey95','grey75','grey55','grey35'),
          groups=c("scrContr", "WTnoCor", "WTposCor", 'WTnegCor' ,'Mutated motifs'),
          newnames=c("scrambled", "WTnoCor", "WTposCor","WTnegCor","Mutated motifs"),
          out_f="plots/figures/FigureS6B.pdf")

FigureS6C(MPRA_data='data/MPRA/association_combined_both_filtered.tsv',
          MPRA_design='data/MPRA/MPRA_labels.tsv',
          cutoff=1,colours=c("white",'grey95','grey75','grey55','grey35'),
          groups=c("scrContr", "WTnoCor", "WTposCor", 'WTnegCor' ,'Mutated motifs'),
          newnames=c("scrambled", "WTnoCor", "WTposCor","WTnegCor","Mutated motifs"),
          out_f="plots/figures/FigureS6C.pdf")

FigureS6F(MPRA_data='data/MPRA/bulk_dna_counts.tsv',
          MPRA_design='data/MPRA/MPRA_labels.tsv',
          groups=c("scrContr", "WTnoCor", "WTposCor"),
          colours=c("white",'grey95','grey75','grey55','grey35'),
          newnames=c("scrambled", "WTnoCor", "WTposCor"),
          out_f="plots/figures/FigureS6F.pdf")

FigureS6G(MPRA_data='data/MPRA/bulk_dna_counts.tsv',
          groups=c("scrContr", "WTposCor", "WTnoCor"),
          out_f="plots/figures/FigureS6G.pdf")

FigureS6HIJ(dnaCounts='data/MPRA/bulk_dna_counts.tsv',
            rnaCounts='data/MPRA/bulk_rna_counts.tsv',
            removel=T,
            whattoplot='ratio',
            out_f="plots/figures/FigureS6H.png") 

FigureS6HIJ(dnaCounts='data/MPRA/bulk_dna_counts.tsv',
            rnaCounts='data/MPRA/bulk_rna_counts.tsv',
            removel=T,
            whattoplot='rna', 
            out_f="plots/figures/FigureS6I.png") 

FigureS6HIJ(dnaCounts='data/MPRA/bulk_dna_counts.tsv',
            rnaCounts='data/MPRA/bulk_rna_counts.tsv',
            removel=T,
            whattoplot='dna',
            out_f="plots/figures/FigureS6J.png") 

FigureS6K(v_enhancer="data/MPRA/vista_forebrain_enhancers.bed",
          MPRA_data="data/MPRA/bulk_MPRAnalyse_out.txt",
          group='WTposCor',
          out_f="plots/figures/FigureS6K.pdf")

FigureS6L(MPRA_data="data/MPRA/bulk_MPRAnalyse_out.txt",
          annotation='results/MPRA_anno.tsv',
          group=c('WTposCor','scrContr'),
          out_f="plots/figures/FigureS6L.pdf")

FigureS6O(data='data/MPRA/immunoMPRA_qPCR_results.txt',
          out_f="plots/figures/FigureS6O")

FigureS6P(MPRA_data_f='data/MPRA/MPRAflow_out/',
          MPRA_design='data/MPRA/MPRA_labels.tsv',
          cutoff=1,
          out_f="plots/figures/FigureS6P.pdf")

FigureS6Q(MPRA_data_f='data/MPRA/MPRAflow_out/',
          out_f="plots/figures/FigureS6Q.pdf")

FigureS6R(MPRA_data_f='data/MPRA/MPRAflow_out/',
          removel=T,
          out_f="plots/figures/FigureS6R.pdf")

FigureS6ST(MPRA_data='results/immunoMPRA_res.tsv',
           annotation='results/MPRA_anno.tsv',
           clusters=c('NSC','IPC','PN1','PN2'),groups=c('scrContr','WTposCor'),
           out_f="plots/figures/FigureS6S.pdf")

FigureS6ST(MPRA_data='results/immunoMPRA_res.tsv',
           annotation='results/MPRA_anno.tsv',
           clusters=c('NSC','IPC','PN1','PN2'),groups=c('scrContr','WTnegCor'),
           out_f="plots/figures/FigureS6T.pdf")
