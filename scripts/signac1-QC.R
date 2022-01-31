library(Signac)
library(Seurat)
library(reticulate)
library(ggplot2)
library(LSD)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
#library(future)
#plan("multiprocess", workers = 8)

set.seed(1234)

source('scripts/config.R')
source('scripts/plot_functions.R')
source('scripts/aux_functions.R')

input_dirs <- paste0(atac_sample_dir,sample_names)
sample_list <- list()
cells_vector <- c()

peak_region_fragments_minfilter <- 5000
peak_region_fragments_maxfilter <-  100000
unique_frags_filter <- 5000
pct_reads_in_peaks_filter <- 15
blacklist_ratio_filter <- 0.02
nucleosome_signal_filter <- 10
minFrags <- 1000
filterFrags <- 8000
max_filterFrags <- 120000
filterTSS <- 9
max_filterTSS <- 25
name <- "Cortex"
by <- 'RG'

doQC <- TRUE
doAssays <- TRUE
seurat_file <- 'data/merged_scATAC_filtered.rds'

if(file.exists(seurat_file)){
  sample_list <- readRDS(seurat_file)
}

#-----------------
# Preparing data
#-----------------

for (i in 1:length(input_dirs)){
  if (doQC){
    h5_file <- paste0(input_dirs[i],"/filtered_peak_bc_matrix.h5")
    fragment.path <- paste0(input_dirs[i],"/fragments.tsv.gz")
    barcodes <- paste0(input_dirs[i],"/singlecell.csv")
    counts <- Read10X_h5(filename = h5_file)
    metadata <- read.csv(
      file = barcodes,
      header = TRUE,
      row.names = 1
    )
    
    cortex.atac <- CreateSeuratObject(
      counts = counts,
      assay = 'peaks',
      project = sample_names[i],
      min.cells = 1,
      meta.data = metadata,
      names.field = 2,
      names.delim = '-'
    )
    
    cortex.atac$orig.ident <- as.factor(i)
    Idents(cortex.atac) <- as.factor(i)
    
    cortex.atac <- SetFragments(
      object = cortex.atac,
      file = fragment.path
    )
    
    cortex.atac <- NucleosomeSignal(object = cortex.atac)
    cortex.atac$nucleosome_signal[is.infinite(cortex.atac$nucleosome_signal)] <- 0
    cortex.atac$pct_reads_in_peaks <- cortex.atac$peak_region_fragments / cortex.atac$nCount_peaks * 100
    cortex.atac$blacklist_ratio <- cortex.atac$blacklist_region_fragments / cortex.atac$peak_region_fragments
    
    ### Calculate TSS 
    gene.ranges <- genes(EnsDb.Mmusculus.v79)
    gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
    tss.ranges <- GRanges(
      seqnames = seqnames(gene.ranges),
      ranges = IRanges(start = start(gene.ranges), width = 2),
      strand = strand(gene.ranges)
    )
    seqlevelsStyle(tss.ranges) <- 'UCSC'
    tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
    tss.ranges <- tss.ranges[!(seqnames(tss.ranges)%in%c('chrX','chrY','chrM')),]
    seqlevels(tss.ranges) <- seqlevels(tss.ranges)[1:19]

    pdf(paste0('plots/scATAC/QC/',sample_names[i],'_QCplot1.pdf'),width=12,height=12)
    print(VlnPlot(
      object = cortex.atac,
      features = c('passed_filters', 'blacklist_ratio', 'nucleosome_signal'),
      pt.size = 0.001,ncol = 4) + NoLegend())
    dev.off()
    
    plot2_a <- VlnPlot(
      object = cortex.atac,
      features = 'peak_region_fragments',
      pt.size = 0.001, log = TRUE) + NoLegend()
    plot2_b <- FeatureScatter(cortex.atac,"peak_region_fragments",'nucleosome_signal', pt.size = 0.01)+ NoLegend()
    plot2_c <- FeatureScatter(cortex.atac,"peak_region_fragments",'blacklist_ratio', pt.size = 0.01)+ NoLegend()
  #  plot2_d <- FeatureScatter(cortex.atac,"peak_region_fragments",'TSS.enrichment', pt.size = 0.01) + NoLegend()
    plot2 <- CombinePlots(plots = list(plot2_a,plot2_b,plot2_c), ncol = 2)+ NoLegend()
    
    pdf(paste0('plots/scATAC/QC/',sample_names[i],'_QCplot2.pdf'),width=12,height=12)
    print(plot2)
    dev.off()
    
    cortex.atac$nucleosome_group <- ifelse(cortex.atac$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
    cortex.atac$nucleosome_group[is.na(cortex.atac$nucleosome_group)] <- 'NS2'
    
    pdf(paste0('plots/scATAC/QC/',sample_names[i],'_FragPlot.pdf'),width=12,height=12)
    print(PeriodPlot(object = cortex.atac, group.by = 'nucleosome_group', region = 'chr1-1-249250621'))
    dev.off()
    
    cortex.atac <- subset(cortex.atac,subset=
                            passed_filters > unique_frags_filter&
                            peak_region_fragments > peak_region_fragments_minfilter&
                            peak_region_fragments < peak_region_fragments_maxfilter&
                            pct_reads_in_peaks > pct_reads_in_peaks_filter&
                            blacklist_ratio < blacklist_ratio_filter&
                            nucleosome_signal< nucleosome_signal_filter)
    cells_keep <- colnames(cortex.atac)
    cells_keep <- gsub('1',i,cells_keep)
    cells_vector <- c(cells_vector,cells_keep)
    
    #-----------------
    # Reading Fragment Files
    #-----------------
    message("Reading in fragment files...")
    fragments <- data.frame(readr::read_tsv(fragment.path, col_names=FALSE))
    
    fragments <- GRanges(
      seqnames = fragments[,1], 
      IRanges(fragments[,2]+1, fragments[,3]), 
      RG = fragments[,4], 
      N = fragments[,5]
    )
    
    message("Filtering Lowly Represented Cells...")
    tabRG <- table(fragments$RG)
    keep <- names(tabRG)[which(tabRG >= minFrags)]
    fragments <- fragments[fragments$RG %in% keep,]
    fragments <- sort(sortSeqlevels(fragments))
    
    #-----------------
    # TSS Profile
    #-----------------
    
    gene.ranges <- genes(EnsDb.Mmusculus.v79)
    gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
    tss.ranges <- GRanges(
      seqnames = seqnames(gene.ranges),
      ranges = IRanges(start = start(gene.ranges), width = 2),
      strand = strand(gene.ranges)
    )
    seqlevelsStyle(tss.ranges) <- 'UCSC'
    tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
    tss.ranges <- tss.ranges[!(seqnames(tss.ranges)%in%c('chrX','chrY','chrM')),]
    seqlevels(tss.ranges) <- seqlevels(tss.ranges)[1:19]
    tss.ranges <- tss.ranges %>% resize(., width = 1, fix = "start") %>% unique
    
    tssProfile <- insertionProfileSingles(feature = tss.ranges, fragments = fragments, 
                                          getInsertions = TRUE, batchSize = 1000)
    tssSingles <- tssProfile$dfall
    tssSingles$uniqueFrags <- 0
    tssSingles[names(tabRG),"uniqueFrags"] <- tabRG
    tssSingles$cellCall <- 0
    tssSingles$cellCall[tssSingles$uniqueFrags >= filterFrags &tssSingles$uniqueFrags < max_filterFrags  & tssSingles$enrichment >= filterTSS& tssSingles$enrichment < max_filterTSS] <- 1
    
    #-----------------
    # Plot Stats
    #-----------------
    tssSingles <- tssSingles[complete.cases(tssSingles),]
    nPass  <- sum(tssSingles$cellCall==1)
    nTotal <- sum(tssSingles$uniqueFrags >= 100&tssSingles$enrichment > 1)
    
    row.names(tssSingles) <- gsub('1',i+4,row.names(tssSingles))
    file_f <- paste0("results/scATAC/",sample_names[i],"_FilterCells.txt")
    write.table(tssSingles,file_f)
    
    pdf(paste0('plots/scATAC/QC/',sample_names[i],'_TSSenrichment.pdf'),width=8,height=8)
    tssSingles <- tssSingles[tssSingles$uniqueFrags > 1000&tssSingles$enrichment > 1,]
    heatscatter(log10(tssSingles$uniqueFrags),tssSingles$enrichment,xlab='Unique Fragments',ylab='TSS enrichment',axes=F,main=sample_names[i],xlim=c(3,5.5),ylim=c(0,30))
    abline(h = filterTSS,v = log10(filterFrags),col='black',lty=2)
    abline(h = max_filterTSS,v = log10(max_filterFrags),col='grey',lty=2)
    magaxis(unlog='x',grid=F,tcl=(-0.5))
    dev.off() 
    
  } else {
    cortex.atac <- sample_list[[sample_names[i]]]
  }
} 

write.table(cells_vector,file='results/cells_pass_noTSS.txt')


