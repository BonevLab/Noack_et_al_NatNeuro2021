library(Signac)
library(Seurat)
library(reticulate)
library(ggplot2)
library(LSD)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)
plan("multiprocess", workers = 8)

set.seed(1234)

source('scripts/config.R')
source('scripts/plot_functions.R')
source('scripts/aux_functions.R')

doAssays <- TRUE
seurat_file <- 'data/mergedNorm_scATAC_filtered.RDS'

filterFrags <- 10000
max_filterFrags <- 120000
filterTSS <- 8
max_filterTSS <- 25

cells_vector <- as.vector(read.table('results/cells_pass_noTSS.txt')$x)
cells_vector <- gsub('1','5',cells_vector)
cells_vector <- gsub('2','6',cells_vector)
files  <- list.files(path = 'results/scATAC/',pattern = 'FilterCells',full.names = T)
tables <- lapply(files, read.table, header = TRUE)
combined.df <- do.call(rbind , tables)
combined.df$cellCall <- 0
combined.df$cellCall[combined.df$uniqueFrags >= filterFrags &combined.df$uniqueFrags < max_filterFrags  & combined.df$enrichment >= filterTSS& combined.df$enrichment < max_filterTSS] <- 1

passed_cells <- cells_vector[cells_vector %in% row.names(combined.df)[combined.df$cellCall==1]]

h5_file <- "/home/hpc/bonev/projects/SC/analysis/scATAC/merged_norm/outs/filtered_peak_bc_matrix.h5"
fragment.path <- '/home/hpc/bonev/projects/SC/analysis/scATAC/merged_norm/outs/fragments.tsv.gz'
barcodes <- "/home/hpc/bonev/projects/SC/analysis/scATAC/merged_norm/outs/singlecell.csv"
fragment.path.filtered <- "data/fragmentsNorm_filtered.tsv"

counts <- Read10X_h5(filename = h5_file)
metadata <- read.csv(
  file = barcodes,
  header = TRUE,
  row.names = 1
)

cortex.atac <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata,
  names.field = 2,
  names.delim = '-'
)

cortex.atac <- subset(cortex.atac,cells=passed_cells)

cortex.atac <- SetFragments(
  object = cortex.atac,
  file = fragment.path
)

if (!file.exists(paste0(fragment.path.filtered, '.bgz'))){
  FilterFragments(
    fragment.path = fragment.path,
    cells = colnames(cortex.atac),
    output.path = fragment.path.filtered
  )
}
cortex.atac <- SetFragments(object = cortex.atac, file = paste0(fragment.path.filtered, '.bgz'))

#levels(cortex.atac) <- seq(1:2)
reidents <- sample_names
names(reidents) <- as.character(seq(1:length(sample_names))+4)
cortex.atac <- RenameIdents(cortex.atac,reidents)

### Perform final QC #######
cortex.atac <- NucleosomeSignal(object = cortex.atac,region = 'chr1-10000000-20000000')
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
cortex.atac <- TSSEnrichment(object = cortex.atac, tss.positions = tss.ranges[1:5000])
##########

pdf(paste0('plots/scATAC/QC/','AggrNorm','_TSSvsFragments.pdf'),width=16,height=8)
par(mfrow=c(1,2))
heatscatter(log10(cortex.atac$passed_filters),(cortex.atac$TSS_fragments/cortex.atac$passed_filters*100),xlab='uniqueFrags',ylab='TSS ratio',main='',cor = F,ylim=c(0,60),xlim=c(3,5.5),axes=F)
magaxis(unlog='x',grid=F,tcl=(-0.5))
heatscatter(log10(cortex.atac$passed_filters),cortex.atac$TSS.enrichment,xlab='uniqueFrags',ylab='TSS enrichment',main='',cor = F,ylim=c(0,15),xlim=c(3,5.5),axes=F)
magaxis(unlog='x',grid=F,tcl=(-0.5))
dev.off()

pdf(paste0('plots/scATAC/QC/','AggrNorm','_TSSenrichment.pdf'),width=16,height=8)
cortex.atac$high.tss <- ifelse(cortex.atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(cortex.atac, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.off()
pdf(paste0('plots/scATAC/QC/','AggrNorm','_QCplot1.pdf'),width=18,height=12)
print(VlnPlot(
  object = cortex.atac,
  features = c('passed_filters', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.001,ncol = 4) + NoLegend())
dev.off()
cortex.atac$nucleosome_group <- ifelse(cortex.atac$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
cortex.atac$nucleosome_group[is.na(cortex.atac$nucleosome_group)] <- 'NS2'

pdf(paste0('plots/scATAC/QC/','AggrNorm','_FragPlot.pdf'),width=12,height=12)
print(PeriodPlot(object = cortex.atac, group.by = 'nucleosome_group', region = 'chr1-10000000-20000000'))
dev.off()

if(doAssays){
  bin_matrix <- GenomeBinMatrix(
    fragments = fragment.path,
    genome = seqlengths(BSgenome.Mmusculus.UCSC.mm10),
    binsize = 5000,
    cells = colnames(cortex.atac),
    chunk = 50
  )
  saveRDS(bin_matrix,file='data/bins5000_norm.RDS')
  
  # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
  gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
  seqlevelsStyle(gene.coords) <- 'UCSC'
  genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
  genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
  
  # create a gene by cell matrix
  gene.activities <- FeatureMatrix(
    fragments = fragment.path,
    features = genebodyandpromoter.coords,
    cells = colnames(cortex.atac),
    chunk = 50
  )
  # convert rownames from chromsomal coordinates into gene names
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]
  rownames(gene.activities)[rownames(gene.activities)==''] <- 'Unknown'
  saveRDS(gene.activities,file='data/gene_activites_norm.RDS')
  
  # extract promoter coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
  prom.coords <- promoters(genes(EnsDb.Mmusculus.v79, filter = ~ tx_biotype == "protein_coding"))
  seqlevelsStyle(prom.coords) <- 'UCSC'
  prom.coords <- keepStandardChromosomes(prom.coords, pruning.mode = 'coarse')
  
  # create a gene by cell matrix
  prom.activities <- FeatureMatrix(
    fragments = fragment.path,
    features = prom.coords,
    cells = colnames(cortex.atac),
    chunk = 50
  )
  prom.key <- prom.coords$gene_name
  names(prom.key) <- GRangesToString(grange = prom.coords)
  rownames(prom.activities) <- prom.key[rownames(prom.activities)]
  rownames(prom.activities)[rownames(prom.activities)==''] <- 'Unknown'
  saveRDS(prom.activities,file='data/prom_activities_norm.RDS')
  
  
  # add the gene activity matrix to the Seurat object as a new assay, and normalize it
  cortex.atac[['GeneBody']] <- CreateAssayObject(counts = gene.activities)
  cortex.atac[['Prom']] <- CreateAssayObject(counts = prom.activities)
  cortex.atac[['Bins5000']] <- CreateAssayObject(counts = bin_matrix)
 # cortex.atac[['peaks']] <- NULL

  
  saveRDS(cortex.atac,seurat_file)
}
