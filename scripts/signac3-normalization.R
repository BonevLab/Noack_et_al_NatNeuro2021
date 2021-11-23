library(Signac)
library(Seurat)
library(reticulate)
library(ggplot2)
library(LSD)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table)
library(harmony)
#library(future)
#plan("multiprocess", workers = 4)

fragment.path <- "data/fragments_filtered.tsv.gz"
doHarmony=TRUE

source('scripts/config.R')
source('scripts/plot_functions.R')
source('scripts/aux_functions.R')
#-----------------------------------

#-----------------------------------------
# NORMALIZATION AND DIMENSIONAL REDUCTION
#-----------------------------------------
cortex.atac <- readRDS('data/merged_scATAC_filtered.RDS')
cortex.atac <- BinarizeCounts(cortex.atac,assay='Bins5000')

#### Select the top 20000 features from the Bins object
nFeatures <- 20000

mat <- cortex.atac@assays$Bins5000@counts
mat <- mat[head(order(Matrix::rowSums(mat),decreasing = TRUE),nFeatures),]
cortex.atac[['TopBins']] <- CreateAssayObject(counts=mat)
rm(mat)

#### Reduce size of the object by setting Bins assays to NULL
#cortex.atac@assays$Bins5000 <- NULL

cortex.atac <- NormalizeData(
  object = cortex.atac,
  assay = 'GeneBody',
  normalization.method = 'LogNormalize',
  scale.factor = median(cortex.atac$nCount_GeneBody)
)

cortex.atac <- NormalizeData(
  object = cortex.atac,
  assay = 'Prom',
  normalization.method = 'LogNormalize',
  scale.factor = median(cortex.atac$nCount_Prom)
)

assay_types <- c('TopBins')
TFIDF_methods <- c(3)
cutoffs <- c('q0')
max_dims <- c(20,25,30)
resolutions <- c(0.6,0.8,1)

#if (do)
for (assay_type in assay_types){
  
  DefaultAssay(cortex.atac) <- assay_type
  
  for (TFIDF_method in TFIDF_methods){
    
    cortex.atac <- RunTFIDF(cortex.atac,method = TFIDF_method,scale.factor = 10000,verbose = T)
    
    for (cutoff in cutoffs){
      
      cortex.atac <- FindTopFeatures(cortex.atac,min.cutoff = cutoff)
      cortex.atac <- RunSVD(
        object = cortex.atac,
        assay = assay_type,
        reduction.key = 'LSI_',
        reduction.name = 'lsi'
      )
      cortex.atac <- RunHarmony(
        object = cortex.atac,
        group.by.vars = 'orig.ident',
        reduction = 'lsi',
        assay.use = assay_type,
        project.dim = FALSE
      )
      for (max_dim in max_dims){
        cortex.atac <- RunUMAP(object = cortex.atac, reduction = ifelse(doHarmony,'harmony','lsi'), dims = 1:max_dim,min.dist = 0.5,spread = 1.5,n.epochs = 2000)
        cortex.atac <- FindNeighbors(object = cortex.atac, reduction = ifelse(doHarmony,'harmony','lsi'), dims = 1:max_dim,force.recalc = T)
        
        for (resolution in resolutions){
          cortex.atac <- FindClusters(object = cortex.atac, verbose = FALSE,resolution = resolution,n.start = 50,n.iter = 50)
          p1 <- DimPlot(object = cortex.atac, label = TRUE,group.by = 'orig.ident') + ggtitle(paste0('N Var Features: ',length(cortex.atac@assays[[assay_type]]@var.features)))
          p2 <- DimPlot(object = cortex.atac, label = TRUE) + NoLegend()
          
          name <- paste0('TFDF_',TFIDF_method,'_assay_',assay_type,'_cutoff_',cutoff,'_dims_',max_dim,'_res',resolution,'_features',length(cortex.atac@assays[[assay_type]]@var.features))
          pdf(paste0('plots/scATAC/',name,'.pdf'),height=8,width=16)
          print(CombinePlots(plots = list(p1, p2))) 
          dev.off()
          DefaultAssay(cortex.atac) <- 'Prom'
          for (set in c('diff_markers','layer_markers','temp_markers','other_markers','regional_markers')){
            pdf(paste0('plots/scATAC/',name,'_',set,'.pdf'),height=10,width=10)
            print(FeaturePlot(cortex.atac, features = get(set),cols = gene_colors(3),pt.size = 0.1,min.cutoff = 'q5',max.cutoff = 'q95',label = T))
            dev.off()
          }
        }
        DefaultAssay(cortex.atac) <- assay_type
      }
    }
  }  
}

DefaultAssay(cortex.atac) <- assay_type

#-------------------------------------------------------------------------------------------------
# Get Cluster Beds
#-------------------------------------------------------------------------------------------------
dirClusters <- "results/LSI-Cluster-Beds_Insertions/"
dir.create(dirClusters)
dirPeaks <- "results/LSI-Cluster-Peaks_Insertions/"
dir.create(dirPeaks)
method <- "q"
cutoff <- 0.05
shift <- -75
extsize <- 150
genome_size <- 'mm'

all_cells <- cortex.atac$seurat_clusters
if (!file.exists('data/all_fragments.RDS')){
  all_reads <- fread(
    file = fragment.path,header = F,
    col.names = c('chr', 'start', 'end', 'cell', 'count'),showProgress = T
  )
  all_reads <- all_reads[with(data = all_reads, expr = order(chr, start)), ]
  saveRDS(all_reads,file='data/all_fragments.RDS')
} else {
  all_reads <- readRDS('data/all_fragments.RDS')
}
if (!file.exists('data/all_insertions.RDS')){
  insertions <- all_reads
  insertions2 <- all_reads
  insertions2$start <- insertions2$end-1
  insertions$end <- insertions$start+1 
  insertions <- rbind(insertions,insertions2)
  insertions <- insertions[with(data = insertions, expr = order(chr, start)), ]
  saveRDS(insertions,file='data/all_insertions.RDS')
} else {
  insertions <- readRDS('data/all_insertions.RDS')
}


cat('source activate macs2;\n',file='scripts/macs.sh',append=F)
for (i in 0:length(levels(all_cells))){
  cluster_cells <- names(all_cells[all_cells==i]) 
  ClusterFragments(reads=insertions,
                  cells=cluster_cells,
                  output.path=paste0(dirClusters,'Cluster',i,'.bed'),
                  assume.sorted = TRUE, verbose = TRUE
                  )
  clusterBedj <- paste0('Cluster',i,".bed")
  cmdPeaks <- sprintf(
    "macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all -B --SPMR", 
    genome_size, 
    paste0('Cluster',i), 
    paste0(dirClusters,clusterBedj), 
    dirPeaks
  )
  if (!is.null(shift) & !is.null(extsize)) {
    cmdPeaks <- sprintf("%s --shift %s --extsize %s", cmdPeaks, shift, extsize)
  }
  if (tolower(method) == "p") {
    cmdPeaks <- sprintf("%s -p %s", cmdPeaks, cutoff)
  }else {
    cmdPeaks <- sprintf("%s -q %s", cmdPeaks, cutoff)
  }
#  message("Running Macs2...")
 # message(cmdPeaks)
  cmdBW <- sprintf(
    "\nsort -k1,1 -k2,2n %s > %s; bedGraphToBigWig %s %s %s;", 
    paste0(dirPeaks,'/Cluster',i,'_treat_pileup.bdg'), 
    paste0(dirPeaks,'/Cluster',i,'.bdg'),
    paste0(dirPeaks,'/Cluster',i,'.bdg'),
    '/home/hpc/bonev/annotations/mm10/mm10.chrom.sizes',
    paste0(dirPeaks,'/Cluster',i,'.bw \n')
  ) 
  cat(cmdPeaks,file='scripts/macs.sh',append=T)
  cat(cmdBW,file='scripts/macs.sh',append=T)
  }
  
  
#-------------------------------------------------------------------------------------------------
# Make Non-Overlapping Peak Set
#-------------------------------------------------------------------------------------------------
df <- data.frame(
  samples = gsub("\\_summits.bed","",list.files(dirPeaks, pattern = "\\_summits.bed", full.names = FALSE)),
  groups = "scATAC",
  summits = list.files(dirPeaks, pattern = "\\_summits.bed", full.names = TRUE)
)

unionPeaks <- extendedPeakSet(
  df = df,
  BSgenome = BSgenome.Mmusculus.UCSC.mm10, 
  extend = 250,
  blacklist = "data/mm10.blacklist.bed",
  nSummits = 200000
)
unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% paste0("chr",c(1:19))]
unionPeaks <- keepSeqlevels(unionPeaks, paste0("chr",c(1:19)))
  
saveRDS(unionPeaks,file='data/unionPeaks.RDS')
rtracklayer::export.bed(unionPeaks, "data/unionPeaks.bed")



