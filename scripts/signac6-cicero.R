library(Signac)
library(Seurat)
library(monocle3)
library(cicero)
library(reticulate)
library(ggplot2)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(org.Mm.eg.db)
library(data.table)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(seriation)
library(scales)


source('scripts/config.R')
source('scripts/aux_functions.R')
fragment.path <- "data/fragments.filtered.gz"
mm10_txdb <- loadDb("/home/hpc/bonev/annotations/mm10/mm10_txdb.sqlite")


##### Initializa the seurat object ########

cortex.atac <- readRDS('data/merged_scATAC_integrated.RDS')
DefaultAssay(cortex.atac) <- 'MACS2peaks'


mdata <- cortex.atac@meta.data
tssWindow <- 2500
flank <- 250*10^3
corCutOff <- 0.35
bsgenome <- BSgenome.Mmusculus.UCSC.mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
orgdb <- org.Mm.eg.db

#Reduced Dimensions
dimred <- cortex.atac@reductions$harmony@cell.embeddings
se <- cortex.atac@assays$MACS2peaks@counts
row.names(se) <- gsub('-','_',row.names(se))
#Get ChromSizes
chromSizes <- seqlengths(bsgenome)[paste0("chr",c(1:19,"X"))]
genome <- data.frame(names(chromSizes),chromSizes)
rownames(genome) <- NULL

obj <- makeCDS(cortex.atac, binarize = TRUE)
obj <- monocle3::detect_genes(obj)
obj <- obj[Matrix::rowSums(exprs(obj)) != 0,] 
obj <- estimate_size_factors(obj)
ciceroObj <- make_cicero_cds(obj, k = 50, reduced_coordinates = dimred)
saveRDS(ciceroObj,file='results/ciceroObj.RDS')


### Run Cicero default #############

mm10_genome <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
mm10_genome <- as.data.frame(mm10_genome[1:21])
mm10_genome$V1 <- row.names(mm10_genome)
row.names(mm10_genome) <- seq(1:nrow(mm10_genome))
mm10_genome <- mm10_genome[,c(2,1)]
colnames(mm10_genome) <- c('V1','V2')
conns <- run_cicero(obj, mm10_genome)
saveRDS(conns,file='results/cicero_conns.RDS')


#Compute Correlations
message("Computing grouped correlations...")
gr <- featureToGR(row.names(ciceroObj))
o <- suppressWarnings(as.matrix( findOverlaps(resize( resize(gr,1,"center"), 2*flank + 1, "center"), resize(gr,1,"center"), ignore.strand=TRUE) ))
o <- data.table::as.data.table(data.frame(i = matrixStats::rowMins(o), j = matrixStats::rowMaxs(o)))
o <- data.frame(o[!duplicated(o),])
o <- o[o[,1]!=o[,2],]
o$cor <- rowCorCpp(o[,1], o[,2], as.matrix(ciceroObj@assays$data$counts), as.matrix(ciceroObj@assays$data$counts))
connections <- data.frame(
  Peak1 = row.names(ciceroObj)[o[,1]], 
  Peak2 = row.names(ciceroObj)[o[,2]], 
  coaccess = o[,3]
)


#Annotate CDS
message("Annotating Cell Data Set...")
genes <- getTxDbGenes(txdb=mm10_txdb,orgdb=orgdb)
names(genes) <- genes$symbol
genes <- resize(genes, 1, "start") %>% resize(tssWindow * 2 + 1, "center")
geneDF <- data.frame(chromosome=seqnames(genes),start=start(genes),end=end(genes), gene=genes$symbol)
obj <- annotate_cds_by_site(obj, geneDF)

#Prepare for Co-Accessibility
nSites <- Matrix::colSums(obj@assays$data$counts)
names(nSites) <- row.names(pData(obj))

#Cicero with Correlations
message("Calculating normalized gene activities...")
ciceroGA <- normalize_gene_activities(build_gene_activity_matrix(obj, connections, coaccess_cutoff = corCutOff), nSites)
ciceroGA_orig <- normalize_gene_activities(build_gene_activity_matrix(obj, conns, coaccess_cutoff = corCutOff), nSites)

seCicero <- SummarizedExperiment(
  assays = SimpleList(gA = ciceroGA),
  rowRanges = genes[rownames(ciceroGA),],
  colData = mdata
)
seCicero_orig <- SummarizedExperiment(
  assays = SimpleList(gA = ciceroGA_orig),
  rowRanges = genes[rownames(ciceroGA_orig),],
  colData = mdata
)

seCiceroLog <- SummarizedExperiment(
  assays = SimpleList(logGA = log2(10^6 * ciceroGA + 1)),
  rowRanges = genes[rownames(ciceroGA),],
  colData = mdata
)

seCiceroLog_orig <- SummarizedExperiment(
  assays = SimpleList(logGA = log2(10^6 * ciceroGA_orig + 1)),
  rowRanges = genes[rownames(ciceroGA_orig),],
  colData = mdata
)
cortex.atac[['cicero_GA']] <- CreateAssayObject(counts = assay(seCicero))
cortex.atac[['ciceroOrig_logGA']] <- CreateAssayObject(counts = assay(seCiceroLog_orig))

plotGene <- 'Hes5'
DefaultAssay(cortex.atac) <- 'Prom'
p1 <- FeaturePlot(
  object = cortex.atac,
  features = plotGene,
  min.cutoff = 'q5',
  pt.size = 0.01,
  max.cutoff = 'q95',
  cols = gene_colors(3),
  ncol = 3
)

DefaultAssay(cortex.atac) <- 'cicero_GA'
p2 <- FeaturePlot(
  object = cortex.atac,
  features = plotGene,
  min.cutoff = 'q5',
  pt.size = 0.01,
  max.cutoff = 'q95',
  cols = gene_colors(3),
  ncol = 3
)


CombinePlots(list(p1,p2)) + theme_cowplot()
saveRDS(cortex.atac,file = 'data/merged_scATAC_integrated_cicero.RDS')



####### Construct cicero object ##############
Dim = "2D"
input.dir <- getwd()
nPC <- 20
cluster.res <- 1
which_cluster <- c("E12_AP_1","E12_AP_2","early_AP","intermediate_AP","late_AP")
which_assay <- 'chromvar'


if (!is.null(which_cluster)){
  cortex.atac <-subset(cortex.atac, idents=levels(Idents(cortex.atac))[(levels(Idents(cortex.atac))%in%which_cluster)])
}

### Re-dimension reduction for 3D rendering #########

if (Dim = "3D") {
  
  print ("Running UMAP 3D")
  
  cortex.atac <- RunUMAP(object = cortex.atac, reduction = "pca", dims = 1:nPC, n.components = 3)
  
  print("Clustering 3D")
  
  cortex.atac <- FindNeighbors(object=cortex.atac, dims=1:nPC)
  cortex.atac <- FindClusters(object=cortex.atac, resolution=cluster.res)
  
}

DefaultAssay(cortex.atac) <- which_assay

### Building the necessary parts for a basic cds ###########

# part one, gene annotations

#gene_annotation <- as.data.frame(rownames(cortex.atac@reductions[["lsi"]]@feature.loadings), row.names = rownames(cortex.atac@reductions[["lsi"]]@feature.loadings))
gene_annotation <- as.data.frame(rownames(cortex.atac), row.names = rownames(cortex.atac))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

#cell_metadata <- as.data.frame(cortex.atac@assays[['MACS2peaks']]@counts@Dimnames[[2]], row.names = cortex.atac@assays[['MACS2peaks']]@counts@Dimnames[[2]])
cell_metadata <- as.data.frame(colnames(cortex.atac), row.names = colnames(cortex.atac))
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix

New_matrix <- cortex.atac@assays[[which_assay]]@data
#New_matrix <- New_matrix[rownames(cortex.atac@reductions[["lsi"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object ######

cds_from_cortex.atac <- new_cell_data_set(expression_matrix,
                                          cell_metadata = cell_metadata,
                                          gene_metadata = gene_annotation)


### Construct and assign the made up partition ##########

recreate.partition <- c(rep(1, length(cds_from_cortex.atac@colData@rownames)))
names(recreate.partition) <- cds_from_cortex.atac@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_cortex.atac@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


### Assign the cluster info

list_cluster <- Idents(object = cortex.atac)
names(list_cluster) <- colnames(cortex.atac)

cds_from_cortex.atac@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster


### Could be a space-holder, but essentially fills out louvain parameters

cds_from_cortex.atac@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds_from_cortex.atac@reducedDims@listData[["UMAP"]] <-cortex.atac@reductions[["umap"]]@cell.embeddings


### Assign feature loading for downstream module analysis

cds_from_cortex.atac@preprocess_aux$gene_loadings <- cortex.atac@reductions[["harmony"]]@feature.loadings


### Learn graph, this step usually takes a significant period of time for larger samples

lineage.table <- Idents(cortex.atac)
indices <- match(pData(cds_from_cortex.atac)$barcode, names(lineage.table))
lineage.table.refined <- lineage.table[indices]
colData(cds_from_cortex.atac)$celltype <- lineage.table.refined

plot_cells(cds_from_cortex.atac ,
           color_cells_by = "celltype",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,label_roots = T,
           label_branch_points=TRUE)


print("Learning graph, which can take a while depends on the sample")
cds_from_cortex.atac <- learn_graph(cds_from_cortex.atac,use_partition = T,learn_graph_control=list(minimal_branch_len=27,geodesic_distance_ratio=1/3))

colData(cds_from_cortex.atac)$celltype <- Idents(object = cortex.atac)
root_cell_list <- grep("E12_AP_2", colData(cds_from_cortex.atac)$celltype)
root_cell_list <- counts(cds_from_cortex.atac)@Dimnames[[2]][root_cell_list]
cds_from_cortex.atac <- order_cells(cds_from_cortex.atac, root_cells = root_cell_list[1])

pdf('plots/scATAC/pseudotime_cicero.pdf',height=8,width=8)
plot_cells(cds_from_cortex.atac ,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=T,
           label_leaves=FALSE,label_roots = F,show_trajectory_graph = F,
           label_branch_points=F)
dev.off()

#Annotate CDS
if(which_assay!='MACS2Peaks'){
  tssWindow <- 10
  message("Annotating Cell Data Set...")
  genes <- getTxDbGenes(txdb=mm10_txdb,orgdb=orgdb)
  names(genes) <- genes$symbol
  genes <- resize(genes, 1, "start") %>% resize(tssWindow * 2 + 1, "center")
  geneDF <- data.frame(chromosome=seqnames(genes),start=start(genes),end=end(genes), gene=genes$symbol)
  fData(cds_from_cortex.atac) <- NULL
  cds_from_cortex.atac <- annotate_cds_by_site(cds_from_cortex.atac, geneDF,maxgap = 0)
}

marker_genes <- c('Fabp7','Pax6','Sox2','Hes5','Eomes','Btg2','Dcx','Mapt','Tubb3','Camk2a','Neurod6','Sox5','Bcl11b','Fezf2','Rorb','Cux1','Cux2','Gad1','Gad2','Sst','Gfap','Olig1','Olig2','Tle4')
time_genes <- c('Hmga2','Cabp1','Hes1','Ccnd2','Suz12','Eed','Fn1','Ddah1','Atg12','Mfge8','Slc1a3')
tf_genes <- c('Id4','Id1','Tcf4','Hey1','Hey2','Zeb1','Meis2','Smarca2','Tcf12','Ldb1','Ldb2','Lhx2','Lmo3','Lmo4','Isl1','Nr2f1')
mitotic_genes <- c('Cdk1','Ube2c','Top2a','Hist1h4e','Hist1h4c','Pcna')

var_features <- row.names(cortex.atac)
cds_subset <- cds_from_cortex.atac[rowData(cds_from_cortex.atac)$gene_short_name %in% var_features,]
mat <-plot_genes_in_pseudotime2(cds_subset,label_by_short_name=FALSE)
mat <- mat[!is.na(rowSums(mat)),]
mat2 <- bin.matrix.rows(mat,bin.size = 10)
indx <- max.col(mat2)
mat <- mat[order(indx),]
mat2 <- mat2[order(indx),]


#ser_dist <- seriate(dist(mat),method='OLO')
ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%(unique(c(marker_genes,time_genes)))), labels = row.names(mat)[which(row.names(mat)%in%unique(c(marker_genes,time_genes)))],labels_gp = gpar(fontsize = 16)))
col.list2 <- inferno(length(unique(cortex.atac$orig.ident)))
names(col.list2) <- unique(cortex.atac$orig.ident)
idents_indx2 <- cortex.atac$orig.ident
idents_indx2 <- idents_indx2[match(colnames(mat),names(idents_indx2))]
col_names2 <- idents_indx2
dup_indices2 <- !(duplicated(col_names2))
mark_pos2 <- round(seq(1:length(col_names2))[dup_indices2] + table(col_names2)/2)
col.list <- hue_pal()(length(levels(Idents(cortex.atac))))
names(col.list) <- levels(Idents(cortex.atac))
idents_indx <- Idents(cortex.atac)
idents_indx <- idents_indx[match(colnames(mat),names(idents_indx))]
col_names <- idents_indx
dup_indices <- !(duplicated(col_names))
mark_pos <- round(seq(1:length(col_names))[dup_indices] + table(col_names)/2)

mat[mat<(-5)] <- -4
mat[mat>5] <- 4
ha = HeatmapAnnotation(foo = anno_mark(at = mark_pos, labels = col_names[dup_indices],labels_gp = gpar(fontsize = 16),labels_rot = 45),cluster = idents_indx ,col = list(cluster=col.list),show_legend = F)
ha2 = HeatmapAnnotation(cluster = idents_indx2 ,col = list(cluster=col.list2),show_legend = F)
pdf('plots/scATAC/APpseudotime_hm.pdf',height=18,width=12)
#DimPlot(cortex.atac,label = T)
Heatmap(as.matrix(mat), name = "Pseudotime Expression",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, right_annotation = ra,top_annotation = ha2,col = viridis(10))
dev.off()













##### Compare cicero accessibility connections with Hi-C
require(misha)
require(tidyr)
require(dplyr)
gsetroot('/home/hpc/bonev/trackdb/mm10')
tss <- gintervals.load("glp.intervals.ucscCanTSS")
peaks1 <- separate(as.data.frame(connections),col = 'Peak1', into = c("chrom1", "start","end1"),sep = "_")
peaks2 <- separate(as.data.frame(connections),col = 'Peak2', into = c("chrom2", "start2","end2"),sep = "_")
head(peaks1)
peaks <- as.data.frame(cbind(peaks1[,1:3],peaks2[,2:4],peaks1[,5],as.character(connections$Peak1),as.character(connections$Peak2)))
colnames(peaks) <- c('chrom1','start1','end1','chrom2','start2','end2','coaccess','Peak1_ID','Peak2_ID')
peaks1 <- gintervals(as.character(peaks[,1]),as.numeric(peaks[,2]),as.numeric(peaks[,3]))
peaks1 <- gintervals.canonic(unique(peaks1))
peaks1_overlap <- gintervals.neighbors(peaks1,tss)
peaks1_overlap <- peaks1_overlap[abs(peaks1_overlap$dist)<10,]
peak_IDs <- as.vector(unite(peaks1_overlap[,1:3],col='peaks_ID',sep='_'))
peaks1_overlap$peak_ID <- peak_IDs$peaks_ID
tss_connections <- peaks[peaks$Peak1_ID%in%peaks1_overlap$peak_ID|peaks$Peak2_ID%in%peaks1_overlap$peak_ID,]



nPeaks <- 1000
interval_window <- 2000
min_dist <- 5e3
max_dist <- 2e6
expand=c(-1500,1500)
domains <- gintervals.load("hic.ncx_Hes5.ins_250_domains_expanded")
tracks <- c('hic.ES.score_k200','hic.ncx_Hes5.score_k200','hic.ncx_Dcx.score_k200')

options(gmax.data.size=1e7)
options(gmultitasking=FALSE)
for (i in 1:length(tracks)){
  gvtrack.create(paste0('v_',tracks[i]),tracks[i],'max')
  gvtrack.iterator.2d(paste0('v_',tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
}

sample_connections <- tss_connections[sample(row.names(tss_connections),nPeaks),]
sample_connections$peakID <- paste0(sample_connections$Peak1_ID,'_',sample_connections$Peak2_ID)
grid <- gintervals.2d(as.character(sample_connections[,1]),as.numeric(sample_connections[,2]),as.numeric(sample_connections[,3]),as.character(sample_connections[,4]),as.numeric(sample_connections[,5]),as.numeric(sample_connections[,6]))
dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
intra_grid <- gintervals.intersect(grid,dom)
inter_2d = construct.grid2(domains,domains,min_dist,max_dist)			
inter_grid <- gintervals.intersect(grid,inter_2d)

intra_scores <- gextract(paste0('v_',tracks),intervals = intra_grid,iterator = intra_grid,band = -c(4e6,5e3))
inter_scores <- gextract(paste0('v_',tracks),intervals = inter_grid,iterator = inter_grid,band = -c(4e6,5e3))

intra_scores$peakID_1 <- paste0(intra_scores[,1],'_',intra_scores[,2],'_',intra_scores[,3])
intra_scores$peakID_2 <- paste0(intra_scores[,4],'_',intra_scores[,5],'_',intra_scores[,6])
intra_scores$peakID <- paste0(intra_scores$peakID_1,'_',intra_scores$peakID_2)

merged_scores <- merge(intra_scores,sample_connections,by='peakID')

