library(monocle3)
library(Seurat)
library(htmlwidgets)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(seriation)
library(scales)
library(doMC)

doMC::registerDoMC(cores=12) 

source('scripts/config.R')
source('scripts/aux_functions.R')

which_cluster <- c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3')
output.dir <- "results/Monocle3_int/"
dir.create(file.path(output.dir), showWarnings = FALSE)
useHarmony=TRUE

### Reading in Seurat object

print("Reading in co-embedded Seurat objects")

cortex <- readRDS('data/coembed_scRNA_scATAC.RDS')
levels(cortex) <- c('NSC','NSC_M','IPC','IPC_M','PN1','PN2','PN3','CR','IN','MG','Mural',"MG+Mural")

##### If using specifing cluster subset the seurat object #####
if (which_cluster!='ALL'){
  cortex <- subset(cortex, idents=levels(Idents(cortex))[(levels(Idents(cortex))%in%which_cluster)])
}

### Building the necessary parts for a basic cds

# part one, gene annotations

gene_annotation <- as.data.frame(rownames(cortex@reductions[["pca"]]@feature.loadings), row.names = rownames(cortex@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- as.data.frame(cortex@assays[["RNA"]]@counts@Dimnames[[2]], row.names = cortex@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix

New_matrix <- cortex@assays[["RNA"]]@data
New_matrix <- New_matrix[rownames(cortex@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

cds_from_cortex <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

cds_from_cortex <- preprocess_cds(cds_from_cortex, num_dim = 20)

### Construct and assign the made up partition

recreate.partition <- c(rep(1, length(cds_from_cortex@colData@rownames)))
names(recreate.partition) <- cds_from_cortex@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_cortex@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


### Assign the cluster info

list_cluster <- Idents(object = cortex)
names(list_cluster) <- cortex@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_cortex@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_from_cortex@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds_from_cortex@reducedDims@listData[["UMAP"]] <-cortex@reductions[["umap"]]@cell.embeddings

### Assign feature loading for downstream module analysis

cds_from_cortex@preprocess_aux$gene_loadings <- cortex@reductions[["pca"]]@feature.loadings

### Reading in cell type information

print("Reading in cell type information")

colData(cds_from_cortex)$celltype <- Idents(object = cortex)
### Learn graph, this step usually takes a significant period of time for larger samples

for (ncenter in seq(550,650,by=10)){
  print("Learning graph, which can take a while depends on the sample")
  
  cds_from_cortex <- learn_graph(cds_from_cortex, use_partition = T, learn_graph_control=list(minimal_branch_len=20,ncenter=600,geodesic_distance_ratio=0.275),close_loop=T,verbose=T)
  
  
  ### Plot cluster info with trajectory
  
  print("Plotting clusters")
  
  pdf(sprintf("%s/clusters.with.trajectory.%s.pdf", output.dir,''), width = 4, height = 4)
  clus <- plot_cells(cds_from_cortex, 
                     color_cells_by = 'cluster',
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE)
  print(clus)
  dev.off()
}
### Plot cell type info with trajectory

print("Plotting cell type info")

pdf(sprintf("%s/celltype.with.trajectory.pdf", output.dir), width = 4, height =4)
ctype <- plot_cells(cds_from_cortex, 
                    color_cells_by = 'celltype',
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)
print(ctype)
dev.off()


### You can also plot out different subclone information, which is not shown here

root_cell_list <- grep("NSC", colData(cds_from_cortex)$celltype)
root_cell_list <- counts(cds_from_cortex)@Dimnames[[2]][root_cell_list]
cds_from_cortex <- order_cells(cds_from_cortex, root_cells = root_cell_list[173])

### Now we plot pseudotime

print("Plotting pseudotime")

pdf(sprintf("%s/pseudotime.with.trajectory.%s.pdf", output.dir, Dim), width = 12, height = 8)
ptime <- plot_cells(cds_from_cortex,show_trajectory_graph = T, 
                    color_cells_by = 'pseudotime',
                    label_groups_by_cluster=T,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)
print(ptime)
dev.off()

pd <- pseudotime(cds_from_cortex)
pgraph <- principal_graph(cds_from_cortex)$UMAP
cortex <- readRDS('data/coembed_scRNA_scATAC.RDS')
cortex <- AddMetaData(cortex,metadata = pd,col.name = 'pseudotime')

#### Extract principal graph coordinates #####

ica_space_df <- t(cds_from_cortex@principal_graph_aux$UMAP$dp_mst) %>% as.data.frame() %>% dplyr::mutate(sample_name = rownames(.),sample_state = rownames(.))
edge_df <- pgraph %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="UMAP_1",
                       source_prin_graph_dim_2="UMAP_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="UMAP_1",
                       target_prin_graph_dim_2="UMAP_2"),
                   by = "target")

saveRDS(edge_df,paste0(output.dir,'ptime_edge_df.RDS'))
###############################################

saveRDS(cortex,'data/coembed_scRNA_scATAC.RDS')
saveRDS(cds_from_cortex,file = paste0(output.dir,'cds_pd.RDS'))
### Now we can subset the cells and do branch expression/pseudotime analysis
### This first part is to identify genes that are differentially expressed on the different paths through the trajectory
#cds_from_cortex_res <- graph_test(cds_from_cortex, neighbor_graph="principal_graph", cores=4)
if (learn_graph_f){
  cds_from_cortex_res <- read.table('results/Monocle3_clusterALL.2D/graph_test.tsv')
  res_ids <- row.names(subset(cds_from_cortex_res, q_value < 0.05))
  gene_module_df <- find_gene_modules(cds_from_cortex[res_ids,], resolution=c(10^seq(-6,-1)))
  write.table(gene_module_df,file='results/Monocle3_clusterALL.2D/gene_module_df.tsv',quote=F,sep='\t',col.names=T,row.names=T)
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds_from_cortex)), 
                                  cell_group=colData(cds_from_cortex)$celltype)
  agg_mat <- aggregate_gene_expression(cds_from_cortex, gene_module_df, cell_group_df)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  pdf('results/Monocle3_clusterALL.2D/geneModules.pdf',height=12,width=8)
  pheatmap::pheatmap(agg_mat,
                     scale="column", clustering_method="ward.D2")
  dev.off()
  sel_genes <- c('Hes1','Ascl1','Neurog1','Neurog2','Eomes','Dcx','Mapt')
  pdf('results/Monocle3_clusterALL.2D/selGenes_PT.pdf',height=10,width=12)
  plot_genes_in_pseudotime(cds_from_cortex[rowData(cds_from_cortex)$gene_short_name %in% sel_genes,],
                           color_cells_by="celltype",
                           min_expr=0.5)
  dev.off()
}
####### Testing Heatmap as function of pseudotime ###########
marker_genes <- c('Sox2','Hes1','Hes5','Fezf2','Sox5','Eomes','Neurog2','Tle4')
time_genes <- c('Bcl11b','Mapt','Pax6','Dcx','Fabp7')
if (which_cluster!='ALL'){
  cortex <- subset(cortex, idents=levels(Idents(cortex))[(levels(Idents(cortex))%in%which_cluster)])
}


###### Generate Models ###################
p2glinks <- readRDS('results/P2G-Links.RDS')
mat_f=c(p2glinks$posCor,p2glinks$negCor)

#### Generate models for P2G links ######
monocle_regression=F

df <- GetAssayData(cortex,slot = 'data',assay='MACS2peaks')
df <- df[row.names(df)%in%gsub('_','-',mat_f$peakName),colSums(df)!=0]
pd <- object$pseudotime
pd <- pd[names(pd)%in%colnames(df)]
df <- df[,order(pd)]
pd <- pd[order(pd)]
if (monocle_regression){
  cell_metadata <- as.data.frame(colnames(df), row.names = colnames(df))
  colnames(cell_metadata) <- "barcode"
  gene_annotation <- as.data.frame(rownames(df), row.names = rownames(df))
  colnames(gene_annotation) <- "gene_short_name"
  cds_subset <- new_cell_data_set(df,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
  mat <- features_in_pseudotime(cds_subset,pd)
} else {
  df <- log2(edgeR::cpm(df)+1)   ### Log normalize  ####
  mat <- adply(df, 1, function(x){
    return(rescale(predict(mgcv::gam(formula = x ~ s(as.numeric(pd), bs = "cs"),method = "REML")),to=c(0,1)))
  }, .parallel = TRUE,.id=NULL)
  row.names(mat) <- row.names(df)
}
mat <- mat[!is.na(rowSums(mat)),]
mat <- mat[,!is.na(colSums(mat))]
mat2 <- bin.matrix.rows(mat,bin.size = 25)
indx <- max.col(mat2)
mat <- mat[order(indx),]
mat2 <- mat2[order(indx),]
saveRDS(mat,file=paste0(output.dir,'pd_P2G_GAMmat.RDS'))
write.table(mat2[,1],paste0(output.dir,'pd_P2G_GAMmat2.tsv'),quote=F,col.names=F,row.names=T,sep='\t')

###########################################

#### Generate models for RNA ######
monocle_regression=F

df <- GetAssayData(cortex,slot = 'data',assay='RNA')
df <- df[row.names(df)%in%cortex@assays$RNA@var.features,]
pd <- cortex$pseudotime
pd <- pd[names(pd)%in%colnames(df)]
df <- df[,order(pd)]
pd <- pd[order(pd)]
df <- df[rowSums(df)>20,]   ### Only for normal RNA object
if (monocle_regression){
  cell_metadata <- as.data.frame(colnames(df), row.names = colnames(df))
  colnames(cell_metadata) <- "barcode"
  gene_annotation <- as.data.frame(rownames(df), row.names = rownames(df))
  colnames(gene_annotation) <- "gene_short_name"
  cds_subset <- new_cell_data_set(df,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
  mat <- features_in_pseudotime(cds_subset,pd)
} else {
  df <- log2(edgeR::cpm(df)+1)   ### Log normalize  ####
  mat <- adply(df, 1, function(x){
    return(rescale(predict(mgcv::gam(formula = x ~ s(as.numeric(pd), bs = "cs"),method = "REML")),to=c(0,1)))
  }, .parallel = TRUE,.id=NULL)
  row.names(mat) <- row.names(df)
}
mat <- mat[!is.na(rowSums(mat)),]
mat <- mat[,!is.na(colSums(mat))]
mat2 <- bin.matrix.rows(mat,bin.size = 25)
indx <- max.col(mat2)
mat <- mat[order(indx),]
mat2 <- mat2[order(indx),]
saveRDS(mat,file=paste0(output.dir,'pd_genes_GAMmat.RDS'))
write.table(mat2[,1],paste0(output.dir,'pd_genes_GAMmat2.tsv'),quote=F,col.names=F,row.names=T,sep='\t')


#### Generate models for Promoter Accessibility/Gene Body ######
monocle_regression=F

df <- GetAssayData(cortex,slot = 'data',assay='Prom') 
#df <- GetAssayData(cortex,slot = 'data',assay='GeneBody')
df <- df[row.names(df)%in%mat_f$gene_name,colSums(df)!=0]
pd <- cortex$pseudotime
pd <- pd[names(pd)%in%colnames(df)]
df <- df[,order(pd)]
pd <- pd[order(pd)]
df <- df[rowSums(df)>10,]   ### Only for normal RNA object
#df <- df[rowSums(df)>10,]
if (monocle_regression){
  cell_metadata <- as.data.frame(colnames(df), row.names = colnames(df))
  colnames(cell_metadata) <- "barcode"
  gene_annotation <- as.data.frame(rownames(df), row.names = rownames(df))
  colnames(gene_annotation) <- "gene_short_name"
  cds_subset <- new_cell_data_set(df,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
  mat <- features_in_pseudotime(cds_subset,pd)
} else {
  df <- log2(edgeR::cpm(df)+1)   ### Log normalize  ####
  mat <- adply(df, 1, function(x){
    return(rescale(predict(mgcv::gam(formula = x ~ s(as.numeric(pd), bs = "cs"),method = "REML")),to=c(0,1)))
  }, .parallel = TRUE,.id=NULL)
  row.names(mat) <- row.names(df)
}
mat <- mat[!is.na(rowSums(mat)),]
mat <- mat[,!is.na(colSums(mat))]
mat2 <- bin.matrix.rows(mat,bin.size = 25)
indx <- max.col(mat2)
mat <- mat[order(indx),]
mat2 <- mat2[order(indx),]
saveRDS(mat,file=paste0(output.dir,'pd_PromAccess_GAMmat.RDS'))
write.table(mat2[,1],paste0(output.dir,'pd_PromAccess_GAMmat2.tsv'),quote=F,col.names=F,row.names=T,sep='\t')

#### Generate models for Motifs ######

df <- GetAssayData(cortex,slot = 'data',assay='chromvar')
df <- df[,colSums(df)!=0]
df <- adply(df,1,function(x){
  return(rescale(x,to=c(0,10)))
})
row.names(df) <- df$X1
df <- df[,-1]
cell_metadata <- as.data.frame(colnames(df), row.names = colnames(df))
colnames(cell_metadata) <- "barcode"
gene_annotation <- as.data.frame(rownames(df), row.names = rownames(df))
colnames(gene_annotation) <- "gene_short_name"
cds_subset <- new_cell_data_set(as.matrix(df),cell_metadata = cell_metadata,gene_metadata = gene_annotation)
pd <- cortex$pseudotime
pd <- pd[names(pd)%in%colnames(df)]

mat <- features_in_pseudotime(cds_subset,pd)
mat <- mat[!is.na(rowSums(mat)),]
mat <- mat[,!is.na(colSums(mat))]
mat2 <- bin.matrix.rows(mat,bin.size = 25)
indx <- max.col(mat2)
mat <- mat[order(indx),]
mat2 <- mat2[order(indx),]
saveRDS(mat,file=paste0(output.dir,'pd_motifs_mat.RDS'))
write.table(mat2[,1],paste0(output.dir,'pd_motifs_mat2.tsv'),quote=F,col.names=F,row.names=T,sep='\t')

#### Generate models for motif Genes ######
# motif.names <- unlist(GetMotifObject(atac.object)@motif.names)
# motif.names <-  capitalize(tolower(motif.names))
# motif.matrix <- GetMotifData(atac.object)
# colnames(motif.matrix) <- as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])
# motifFactors <- sapply(strsplit(colnames(motif.matrix), "\\(|\\:|\\."), function(x) x[[1]])
# 
# df <- GetAssayData(cortex,slot = 'data',assay='RNA')
# df <- df[row.names(df)%in%motifFactors,]
# pd <- cortex$pseudotime
# pd <- pd[names(pd)%in%colnames(df)]
# df <- df[,order(pd)]
# pd <- pd[order(pd)]
# df <- df[rowSums(df)>1,] 
# if (monocle_regression){
#   cell_metadata <- as.data.frame(colnames(df), row.names = colnames(df))
#   colnames(cell_metadata) <- "barcode"
#   gene_annotation <- as.data.frame(rownames(df), row.names = rownames(df))
#   colnames(gene_annotation) <- "gene_short_name"
#   cds_subset <- new_cell_data_set(df,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
#   mat <- features_in_pseudotime(cds_subset,pd)
# } else {
#   df <- log2(edgeR::cpm(df)+1)   ### Log normalize  ####
#   mat <- adply(df, 1, function(x){
#     return(scales::rescale(predict(mgcv::gam(formula = x ~ s(as.numeric(pd), bs = "cs"),method = "REML")),to=c(0,1)))
#   }, .parallel = TRUE,.id=NULL)
#   row.names(mat) <- row.names(df)
# }

