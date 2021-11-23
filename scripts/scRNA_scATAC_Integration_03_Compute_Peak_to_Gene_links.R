#Clustering and scATAC-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
library(plyr)
library(dplyr)
require(ComplexHeatmap)
library(circlize)
library(LSD)
set.seed(1)
####################################################
#Functions
####################################################
source('scripts/config.R')
source('scripts/aux_functions.R')

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


getGeneGTF <- function(file){
  #Import
  message("Reading in GTF...")
  importGTF <- rtracklayer::import(file)
  #Exon Info
  message("Computing Effective Exon Lengths...")
  exonGTF <- importGTF[importGTF$type=="exon",]
  exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
  exonReduced <- unlist(exonList, use.names=TRUE)
  mcols(exonReduced)$gene_id <- names(exonReduced)
  mcols(exonReduced)$widths <- width(exonReduced)
  exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
  exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>% 
    unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
  #Gene Info
  message("Constructing gene GTF...")
  geneGTF1 <- importGTF[importGTF$type=="gene",]
  geneGTF2 <- GRanges(
      seqnames=paste0("chr",seqnames(geneGTF1)),
      ranges=ranges(geneGTF1),
      strand=strand(geneGTF1),
      gene_name=geneGTF1$gene_name,
      gene_id=geneGTF1$gene_id
    ) %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
  mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
  return(geneGTF2)
}

Rcpp::sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX)-1 > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY)-1 > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);
    } 

    return(cor);

  }'
)

nSample <- function(x, n, type="v"){
  if(type=="v"){
    if(length(x) > n){
      s <- x[sample(seq_along(x),n)]
    }else{
      s <- x
    }
  }else if(type=="c"){
    if(ncol(x) > n){
      s <- x[,sample(seq_len(ncol(x)),n),drop=F]
    }else{
      s <- x
    }
  }else if(type=="r"){
    if(nrow(x) > n){
      s <- x[sample(seq_len(nrow(x)),n),,drop=F]
    }else{
      s <- x
    }
  }else{
    stop(paste0("type ",type," unrecognized..."))
  }
  return(s)
}

getNullCorrelations <- function(seA, seB, o, n){

  set.seed(1)
  o$seq <- seqnames(seA)[o$A]

  nullCor <- lapply(seq_along(unique(o$seq)), function(i){

    #Get chr from olist
    chri <- unique(o$seq)[i]
    cat(paste0(chri), "\n")

    #Randomly get n seA
    transAidx <- nSample( which(as.character(seqnames(seA)) != chri), n, "v")

    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$B))

    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])

    seSubA <- seA[idxA]
    seSubB <- seB[idxB]

    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)

    colnames(grid) <- c("A", "B")
    out <- rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    out <- na.omit(out)

    return(out)

  }) %>% SimpleList

  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)

  return(list(summaryDF, unlist(nullCor)))

}

getQuantiles <- function(v, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}


####################################################
#Input Data
####################################################

#Input Files
scATAC_file <- "results/scATAC-Merged-KNN-SVD.RDS"
scRNA_file <- "results/scRNA-Merged-KNN-SVD.RDS"
gtf_file <- '/home/hpc/bonev/annotations/mm10/cellranger_rna/genes/genes.gtf'

#Params
fixA <- "center"
fixB <- "start"
associationWindow <- 2 * 500*10^3 + 1 #+-500 Kb
corCutOff <- 0.35 #Pearson Correlation Cutoff
fdrCutOff <- 0.1 #FDR Cutoff
distCutOff <- 5000 #Min Dist to TSS

#Input Summarized Experiments Log2 Normalize
seA <- readRDS(scATAC_file) #Aggregated scATAC Summarized Experiment
seB <- readRDS(scRNA_file) #Aggregated scRNA Summarized Experiment
assay(seA) <- log2(edgeR::cpm(assay(seA))/100+1)
assay(seB) <- log2(edgeR::cpm(assay(seB))/100+1)

#Resize B to association Window
seBWindow <- resize(rowRanges(seB), width = 1, fix = fixB) %>%
  {suppressWarnings(resize(., width = associationWindow, fix = "center"))} %>% trim(.)

#Keep only seA within association window
seA <- seA[unique(queryHits(findOverlaps(resize(seA,1,fixA), seBWindow, ignore.strand = TRUE)))]

#Getting distances
message("Getting Distances...")
o <- findOverlaps(seBWindow, resize(rowRanges(seA),1,fixA), ignore.strand = TRUE)

#Get Distance from Fixed point A B correct for minus stranded
mcols(o)$distance <- start(resize(rowRanges(seA),1,fixA))[subjectHits(o)] - start(resize(rowRanges(seB),1,fixB))[queryHits(o)]
mcols(o)$distance[which(as.character(strand(rowRanges(seB)))[queryHits(o)]=="-")] <- -1*mcols(o)$distance[which(as.character(strand(rowRanges(seB)))[queryHits(o)]=="-")]

#Add other info
o <- DataFrame(o)
colnames(o) <- c("B","A","distance")
o <- o[,c("A","B","distance")]

nullCor <- getNullCorrelations(seA, seB, o, 1000)
o$Correlation <- rowCorCpp(as.integer(o[,1]), as.integer(o[,2]), assay(seA), assay(seB))
o$VarAssayA <- getQuantiles(matrixStats::rowVars(assay(seA)))[o$A]
o$VarAssayB <- getQuantiles(matrixStats::rowVars(assay(seB)))[o$B]
o$Pval <- 2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
o$FDR <- p.adjust(o$Pval, method = "fdr")
o <- o[!is.na(o$Correlation),]
#Get GTF
gtf <- getGeneGTF(gtf_file)
tssRNA <- resize(gtf, 1, "start")
### Correct TSS for ucsc known gene coordinates
tss <- gintervals.load("glp.intervals.ucscCanTSS")
tss$strand='*'
tss <- tss[match(tssRNA$gene_name,as.character(tss$geneName)),]
ranges(tssRNA)[!is.na(tss$start)] <- IRanges(start=tss$start[!is.na(tss$start)],end=tss$start[!is.na(tss$start)])  
######################################
strand(tssRNA) <- "*"
peakLinks <- rowRanges(seA)[o[,1]]
geneLinks <- rowRanges(seB) %>% resize(1, "start") %>% {.[o[,2]]}
mcolsLinks <- data.frame(geneLinks)[,c("seqnames","start","strand","gene_name","gene_id","exonLength")]
colnames(mcolsLinks) <- c("gene_chr","gene_start","gene_strand","gene_name","gene_id","exonLength")
mcolsLinks <- cbind(mcolsLinks, data.frame(o))
mcolsLinks$nearestGene <- tssRNA$gene_name[subjectHits(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))]
mcolsLinks$nearestGeneDistance <- mcols(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))$distance
mcols(peakLinks) <- mcolsLinks
peakLinks$peakName <- paste(seqnames(peakLinks), start(peakLinks), end(peakLinks), sep = "_")
peakLinks$sigCorrelation <- peakLinks$Correlation >= corCutOff & peakLinks$FDR <= fdrCutOff & abs(peakLinks$distance) >= distCutOff
linksSig <- peakLinks[which(peakLinks$sigCorrelation)]

outMatch <- list(
  seA = seA[unique(mcols(linksSig)$peakName),], 
  seB = seB[unique(mcols(linksSig)$gene_name),], 
  posCor = linksSig,
  linksAll = peakLinks
  )

saveRDS(outMatch, "results/P2G-Links.RDS")

##### Compare connections with Hi-C
require(misha)
require(tidyr)
require(plyr)
require(dplyr)
gsetroot('/home/hpc/bonev/trackdb/mm10')
tss <- as.data.frame(resize(tssRNA, 2, "start"))

nPeaks <- 1000
interval_window <- 2500
min_dist <- 5e3
max_dist <- 5e5
expand=c(-5000,5000)
domains <- gintervals.load("hic.E14_NSC.ins_250_domains_expanded")
tracks <- c('hic.E14_NSC.score_k100','hic.E14_IPC.score_k100','hic.E14_PN.score_k100')
methyl_tracks <- c("methylation.E14_NSC_10x","methylation.E14_IPC_10x","methylation.E14_PN_10x")
iue_tracks <- c("hic.GFP_IUE24h.score_k100","hic.NGN2_IUE24h.score_k100","hic.others.GFP_IUE_rep1_score_k100","hic.others.GFP_IUE_rep2_score_k100","hic.others.GFP_IUE_rep3_score_k100","hic.others.NGN2_IUE_rep1_score_k100","hic.others.NGN2_IUE_rep2_score_k100","hic.others.NGN2_IUE_rep3_score_k100")
iue_cpg_tracks <- c("methylation.GFP_IUE24h_CpG_cov10x","methylation.NGN2_IUE24h_CpG_cov10x")
iue_gpc_tracks <- c("methylation.GFP_IUE24h_GpC_cov10x","methylation.NGN2_IUE24h_GpC_cov10x")

options(gmax.data.size=5e7)
options(gmultitasking=F)
for (i in 1:length(tracks)){
  #gvtrack.rm(paste0('v_',tracks[i]))
  gvtrack.create(paste0('v_',tracks[i]),tracks[i],'max')
  gvtrack.iterator.2d(paste0('v_',tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
}

for (i in 1:length(methyl_tracks)){
  gvtrack.create(paste0('v_',methyl_tracks[i]),methyl_tracks[i],'avg')
  gvtrack.iterator(paste0('v_',methyl_tracks[i]), sshift=-100, eshift=100)
}

########################

outMatch <- readRDS("results/P2G-Links.RDS")
linksSig = outMatch$posCor

grid <- gintervals.2d(chroms1 = seqnames(linksSig),starts1 = (start(linksSig)+end(linksSig))/2,ends1 = (start(linksSig)+end(linksSig))/2+1,chroms2 = seqnames(linksSig),starts2 = linksSig$gene_start,ends2 = linksSig$gene_start+1)
grid <- gintervals.canonic(grid)
grid$IDs <- as.vector(unite(grid[,c(1,2,5)],col='peaks_ID',sep='_'))
linksSig_IDs <- as.data.frame(cbind(as.vector(seqnames(linksSig)),as.numeric((start(linksSig)+end(linksSig))/2),as.numeric(linksSig$gene_start)))
linksSig$IDs <- unite(linksSig_IDs,col='peaks_ID',sep='_')[,1]

dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
intra_grid <- gintervals.intersect(grid,dom)
intra_grid1 <- intra_grid[intra_grid$start2>intra_grid$start1,]
intra_grid2 <- intra_grid[intra_grid$start1>intra_grid$start2,]
intra_grid2 <- intra_grid2[,c(4:6,1:3)]
colnames(intra_grid2) <- colnames(intra_grid1)
intra_indx = as.vector(unite(intra_grid[,c(1,2,5)],col='peaks_ID',sep='_'))[,1]		
inter_grid <- grid[!(grid$IDs[,1]%in%intra_indx),-7]
inter_grid1 <- inter_grid[inter_grid$start2>inter_grid$start1,]
inter_grid2 <- inter_grid[inter_grid$start1>inter_grid$start2,]
inter_grid2 <- inter_grid2[,c(4:6,1:3)]
colnames(inter_grid2) <- colnames(inter_grid1)

intra_scores1<- gextract(paste0('v_',tracks),intervals = intra_grid1,iterator = intra_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
intra_scores2 <- gextract(paste0('v_',tracks),intervals = intra_grid2,iterator = intra_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))
inter_scores1 <- gextract(paste0('v_',tracks),intervals = inter_grid1,iterator = inter_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
inter_scores2 <- gextract(paste0('v_',tracks),intervals = inter_grid2,iterator = inter_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))

intra_scores2 <- intra_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
colnames(intra_scores2)[1:6] <- colnames(intra_scores1)[1:6]
inter_scores2 <- inter_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
colnames(inter_scores2)[1:6] <- colnames(inter_scores1)[1:6]
intra_scores <- rbind(intra_scores1,intra_scores2)
inter_scores <- rbind(inter_scores1,inter_scores2)
intra_scores$peakID <- as.vector(unite(intra_scores[,c(1,2,5)],col='peaks_ID',sep='_')$peaks_ID)
inter_scores$peakID <- as.vector(unite(inter_scores[,c(1,2,5)],col='peaks_ID',sep='_')$peaks_ID)
intra_scores$domain <- 'intraTAD'
inter_scores$domain <- 'interTAD'
row.names(intra_scores) <- paste0(row.names(intra_scores),'_intra')
row.names(inter_scores) <- paste0(row.names(inter_scores),'_inter')
grid <- rbind(intra_scores,inter_scores)
grid <- grid[!is.na(grid$peakID),]
linksSig <- linksSig[match(grid$peakID,linksSig$IDs)]
linksSig$NSCscore <- grid$v_hic.E14_NSC.score_k100
linksSig$IPCscore <- grid$v_hic.E14_IPC.score_k100
linksSig$PNscore <- grid$v_hic.E14_PN.score_k100

linksSig$GFP_score <- grid$v_hic.GFP_IUE24h.score_k100
linksSig$NGN2_score <- grid$v_hic.NGN2_IUE24h.score_k100
linksSig$GFP_rep1_score <- grid$v_hic.others.GFP_IUE_rep1_score_k100
linksSig$GFP_rep2_score <- grid$v_hic.others.GFP_IUE_rep2_score_k100
linksSig$GFP_rep3_score <- grid$v_hic.others.GFP_IUE_rep3_score_k100
linksSig$NGN2_rep1_score <- grid$v_hic.others.NGN2_IUE_rep1_score_k100
linksSig$NGN2_rep2_score <- grid$v_hic.others.NGN2_IUE_rep2_score_k100
linksSig$NGN2_rep3_score <- grid$v_hic.others.NGN2_IUE_rep3_score_k100


linksSig$domain <- grid$domain
linksSig$deltaHiC <- rowMaxs(as.matrix(mcols(linksSig[,c('NSCscore','IPCscore','PNscore')])))-rowMins(as.matrix(mcols(linksSig[,c('NSCscore','IPCscore','PNscore')])))
#linksSig <- linksSig[order(linksSig$FDR,linksSig$gene_name,100-abs(linksSig$deltaHiC),decreasing=F),]

distal_grid <- gintervals(chroms = grid$chrom1,starts = grid$start1,ends = grid$end1)
distal_grid <- gintervals.canonic(distal_grid)
prom_grid <- gintervals(chroms = grid$chrom2,starts = grid$start2,ends = grid$end2)
prom_grid <- gintervals.canonic(prom_grid)

distal_meth <- gextract(paste0('v_',methyl_tracks),intervals = distal_grid,iterator = distal_grid)
colnames(distal_meth) <- gsub('v_methylation','distal',colnames(distal_meth))
distal_meth$intervalID <- paste0(distal_meth$chrom,':',distal_meth$start,'-',distal_meth$end)
prom_meth <- gextract(paste0('v_',methyl_tracks),intervals = prom_grid,iterator = prom_grid)
colnames(prom_meth) <- gsub('v_methylation','prom',colnames(prom_meth))
prom_meth$intervalID <- paste0(prom_meth$chrom,':',prom_meth$start,'-',prom_meth$end)

distalIDs <- paste0(seqnames(linksSig),':',(start(linksSig)+end(linksSig))/2,'-',(start(linksSig)+end(linksSig))/2+1)
promIDs <- paste0(seqnames(linksSig),':',linksSig$gene_start,'-',ends2 = linksSig$gene_start+1)
comb_meth <- cbind(distal_meth[match(distalIDs,distal_meth$intervalID),4:(ncol(distal_meth)-1)],prom_meth[match(promIDs,prom_meth$intervalID),4:(ncol(prom_meth)-1)])
mcols(linksSig) <- cbind(mcols(linksSig),comb_meth)

outMatch$posCor <- linksSig
saveRDS(outMatch, "results/P2G-Links.RDS")

######################################################
###Extract Hi-C connectivity for anticorrelated peaks ####
##########################################################

outMatch <- readRDS("results/P2G-Links.RDS")
peakLinks = outMatch$linksAll
peakLinks$sigCorrelation <- peakLinks$Correlation <= (corCutOff*(-1)) & peakLinks$FDR <= fdrCutOff & abs(peakLinks$distance) >= distCutOff
linksSig <- peakLinks[which(peakLinks$sigCorrelation)]

##### Compare cicero accessibility connections with Hi-C
grid <- gintervals.2d(chroms1 = seqnames(linksSig),starts1 = (start(linksSig)+end(linksSig))/2,ends1 = (start(linksSig)+end(linksSig))/2+1,chroms2 = seqnames(linksSig),starts2 = linksSig$gene_start,ends2 = linksSig$gene_start+1)
grid <- gintervals.canonic(grid)
grid$IDs <- as.vector(unite(grid[,c(1,2,5)],col='peaks_ID',sep='_'))
linksSig_IDs <- as.data.frame(cbind(as.vector(seqnames(linksSig)),as.numeric((start(linksSig)+end(linksSig))/2),as.numeric(linksSig$gene_start)))
linksSig$IDs <- unite(linksSig_IDs,col='peaks_ID',sep='_')[,1]

dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
intra_grid <- gintervals.intersect(grid,dom)
intra_grid1 <- intra_grid[intra_grid$start2>intra_grid$start1,]
intra_grid2 <- intra_grid[intra_grid$start1>intra_grid$start2,]
intra_grid2 <- intra_grid2[,c(4:6,1:3)]
colnames(intra_grid2) <- colnames(intra_grid1)
intra_indx = as.vector(unite(intra_grid[,c(1,2,5)],col='peaks_ID',sep='_'))[,1]		
inter_grid <- grid[!(grid$IDs[,1]%in%intra_indx),-7]
inter_grid1 <- inter_grid[inter_grid$start2>inter_grid$start1,]
inter_grid2 <- inter_grid[inter_grid$start1>inter_grid$start2,]
inter_grid2 <- inter_grid2[,c(4:6,1:3)]
colnames(inter_grid2) <- colnames(inter_grid1)

intra_scores1<- gextract(paste0('v_',tracks),intervals = intra_grid1,iterator = intra_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
intra_scores2 <- gextract(paste0('v_',tracks),intervals = intra_grid2,iterator = intra_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))
inter_scores1 <- gextract(paste0('v_',tracks),intervals = inter_grid1,iterator = inter_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
inter_scores2 <- gextract(paste0('v_',tracks),intervals = inter_grid2,iterator = inter_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))

intra_scores2 <- intra_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
colnames(intra_scores2)[1:6] <- colnames(intra_scores1)[1:6]
inter_scores2 <- inter_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
colnames(inter_scores2)[1:6] <- colnames(inter_scores1)[1:6]
intra_scores <- rbind(intra_scores1,intra_scores2)
inter_scores <- rbind(inter_scores1,inter_scores2)
intra_scores$peakID <- as.vector(unite(intra_scores[,c(1,2,5)],col='peaks_ID',sep='_')$peaks_ID)
inter_scores$peakID <- as.vector(unite(inter_scores[,c(1,2,5)],col='peaks_ID',sep='_')$peaks_ID)
intra_scores$domain <- 'intraTAD'
inter_scores$domain <- 'interTAD'
row.names(intra_scores) <- paste0(row.names(intra_scores),'_intra')
row.names(inter_scores) <- paste0(row.names(inter_scores),'_inter')
grid <- rbind(intra_scores,inter_scores)
grid <- grid[!is.na(grid$peakID),]
linksSig <- linksSig[match(grid$peakID,linksSig$IDs)]
linksSig$NSCscore <- grid$v_hic.E14_NSC.score_k100
linksSig$IPCscore <- grid$v_hic.E14_IPC.score_k100
linksSig$PNscore <- grid$v_hic.E14_PN.score_k100

linksSig$GFP_score <- grid$v_hic.GFP_IUE24h.score_k100
linksSig$NGN2_score <- grid$v_hic.NGN2_IUE24h.score_k100
linksSig$GFP_rep1_score <- grid$v_hic.others.GFP_IUE_rep1_score_k100
linksSig$GFP_rep2_score <- grid$v_hic.others.GFP_IUE_rep2_score_k100
linksSig$GFP_rep3_score <- grid$v_hic.others.GFP_IUE_rep3_score_k100
linksSig$NGN2_rep1_score <- grid$v_hic.others.NGN2_IUE_rep1_score_k100
linksSig$NGN2_rep2_score <- grid$v_hic.others.NGN2_IUE_rep2_score_k100
linksSig$NGN2_rep3_score <- grid$v_hic.others.NGN2_IUE_rep3_score_k100


linksSig$domain <- grid$domain
linksSig$deltaHiC <- rowMaxs(as.matrix(mcols(linksSig[,c('NSCscore','IPCscore','PNscore')])))-rowMins(as.matrix(mcols(linksSig[,c('NSCscore','IPCscore','PNscore')])))
#linksSig <- linksSig[order(linksSig$FDR,linksSig$gene_name,100-abs(linksSig$deltaHiC),decreasing=F),]

distal_grid <- gintervals(chroms = grid$chrom1,starts = grid$start1,ends = grid$end1)
distal_grid <- gintervals.canonic(distal_grid)
prom_grid <- gintervals(chroms = grid$chrom2,starts = grid$start2,ends = grid$end2)
prom_grid <- gintervals.canonic(prom_grid)

distal_meth <- gextract(paste0('v_',methyl_tracks),intervals = distal_grid,iterator = distal_grid)
colnames(distal_meth) <- gsub('v_methylation','distal',colnames(distal_meth))
distal_meth$intervalID <- paste0(distal_meth$chrom,':',distal_meth$start,'-',distal_meth$end)
prom_meth <- gextract(paste0('v_',methyl_tracks),intervals = prom_grid,iterator = prom_grid)
colnames(prom_meth) <- gsub('v_methylation','prom',colnames(prom_meth))
prom_meth$intervalID <- paste0(prom_meth$chrom,':',prom_meth$start,'-',prom_meth$end)

distalIDs <- paste0(seqnames(linksSig),':',(start(linksSig)+end(linksSig))/2,'-',(start(linksSig)+end(linksSig))/2+1)
promIDs <- paste0(seqnames(linksSig),':',linksSig$gene_start,'-',ends2 = linksSig$gene_start+1)
comb_meth <- cbind(distal_meth[match(distalIDs,distal_meth$intervalID),4:(ncol(distal_meth)-1)],prom_meth[match(promIDs,prom_meth$intervalID),4:(ncol(prom_meth)-1)])
mcols(linksSig) <- cbind(mcols(linksSig),comb_meth)

outMatch$negCor <- linksSig
saveRDS(outMatch, "results/P2G-Links.RDS")

######################################################
###Extract Hi-C connectivity for non-correlated peaks ####
##########################################################

outMatch <- readRDS("results/P2G-Links.RDS")
peakLinks = outMatch$linksAll
peakLinks$sigCorrelation <- peakLinks$Correlation >= (corCutOff*(-1)) & peakLinks$Correlation <= corCutOff & abs(peakLinks$distance) >= distCutOff
linksSig <- peakLinks[which(peakLinks$sigCorrelation)]
linksSig <- linksSig[sample(1:length(linksSig),length(outMatch$linksSig)),]
##### Compare cicero accessibility connections with Hi-C
grid <- gintervals.2d(chroms1 = seqnames(linksSig),starts1 = (start(linksSig)+end(linksSig))/2,ends1 = (start(linksSig)+end(linksSig))/2+1,chroms2 = seqnames(linksSig),starts2 = linksSig$gene_start,ends2 = linksSig$gene_start+1)
grid <- gintervals.canonic(grid)
grid$IDs <- as.vector(unite(grid[,c(1,2,5)],col='peaks_ID',sep='_'))
linksSig_IDs <- as.data.frame(cbind(as.vector(seqnames(linksSig)),as.numeric((start(linksSig)+end(linksSig))/2),as.numeric(linksSig$gene_start)))
linksSig$IDs <- unite(linksSig_IDs,col='peaks_ID',sep='_')[,1]

dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
intra_grid <- gintervals.intersect(grid,dom)
intra_grid1 <- intra_grid[intra_grid$start2>intra_grid$start1,]
intra_grid2 <- intra_grid[intra_grid$start1>intra_grid$start2,]
intra_grid2 <- intra_grid2[,c(4:6,1:3)]
colnames(intra_grid2) <- colnames(intra_grid1)
intra_indx = as.vector(unite(intra_grid[,c(1,2,5)],col='peaks_ID',sep='_'))[,1]		
inter_grid <- grid[!(grid$IDs[,1]%in%intra_indx),-7]
inter_grid1 <- inter_grid[inter_grid$start2>inter_grid$start1,]
inter_grid2 <- inter_grid[inter_grid$start1>inter_grid$start2,]
inter_grid2 <- inter_grid2[,c(4:6,1:3)]
colnames(inter_grid2) <- colnames(inter_grid1)

intra_scores1<- gextract(paste0('v_',tracks),intervals = intra_grid1,iterator = intra_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
intra_scores2 <- gextract(paste0('v_',tracks),intervals = intra_grid2,iterator = intra_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))
inter_scores1 <- gextract(paste0('v_',tracks),intervals = inter_grid1,iterator = inter_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
inter_scores2 <- gextract(paste0('v_',tracks),intervals = inter_grid2,iterator = inter_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))

intra_scores2 <- intra_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
colnames(intra_scores2)[1:6] <- colnames(intra_scores1)[1:6]
inter_scores2 <- inter_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
colnames(inter_scores2)[1:6] <- colnames(inter_scores1)[1:6]
intra_scores <- rbind(intra_scores1,intra_scores2)
inter_scores <- rbind(inter_scores1,inter_scores2)
intra_scores$peakID <- as.vector(unite(intra_scores[,c(1,2,5)],col='peaks_ID',sep='_')$peaks_ID)
inter_scores$peakID <- as.vector(unite(inter_scores[,c(1,2,5)],col='peaks_ID',sep='_')$peaks_ID)
intra_scores$domain <- 'intraTAD'
inter_scores$domain <- 'interTAD'
row.names(intra_scores) <- paste0(row.names(intra_scores),'_intra')
row.names(inter_scores) <- paste0(row.names(inter_scores),'_inter')
grid <- rbind(intra_scores,inter_scores)
grid <- grid[!is.na(grid$peakID),]
linksSig <- linksSig[match(grid$peakID,linksSig$IDs)]
linksSig$NSCscore <- grid$v_hic.E14_NSC.score_k100
linksSig$IPCscore <- grid$v_hic.E14_IPC.score_k100
linksSig$PNscore <- grid$v_hic.E14_PN.score_k100

linksSig$GFP_score <- grid$v_hic.GFP_IUE24h.score_k100
linksSig$NGN2_score <- grid$v_hic.NGN2_IUE24h.score_k100
linksSig$GFP_rep1_score <- grid$v_hic.others.GFP_IUE_rep1_score_k100
linksSig$GFP_rep2_score <- grid$v_hic.others.GFP_IUE_rep2_score_k100
linksSig$GFP_rep3_score <- grid$v_hic.others.GFP_IUE_rep3_score_k100
linksSig$NGN2_rep1_score <- grid$v_hic.others.NGN2_IUE_rep1_score_k100
linksSig$NGN2_rep2_score <- grid$v_hic.others.NGN2_IUE_rep2_score_k100
linksSig$NGN2_rep3_score <- grid$v_hic.others.NGN2_IUE_rep3_score_k100


linksSig$domain <- grid$domain
linksSig$deltaHiC <- rowMaxs(as.matrix(mcols(linksSig[,c('NSCscore','IPCscore','PNscore')])))-rowMins(as.matrix(mcols(linksSig[,c('NSCscore','IPCscore','PNscore')])))
#linksSig <- linksSig[order(linksSig$FDR,linksSig$gene_name,100-abs(linksSig$deltaHiC),decreasing=F),]

distal_grid <- gintervals(chroms = grid$chrom1,starts = grid$start1,ends = grid$end1)
distal_grid <- gintervals.canonic(distal_grid)
prom_grid <- gintervals(chroms = grid$chrom2,starts = grid$start2,ends = grid$end2)
prom_grid <- gintervals.canonic(prom_grid)

distal_meth <- gextract(paste0('v_',methyl_tracks),intervals = distal_grid,iterator = distal_grid)
colnames(distal_meth) <- gsub('v_methylation','distal',colnames(distal_meth))
distal_meth$intervalID <- paste0(distal_meth$chrom,':',distal_meth$start,'-',distal_meth$end)
prom_meth <- gextract(paste0('v_',methyl_tracks),intervals = prom_grid,iterator = prom_grid)
colnames(prom_meth) <- gsub('v_methylation','prom',colnames(prom_meth))
prom_meth$intervalID <- paste0(prom_meth$chrom,':',prom_meth$start,'-',prom_meth$end)

distalIDs <- paste0(seqnames(linksSig),':',(start(linksSig)+end(linksSig))/2,'-',(start(linksSig)+end(linksSig))/2+1)
promIDs <- paste0(seqnames(linksSig),':',linksSig$gene_start,'-',ends2 = linksSig$gene_start+1)
comb_meth <- cbind(distal_meth[match(distalIDs,distal_meth$intervalID),4:(ncol(distal_meth)-1)],prom_meth[match(promIDs,prom_meth$intervalID),4:(ncol(prom_meth)-1)])
mcols(linksSig) <- cbind(mcols(linksSig),comb_meth)

outMatch$noCor <- linksSig

saveRDS(outMatch, "results/P2G-Links.RDS")
names(outMatch)[3:6] <- c('posCor','all','negCor','noCor')
links <- outMatch$posCor
peaks <- data.frame(peak_id=links$peakName,gene_short_name=links$gene_name,coaccess=links$Correlation)
write.csv(x = peaks, file = "results/cellOracle/posCor_connections.csv")

links <- outMatch$negCor
peaks <- data.frame(peak_id=links$peakName,gene_short_name=links$gene_name,coaccess=links$Correlation)
write.csv(x = peaks, file = "results/cellOracle/negCor_connections.csv")

links <- outMatch$noCor
peaks <- data.frame(peak_id=links$peakName,gene_short_name=links$gene_name,coaccess=links$Correlation)
write.csv(x = peaks, file = "results/cellOracle/noCor_connections.csv")
