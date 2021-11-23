library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(magrittr)
library(ggplot2)
library(Rcpp)
library(viridis)
library(reshape2)
library(org.Mm.eg.db)
library(tidyr)

#--------------------------------------------
# Functions
#--------------------------------------------

sourceCpp(code='
  #include <Rcpp.h>
          
          using namespace Rcpp;
          using namespace std;
          
          // [[Rcpp::export]]
          IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
          if(x1.size() != y1.size()){
          stop("width must equal size!");
          }
          IntegerVector x = clone(x1);
          IntegerVector y = clone(y1);
          int n = x.size();
          IntegerVector rx = seq(xmin,xmax);
          IntegerVector ry = seq(ymin,ymax);
          IntegerMatrix mat( ry.size() , rx.size() );
          int xi,yi;
          for(int i = 0; i < n; i++){
          xi = (x[i] - xmin);
          yi = (y[i] - ymin);
          if(yi >= 0 && yi < ry.size()){
          if(xi >= 0 && xi < rx.size()){
          mat( yi , xi ) = mat( yi , xi ) + 1; 
          }
          }
          }
          return mat;
          }'
)

insertionProfileSingles <- function(feature, fragments, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){
  
  insertionProfileSingles_helper <- function(feature, fragments, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){
    #Convert To Insertion Sites
    if(getInsertions){
      insertions <- c(
        GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
        GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
      )
      by <- "RG"
    }else{
      insertions <- fragments
    }
    remove(fragments)
    gc()
    
    #center the feature
    center <- unique(resize(feature, width = 1, fix = fix, ignore.strand = FALSE))
    
    #get overlaps between the feature and insertions only up to flank bp
    overlap <- DataFrame(findOverlaps(query = center, subject = insertions, maxgap = flank, ignore.strand = TRUE))
    overlap$strand <- strand(center)[overlap[,1]]
    overlap$name <- mcols(insertions)[overlap[,2],by]
    overlap <- transform(overlap, id=match(name, unique(name)))
    ids <- length(unique(overlap$name))
    
    #distance
    overlap$dist <- NA
    minus <- which(overlap$strand == "-")
    other <- which(overlap$strand != "-")
    overlap$dist[minus] <- start(center[overlap[minus,1]]) - start(insertions[overlap[minus,2]])
    overlap$dist[other] <- start(insertions[overlap[other,2]]) - start(center[overlap[other,1]])
    
    #Insertion Mat
    profile_mat <- tabulate2dCpp(x1 = overlap$id, y1 = overlap$dist, xmin = 1, xmax = ids, ymin = -flank, ymax = flank)
    colnames(profile_mat) <- unique(overlap$name)
    profile <- rowSums(profile_mat)
    
    #normalize
    profile_mat_norm <- apply(profile_mat, 2, function(x) x/max(mean(x[c(1:norm,(flank*2-norm+1):(flank*2+1))]), 0.5)) #Handles low depth cells
    profile_norm <- profile/mean(profile[c(1:norm,(flank*2-norm+1):(flank*2+1))])
    
    #smooth
    profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) zoo::rollmean(x, smooth, fill = 1))
    profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)
    
    #enrichment
    max_finite <- function(x){
      suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
    }
    e_mat <- apply(profile_mat_norm_smooth, 2, function(x) max_finite(x[(flank-range):(flank+range)]))
    names(e_mat) <- colnames(profile_mat_norm_smooth)
    e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])
    
    #Summary
    df_mat <- data.frame(
      enrichment = e_mat,
      insertions = as.vector(table(mcols(insertions)[,by])[names(e_mat)]),
      insertionsWindow = as.vector(table(overlap$name)[names(e_mat)])
    )
    df_sum <- data.frame(bp = (-flank):flank, profile = profile, norm_profile = profile_norm, smooth_norm_profile = profile_norm_smooth, enrichment = e)
    rownames(df_sum) <-  NULL
    
    return(list(df = df_sum, dfall = df_mat, profileMat = profile_mat_norm, profileMatSmooth = profile_mat_norm_smooth))
  }
  
  uniqueTags <- as.character(unique(mcols(fragments)[,by]))
  splitTags <- split(uniqueTags, ceiling(seq_along(uniqueTags)/batchSize))
  
  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  batchTSS <- lapply(seq_along(splitTags), function(x){
    setTxtProgressBar(pb, round(x * 100/length(splitTags), 0))
    profilex <- insertionProfileSingles_helper(
      feature=feature, 
      fragments=fragments[which(mcols(fragments)[,by] %in% splitTags[[x]])], 
      by = by, 
      getInsertions = getInsertions,
      fix = fix, 
      flank = flank, 
      norm = norm, 
      smooth = smooth, 
      range = range
    )
    
    return(profilex)
  })
  df <- lapply(batchTSS, function(x) x$df) %>% Reduce("rbind",.)
  dfall <- lapply(batchTSS, function(x) x$dfall) %>% Reduce("rbind",.)
  profileMat <- lapply(batchTSS, function(x) x$profileMat) %>% Reduce("cbind",.)
  profileMatSmooth <- lapply(batchTSS, function(x) x$profileMatSmooth) %>% Reduce("cbind",.)
  return(list(df = df, dfall = dfall, profileMat = profileMat, profileMatSmooth = profileMatSmooth))
}


countInsertions <- function(query, fragments, by = "RG"){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

extendedPeakSet <- function(df, BSgenome = NULL, extend = 250, blacklist = NULL, nSummits = 100000){
  #Helper Functions
  readSummits <- function(file){
    df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end","name","score"))))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
  }
  nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
    stopifnot(by %in% colnames(mcols(gr)))
    clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
      o <- findOverlaps(gr,r)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      gr <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
    }
    if(verbose){
      message("Converging", appendLF = FALSE)
    }
    i <-  0
    gr_converge <- gr
    while(length(gr_converge) > 0){
      if(verbose){
        message(".", appendLF = FALSE)
      }
      i <-  i + 1
      gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
      gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
      if(i == 1){ #if i=1 then set gr_all to clustered
        gr_all <- gr_selected
      }else{
        gr_all <- c(gr_all, gr_selected)
      } 
    }
    if(verbose){
      message("\nSelected ", length(gr_all), " from ", length(gr))
    }
    gr_all <- sort(sortSeqlevels(gr_all))
    return(gr_all)
  }
  #Check-------
  stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(BSgenome))
  stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- rtracklayer::import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))
  #------------
  #Time to do stuff
  chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  groups <- unique(df$groups)
  groupGRList <- GenomicRanges::GenomicRangesList(lapply(seq_along(groups), function(i){
    df_group = df[which(df$groups==groups[i]),]
    grList <- GenomicRanges::GenomicRangesList(lapply(paste0(df_group$summits), function(x){
      extended_summits <- readSummits(x) %>%
        resize(., width = 2 * extend + 1, fix = "center") %>%     
        subsetByOverlaps(.,chromSizes,type="within") %>%
        subsetByOverlaps(.,blacklist,invert=TRUE) %>%
        nonOverlappingGRanges(., by="score", decreasing=TRUE)
      extended_summits <- extended_summits[order(extended_summits$score,decreasing=TRUE)]
      if(!is.null(nSummits)){
        extended_summits <- head(extended_summits, nSummits)
      }
      mcols(extended_summits)$scoreQuantile <- trunc(rank(mcols(extended_summits)$score))/length(mcols(extended_summits)$score)
      extended_summits
    }))
    #Non Overlapping
    grNonOverlapping <- nonOverlappingGRanges(unlist(grList), by = "scoreQuantile", decreasing = TRUE)
    #Free Up Memory
    remove(grList)
    gc()
    grNonOverlapping
  }))
  grFinal <- nonOverlappingGRanges(unlist(groupGRList), by = "scoreQuantile", decreasing = TRUE)
  grFinal <- sort(sortSeqlevels(grFinal))
  return(grFinal)
}

groupSums <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

ClusterFragments <- function(
  reads,
  reads2=NULL,
  cells,
  output.path,
  assume.sorted = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Retaining ", length(x = cells), " cells")
    message("Reading fragments")
  }
  reads <- reads[reads$cell %in% cells, ]
  if (!is.null(reads2)){
    reads2 <- reads2[reads2$cell %in% cells, ]
    reads <- rbind(reads,reads2)    
  } 
  if (!assume.sorted) {
    if (verbose) {
      message("Sorting fragments")
    }
    reads <- reads[with(data = reads, expr = order(chr, start)), ]
  }
  if (verbose) {
    message("Writing output")
  }
  fwrite(
    x = reads,
    file = output.path,
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE,
    sep = '\t'
  )
  rm(reads)
  invisible(x = gc())
}

countInsertions2 <- function(query, fragments,cells, by = "RG",window){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(cells)))
  colnames(sparseM) <- c(unique(overlapDF$name),cells[!cells%in%unique(overlapDF$name)])
  row.names(sparseM) <- start(resize(query,width = 1))+window/2
  return(sparseM)
}

grToFeature <- function(gr,sep='-'){
  #gr <- separate(as.data.frame(gr),col='gr', c("chrom", "bp1","bp2"), sep)
  peakinfo <- data.frame(
    row.names = paste(seqnames(gr),start(gr),end(gr),sep="_"),
    site_name = paste(seqnames(gr),start(gr),end(gr),sep="_"),
    chr = gsub("chr","",as.character(seqnames(gr))),
    bp1 = start(gr),
    bp2 = end(gr)
  )
  return(peakinfo)
}

featureToGR <- function(feature,pattern="_"){
  featureSplit <- stringr::str_split(paste0(feature), pattern =pattern , n = 3, simplify = TRUE)
  gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
  return(gr)
}

makeCDS <- function(se, binarize = TRUE){
  peakinfo <- separate(as.data.frame(row.names(se)),col=1,into=c("chr", "bp1","bp2"), sep="-",convert=T)
  peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
  row.names(peakinfo) <- peakinfo$site_name
  mat <- GetAssayData(se,slot='counts')
  row.names(mat) <- gsub('-','_',row.names(mat))
  if(binarize){
    mat@x[which(mat@x > 0)] <- 1
  }
  cellinfo <- se@meta.data
  cellinfo$cells <- rownames(cellinfo)
  cellinfo <- cellinfo[cellinfo$cells%in%colnames(mat),]
  cds <-  suppressWarnings(new_cell_data_set(mat,
                                             cell_metadata = cellinfo,
                                             gene_metadata = peakinfo))
  cds <- cds[order(fData(cds)$chr, fData(cds)$bp1),]
  return(cds)
}


sourceCpp(code='
          #include <Rcpp.h>
          
          using namespace Rcpp;
          using namespace std;
          
          // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
          // [[Rcpp::export]]
          Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
          
          if(X.ncol() != Y.ncol()){
          stop("Columns of Matrix X and Y must be equal length!");
          }
          
          if(max(idxX) > X.nrow()){
          stop("Idx X greater than nrow of Matrix X");
          }
          
          if(max(idxY) > Y.nrow()){
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

getTxDbGenes <- function(txdb = NULL, orgdb = NULL, gr = NULL, ignore.strand = TRUE){
  
  if (is.null(genome)) {
    if (is.null(txdb) | is.null(orgdb)) {
      stop("If no provided genome then you need txdb and orgdb!")
    }
  }
  
  if (is.null(gr)) {
    genes <- GenomicFeatures::genes(txdb)
  }else {
    genes <- suppressWarnings(subsetByOverlaps(GenomicFeatures::genes(txdb), gr, ignore.strand = ignore.strand))
  }
  
  if (length(genes) > 1) {
    mcols(genes)$symbol <- suppressMessages(mapIds(orgdb, 
                                                   keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID", 
                                                   multiVals = "first"))
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
    names(genes) <- NULL
    out <- genes
  }else {
    out <- GRanges(seqnames(gr), ranges = IRanges(0, 0), gene_id = 0, symbol = "none")[-1]
  }
  
  return(out)
  
}

construct.grid = function(interv1,interv2,min_dist,max_dist){
  return(ddply(interv1, .(chrom), function(i1) {
    i2 = interv2[as.character(interv2$chrom) == as.character(i1$chrom[1]),]
    if (nrow(i2) ==0) {
      return(c())
    }
    g = expand.grid(1:nrow(i1), 1:nrow(i2))
    g = g[i2$start[g$Var2]-i1$start[g$Var1] > min_dist & i2$start[g$Var2]-i1$start[g$Var1] < max_dist,]
    grid = cbind(i1[g$Var1,c("chrom", "start", "end")], i2[g$Var2,c("chrom", "start", "end")])
    colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    return(grid)
  })[,-1] )
}

construct.grid2 = function(interv1,interv2,min_dist,max_dist){
  return(ddply(interv1, .(chrom), function(i1) {
    i2 = interv2[as.character(interv2$chrom) == as.character(i1$chrom[1]),]
    if (nrow(i2) ==0) {
      return(c())
    }
    g = expand.grid(1:nrow(i1), 1:nrow(i2))
    g = g[i2$start[g$Var2]-i1$start[g$Var1] > min_dist &
            i2$start[g$Var2]-i1$start[g$Var1] < max_dist,]
    grid = cbind(i1[g$Var1,c("chrom", "start", "end")], i2[g$Var2,c("chrom", "start", "end")])
    colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    return(grid)
  })[,-1] )
}


groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

groupSds <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE) {
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gs <- lapply(unique(groups), function(x) {
    if (sparse) {
      matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
    }
    else {
      matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gs) <- unique(groups)
  return(gs)
}

createPseudoBulk <- function(mat, groups, labels, minCells = 100, maxCells = 500, minReps = 3, ceiling = 1, prior.count = 3, nSim = 1000, distMethod = "vars", seed = 1){
  
  calcDiff <- function(mat, method = "vars"){
    if(tolower(method)=="vars"){
      sum(matrixStats::rowVars(mat))/ncol(mat)
    }else if(tolower(method)=="euclidean"){
      sum(dist(t(mat))/ncol(mat))
    }else{
      stop("Error method not found!")
    }
  }
  
  sumCells <- function(mat, groups, maxCells = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
      idx <- which(groups == x)
      if(!is.null(maxCells)){
        idx <- sample(idx, size = min(maxCells, length(idx)), replace = FALSE)
      }
      if (sparse) {
        Matrix::rowSums(mat[, idx, drop = FALSE], na.rm = na.rm)
      }
      else{
        rowSums(mat[, idx, drop = FALSE], na.rm = na.rm)
      }
    }) %>% Reduce("cbind", .) %>% as.matrix
    colnames(gm) <- unique(groups)
    return(gm)
  }
  
  message(paste0("Setting Seed = ",seed))
  set.seed(seed)
  
  names(labels) <- colnames(mat)
  groupList <- split(labels,groups)
  
  if(minReps <= 1){
    stop("Minimum 2 replicates!")
  }
  
  if(is.numeric(ceiling)){
    message(paste0("Setting ceiling of input matrix to ", ceiling, "!"))
    mat@x[mat@x>ceiling]<-ceiling
  }
  
  #--------------------------------------------------------------------------
  # Constructing Bulk Pseudo ATAC Matrix v1.0
  #--------------------------------------------------------------------------
  
  pseudoAll <- lapply(seq_along(groupList), function(x){
    
    message(sprintf("####################################\n  Groups %s of %s : %s\n####################################",x,length(groupList), names(groupList)[x]))
    start <- Sys.time()
    groupx <- groupList[[x]]
    matx <- mat[,names(groupx)]
    bioReps <- names(table(groupx))[which(table(groupx)>minCells)]
    nToSim <- minReps - length(bioReps)
    
    if(length(bioReps) >= minReps){
      
      #--------------------------------------------------------------------------
      # If there is enough biological samples passing the minimal cells, great
      # good to go! Just merge into true pseudo bulk replicates!
      #--------------------------------------------------------------------------
      
      message(sprintf("Found %s Bio Reps which is more or equal than required (%s)", length(bioReps), minReps))
      nBio <- length(bioReps)
      groupBio <- groupx[which(groupx %in% bioReps)]
      pseudoBio <- sumCells(matx[,names(groupBio)],groups=groupBio,sparse=TRUE,na.rm=TRUE, maxCells = maxCells)
      nBio <- table(groupBio)[bioReps]
      if(!is.null(maxCells)){
        nBio[nBio > maxCells] <- maxCells
      }
      
      colnames(pseudoBio) <- lapply(seq_along(bioReps),function(k){
        paste0(names(groupList)[x],"._.BRep_",bioReps[k],".",nBio[which(names(nBio)==bioReps[k])],".FALSE")
      })
      pseudoMat <- pseudoBio
      
    }else if(length(bioReps) > 0 & ((length(groupx[groupx %ni% bioReps]) + 1) / min(nToSim, 2)) > minCells){
      
      #--------------------------------------------------------------------------
      # If there is at least 1 biological sample with the minimum cells but not
      # as many as required, we will make pseudo replicates with the true replicate
      # to attempt to capture real biological varation
      #--------------------------------------------------------------------------
      
      message("PSA : To ensure minimum replicates, simulation must be performed!")
      groupBio <- groupx[which(groupx %in% bioReps)]
      
      nBio <- table(groupBio)
      if(length(bioReps) == 1){
        pseudoBio <- Matrix::rowSums(matx[,names(groupBio)],na.rm=TRUE)
      }else{
        pseudoBio <- sumCells(matx[,names(groupBio)],groups=groupBio,sparse=TRUE,na.rm=TRUE, maxCells = maxCells)
      }
      pseudoBioLog <- edgeR::cpm(pseudoBio, log = TRUE, prior.count = prior.count)
      
      #Determine how sampling is to be performed, ideally we could have the minimum cells * non bio cells of cells!
      nSplit <- floor(length(groupx[groupx %ni% bioReps])/(nToSim))
      if(!is.null(maxCells)){
        nSplit <- min(nSplit, maxCells)
      }
      if(nSplit < minCells){
        message("Splitting cluster into overlapping cells using sampling with replacement BE CAREFUL OF LOW VARIANCE!!!!")
        replacement <- TRUE
        nSplit <- minCells
      }else{
        replacement <- FALSE
      }
      
      for(i in seq_len(nSim)){
        
        #Figure out how to split Matrix!
        randOrder <- sample(seq_len(ncol(matx)), ncol(matx))
        if(replacement){
          splitList <- lapply(seq_len(nToSim), function(x) sample(seq_len(ncol(matx)), size = nSplit, replace = replacement))
        }else{
          splitList <- split(randOrder, ceiling(seq_along(randOrder)/nSplit))[seq_len(nToSim)]
        }
        
        if(i == 1){
          
          pseudoMin <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)
          pseudoMinLog <- edgeR::cpm(pseudoMin, log = TRUE, prior.count = prior.count)
          diffMax <- calcDiff(cbind(pseudoBioLog, pseudoMinLog), method = distMethod)
          
        }else{
          
          pseudoI <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)
          
          pseudoILog <- edgeR::cpm(pseudoI, log = TRUE, prior.count = prior.count)
          diffI <- calcDiff(cbind(pseudoBioLog, pseudoILog), method = distMethod)
          message(sprintf("Trial %s, Current Distance = %s , Max Distance = %s", i, diffI, diffMax))
          
          if(diffI > diffMax){
            message("Found new maxima pseudo to be conservative...")
            pseudoMin <- pseudoI
            pseudoMinLog <- pseudoILog
            diffMax <- diffI
          }
          
        }
        
      }
      
      pseudoBio <- as.matrix(pseudoBio)
      pseudoMin <- as.matrix(pseudoMin)
      colnames(pseudoBio) <- lapply(seq_along(bioReps),function(k){
        paste0(names(groupList)[x],"._.BRep_",bioReps[k],".",nBio[which(names(nBio)==bioReps[k])],".FALSE")
      })
      colnames(pseudoMin) <- paste0(names(groupList)[x],"._.Rep",seq_len(nToSim),".",nSplit,".",replacement)
      pseudoMat <- cbind(pseudoBio, pseudoMin)
      
    }else{
      
      #--------------------------------------------------------------------------
      # If there is not at least 1 sample with enough cells we will bootstrap replicates.
      # This is not preferred as we will have a large underestimate of true biological
      # variation.
      #--------------------------------------------------------------------------
      
      message("PSA : No representation by at least one separate rep, please be cautious with result!")
      
      #Determine how sampling is to be performed, ideally we could have the minimum cells * non bio cells of cells!
      nToSim <- minReps
      nSplit <- floor(length(groupx[groupx %ni% bioReps])/(nToSim))
      if(floor(1.5 * minCells) > length(groupx)){
        nToSim <- 2
        nSplit <- floor(2 / 3 * length(groupx))
        message(sprintf("Warning! Group size (%s) is smaller than the 3/2 * minimal number of cells (%s)",length(groupx),floor(1.5 * minCells)))
        message(sprintf("To deal with this, we will sample %s replicates at 2/3 the group size (%s of %s)", nToSim, nSplit, length(groupx)))
      }
      if(!is.null(maxCells)){
        nSplit <- min(nSplit, maxCells)
      }
      if(nSplit < minCells){
        message("Splitting cluster into overlapping cells using sampling with replacement BE CAREFUL OF LOW VARIANCE!!!!")
        replacement <- TRUE
        nSplit <- minCells
      }else{
        replacement <- FALSE
      }
      
      for(i in seq_len(nSim)){
        
        #Figure out how to split Matrix!
        randOrder <- sample(seq_len(ncol(matx)), ncol(matx))
        if(replacement){
          splitList <- lapply(seq_len(nToSim), function(x) sample(seq_len(ncol(matx)), size = minCells, replace = replacement))
        }else{
          splitList <- split(randOrder, ceiling(seq_along(randOrder)/nSplit))[seq_len(nToSim)]
        }
        
        if(i == 1){
          
          pseudoMin <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)
          pseudoMinLog <- edgeR::cpm(pseudoMin, log = TRUE, prior.count = prior.count)
          diffMax <- calcDiff(pseudoMinLog, method = distMethod)
          
        }else{
          
          pseudoI <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)
          pseudoILog <- edgeR::cpm(pseudoI, log = TRUE, prior.count = prior.count)
          diffI <- calcDiff(pseudoILog, method = distMethod)
          
          message(sprintf("Trial %s, Current Distance = %s , Max Distance = %s", i, diffI, diffMax))
          if(diffI > diffMax){
            message("Found new maxima pseudo to be conservative...")
            pseudoMin <- pseudoI
            pseudoMinLog <- pseudoILog
            diffMax <- diffI
          }
          
        }
      }
      
      pseudoMin <- as.matrix(pseudoMin)
      colnames(pseudoMin) <- paste0(names(groupList)[x],"._.Rep_",seq_len(nToSim),".",nSplit,".",replacement)
      pseudoMat <- pseudoMin
      
    }
    
    print(Sys.time() - start)
    
    pseudoMat 
    
  }) %>% Reduce("cbind",.)
  
  out <- list(pseudoMat = pseudoAll, groupMat = sumCells(mat, groups=groups, sparse=TRUE, na.rm=TRUE))
  
  return(out)
  
}


uniqueFeatures <- function(
  mat, groups, padj = 0.01, minSdRatio = 0.001, 
  minLFC = 0.25, zCutoff = 1, breakPt = "last",
  padjMethod = "fdr", clusterCols = FALSE, 
  sparse = FALSE, twoWay = FALSE, groupMin = 10, 
  minGroupSize = 1, maxGroupSize = NULL){
  
  #----------------------
  # Functions
  #----------------------
  binarizeMatrix <- function(matMean, matSd, cutoff, method, minSdRatio, minLFC){
    
    binarizeVector <- function(vMean, vSd, cutoff, method, minSdRatio, minLFC){
      
      #Method
      if(method == "mean"){
        vTest = vMean
      }else{
        vTest <- vMean - cutoff * vSd    
      }
      
      #Order lowest to highest
      idx <- order(vTest)
      vTest <- vTest[idx]
      vMean <- vMean[idx]
      vSd <- vSd[idx]
      
      #Set which are too low for evaluation ie low Sd that is probably due to 0s or weird offsets
      sdToLow <- which(vSd < minSdRatio*vMean | vSd < 10^-5)
      vBinarySd <- rep(1,length(vSd))
      if(length(sdToLow) > 0){
        vBinarySd[sdToLow] <- 0
      }
      
      #Initialize
      vMeanCutSd <- vMean + cutoff*vSd;
      maxValue <- vMeanCutSd[1]
      maxBSd <- vBinarySd[1]
      maxMean <- vMean[1]
      
      #Create out vector and initialize breakpoint
      n <- length(vMean)
      out <- rep(0, n)
      breakPoint <- 0
      out[1] <- breakPoint
      
      #Evaluate
      for(i in seq(2,n)){
        
        #Check if break point assuming log space
        if((vTest[i] - maxValue) > 0 & (vMean[i] - maxMean) >= minLFC & maxBSd != 0){
          breakPoint <- breakPoint + 1
        }
        
        #Set current value of break point
        out[i] <- breakPoint
        
        #Keep Max value observed
        if(vMeanCutSd[i] > maxValue){
          maxValue <- vMeanCutSd[i]
          maxMean <- vMean[i]
          maxBSd <- vBinarySd[i]
        }
      }
      
      out <- out[order(idx)]
      return(out)
      
    }
    
    #Create binary matrix
    bMat <- matrix(NA,nrow=nrow(matMean),ncol=ncol(matMean))
    for(i in seq_len(nrow(matMean))){
      if(i%%5000==0){message(sprintf("%s of %s (percent = %s)", i, nrow(bMat), round(100*i/nrow(bMat),2)))}
      bMat[i,] <- binarizeVector(matMean[i,],matSd[i,],cutoff, method, minSdRatio, minLFC) 
    }
    
    #Add names
    colnames(bMat) <- colnames(matMean)
    rownames(bMat) <- rownames(matMean)
    return(bMat)
    
  }
  
  idxRow <- seq_len(nrow(mat))
  #-----------------------------------------------
  #Within Group Statistics
  #-----------------------------------------------
  message("Getting Within Group Stats...")
  intraMean <- groupMeans(mat, groups = groups, sparse = sparse, na.rm = TRUE)
  rownames(intraMean) <- rownames(mat)
  
  intraSd <- groupSds(mat, groups = groups, sparse = sparse, na.rm = TRUE)
  rownames(intraSd) <- rownames(mat)
  
  #-----------------------------------------------
  #Binarize Rows of Matrix
  #-----------------------------------------------
  message("Binarizing Features...")
  if (twoWay) {
    binarizedMat <- binarizeMatrix(intraMean, intraSd, zCutoff, "meanSd", minSdRatio, minLFC)
  }else {
    binarizedMat <- binarizeMatrix(intraMean, intraSd, zCutoff, "mean", minSdRatio, minLFC)
  }
  colnames(binarizedMat) <- colnames(intraMean)
  rownames(binarizedMat) <- rownames(intraMean)
  for(i in seq_len(nrow(binarizedMat))){
    bvi <- binarizedMat[i,]
    if(tolower(breakPt) == "last"){
      bvi[bvi!=max(bvi)] <- 0
      bvi[bvi>0] <- 1
    }else{
      bvi[bvi<1] <- 0
      bvi[bvi>0] <- 1
    }
    binarizedMat[i,] <- bvi
  }
  message(sprintf("Successful Binarization of %s Features...", sum(rowSums(binarizedMat) > 0)))
  
  #-----------------------------------------------
  #Get Test Statistics
  #-----------------------------------------------
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  groupSplit <- split(colnames(mat), groups)
  bpval <- unlist(lapply(seq_len(nrow(binarizedMat)), function(i){
    setTxtProgressBar(pb,round(i*100/nrow(binarizedMat),0))
    if(any(binarizedMat[i,]>0)){
      cu <- as.character(unlist(groupSplit[names(which(binarizedMat[i,]==max(binarizedMat[i,])))]))
      cd <- as.character(unlist(groupSplit[names(which(binarizedMat[i,]!=max(binarizedMat[i,])))]))
      mati <- mat[i,,drop=F]
      pval <- t.test(mati[,cu],mati[,cd])$p.value
    }else{
      pval <- 1
    }
    pval
  }))
  bpadj <- p.adjust(bpval, method=padjMethod) #sum(binarizedMat + noSdMat < 1))
  
  message(sprintf("\nFiltering by signficance %s of %s...", sum(bpadj < padj), sum(bpval < padj)))
  idxKeep <- which(bpadj < padj)
  bpadj <- bpadj[idxKeep]
  idxRow <- idxRow[idxKeep]
  binarizedMat <- binarizedMat[idxKeep, ]
  mat <- mat[idxKeep, ]
  intraMean <- intraMean[idxKeep, ]
  intraSd <- intraSd[idxKeep, ]
  
  #-----------------------------------------------
  #Determine which rows are above min group size and max group size to keep
  #-----------------------------------------------
  if (!is.null(maxGroupSize)) {
    idxKeep <- which(rowSums(binarizedMat) < (maxGroupSize + 1) & rowSums(binarizedMat) > (minGroupSize - 1))
  }else {
    idxKeep <- which(rowSums(binarizedMat) > (minGroupSize - 1))
  }
  message(sprintf("Filtering Features that are within Group Size %s of %s...", length(idxKeep), nrow(binarizedMat)))
  idxRow <- idxRow[idxKeep]
  bpadj <- bpadj[idxKeep]
  binarizedMat <- binarizedMat[idxKeep, ]
  mat <- mat[idxKeep, ]
  intraMean <- intraMean[idxKeep, ]
  intraSd <- intraSd[idxKeep, ]
  
  #-----------------------------------------------
  #Determine Pattern Occurences using data.table
  #-----------------------------------------------
  binarizedTable <- data.table::as.data.table(binarizedMat)[,.N,by=c(colnames(binarizedMat))]
  binarizedTable <- data.frame(binarizedTable[which(binarizedTable$N > groupMin),])[,which(colnames(binarizedTable) %ni% "N")]
  idxKeep <- unlist(lapply(seq_len(nrow(binarizedTable)), function(x){
    idx1  <- which(binarizedTable[x,,drop=TRUE] > 0)
    rs1   <- which(rowSums(binarizedMat[,idx1,drop=FALSE]) == length(idx1))
    rs2   <- which(rowSums(binarizedMat[,-idx1,drop=FALSE]) == 0)
    idxBM <- intersect(rs1,rs2)
    idxBM
  }))
  message(sprintf("Filtering Features Pattern Appearances %s of %s...", length(idxKeep), nrow(binarizedMat)))
  idxRow <- idxRow[idxKeep]
  bpadj <- bpadj[idxKeep]
  binarizedMat <- binarizedMat[idxKeep, ]
  mat <- mat[idxKeep, ]
  intraMean <- intraMean[idxKeep, ]
  intraSd <- intraSd[idxKeep, ]
  
  #-----------------------------------------------
  #Organize the output for maximum interpretation
  #-----------------------------------------------
  message(sprintf("Found %s unique elements!", length(idxKeep)))
  message("Finalizing Output...")
  colClust <- hclust(as.dist(1 - cor(intraMean)))
  colOrder <- unique(groups)[colClust$order]
  if(clusterCols){
    binarizedMat <- binarizedMat[, colClust$order]
  }
  idxOrdered <- do.call("order", c(as.data.frame(binarizedMat)[seq_len(ncol(binarizedMat))], list(decreasing = TRUE)))
  binarizedMat <- binarizedMat[idxOrdered, ]
  idxRow <- idxRow[idxOrdered]
  if(clusterCols){
    intraMean <- intraMean[idxOrdered, colClust$order]
    mat <- mat[idxOrdered, order(match(groups, colOrder))]
  } else {
    intraMean <- intraMean[idxOrdered,]
    mat <- mat[idxOrdered, ]
  }
  bpadj <- bpadj[idxOrdered]
  
  #-----------------------------------------------
  #Time to label each row
  #-----------------------------------------------
  binarizedTable <- data.frame(data.table::as.data.table(binarizedMat)[,.N,by=c(colnames(binarizedMat))])
  rownames(binarizedTable) <- paste0("Feature_",seq_len(nrow(binarizedTable)))
  dfUniuqe <- lapply(seq_len(nrow(binarizedTable)), function(x){
    idx1 <- which(binarizedTable[x,-ncol(binarizedTable),drop=TRUE] > 0)
    rs1 <- which(rowSums(binarizedMat[,idx1,drop=FALSE]) == length(idx1))
    rs2 <- which(rowSums(binarizedMat[,-idx1,drop=FALSE])==0)
    idxBM <- intersect(rs1,rs2)
    data.frame(feature = rownames(binarizedTable)[x], rows = idxBM)
  }) %>% Reduce("rbind",.)
  
  #-----------------------------------------------
  #Return Output
  #-----------------------------------------------
  out <- list(mat = mat, binaryMat = binarizedMat, groupMat = intraMean, binarizedTable = binarizedTable, dfFeature = dfUniuqe, rowOrder = idxRow, padj = bpadj)
  return(out)
  
}  

'%ni%' <- Negate('%in%')

normalize_expr_data <- function(cds,
                                norm_method = c("log", "size_only", "none"),
                                pseudo_count = NULL) {
  norm_method <- match.arg(norm_method)
  
  FM <- SingleCellExperiment::counts(cds)
  
  # If we're going to be using log, and the user hasn't given us a
  # pseudocount set it to 1 by default.
  if (is.null(pseudo_count)){
    if(norm_method == "log")
      pseudo_count <- 1
    else
      pseudo_count <- 0
  }
  
  if (norm_method == "log") {
    # If we are using log, normalize by size factor before log-transforming
    
    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
    
    if (pseudo_count != 1 || is_sparse_matrix(SingleCellExperiment::counts(cds)) == FALSE){
      FM <- FM + pseudo_count
      FM <- log2(FM)
    } else {
      FM@x = log2(FM@x + 1)
    }
    
  } else if (norm_method == "size_only") {
    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
    FM <- FM + pseudo_count
  }
  return (FM)
}

#Nearest Neighbor differential
findNN <- function(query, reference, method = "euclidean"){
  findClosest <- function(x, m, method = "euclidean"){
    if(method=="euclidean"){
      which.min(sqrt(colSums((t(m) - x) * (t(m) - x))))
    }else if(method=="pearson"){
      which.max(cor(t(m),x,method = method)[,1])
    }else if(method=="spearman"){
      which.max(cor(t(m),x,method = method)[,1])
    }
  }
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  mat <- data.frame(matrix(ncol = 4, nrow = nrow(query)))
  colnames(mat) <- c("x", "i", "y", "j")
  for(i in seq_len(nrow(query))){
    setTxtProgressBar(pb,round(i*100/nrow(query),0))
    j <- findClosest(query[i,], reference, method)
    mat[i,] <- c(x = rownames(query)[i], i = i, y = rownames(reference)[j], j = j)
  }
  return(mat)
}

sourceCpp(code='
          #include <Rcpp.h>
          
          using namespace Rcpp;
          using namespace std;
          
          // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
          // [[Rcpp::export]]
          Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
          
          if(X.ncol() != Y.ncol()){
          stop("Columns of Matrix X and Y must be equal length!");
          }
          
          if(max(idxX) > X.nrow()){
          stop("Idx X greater than nrow of Matrix X");
          }
          
          if(max(idxY) > Y.nrow()){
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


plot_genes_in_pseudotime2 <-function(cds_subset,
                                     trend_formula="~ splines::ns(pseudotime, df=4)",
                                     label_by_short_name=TRUE){
  
  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  
  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[,is.finite(colData(cds_subset)$pseudotime)]
  
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  #cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- inner_join(as_tibble(cds_exprs),as_tibble(cds_rowData),by=c('f_id'='gene_short_name'))
  #cds_exprs <- merge(as.data.frame(cds_exprs), cds_colData, by.x = "Cell", by.y = "barcode")
  cds_exprs <- inner_join(cds_exprs,as_tibble(cds_colData),by=c("Cell"="barcode"))
  cds_exprs$adjusted_expression <- cds_exprs$expression
  
  cds_exprs$feature_label <- cds_exprs$f_id

  
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  
  new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula,cores=1)
  
  model_expectation <- model_predictions(model_tbl,
                                         new_data = colData(cds_subset))
  
  colnames(model_expectation) <- colnames(cds_subset)
  model_expectation <-as.data.frame(t(model_expectation))
  model_expectation <- as.data.frame(scale(model_expectation,center = F,scale=colMaxs(as.matrix(model_expectation))))
  model_expectation$pseudotime <- colData(cds_subset)$pseudotime
  model_expectation <- model_expectation[order(model_expectation$pseudotime),]
  #model_expectation <- model_expectation[,-c(which(colSums(model_expectation)>=(as.numeric(nrow(model_expectation))-10)))]
  #idx_order <- colwise(function(x){which(x==1)})(model_expectation[,-ncol(model_expectation)])
  model_expectation <- model_expectation[,-ncol(model_expectation)]
  #model_expectation <- model_expectation[,order(idx_order)]
  return(t(model_expectation))
}

bin.matrix.rows <- function(m, bin.size) {
  splitv <- rep(1:round((ncol(m)/bin.size)), each=bin.size)
  splitv <- append(splitv, rep(NA,ncol(m)-length(splitv)))
  bmat <- t(apply(m,1, function(x) {unlist(lapply(split(x,splitv), mean))}))
  cnames <- round(unlist(lapply(split(as.integer(colnames(m)), splitv),median)))
  colnames(bmat) <- cnames
  bmat
}

FindMotifs_fisher <- function(
  object,
  features,
  background = 40000,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- DefaultAssay(object = object)
  if (is(object = background, class2 = 'numeric')) {
    if (verbose) {
      message("Selecting background regions to match input sequence characteristics")
    }
    background <- MatchRegionStats(
      meta.feature = GetAssayData(object = object, assay = assay, slot = 'meta.features'),
      regions = features,
      n = background,
      verbose = verbose,
      ...
    )
  }
  if (verbose) {
    message('Testing motif enrichment in ', length(x = features), ' regions')
  }
  motif.all <- GetMotifData(object = object, assay = assay, slot = 'data')
  motif.names <- GetMotifData(object = object, assay = assay, slot = 'motif.names')
  query.motifs <- motif.all[features, ]
  background.motifs <- motif.all[background, ]
  query.counts <- colSums(x = query.motifs)
  background.counts <- colSums(x = background.motifs)
  percent.observed <- query.counts / length(x = features) * 100
  percent.background <- background.counts / length(x = background) * 100
  fold.enrichment <- percent.observed / percent.background
  p.list <- c()
  for (i in seq_along(along.with = query.counts)) {
    dat <- data.frame(Background=c(background.counts[[i]],nrow(background.motifs) - background.counts[[i]]),Query=c(query.counts[[i]],length(features)-query.counts[[i]]))
    p.list[[i]] <- fisher.test(dat)$p.value
  }
  results <- data.frame(
    motif = names(x = query.counts),
    observed = query.counts,
    background = background.counts,
    percent.observed = percent.observed,
    percent.background = percent.background,
    fold.enrichment = fold.enrichment,
    pvalue = p.list,
    motif.name = as.vector(x = unlist(x = motif.names[names(x = query.counts)])),
    stringsAsFactors = FALSE
  )
  if (nrow(x = results) == 0) {
    return(results)
  } else {
    return(results[with(data = results, expr = order(pvalue, -fold.enrichment)), ])
  }
}

features_in_pseudotime <-function(cds_subset,pd,
                                   trend_formula="~ splines::ns(pseudotime, df=4)",
                                   label_by_short_name=TRUE){
  
  colData(cds_subset)$pseudotime <- pd
  
  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[,is.finite(colData(cds_subset)$pseudotime)]
  
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  #cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- inner_join(as_tibble(cds_exprs),as_tibble(cds_rowData),by=c('f_id'='gene_short_name'))
  #cds_exprs <- merge(as.data.frame(cds_exprs), cds_colData, by.x = "Cell", by.y = "barcode")
  cds_exprs <- inner_join(cds_exprs,as_tibble(cds_colData),by=c("Cell"="barcode"))
  cds_exprs$adjusted_expression <- cds_exprs$expression
  
  cds_exprs$feature_label <- cds_exprs$f_id
  
  
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  
  new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula,cores=1)
  
  model_expectation <- model_predictions(model_tbl,
                                         new_data = colData(cds_subset))
  
  colnames(model_expectation) <- colnames(cds_subset)
  model_expectation <-as.data.frame(t(model_expectation))
  model_expectation <- as.data.frame(scale(model_expectation,center = F,scale=colMaxs(as.matrix(model_expectation))))
  model_expectation$pseudotime <- colData(cds_subset)$pseudotime
  model_expectation <- model_expectation[order(model_expectation$pseudotime),]
  #model_expectation <- model_expectation[,-c(which(colSums(model_expectation)>=(as.numeric(nrow(model_expectation))-10)))]
  #idx_order <- colwise(function(x){which(x==1)})(model_expectation[,-ncol(model_expectation)])
  model_expectation <- model_expectation[,-ncol(model_expectation)]
  #model_expectation <- model_expectation[,order(idx_order)]
  return(t(model_expectation))
}

generate_linkedGIntervals <- function(mat,clusters=c('NSC','IPC','PN'),mat_name='posCor'){
  mat$PN1 <- 0
  for (select_group in clusters){
    indx <- colnames(mat)[grep(select_group,colnames(mat[,1:7]))]
    subset_peaks <- mat[rowSums(mat[,indx,drop=F])>=1,]
    
    name <- paste0('results/beds/',select_group,'_',mat_name,'_GA.bed')
    bed_f <- DataFrame(subset_peaks[,c('gene_chr','gene_start','gene_start')])
    colnames(bed_f) <- c('chrom','start','end')
    bed_f[subset_peaks$gene_strand=='+','end'] <- bed_f$end[subset_peaks$gene_strand=='+']+1
    bed_f[subset_peaks$gene_strand=='-','start'] <- bed_f$start[subset_peaks$gene_strand=='-']-1
    bed_f1 <- makeGRangesFromDataFrame(bed_f)
    rtracklayer::export.bed(bed_f1,con = name)
    
    bed_f2 <- separate(subset_peaks[,'peakName',drop=F],col='peakName',into=c('chrom','start','end'),sep='_')
    bed_f <- as.data.frame(bed_f)
    g_beds <- gintervals.2d(chroms1 = bed_f[,1],starts1 = bed_f[,2],ends1 = bed_f[,3],chroms2 = bed_f2[,1],starts2 = (as.numeric(bed_f2[,2])+as.numeric(bed_f2[,3]))/2,ends2 = (as.numeric(bed_f2[,2])+as.numeric(bed_f2[,3]))/2+1)
    g_beds_DA <- g_beds[,4:6]
    g_beds_GA <- g_beds[,1:3]
    colnames(g_beds_DA) <-  c('chrom','start','end')
    colnames(g_beds_GA) <-  c('chrom','start','end')
    name_GA <- gsub('\\.bed','',name)
    save(g_beds_GA,file = name_GA)
    
    name <- paste0('results/beds/',select_group,'_',mat_name,'_DA.bed')
    bed_f2 <- makeGRangesFromDataFrame(bed_f2)
    rtracklayer::export.bed(bed_f2,con = name)
    name_DA <- gsub('\\.bed','',name)
    save(g_beds_DA,file = name_DA)  
    
    name <- paste0('plots/HiC/',select_group,'_',mat_name,'.pdf')
    pdf(name,width=8,height=8)
    plot(density(subset_peaks$NSCscore,na.rm=T),col='red',ylim=c(0,0.014))
    lines(density(subset_peaks$IPCscore,na.rm=T),col='green')
    lines(density(subset_peaks$PNscore,na.rm=T),col='blue')  
    dev.off()
    
    name <- paste0('plots/HiC/',select_group,'_',mat_name,'_meth.pdf')
    pdf(name,width=8,height=8)
    plot(density(subset_peaks$distal.E14_NSC_10x,na.rm=T),col='red')
    lines(density(subset_peaks$distal.E14_IPC_10x,na.rm=T),col='green')
    lines(density(subset_peaks$distal.E14_PN_10x,na.rm=T),col='blue')  
    dev.off()
  }
}

rankMotifs_HiCscore <- function(df=p2glinks$posCor,anno_mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T),chip_matrix=NULL,out_f,domain='intraTAD',include_prom=T,padj=1,minLFC=0,minSdRatio=0,zCutoff=0){
  if (!is.null(chip_matrix)){
    motif.matrix <- readRDS(chip_matrix)
    out_f <- gsub('motif','chip',out_f)
  } else {
    motif.names <- unlist(GetMotifObject(atac.object)@motif.names)
    motif.names <-  capitalize(tolower(motif.names))
    motif.matrix <- GetMotifData(atac.object)
    colnames(motif.matrix) <- as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])
  }
  if(!is.null(domain)){
    df <- df[df$domain==domain]
    out_f <- gsub('\\.tsv',paste0('_',domain,'.tsv'),out_f)
  }
  if(include_prom){
    df_prom <- df
    df_prom$Correlation <- 1
    df_prom$distance <- 0
    ranges(df_prom) <- IRanges(start=df$gene_start-2000,end=df$gene_start+200)
    df_prom <- df_prom[!duplicated(df_prom$gene_name)]
    all_features <- StringToGRanges(row.names(atac.object), sep = c("-", "-"))
    prom_overlaps <- findOverlaps(df_prom,all_features)
    #   prom_overlaps <- prom_overlaps[!is.na(prom_overlaps@to)]
    df_prom <- df_prom[prom_overlaps@from]
    ranges(df_prom) <- ranges(all_features[prom_overlaps@to])
    df_prom$peakName <- paste0(seqnames(df_prom),'_',start(df_prom),'_',end(df_prom))
  }
  
  df$label <- paste0(df$gene_name,':',df$distance)
  df_prom$label <- paste0(df_prom$gene_name,':',df_prom$distance)
  anno_mat$label <- paste0(anno_mat$gene_name,':',anno_mat$distance)
  anno_mat <- anno_mat[match(df$label,anno_mat$label),]
  prom_anno_mat <- anno_mat[match(df_prom$gene_name,anno_mat$gene_name),]
  
  df$cluster <- colnames(anno_mat[,1:7])[max.col(anno_mat[,1:7])]
  df_prom$cluster <- colnames(prom_anno_mat[,1:7])[max.col(prom_anno_mat[,1:7])]
  prom.motif.matrix <- motif.matrix[match(gsub('_','-',df_prom$peakName),row.names(motif.matrix)),]
  motif.matrix <- motif.matrix[match(gsub('_','-',df$peakName),row.names(motif.matrix)),]
  sel_df <- cbind(mcols(df),motif.matrix)
  prom_sel_df <- cbind(mcols(df_prom),prom.motif.matrix)
  
  #### Create permutation table #####
  perm_df <- data.frame(cluster=NA,NSCmean=NA,IPCmean=NA,PNmean=NA,FC_IPCvsNSC=NA,FC_PNvsNSC=NA,FC_PNvsIPC=NA)
  for (s in unique(df$cluster)){
    x <- sel_df[sel_df$cluster==s,]
    for (i in 1:1000){
      subset_df <- x[sample(row.names(x),100),c('NSCscore','IPCscore','PNscore')]
      score_mat <- colMeans(as.matrix(subset_df),na.rm=T)
      diff_mat <- colMeans(matrix(c(subset_df$IPCscore-subset_df$NSCscore,subset_df$PNscore-subset_df$NSCscore,subset_df$PNscore-subset_df$IPCscore),byrow = F,nrow=nrow(subset_df),ncol=3),na.rm=T)
      perm_df <- rbind(perm_df,c(s,score_mat,diff_mat))
    }
  }
  perm_df <- perm_df[-1,]
  perm_df[,2:7] <- mapply(perm_df[,2:7],FUN=as.numeric)
  ####################################
  
  
  results_df <- data.frame(motif=NA,cluster=NA,binding=NA,NSCmean=NA,IPCmean=NA,PNmean=NA,FC_IPCvsNSC=NA,FC_PNvsNSC=NA,FC_PNvsIPC=NA,PValue_NSCmean=NA,PValue_IPCmean=NA,PValue_PNmean=NA,PValue_FC_IPCvsNSC=NA,PValue_FC_PNvsNSC=NA,PValue_FC_PNvsIPC=NA)
  prom_motifs <- ddply(as.data.frame(prom_sel_df),.(gene_name),function(x){return(colSums(x[,(ncol(mcols(df))+1):ncol(x)],na.rm=T))})
  results_list <- list()
  for (s in unique(df$cluster)){
    x <- sel_df[sel_df$cluster==s,]
    x_perm <- perm_df[perm_df$cluster==s,-1]
    for (i in (ncol(mcols(df))+1):ncol(sel_df)){
      target_TF <- colnames(x)[i]
      subset_df <- x[x[,i]==1,c('gene_name','NSCscore','IPCscore','PNscore')]
      subset_df <- subset_df[complete.cases(subset_df),]
      if(nrow(subset_df)<3){next}
      subset_df$gene_motif <- prom_motifs[match(subset_df$gene_name,prom_motifs$gene_name),gsub('\\(|\\:|\\)','.',target_TF)]
      if(sum(is.na(subset_df$gene_motif))>0){
        subset_df$gene_motif[is.na(subset_df$gene_motif)] <- 0
      }
      res_single <- c(target_TF,x$cluster[1],'single',colMeans(as.matrix(subset_df[subset_df$gene_motif==0,2:4]),na.rm=T))
      res_double <- c(target_TF,x$cluster[1],'double',colMeans(as.matrix(subset_df[subset_df$gene_motif!=0,2:4]),na.rm=T))
      diff_mat_single <- colMeans(matrix(c(subset_df[subset_df$gene_motif==0,'IPCscore']-subset_df[subset_df$gene_motif==0,'NSCscore'],subset_df[subset_df$gene_motif==0,'PNscore']-subset_df[subset_df$gene_motif==0,'NSCscore'],subset_df[subset_df$gene_motif==0,'PNscore']-subset_df[subset_df$gene_motif==0,'IPCscore']),byrow = F,nrow=nrow(subset_df[subset_df$gene_motif==0,]),ncol=3),na.rm=T)
      diff_mat_double <- colMeans(matrix(c(subset_df[subset_df$gene_motif!=0,'IPCscore']-subset_df[subset_df$gene_motif!=0,'NSCscore'],subset_df[subset_df$gene_motif!=0,'PNscore']-subset_df[subset_df$gene_motif!=0,'NSCscore'],subset_df[subset_df$gene_motif!=0,'PNscore']-subset_df[subset_df$gene_motif!=0,'IPCscore']),byrow = F,nrow=nrow(subset_df[subset_df$gene_motif!=0,]),ncol=3),na.rm=T)
      pvalue_single <- c()
      pvalue_double <- c()
      for (k in 1:6){
        pvalue_single <- c(pvalue_single,sum(x_perm[,k]>=as.numeric(c(res_single[4:6],diff_mat_single)[k]))/1000)
        pvalue_double <- c(pvalue_double,sum(x_perm[,k]>=as.numeric(c(res_double[4:6],diff_mat_double)[k]))/1000)
      }
      results_df <- rbind(results_df,c(res_single,diff_mat_single,pvalue_single),c(res_double,diff_mat_double,pvalue_double))
      results_list[[paste0(s,'_',target_TF,'_single')]] <- as.matrix(subset_df[subset_df$gene_motif==0,2:4])
      results_list[[paste0(s,'_',target_TF,'_double')]] <- as.matrix(subset_df[subset_df$gene_motif!=0,2:4])
    }
  }
  results_df <- results_df[-1,]
  results_df[,4:ncol(results_df)] <- mapply(results_df[,4:ncol(results_df)],FUN=as.numeric)
  
  results_df <- results_df[complete.cases(results_df),]
  
  uf_rna <- suppressMessages(uniqueFeatures(
    edgeR::cpm(assay(seRNA_all),log=TRUE,prior.count=1),
    groups = colData(seRNA_all)$Group,
    padj = padj,
    minSdRatio = minSdRatio,
    minLFC = minLFC,
    zCutoff = zCutoff,
    breakPt = "last",
    groupMin = 0,
    maxGroupSize = 2,
    clusterCols = F
  ))
  
  results_df$gene <- sapply(strsplit(as.character(results_df$motif), "\\(|\\:|\\."), function(x) x[[1]])
  rna_clusters <- uf_rna$binaryMat
  rna_clusters <- rna_clusters[match(results_df$gene,row.names(rna_clusters)),]
  results_df$gene_cluster <- gsub('_M|1|2','',colnames(rna_clusters)[max.col(rna_clusters)])
  
  #results_df$cluster <- gsub('|1|2','',results_df$cluster)
  write.table(results_df,out_f,quote=F,col.names=T,row.names=F,sep='\t')
  write.table(perm_df,gsub('\\.tsv','_perm.tsv',out_f),quote=F,col.names=T,row.names=F,sep='\t')
  saveRDS(results_list,file = gsub('\\.tsv','.RDS',out_f))
}

rankMotifs_meth <- function(df=p2glinks$posCor,anno_mat=read.table('results/P2G_binaryMat_posCor.tsv',header=T),chip_matrix=NULL,out_f,domain='intraTAD',include_prom=F,padj=1,minLFC=0,minSdRatio=0,zCutoff=0){
  if (!is.null(chip_matrix)){
    motif.matrix <- readRDS(chip_matrix)
    out_f <- gsub('motif','chip',out_f)
  } else {
    motif.names <- unlist(GetMotifObject(atac.object)@motif.names)
    motif.names <-  capitalize(tolower(motif.names))
    motif.matrix <- GetMotifData(atac.object)
    colnames(motif.matrix) <- as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])
  }
  if(!is.null(domain)){
    df <- df[df$domain==domain]
    out_f <- gsub('\\.tsv',paste0('_',domain,'.tsv'),out_f)
  }
  if(include_prom){
    df_prom <- df
    df_prom$Correlation <- 1
    df_prom$distance <- 0
    ranges(df_prom) <- IRanges(start=df$gene_start-2000,end=df$gene_start+200)
    df_prom <- df_prom[!duplicated(df_prom$gene_name)]
    all_features <- StringToGRanges(row.names(atac.object), sep = c("-", "-"))
    prom_overlaps <- findOverlaps(df_prom,all_features)
    #   prom_overlaps <- prom_overlaps[!is.na(prom_overlaps@to)]
    df_prom <- df_prom[prom_overlaps@from]
    ranges(df_prom) <- ranges(all_features[prom_overlaps@to])
    df_prom$peakName <- paste0(seqnames(df_prom),'_',start(df_prom),'_',end(df_prom))
  }
  
  df$label <- paste0(df$gene_name,':',df$distance)
  df_prom$label <- paste0(df_prom$gene_name,':',df_prom$distance)
  anno_mat$label <- paste0(anno_mat$gene_name,':',anno_mat$distance)
  anno_mat <- anno_mat[match(df$label,anno_mat$label),]
  prom_anno_mat <- anno_mat[match(df_prom$gene_name,anno_mat$gene_name),]
  
  df$cluster <- colnames(anno_mat[,1:7])[max.col(anno_mat[,1:7])]
  df_prom$cluster <- colnames(prom_anno_mat[,1:7])[max.col(prom_anno_mat[,1:7])]
  prom.motif.matrix <- motif.matrix[match(gsub('_','-',df_prom$peakName),row.names(motif.matrix)),]
  motif.matrix <- motif.matrix[match(gsub('_','-',df$peakName),row.names(motif.matrix)),]
  sel_df <- cbind(mcols(df),motif.matrix)
  prom_sel_df <- cbind(mcols(df_prom),prom.motif.matrix)
  
  #### Create permutation table #####
  perm_df <- data.frame(cluster=NA,NSCmean=NA,IPCmean=NA,PNmean=NA,FC_IPCvsNSC=NA,FC_PNvsNSC=NA,FC_PNvsIPC=NA)
  for (s in unique(df$cluster)){
    x <- sel_df[sel_df$cluster==s,]
    for (i in 1:1000){
      subset_df <- x[sample(row.names(x),100),c('distal.E14_NSC_10x','distal.E14_IPC_10x','distal.E14_PN_10x')]
      colnames(subset_df) <- c('NSCscore','IPCscore','PNscore')
      score_mat <- colMeans(as.matrix(subset_df),na.rm=T)
      diff_mat <- colMeans(matrix(c(subset_df$IPCscore-subset_df$NSCscore,subset_df$PNscore-subset_df$NSCscore,subset_df$PNscore-subset_df$IPCscore),byrow = F,nrow=nrow(subset_df),ncol=3),na.rm=T)
      perm_df <- rbind(perm_df,c(s,score_mat,diff_mat))
    }
  }
  perm_df <- perm_df[-1,]
  perm_df[,2:7] <- mapply(perm_df[,2:7],FUN=as.numeric)
  ####################################
  
  
  results_df <- data.frame(motif=NA,cluster=NA,binding=NA,NSCmean=NA,IPCmean=NA,PNmean=NA,FC_IPCvsNSC=NA,FC_PNvsNSC=NA,FC_PNvsIPC=NA,PValue_NSCmean=NA,PValue_IPCmean=NA,PValue_PNmean=NA,PValue_FC_IPCvsNSC=NA,PValue_FC_PNvsNSC=NA,PValue_FC_PNvsIPC=NA)
  prom_motifs <- ddply(as.data.frame(prom_sel_df),.(gene_name),function(x){return(colSums(x[,(ncol(mcols(df))+1):ncol(x)],na.rm=T))})
  results_list <- list()
  for (s in unique(df$cluster)){
    x <- sel_df[sel_df$cluster==s,]
    x_perm <- perm_df[perm_df$cluster==s,-1]
    for (i in (ncol(mcols(df))+1):ncol(sel_df)){
      target_TF <- colnames(x)[i]
      subset_df <- x[x[,i]==1,c('gene_name','distal.E14_NSC_10x','distal.E14_IPC_10x','distal.E14_PN_10x')]
      colnames(subset_df) <- c('gene_name','NSCscore','IPCscore','PNscore')
      subset_df <- subset_df[complete.cases(subset_df),]
      if(nrow(subset_df)<3){next}
      subset_df$gene_motif <- prom_motifs[match(subset_df$gene_name,prom_motifs$gene_name),gsub('\\(|\\:|\\)','.',target_TF)]
      if(sum(is.na(subset_df$gene_motif))>0){
        subset_df$gene_motif[is.na(subset_df$gene_motif)] <- 0
      }
      res_single <- c(target_TF,x$cluster[1],'single',colMeans(as.matrix(subset_df[subset_df$gene_motif==0,2:4]),na.rm=T))
      res_double <- c(target_TF,x$cluster[1],'double',colMeans(as.matrix(subset_df[subset_df$gene_motif!=0,2:4]),na.rm=T))
      diff_mat_single <- colMeans(matrix(c(subset_df[subset_df$gene_motif==0,'IPCscore']-subset_df[subset_df$gene_motif==0,'NSCscore'],subset_df[subset_df$gene_motif==0,'PNscore']-subset_df[subset_df$gene_motif==0,'NSCscore'],subset_df[subset_df$gene_motif==0,'PNscore']-subset_df[subset_df$gene_motif==0,'IPCscore']),byrow = F,nrow=nrow(subset_df[subset_df$gene_motif==0,]),ncol=3),na.rm=T)
      diff_mat_double <- colMeans(matrix(c(subset_df[subset_df$gene_motif!=0,'IPCscore']-subset_df[subset_df$gene_motif!=0,'NSCscore'],subset_df[subset_df$gene_motif!=0,'PNscore']-subset_df[subset_df$gene_motif!=0,'NSCscore'],subset_df[subset_df$gene_motif!=0,'PNscore']-subset_df[subset_df$gene_motif!=0,'IPCscore']),byrow = F,nrow=nrow(subset_df[subset_df$gene_motif!=0,]),ncol=3),na.rm=T)
      pvalue_single <- c()
      pvalue_double <- c()
      for (k in 1:6){
        pvalue_single <- c(pvalue_single,sum(x_perm[,k]>=as.numeric(c(res_single[4:6],diff_mat_single)[k]))/1000)
        pvalue_double <- c(pvalue_double,sum(x_perm[,k]>=as.numeric(c(res_double[4:6],diff_mat_double)[k]))/1000)
      }
      results_df <- rbind(results_df,c(res_single,diff_mat_single,pvalue_single),c(res_double,diff_mat_double,pvalue_double))
      results_list[[paste0(s,'_',target_TF,'_single')]] <- as.matrix(subset_df[subset_df$gene_motif==0,2:4])
      results_list[[paste0(s,'_',target_TF,'_double')]] <- as.matrix(subset_df[subset_df$gene_motif!=0,2:4])
    }
  }
  results_df <- results_df[-1,]
  results_df[,4:ncol(results_df)] <- mapply(results_df[,4:ncol(results_df)],FUN=as.numeric)
  
  results_df <- results_df[complete.cases(results_df),]
  
  uf_rna <- suppressMessages(uniqueFeatures(
    edgeR::cpm(assay(seRNA_all),log=TRUE,prior.count=1),
    groups = colData(seRNA_all)$Group,
    padj = padj,
    minSdRatio = minSdRatio,
    minLFC = minLFC,
    zCutoff = zCutoff,
    breakPt = "last",
    groupMin = 0,
    maxGroupSize = 2,
    clusterCols = F
  ))
  
  results_df$gene <- sapply(strsplit(as.character(results_df$motif), "\\(|\\:|\\."), function(x) x[[1]])
  rna_clusters <- uf_rna$binaryMat
  rna_clusters <- rna_clusters[match(results_df$gene,row.names(rna_clusters)),]
  results_df$gene_cluster <- gsub('_M|1|2','',colnames(rna_clusters)[max.col(rna_clusters)])
  
  #results_df$cluster <- gsub('|1|2','',results_df$cluster)
  write.table(results_df,out_f,quote=F,col.names=T,row.names=F,sep='\t')
  write.table(perm_df,gsub('\\.tsv','_perm.tsv',out_f),quote=F,col.names=T,row.names=F,sep='\t')
  saveRDS(results_list,file = gsub('\\.tsv','.RDS',out_f))
}



plot_sankey <- function(true_labels, prediction, links = NULL, custom_node_color = NULL, custom_link_color = NULL){
  #------------------------------------------------------------------------------------------
  # function to create Sankey plot (via d3 network) to visualise true cell type vs prediction
  #------------------------------------------------------------------------------------------
  # input: true_labels = vector of celltype labels of the projection dataset
  #        prediction = vector of predicted labels for the projection dataset
  #        links = optional vector how to color the links
  #        custom_node_color = dataframe with 2 columns: node --> levels of labels (e.g. "n1"/ "unassigned")
  #                                                      color --> corresponding color for each label
  #        custom_link_color = dataframe with 2 columns: link --> levels of 'links'-vector
  #                                                      color --> corresponding color for each label
  # output: sankey plot (d3 network)
  library(plyr)
  library(networkD3)
  df <- data.frame(celltype = true_labels,
                   prediction = prediction,
                   value = rep(1,length(true_labels)))
  nodes <- as.data.frame(cbind(c(levels(df$celltype), levels(df$prediction))))
  
  
  if(!is.null(links)){
    df$res <- links
  }
  
  # convert group names into ascending numbers starting from 0 (required for the sankeyNetwork function)
  l <- length(levels(df$celltype))
  levels(df$celltype) <- c(0:(l-1))
  levels(df$prediction) <- c(l : (l+length(levels(df$prediction))-1))
  df$celltype <- as.numeric(as.character(df$celltype))
  df$prediction <- as.numeric(as.character(df$prediction))
  
  
  if(!is.null(custom_link_color)){
    # filter for link labels that are actually present
    link_filter <- dplyr::filter(custom_link_color, link %in% levels(as.factor(links)))  
  }
  
  
  if(!is.null(custom_node_color) | !is.null(custom_link_color)){
    
    
    my_color <- paste0('d3.scaleOrdinal().domain([',
                       paste(shQuote(custom_node_color$node), collapse=", "),',',
                       paste(shQuote(link_filter$link), collapse=", "),
                       ']).range([',
                       paste(shQuote(custom_node_color$color), collapse=", "),',',
                       paste(shQuote(link_filter$color), collapse=", "),
                       '])')
  } else {
    my_color <- networkD3::JS("d3.scaleOrdinal(d3.schemeCategory20);")
  }
  
  networkD3::sankeyNetwork(Links = df, Nodes = nodes, Source = 'celltype', Target = 'prediction',
                           LinkGroup = "res", colourScale = my_color,
                           Value = 'value', NodeID = 'V1',
                           units = '', fontSize = 24, nodeWidth = 30,height=800,width=1000)
  
  
}

save_sankey <- function(file_f,out_f, width = 800, height = 800){
  require(networkD3)
  require(webshot)
  webshot(file_f, file = out_f,
          vwidth = width, vheight = height)
}


GetFragments <- function(
  object,
  assay = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  tools <- slot(object = object, name = 'tools')
  if ('fragments' %in% names(x = tools)) {
    if (assay %in% names(x = tools$fragments)) {
      fragment.path <- tools$fragments[[assay]]
    } else {
      stop('Fragment file not supplied for the requested assay')
    }
  } else {
    stop('Fragment file not set.
         Run SetFragments to set the fragment file path.')
  }
  if (!(all(file.exists(fragment.path, paste0(fragment.path, '.tbi'))))) {
    stop('Requested file does not exist or is not indexed')
  } else {
    return(fragment.path)
  }
  }

MultiRegionCutMatrix <- function(
  object,
  regions,
  assay = NULL,
  cells = NULL,
  verbose = FALSE
) {
  fragment.path <- GetFragments(object = object, assay = assay)
  tabix.file <- TabixFile(file = fragment.path)
  open(con = tabix.file)
  cm.list <- lapply(
    X = seq_along(along.with = regions),
    FUN = function(x) {
      CutMatrix(
        object = object,
        assay = assay,
        tabix.file = tabix.file,
        region = regions[x, ],
        verbose = verbose
      )
    }
  )
  cm <- Reduce(f = `+`, x = cm.list)
  close(con = tabix.file)
  return(cm)
}


CreateRegionPileupMatrix <- function(
  object,
  regions,
  upstream = 1000,
  downstream = 1000,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  # extend upstream and downstream from midpoint
  regions <- Extend(
    x = regions,
    upstream = upstream,
    downstream = downstream,
    from.midpoint = TRUE
  )
  # split into strands
  on_plus <- strand(x = regions) == "+" | strand(x = regions) == "*"
  plus.strand <- regions[on_plus, ]
  minus.strand <- regions[!on_plus, ]
  
  # get cut matrices for each strand
  if (verbose) {
    message("Finding + strand cut sites")
  }
  cut.matrix.plus <- MultiRegionCutMatrix(
    regions = plus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  if (verbose) {
    message("Finding - strand cut sites")
  }
  cut.matrix.minus <- MultiRegionCutMatrix(
    regions = minus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  
  # reverse minus strand and add together
  full.matrix <- cut.matrix.plus + cut.matrix.minus[, rev(
    x = colnames(x = cut.matrix.minus)
  )]
  colnames(full.matrix) <- -upstream:downstream
  return(full.matrix)
}


GetGroups <- function(
  object,
  group.by,
  idents
) {
  if (is.null(x = group.by)) {
    obj.groups <- Idents(object = object)
  } else {
    obj.md <- object[[group.by]]
    obj.groups <- obj.md[, 1]
    names(obj.groups) <- rownames(x = obj.md)
  }
  if (!is.null(idents)) {
    obj.groups <- obj.groups[obj.groups %in% idents]
  }
  return(obj.groups)
}



ApplyMatrixByGroup <- function(
  mat,
  groups,
  fun,
  normalize = TRUE,
  group.scale.factors = NULL,
  scale.factor = NULL
) {
  if (normalize) {
    if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
      stop("If normalizing counts, supply group scale factors")
    }
  }
  results <- list()
  all.groups <- unique(x = groups)
  for (i in seq_along(along.with = all.groups)) {
    pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    results[[i]] <- data.frame(
      group = all.groups[[i]],
      count = totals,
      position = as.numeric(colnames(x = mat)),
      stringsAsFactors = FALSE
    )
  }
  coverages <- as.data.frame(
    x = do.call(what = rbind, args = results), stringsAsFactors = FALSE
  )
  if (normalize) {
    scale.factor <- SetIfNull(
      x = scale.factor, y = median(x = group.scale.factors)
    )
    coverages$norm.value <- coverages$count /
      group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
}


SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

AddToMisc <- function(
  object,
  new.data,
  save.as,
  assay = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  misc.slot <- SetIfNull(x = Misc(object = object[[assay]]), y = list())
  if (!inherits(x = misc.slot, what = 'list')) {
    warning("Misc slot already occupied")
  } else{
    misc.slot[[save.as]] <- new.data
    object[[assay]]@misc <- misc.slot
  }
  return(object)
}

enrichGO_wrapper <- function(genes,qvalue.cutoff=0.05,orgDB=org.Mm.eg.db,organism='mmu'){
  gene.df <- bitr(genes, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = orgDB)
  ego <- enrichGO(gene         = unique(gene.df$ENTREZID),
                  OrgDb         = orgDB,
                  keyType       = 'ENTREZID',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = qvalue.cutoff)
  p <- simplify(ego)
  return(p)
}

extract_MotifPos <- function(df,pwm,chip_f=NULL,motif.name='Neurog2(var.2)'){
  motif.names <- capitalize(tolower(as.vector(name(pwm))))
  if(!is.null(chip_f)){
    chip <- read.table(chip_f)
    colnames(chip) <- c('chrom','start','end')
    df <- makeGRangesFromDataFrame(chip)
  }
  motif_pos <- unlist(matchMotifs(pwm[[which(motif.names==motif.name)]], df, genome = "mm10", out = "positions") )
  motif_6 <- data.frame(seqnames=seqnames(motif_pos),starts=start(motif_pos),ends=end(motif_pos),names=c(rep(".", length(motif_pos))),scores=elementMetadata(motif_pos[,1]),strands=strand(motif_pos))
  return(motif_6)
}
