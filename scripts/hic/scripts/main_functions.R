##########################################################################
#'  Classify domains based on trans contacts into k number of compartments
#'
#' \code{analyzeCompartments}
#'
#' This function calculates the trans contacts between domains and clusters them into k number of compartments as described by Nagano et al., Cell 2017. Output file is an interval track with each domain associated with its compartment identity.
#'
#' @param cells Which condition to work on. Vector of conditions is also accepted.
#' @param domains Either pre-computed domains in gintervals format or NULL if domains should be computed based on the insulation score.
#' @param hic_tracks Which hic tracks to use. Default to all tracks set in the config. Only tracks containing the condition(s) name will be retained.
#' @param ins_tracks Which insulation tracks to use. Default to insulation tracks set in the config. Only tracks containing the condition(s) name will be retained.
#' @param ins_dom_thresh Threshold of the regions with highest insulation to consider as domain boundaries (e.g. 0.1 represents the regions with the top 10 percent of insulation)
#' @param min_domainSize Domains with a size below this limit are discarded.
#' @param k Number of compartments to use for k-means based clustering of the trans contacts between domains.
#' @param insulation Insulation window used to generate the insulation tracks.
#' @param path Folder where to store the resulting data
#'
#' @return R list containing compartment classification per cell.
#'
#' @examples
#'
#' analyzeCompartments(cells='eGFP',ins_dom_thresh=0.1,insulation=100)    #One sample using insulation window of 100kb. Default number of compartments - 2.
#' analyzeCompartments(cells='eGFP',ins_dom_thresh=0.1,insulation=100)  #Two samples, all chromosomes
#'
#' @export
##########################################################################################################
analyzeCompartments <- function(cells=cells,domains=NULL,hic_tracks=all_tracks,insTracks=ins_tracks,ins_dom_thresh=0.1,min_domainSize=5e4,k=2,insulation=250,path=NULL){
	res <- list()
	trim_below=-2
	for (cell in cells){
		tracks <- hic_tracks[grep(cell,hic_tracks)]
		options(gmax.data.size=5e07)
		if (is.null(domains)){
  		d_name <- paste0('hic.',cell,'.ins',insulation,'_k',k,'_domains')
  		ins_track <- ins_tracks[grep(cell,insTracks)]
  		d <- get_domains(ins_track,min_domainSize=min_domainSize,ins_dom_thresh=ins_dom_thresh)
		} else {
		  d_name <- paste0(domains,'_k',k,'_domains')
		  d <- gintervals.load(domains)[,1:3]
		}
  	d$len <- d$end - d$start
		cat("Creating intervals ...\n")
		d_inds = as.data.frame(expand.grid(1:nrow(d), 1:nrow(d)))
		d_ints = data.frame(chrom1=d[d_inds[,1],'chrom'], start1=d[d_inds[,1],'start'], end1=d[d_inds[,1],'end'], chrom2=d[d_inds[,2],'chrom'], start2=d[d_inds[,2],'start'], end2=d[d_inds[,2],'end'])
		ints = gintervals.canonic(d_ints)
		trans_ints <- ints[ as.character(ints$chrom1) != as.character(ints$chrom2),]
		d_cov <- ints
		d_cov[,7:8] <- NA
		colnames(d_cov)[7:8] <- c('cov','intervalID')
		d_cov$intervalID <- seq(1:nrow(d_cov))
		ps_tn = sprintf("%s_area", tracks)
		for (t in 1:length(tracks)){
			gvtrack.create(ps_tn[t], tracks[t], "area")
		}
		cat(sprintf("Extracting from %s ...\n", cell))
		trans_cov = gextract(paste0(ps_tn,collapse='+'), intervals=trans_ints, iterator=trans_ints,colnames='cov')
		trans_cov$intervalID <- as.numeric(row.names(trans_ints))
		d_cov[d_cov$intervalID %in% trans_cov$intervalID, 'cov'] <- trans_cov$cov
		m = matrix(nrow=nrow(d), ncol=nrow(d), data=0)
		d_reord_inds = d_inds[attr(ints, 'mapping'),]
		m[cbind(d_reord_inds[,1], d_reord_inds[,2])] = d_cov$cov
		s = matrix(nrow=nrow(m), ncol=ncol(m), data=0)
		s[ cbind(d_inds[,1], d_inds[,2]) ] = d[ d_inds[,1], 'len'] * d[ d_inds[,2], 'len']
		s = (s * sum(m+1e10, na.rm=T)) / sum(s[!is.na(m)])
		m = log2( (m+1e10) / s)
		if (!is.na(trim_below)) {
			m[ m < trim_below ] = trim_below
		}
		message("clustering...")
		cx = hc = km = NULL
		km = TGLKMeans_wrapper(m, paste0(main_f,'temp/',cell,'_comp.tab'), k)
		d$len <- d$end-d$start
		d$cluster = km$cluster
		if (!gintervals.exists(d_name)){gintervals.save(d_name,intervals=d)}
		res[[cell]] <- d
		if(!is.null(path)){write.table(d,paste0(path,'/',cell,'_k',k,'_ins',insulation,'.tsv'),quote=F,col.names=T,row.names=F,sep='\t')}
	}
	return(res)
}

###################################################################################################
#'  Extract Binned observed data from a hic track
#'
#' \code{extractBinned}
#'
#' This function calculates a hic matrix with a given bin size per replicate or per condition (pooled replicates) and saves the data to a file for subsequent plotting.
#' Output file is a list of list with chromosomes and conditions as indices. KR matrix balancing is done automatically on the full chromosome and stored separately.
#'
#' @param trackdb misha database to use. Defaults to the one set in the config file.
#' @param tracks Which tracks to consider when plotting. Only tracks containing the selected cell parameter will be used. Default to all tracks set in the config.
#' @param cells Which condition to work on. Vector of conditions is also accepted.
#' @param chrs Which chromosomes to extract. Should be a valid chromosome name. Could also contain a vector of chrom names (recommended to run all chromosomes at the same time).
#' @param path Where to store the extracted data. Default is "data/extractedBins"
#' @param binSize Resolution of the hic matrix to extract
#' @param file_f Full path to the output data file.
#'
#' @examples
#' path=paste0(main_f,'data/extractedBins/')
#' file_f=paste0(path,eGFP_100kb_chr2L')
#' extractBinned(cells='eGFP',chrs='chr2L',binSize=1e5,path=path,file_f=file_f)    #One sample, one chromosome
#' extractBinned(cells=c('eGFP','Baf'),chrs=gintervals.all()[,1],path=path,binSize=1e5,file_f=paste0(path,eGFP_100kb_allChr'))  #Two samples, all chromosomes
#'
#' @export
##########################################################################################################
extractBinned <- function(trackdb=trackdb,tracks=all_tracks,cells=cells,chrs=chrs,path=path,binSize=binSize,file_f=file_f){
	for (t in 1:length(tracks)){
		gvtrack.create(paste0('v_',tracks[t]), tracks[t], "area")
	}
	shuffled_tracks <- paste0(tracks,'_shuffle')
	for (t in 1:length(tracks)){
		gvtrack.create(paste0('v_',shuffled_tracks[t]), shuffled_tracks[t], "area")
	}
	if (file.exists(file_f)){
		df_list <- get(load(file_f))
	} else {
		df_list <- list()
	}
	for (cell in cells){
		temp_list <- list()
		for (chr in chrs){
			iter2d <- giterator.intervals(intervals=gintervals.2d(chr), iterator=c(binSize,binSize), band=c(-gintervals(chr)$end,0))
			df_obs <- gextract(paste0('v_',tracks[grep(cell,tracks)],collapse='+'),intervals=iter2d,iterator=iter2d)
			df_obs <- df_obs[,-ncol(df_obs)]
			df_obs$bin1 <- paste0(gsub('chr','',df_obs$chrom1),'_',df_obs[,2]/binSize)
			df_obs$bin2 <- paste0(gsub('chr','',df_obs$chrom1),'_',df_obs[,5]/binSize)
			df_obs <- df_obs[,7:ncol(df_obs)]
			colnames(df_obs)[1] <- 'count'
			is.na(df_obs) <- sapply(df_obs, is.infinite)
			df_exp <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)],collapse='+'),')/2'),intervals=iter2d,iterator=iter2d)
			df_exp <- df_exp[,-ncol(df_exp)]
			df_exp$bin1 <- paste0(gsub('chr','',df_exp$chrom1),'_',df_exp[,2]/binSize)
			df_exp$bin2 <- paste0(gsub('chr','',df_exp$chrom1),'_',df_exp[,5]/binSize)
			df_exp <- df_exp[,7:ncol(df_exp)]
			colnames(df_exp)[1] <- 'count'
			is.na(df_exp) <- sapply(df_exp, is.infinite)
			df_cast_obs <- dcast(df_obs,factor(bin1,levels=unique(bin1))~factor(bin2,levels=unique(bin2)),value.var='count',)
			colnames(df_cast_obs)[1] <- 'bin'
			row.names(df_cast_obs) <- df_cast_obs$bin
			df_cast_obs <- df_cast_obs[,-1]
			df_cast_obs[lower.tri(df_cast_obs)] = t(df_cast_obs)[lower.tri(df_cast_obs)]
			df_cast_exp <- dcast(df_exp,factor(bin1,levels=unique(bin1))~factor(bin2,levels=unique(bin2)),value.var='count',)
			colnames(df_cast_exp)[1] <- 'bin'
			row.names(df_cast_exp) <- df_cast_exp$bin
			df_cast_exp <- df_cast_exp[,-1]
			df_cast_exp[lower.tri(df_cast_exp)] = t(df_cast_exp)[lower.tri(df_cast_exp)]
			df_norm <- log2(df_cast_obs/df_cast_exp)
			temp_list[[as.character(chr)]]$obs <- df_cast_obs
			temp_list[[as.character(chr)]]$obs_KR <- KRnorm(as.matrix(df_cast_obs))
			temp_list[[as.character(chr)]]$exp <- df_cast_exp
			#temp_list[[as.character(chr)]]$exp_KR <- KRnorm(as.matrix(df_cast_exp))
			temp_list[[as.character(chr)]]$norm <- df_norm
			#temp_list[[as.character(chr)]]$norm_KR <- KRnorm(as.matrix(df_norm))
		}
		df_list[[cell]] <- temp_list
	}
	save(df_list,file=file_f)
}

# calculate_eigen <-function(trackdb=trackdb,tracks=all_tracks,cells=cells,chrs=chrs,path_f=path,bin=binSize,gen=genome,cworld=cworld_path,chrom_f=chroms_sizes_f){
# 	for (t in 1:length(tracks)){
# 		gvtrack.create(paste0('v_',tracks[t]), tracks[t], "area")
# 	}
# 	for (cell in cells){
# 		for (chr in chrs){
# 			bins_2d <- giterator.intervals(intervals=gintervals.2d(chr), iterator=c(bin,bin))
# 			for (cell in cells){
# 				message(cell,'_chr',chr)
# 				if(!file.exists(paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'))|(file.info(paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'))$size==0)){
# 					dir.create(paste0(path_f,cell,'/',bin/1000,'kb'),showWarnings = FALSE, recursive=T)
# 					df <- gextract(paste0('v_',tracks[grep(cell,tracks)],collapse='+'),intervals=bins_2d,iterator=bins_2d)
# 					df <- df[,-ncol(df)]
# 					df$bin1 <- paste0(df[,2]/bin,'|',gen,'|chr',chr,':',df[,2],'-',df[,3])
# 					df$bin2 <- paste0(df[,5]/bin,'|',gen,'|chr',chr,':',df[,5],'-',df[,6])
# 					df <- df[,7:ncol(df)]
# 					colnames(df)[1] <- 'count'
# 					df_cast <- dcast(df,factor(bin1,levels=unique(bin1))~factor(bin2,levels=unique(bin2)),value.var='count',)
# 					colnames(df_cast)[1] <- 'bin'
# 					write.table(df_cast,paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'),col.names=T,row.names=F,quote=F,sep='\t')
# 					command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ;'),'perl -I ',cworld,'/lib/ ',cworld,'/scripts/perl/matrix2compartment.pl -i ',paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'),' --et --minDist 2000000 -o ',path_f,cell,'/',bin/1000,'kb/chr',chr,'_eigen')
# 					system(command_f,wait = TRUE)
# 				} else if (!file.exists(paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_eigen.zScore.eigen1.bedGraph'))){
# 					command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ;'),'perl -I ',cworld,'/lib/ ',cworld,'/scripts/perl/matrix2compartment.pl -i ',paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'),' --et --minDist 2000000 -o ',path_f,cell,'/',bin/1000,'kb/chr',chr,'_eigen')
# 					system(command_f,wait = TRUE)
# 				} else {
# 					message('Eigen track exist. Skipping...')
# 				}
# 			}
# 		}
# 		command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ;'),'cat *eigen1.bedGraph > combined.bedGraph')
# 		system(command_f,wait = TRUE)
# 		command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ; grep -v bedGraph combined.bedGraph > combined_pol.bedGraph; mv combined_pol.bedGraph combined.bedGraph; sort -k1,1 -k2,2n combined.bedGraph > combined_pol.bedGraph; mv combined_pol.bedGraph combined.bedGraph; '),'bedGraphToBigWig combined.bedGraph ',chrom_f,' ',paste0(path_f,cell,'_',bin/1000,'kb','.bw'))
# 		system(command_f,wait = TRUE)
# 		if(!gtrack.exists(paste0('eigen.',cell,'_',bin/1000,'kb'))){gtrack.import(track=paste0('eigen.',cell,'_',bin/1000,'kb'),description=paste0('eigen vector based on ',bin/1000,'kb'),file=paste0(path_f,cell,'/',bin/1000,'kb/combined.bedGraph'),binsize=bin)}
# 	}
# }

###############################################################################################################
#'  Use eigenvector values or chipseq track to rank bins and calculate hic contact enrichment between bin pairs
#'
#' \code{rank_matrix}
#'
#' This function ranks genomic bins based on their eigenvector value (or chip enrichment) into a predefined number of ranks and calculates the contact enrichment (observed/expected) for each bin-pair.
#' The result is output as a file,which is used for plotting by \code{plot_rankMatrix}
#'
#' @param trackdb misha database to use. Defaults to the one set in the config file.
#' @param tracks Which hic tracks to consider when plotting. Only tracks containing the selected cell parameter will be used. Default to all tracks set in the config.
#' @param eigen_tracks Which linear tracks to use to rank genomic bins. Only tracks containing the string with cell parameter will be considered. Default is eigen tracks.
#' @param cells Which condition to work on. Vector of conditions is also accepted.
#' @param chrs Which chromosomes to extract. Should be a valid chromosome name. Default is all chromsomes.
#' @param path Where to store the extracted data. Default is "data/extractedBins"
#' @param binSize Size of the genomic bins to use.
#' @param min_dist Only bin pairs separated by at least this distance are considered. Necessary to avoid working too close to the diagonal.
#' @param ranks How many ranks to use. Recommended 50 for sparse data and 100 for high-resolution.
#' @param file_f Full path to the output data file.
#'
#' @examples
#'
#' path=paste0(main_f,'/data/extractedBins/')
#' file_name <- paste0(path,'rankMatrix_100kb_ranks50')
#' rank_matrix(eigen_tracks='eigen.eGFP_50kb',cells='eGFP',path=path,binSize=2e5,min_dist=1e7,ranks=50,file_f=file_name)   #One sample, all chromosomes.
#' rank_matrix(eigen_tracks=c('eigen.eGFP_50kb','eigen.Baf_50kb'),cells=('eGFP','Baf'),chrs='2L',path=path,binSize=1e5,min_dist=1e7,ranks=50,file_f=file_name)   #Two samples, one chromosome.
#'
#' @export
##########################################################################################################
rank_matrix <- function(trackdb=trackdb,tracks=all_tracks,eigen_tracks=eigen_tracks,cells=cells,chrs=gintervals.all()$chrom,path=paste0(main_f,'/data/extractedBins/'),binSize=1e5,min_dist=1e7,ranks=50,file_f=paste0(main_f,'/data/extractedBins/test')){
	options(gmax.data.size=5e7)
  for (t in 1:length(eigen_tracks)){
		gvtrack.create(paste0('v_',eigen_tracks[t]), eigen_tracks[t], 'avg')
	}
	for (t in 1:length(tracks)){
		gvtrack.create(paste0('v_',tracks[t]), tracks[t], "area")
	}
	shuffled_tracks <- paste0(tracks,'_shuffle')
	for (t in 1:length(tracks)){
		gvtrack.create(paste0('v_',shuffled_tracks[t]), shuffled_tracks[t], "area")
	}
	temp_list2 <- list()
	for (cell in cells){
		eigen <- gextract(paste0('v_',eigen_tracks[grep(cell,eigen_tracks)]),intervals=gintervals(chrs),iterator=binSize)
		eigen$bin <- paste0(gsub('chr','',eigen$chrom),'_',eigen[,2]/binSize)
		intervals1 <- eigen[,1:3]
		intervals2 <- eigen[,1:3]
		eigen <- eigen[order(eigen[,4]),]
		eigen$rank <- ntile(eigen[,4], ranks)
		iter2d <- giterator.intervals(intervals=gintervals.2d(chrs), iterator=c(binSize,binSize), band=c(-max(gintervals(chrs)$end),-min_dist))
		message('Working with grid: ',nrow(iter2d))

		if (nrow(iter2d)>10000000){
				chunk <- 10000000
				n <- nrow(iter2d)
				r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
				d <- split(iter2d,r)
				d1 <- lapply(d, function(x) gextract(paste0('v_',tracks[grep(cell,tracks)],collapse='+'),intervals=x,iterator=x))
				df_obs <- do.call("rbind", d1)
				row.names(df_obs) <- seq(1,nrow(df_obs))
				rm(d)
				rm(d1)
				rm(r)
			} else {
				df_obs <- gextract(paste0('v_',tracks[grep(cell,tracks)],collapse='+'),intervals=iter2d,iterator=iter2d)
			}
		df_obs <- df_obs[,-ncol(df_obs)]
		df_obs$bin1 <- paste0(gsub('chr','',df_obs$chrom1),'_',df_obs[,2]/binSize)
		df_obs$bin2 <- paste0(gsub('chr','',df_obs$chrom1),'_',df_obs[,5]/binSize)
		df_obs <- df_obs[,7:ncol(df_obs)]
		colnames(df_obs)[1] <- 'count'
		is.na(df_obs) <- sapply(df_obs, is.infinite)

		if (nrow(iter2d)>10000000){
				chunk <- 10000000
				n <- nrow(iter2d)
				r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
				d <- split(iter2d,r)
				d1 <- lapply(d, function(x) gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)],collapse='+'),')/2'),intervals=x,iterator=x))
				df_exp <- do.call("rbind", d1)
				row.names(df_exp) <- seq(1,nrow(df_exp))
				rm(d)
				rm(d1)
				rm(r)
			} else {
				df_exp <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)],collapse='+'),')/2'),intervals=iter2d,iterator=iter2d)
			}
		df_exp <- df_exp[,-ncol(df_exp)]
		df_exp$bin1 <- paste0(gsub('chr','',df_exp$chrom1),'_',df_exp[,2]/binSize)
		df_exp$bin2 <- paste0(gsub('chr','',df_exp$chrom1),'_',df_exp[,5]/binSize)
		df_exp <- df_exp[,7:ncol(df_exp)]
		colnames(df_exp)[1] <- 'count'
		is.na(df_exp) <- sapply(df_exp, is.infinite)

		rank_matrix <- matrix(nrow=ranks,ncol=ranks)
		obs_matrix <- matrix(nrow=ranks,ncol=ranks)
		exp_matrix <- matrix(nrow=ranks,ncol=ranks)
		for (j in 1:ranks){
			for (k in 1:ranks){
				bins1 <- eigen[eigen$rank==j,'bin']
				bins2 <- eigen[eigen$rank==k,'bin']
				sum_obs <- sum(df_obs[(df_obs$bin1 %in% bins1)&(df_obs$bin2 %in% bins2),'count'],na.rm=T)
				sum_exp <- sum(df_exp[(df_exp$bin1 %in% bins1)&(df_exp$bin2 %in% bins2),'count'],na.rm=T)
				obs_matrix[j,k] <- sum_obs
				exp_matrix[j,k] <- sum_exp
				rank_matrix[j,k] <- log2(sum_obs/sum_exp)
			}
		}
		temp_list2[[cell]]$obs <- df_obs
		temp_list2[[cell]]$exp <- df_exp
		temp_list2[[cell]]$rank_matrix <- rank_matrix
		temp_list2[[cell]]$obs_matrix <- obs_matrix
		temp_list2[[cell]]$exp_matrix <- exp_matrix
	}
	temp_list2$binSize <- binSize
	temp_list2$ranks <- ranks
	temp_list2$cells <- cells
	temp_list2$chrs <- chrs
	temp_list2$eigen_tracks <- eigen_tracks
# 	file_f <- paste0(path,'rankMatrix_',binSize/1000,'kb_ranks',ranks)
	save(temp_list2,file=file_f)
}

#####################################################################################################################################################
#'  Calculate compartment A-score as defined in Nagano et al., Nature 2017. Just a wrapper for the functions in https://bitbucket.org/tanaylab/schic2
#'
#' \code{calc_Ascore}
#'
#' This function calculates a score for each domain which represents the percentage of all trans interactions that the domain makes with A-type domains.
#' The result is output as a R list,which is saved in the subdirectory '/table'
#'
#' @param sch_data_dir misha database to use. Defaults to the one set in the config file.
#' @param cells Which condition to work on. Vector of conditions is also accepted.
#' @param tracks Which tracks to consider when plotting. Only tracks containing the selected cell parameter will be used. Default to all tracks set in the config.
#' @param insulation Insulation window on which compartments have been previously identified.
#' @param chrom_f File with chromosome file sizes. Read from the config.
#'
#' @examples
#'
#' calc_Ascore(sch_data_dir=paste0(main_f,'analysis/compartments/'),cells='eGFP',insulation=100)   #Two samples, one chromosome.
#'
#' @export
##########################################################################################################
calc_Ascore <- function(sch_data_dir=path,cells=cells,tracks=all_tracks,insulation=250,chrom_f=chrom_sizes_f){
	gdir.create('ascore', showWarnings=F)
	dir.create(paste0(main_f,sch_data_dir), showWarnings = FALSE,recursive=T)
	sch_table_dir <<- sprintf("%s/tables", sch_data_dir)
	dir.create(sch_table_dir, showWarnings = FALSE,recursive=T)
	sch_fig_dir <<-  sprintf("%s/figs", sch_data_dir)
	dir.create(sch_fig_dir, showWarnings = FALSE,recursive=T)
	sch_rdata_dir <<-  sprintf("%s/rdata", sch_data_dir)
	dir.create(sch_rdata_dir, showWarnings = FALSE,recursive=T)
	temp_list <- list()
	for (cell in cells){
		d <- gintervals.load(paste0('hic.',cell,'.ins',insulation,'_k2_domains'))
		ab_tags = c('B', 'A')
		d$ab = ab_tags[d$cluster]
		tracks <- all_tracks[grep(cell,all_tracks)]
		nms = paste0(rep(cell,length(tracks)),'_rep',1:length(tracks))
		d = d[complete.cases(d), ]
		m = matrix(0, nrow(d), length(tracks))
		rownames(m) = paste0(as.character(d$chrom),'_',d$start)
		colnames(m) = nms
		for (i in 1:length(tracks)) {
			rv = calc_domains_a_score(tn=tracks[i], d=d, use_trans=T, rebuild=F)
		#	rv <- rv[match(row.names(m),row.names(rv)),]
			m[rownames(rv), nms[i]] = rv$trans_A_score
		  }
		out_f <- d[,1:4]
		out_f[,4] <- rowMeans(m,na.rm=T)
		out_nm <- paste0(sch_table_dir,'/',cell,'_Ascore.bedGraph')
		out_f <- out_f[complete.cases(out_f),]
		write.table(out_f,out_nm,quote=F,col.names=F,row.names=F,sep='\t')
		command_f <- paste0('sort -k1,1 -k2,2n ',out_nm,' > ',paste0(out_nm,'.sorted'),'; mv ',paste0(out_nm,'.sorted'),' ',out_nm,'; bedGraphToBigWig ',out_nm,' ',chrom_f,' ',gsub('bedGraph','bw',out_nm))
	  	system(command_f,wait = TRUE)
		if(!gtrack.exists(paste0('ascore.',cell,'_',insulation,'kb'))){gtrack.import(track=paste0('ascore.',cell,'_',insulation,'kb'),description=paste0('Ascore based on ',insulation,'kb'),file=gsub('bedGraph','bw',out_nm),binsize=0)}
	  temp_list[[cell]] <- m
	}
	save(temp_list,file=paste0(sch_table_dir,'/Ascore_combined'))
}

######################################################################################################################################
#'  Use either domain type (Nagano et al., 2017) or eigenvector values to determine the locations and change in compartment boundaries
#'
#' \code{comp_boundaries}
#'
#' This function determines the genomic locations which act as a borders between domains (or regions) with different compartment identity. It operates in two major modes controlled by the paramater type - per domain or based on the eigenvector value in genomic bins.
#' The result is a table containing the compartment borders per cell and a table listing the stats of compartment boundaries.
#'
#' @param path Where to store the extracted data. Default is "data/extractedBins"
#' @param type One of "domain" or "eigen"
#' @param cells Which condition to work on. Vector of conditions is also accepted.
#' @param window_f Insulation window on which compartments have been previously identified when in "domain" or eigenvector resolution.
#' @param write_borders Boolean controlling whether to write the compartment borders per condition in a file.
#'
#' @return Returns a table with the results.
#'
#' @examples
#'
#' res <- comp_boundaries(type='domain',cells=cells,window=100,write_borders=T)    # Output compartment borders based on domains identified using 100kb window.
#' res <- comp_boundaries(type='eigen',cells=cells,write_borders=T)    # Output compartment borders based on eigenvector values.
#'
#' @export
##########################################################################################################
comp_boundaries <- function(path_f=paste0(main_f,'/analysis/compartments/'),type='amos',cells=cells,window_f=250,write_borders=F){
	res_m <- as.data.frame(matrix(data = NA, nrow = length(cells), ncol = 2))
	colnames(res_m) <- c('B-A','A-B')
	row.names(res_m) <- cells
	for (cell in cells){
		if(type!='eigen'){
			d <- gintervals.load(paste0('hic.',cell,'.ins',window_f,'_k2_domains'))
			d <- d[d$chrom!='chrY',]
			d <- d[d$chrom!='chrM',]
			d <- d[grep('Het',d$chrom,invert=T),]
			d <- d[grep('extra',d$chrom,invert=T),]
			d <- d[grep('chrU',d$chrom,invert=T),]
			comp <- ddply(d,.(chrom),function(x){
				int <- x[1,]
				for (i in 2:nrow(x)){
					if(x[i,5]!=x[(i-1),5]){
						int <- rbind(int,x[i,])
					} else {
						int[nrow(int),3] <- x[i,3]
						}
				}
				return(int[,-c(1,4)])
			})

			df <- gintervals.diff(gintervals.all(),comp)
			df <- df[df$chrom!='chrM',]
			df <- df[df$chrom!='chrY',]
			df <- df[grep('Het',df$chrom,invert=T),]
			df <- df[grep('extra',df$chrom,invert=T),]
			df <- df[grep('chrU',df$chrom,invert=T),]
			df$start=df$start-2
			df$start[df$start<=0] <- 1
			df2 <- gintervals.neighbors(df,comp,maxdist=1)
			if (write_borders){
				temp <- df2
				temp[,2] <- temp[,2]+2
				temp[temp[,7]==1,7] <- 'B'
				temp[temp[,7]==2,7] <- 'A'
				write.table(temp[,c(1:3,7)],paste0(path_f,cell,'_compBorders_',type,'.bed'),quote=F,col.names=F,row.names=F,sep='\t')
			}
			res_m[cell,1] <- nrow(df2[df2[,7]==1,])/nrow(d)*100
			res_m[cell,2] <- nrow(df2[df2[,7]==2,])/nrow(d)*100

		} else {
			compA <- gscreen(paste0(paste0('eigen.',cell,'_',window_f,'kb'),'>0'))
			compA$length=compA$end-compA$start
			compA$type=2
	#		message(cell,' comp A N:',nrow(compA))
			compB <- gscreen(paste0(paste0('eigen.',cell,'_',window_f,'kb'),'<0'))
			compB$length=compB$end-compB$start
			compB$type=1
	#		message(cell,' comp B N:',nrow(compB))

			comp <- rbind(compA,compB)
			comp[comp$start!=0,2] <- comp$start[comp$start!=0]+1
			comp <- comp[order(comp$chrom,comp$start,comp$end),]
			df <-gintervals.diff(gintervals.all(),comp)
			df <- df[df$chrom!='chrM',]
			df <- df[df$chrom!='chrY',]
			df <- df[grep('Het',df$chrom,invert=T),]
			df <- df[grep('extra',df$chrom,invert=T),]
			df <- df[grep('chrU',df$chrom,invert=T),]
			df$start=df$start-2
			df$start[df$start<=0] <- 1
			df2 <- gintervals.neighbors(df,comp,maxdist=1)
			if (write_borders){
				temp <- df2
				temp[,2] <- temp[,2]+2
				temp[temp[,8]==1,8] <- 'B'
				temp[temp[,8]==2,8] <- 'A'
				write.table(temp[,c(1:3,8)],paste0(path_f,cell,'_compBorders_',type,'.bed'),quote=F,col.names=F,row.names=F,sep='\t')
			}
			res_m[cell,2] <- nrow(df2[df2$type==1,])
			res_m[cell,1] <- nrow(df2[df2$type==2,])
		}
	}
	return(res_m)
}


##############################################################################################################################################
#'  Use either eigenvector value or Ascore (Nagano et al, 2017) to calculate a correlation between compartment strength and chipseq enirchment
#'
#' \code{comp_correlation}
#'
#' This function extracts the value of the eigenvector or Ascore per defined bin and calculates the correlation between the compartment strength and the signal enrichment from a linear tracks (such as ChIPseq or RNAseq)
#' result is displayed directly
#'
#' @param cells Which condition to work on. Vector of conditions is also accepted.
#' @param e_tracks Either eigen or Ascore tracks, default is eigen tracks defined in the configuration file.
#' @param tracks Which linear tracks to use as an input for the correlation.
#' @param binSize Bin size to calculate correlation on.
#' @param type One of "chip", "rna" or "repli" for replication timing tracks.
#' @param cutoff Only consider the regions with value in the desired track above/under this quantile (above when >50 percent and below when <50 percent)
#' @param method Method to calculate correlation.
#'
#' @examples
#'
#' comp_correlation(cells=cells,e_tracks=eigen_tracks,binSize=1e5,tracks=c('chipseq_RPM.D0_H3K4me3','chipseq_RPM.D0_H3K9me3'),type='chip')    # Correlation between eigenvector value and chipseq signal enrichment in 100kb bins.
#' comp_correlation(cells=cells,e_tracks=gtrack.ls('ascore'),binSize=1e5,tracks=gtrack.ls('rnaseq','D0'),type='rna')    # Correlation between eigenvector value and .
#'
#' @export
##########################################################################################################
comp_correlation <- function(cells=cells,e_tracks=eigen_tracks,tracks=chipseq_tracks,binSize=1e5,type='chip',cutoff=NULL,method='spearman'){
	for (t in 1:length(e_tracks)){
		gvtrack.create(paste0('v_',e_tracks[t]), e_tracks[t], 'avg')
	}
	if (type=='chip'){
  	for (t in 1:length(tracks)){
  		gvtrack.create(paste0('v_',tracks[t]),tracks[t], 'global.percentile')
  	}
	} else if (type=='repli'){
		for (t in 1:length(tracks)){
			gvtrack.create(paste0('v_',tracks[t]), tracks[t], 'avg')
		}
	}
	if (type=='rna'){
	  rev_tracks <- tracks[grep('rev',tracks,ignore.case =T)]
		for_tracks <- tracks[grep('for',tracks,ignore.case =T)]
		for (t in 1:length(for_tracks)){
			gvtrack.create(paste0('v_',for_tracks[t]), for_tracks[t], 'global.percentile')
			gvtrack.create(paste0('v_',rev_tracks[t]), rev_tracks[t], 'global.percentile')
		}
	}
	cor_matrix <- data.frame(matrix(NA, nrow = length(cells), ncol = length(tracks)))
	row.names(cor_matrix) <- cells
	colnames(cor_matrix) <- tracks
	for (cell in cells){
		if (type=='chip'){
			df <- gextract(e_tracks[grep(cell,e_tracks)],paste0('-log2(1-v_',tracks,')'),intervals=gintervals.all(),iterator=binSize)
		} else if (type=='repli'){
			df <- gextract(e_tracks[grep(cell,e_tracks)],paste0('v_',tracks),intervals=gintervals.all(),iterator=binSize)
		} else if (type=='rna'){
#		track_expression <- paste0('-log2(1-((',paste0(paste0('v_',for_tracks),collapse='+'),')/',length(for_tracks),'+(',paste0(paste0('v_',rev_tracks),collapse='+'),')/',length(rev_tracks),')/2)')
		track_expression <- c(paste0('v_',for_tracks),paste0('v_',rev_tracks))
		  print(track_expression)
			df <- gextract(e_tracks[grep(cell,e_tracks)],track_expression,intervals=gintervals.all(),iterator=binSize)
		}
	  if (type=='rna') {f_tracks <- for_tracks} else {f_tracks <- tracks}
	  for (track in f_tracks){
			if (!is.null(cutoff)){
			  if(cutoff>0.5){idx <- df[,grep(track,colnames(df))]>quantile(df[,grep(track,colnames(df))],cutoff,na.rm=T)} else {idx <- df[,grep(track,colnames(df))]<quantile(df[,grep(track,colnames(df))],cutoff,na.rm=T)}
			  df <- df[idx,]
			}
		  if (type!='rna') {cor_matrix[grep(cell,row.names(cor_matrix)),grep(track,colnames(cor_matrix))] <- cor(df[,grep(track,colnames(df))],df[,4],use='complete.obs',method=method)}
		  else {}
		}
	}
	return(cor_matrix)
	for (t in 1:length(tracks)){
		gvtrack.rm(paste0('v_',tracks[t]))
		gvtrack.rm(paste0('v_',rev_tracks[t]))
	}
}

##############################################################################################################################################
#' Calculate the correlation between several hic samples within a certain distance band
#'
#' \code{corrSamples}
#'
#' This function will extract the contacts within a certain distance band using predetermined bin sizes and return a list containing the bins above a cutoff value which can be used for pearson or spearman correlation.
#'
#' @param tracks Which tracks to work on. Defaults to all tracks listed in the config file
#' @param binSizes Bin sizes to calculate correlation on. One or vector of values are accepted.
#' @param minDist Minimum distance which to consider when constructing intervals.Recommended at least 10kb.
#' @param maxDist Maximum distance which to consider when constructing intervals. Defaults to the biggest chromosome.
#' @param cutoff Any bins with fewer number of contacts in at least one sample will be discarded. Recommended value - 5.
#' @param cor_method Which correlation method to use. One out of 'pearson' or 'spearman'.
#' @param cells_vector Vector of conditions, repeat each as many tines as there are replicates.
#' @param path Folder where to store the generated heatmap.
#'
#' @return R list containing the extracted values per binSize before applying cutoff value
#'
#' @examples
#'
#' corrSamples(tracks=all_tracks,binSizes=1e5,minDist=1e4,maxDist=max(gintervals.all()$end),cutoff=10,cells_vector=rep(cells,each=2),path=paste0(main_f,'analysis/compartments/'))    # Pearson correlation with a binsize of 100kb across all bins from 1kb to max considering only bins where there is at least 10 reads in any one condition.
#'
#' @export
##########################################################################################################
corrSamples <- function(tracks=all_tracks,binSizes=c(1e5,2.5e5,5e5),minDist=1e4,maxDist=max(gintervals.all()$end),cutoff=5,cor_method='pearson',cells_vector=rep(cells[1],2),path=paste0(main_f,'/analysis/')){
  options(gmax.data.size=2e8)
  for(track in tracks){gvtrack.create(paste0("v_", track) , track, "weighted.sum")}
	cis <- gintervals.2d.all()
	cis <- cis[cis$chrom1==cis$chrom2,]
	temp_list <- list()
	for (binSize in binSizes)
	{
		iter2d <- giterator.intervals(intervals=cis, iterator=c(binSize,binSize))
		if (nrow(iter2d)>5000000){
			chunk <- 5000000
			n <- nrow(iter2d)
			r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
			d <- split(iter2d,r)
			d1 <- lapply(d, function(x) gextract(paste0('v_',tracks), x, iterator=x, band=-c(maxDist,minDist)))
			result <- do.call("rbind", d1)
			row.names(result) <- seq(1,nrow(result))
			rm(d)
			rm(d1)
			rm(r)
		} else {
			result <- gextract(paste0('v_',tracks), iter2d, iterator=iter2d, band=-c(maxDist,minDist))
		}
		temp_result <- result[rowMaxs(as.matrix(result[,7:(ncol(result)-1)]),na.rm=T)>cutoff,]
		temp_list[[binSize]] <- result[,-ncol(result)]
		result <- temp_result
		res <- cor(result[,7:(ncol(result)-1)],use='complete.obs',method=cor_method)
#		colnames(res) <- sapply(strsplit(tracks,'\\.'),function(x){x[3]})
#		row.names(res) <- sapply(strsplit(tracks,'\\.'),function(x){x[3]})
		distance <- dist(1-res)
		ncells <- length(unique(cells_vector))
		cols <- rainbow(ncells)
		col<- colorRampPalette(c("blue", "white", "red"))(100)
		write.table(res,paste0(path,'cor_bin',binSize/1000,'_minDist_',minDist,'_maxDist',maxDist,'.tsv'),quote=F,col.names=T,row.names=T,sep='\t')
		res <- round(res, 2)
		annotation_df <- data.frame(cell_Type=factor(cells_vector))
		colnames(annotation_df) <- c('cell_type')
		row.names(annotation_df) <- colnames(res)
		pheatmap(res,color=col,breaks=seq(0.7,1,length=101),clustering_distance_rows=distance,clustering_distance_cols=distance,clustering_method='ward.D', display_numbers = TRUE,fontsize_number = 10,border_color = "black",annotation_col=annotation_df,filename=paste0(path,'bin',binSize/1000,'_minDist_',minDist,'_maxDist',maxDist,'_heatmap.pdf'),height=10,width=10)

	}
	file_f <- paste0(path,'bin',binSize/1000,'_minDist_',minDist,'_maxDist',maxDist)
	save(temp_list,file=file_f)
# 	return(temp_list)
}

##############################################################################################################################################
#' Identifies regions where there is a compartment switch and outputs clusters based in the eigenvector value and the transcriptional output in each cluster.
#'
#' \code{comp_expression}
#'
#' This function will extract eigenvector values and determine regions where there is an change from positive to negative (or vice versa) - i.e. compartment switch. If a cutoff value is specified (percentage quantile from 1-100) it will consider only those regions where the difference between maximum and minimum eigenvalue is higher than the quantile specified. It will then cluster the eigenvalues using kmeans and calculate the average expression values (fpkm) across genes within the resulting clusters.
#'
#' @param cells Which condition to work on. Vector of conditions is also accepted. Defaults to conditions listed in the config file.
#' @param tracks Which eigen tracks to work on. Defaults to all eigen tracks listed in the config file. Only eigen tracks for the selected conditions will be further considered.
#' @param binSize Bin size to extract eigenvalue. Should be equal the binsize with which the eigenvectors were created.
#' @param beanplot_f Boolean whether to use beanplots or boxplots for expression plotting.
#' @param path_f Where to store the data.
#' @param cutoff Any bins with fewer number of contacts in at least one sample will be discarded. Recommended value: 5.
#' @param cell_colors Vector of colors to use when plotting the correlation heatmap. Should be same length as the number of conditions.
#' @param tss gintervals holding the information of gene promoters. Defined in the config file.
#' @param expression_f Expression file containing gene expression information. Format should be a tab-separated file where each row is a separate gene and the columns are individual replicates. Column names should match the same names defined in the config file.
#'
#' @return A vector with p-value of difference in expression between first, 2nd and last condition per cluster.
#'
#' @examples
#'
#' comp_expression(cells=cells,K=6,beanplot_f=T,cutoff=NULL,cell_colors=rainbow(length(cells)),path_f=paste0(main_f,'analysis/compartments/'),expression_f="/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt")
#'
#' @export
##########################################################################################################
comp_expression <- function(cells=cells,tracks=eigen_tracks,binSize=1e5,beanplot_f=T,path_f=path,K,cutoff=NULL,cell_colors=rainbow(length(cells)),tss=gintervals.load(tss_f),expression_f="/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt"){
	require(beanplot)
	tracks <- as.vector(unlist(sapply(cells,function(x){tracks[grep(x,tracks)]})))
	df <- gextract(tracks,gintervals.all(),iterator=binSize)
	df <- df[complete.cases(df),]

	df <- df[,-ncol(df)]
	idx <- rowSums(df[,4:ncol(df)]>0)
	idx <- (idx!=0)&(idx!=length(tracks))
	eigen_df <- df[idx,]
	if (!is.null(cutoff)) {
		logFC <- abs(rowMaxs(as.matrix(eigen_df[,4:ncol(eigen_df)]),na.rm=T)-rowMins(as.matrix(eigen_df[,4:ncol(eigen_df)]),na.rm=T))
		eigen_df <- eigen_df[logFC>quantile(logFC,cutoff),]
	}

	eigen_clust <- cluster_and_plot(eigen_df[,4:ncol(eigen_df)], paste0(path_f,"t"), paste0(path_f,"eigenCluster_",K,".png"), K, image_to_file=TRUE, zlim=c(-max(eigen_df[,4:ncol(eigen_df)]), max(eigen_df[,4:ncol(eigen_df)])), fig.height=3200, fig.width=1600, shades=colorRampPalette(c("blue","white","red")), shade_count=100)

	peaks_df <- eigen_df[,1:3]
	peaks_df$cluster <- eigen_clust$clust

	peaks_df_inters <- gintervals.neighbors(gintervals.load(tss_f),peaks_df)
	tss_peaks <- peaks_df_inters[abs(peaks_df_inters$dist)==0,]

	#fpkm_f <- read.table(expression_f,header=T)
	#clust_colors=rainbow(K)

	#fpkm <- merge(fpkm_f,tss_peaks[,c('geneName','cluster')],by.x='gene_name',by.y='geneName')
	#write.table(fpkm,paste0(path_f,"eigenFPKM_K",K,".txt"),row.names=F,col.names=T,quote=F,sep='\t')
	#fpkm_df <- fpkm[,-1]
	#fpkm_df <- unique(fpkm_df)

	#fpkm_means <- as.data.frame(matrix(NA,ncol=length(cells),nrow=nrow(fpkm_df)))
	#colnames(fpkm_means) <- cells
	#for (i in 1:length(cells)){
	#	fpkm_means[,i] <- rowMeans(fpkm_df[,grep(cells[i],colnames(fpkm_df))],na.rm=T)
	#	}
#	fpkm_means$cluster=fpkm_df$cluster

	#stat_test <- as.data.frame(matrix(NA,nrow=K,ncol=2))
	#colnames(stat_test) <- c(paste0(cells[1],'vs',cells[2]),paste0(cells[1],'vs',cells[length(cells)]))

	pdf(paste0(path_f,"eigenCluster_expression",K,".pdf"),width=4,height=K*4)
	par(mfrow=c(K,1))
	for (i in c(rev(seq(1,K,by=1)))){
		x <- fpkm_means[fpkm_means$cluster==i,]
		x <- x[,-ncol(x)]
		stat_test[i,1] <- wilcox.test(x[,1],x[,2],paired=T)$p.value
		stat_test[i,2] <- wilcox.test(x[,1],x[,ncol(x)],paired=T)$p.value
		if (beanplot_f){
			beanplot(x+0.01, what=c(1,1,1,0),log='y', side='no', col=as.list(cell_colors),names=cells,main=paste0('Cluster: ',i,' Number of genes: ',nrow(x)))
		} else {
			boxplot(log10(x+1),col=cell_colors,names=cells,notch=FALSE,outline=FALSE,ylab='log10(FPKM+1)',main=paste0('Cluster: ',i,' Number of genes: ',nrow(x)))
			points(1:4,colMeans(log10(x+1),na.rm=T),pch='x',col='darkgrey',lwd=3)
		}
	}
	dev.off()
	return(stat_test)
}

##############################################################################################################################################
#' Calculates and plots the cis decay profiles for selected tracks.
#'
#' \code{run_cisDecay}
#'
#' This function will calculate the distribution of contacts across range of distances. Works in two major modes: log10-log10 or normalized frequency in log2 base intervals. In addition when run in log10 mode it will also output the linear coefficient sigma over a defined distance interval which can be used to estimate chromatin compaction.
#'
#' @param tracks Which hic tracks to work on. Defaults to all tracks listed in the config file. Only tracks for the selected conditions will be further considered.
#' @param cells Which condition to work on. Vector of conditions is also accepted. Defaults to conditions listed in the config file.
#' @param path Where to store the data.
#' @param log_base One of 2 or 10 indicating whether to extract data with log2 or log10 genomic bins.
#' @param norm.by.binsize Boolean value whether to account for the size of the bin when run in log10 mode. Default - TRUE.
#' @param norm.by.total.obs Boolean value whether to normalize for the number of total contacts when run in log10 mode in order to compare across samples. Default - TRUE.
#' @param file_f Name of the file which stores the data. It will be saved in the directory path defined above.
#' @param minDist Minimum distance to consider when counting contacts. Default value - 1000.
#' @param maxDist Maximum distance to consider when counting contacts. Default value - 2e8. It is recommended to extract the data using the length of the biggest chromosome as a guide for this value. It can be later adjusted for plotting.
#' @param coeff_range Vector of distances where to calculate the linear coefficient sigma as described above.
#' @param colors Vector of colors to use when plotting. Should be same length as cellsToPlot.
#' @param out_f Name of the output png file to be generated.
#' @param cellsToPlot Which conditions to plot.Should be the same or a subset of the conditions defined in cells.
#' @param alpha If dispersion colors are not defined, it will use create them based on the cell colors and the transparency defined with this paramater.
#' @param disp_colors Vector of colors to use when plotting the dispersion values. Should be same length as cellsToPlot.
#'
#' @examples
#'
#' run_cisDecay(tracks=all_tracks,cells=cells,path=paste0(main_f,'analysis/cis_decay/'),log_base=2,file_f='example',maxDist=1e8,colors=c('black','red','blue'),alpha=0.5,disp_colors=c('gray','pink','lightblue'))
#'
#' @export
##########################################################################################################
run_cisDecay <-function(tracks=all_tracks,cells=cells,path=paste0(main_f,'analysis/cis_decay/'),log_base=2,norm.by.binsize=T,norm.by.total.obs=T,file_f='test',minDist=1000,maxDist=2e8,coeff_range=c(1e5,2e6),colors=rainbow(length(cells)),out_f=paste0(main_f,'analysis/cis_decay/cisDecay_log',log_base,'.pdf'),cellsToPlot=cells,alpha=0.7,disp_colors=NULL){
	file_f <- paste0(path,file_f,'_log',log_base)
	tracks <- as.vector(unlist(sapply(cells,function(x){tracks[grep(x,tracks)]})))
	if (length(colors) != length(cells)) {stop("Number of colors is not equal to number of cells")}
	if (!file.exists(file_f)){
		start_n <- round(log(minDist) / log(log_base),2)
		end_n <- round(log(maxDist) / log(log_base),2)
		dist_f = log_base^seq(start_n, end_n, by=0.1/log_base)
		dist = seq(start_n, end_n, by=0.1/log_base)[-1]
		a = adply(tracks, 1, function(x) { g = gcis_decay(x, dist_f, gintervals.all(), gintervals.all(), gintervals.2d.all()); return(as.vector(g[,1])) })
		a$type <- as.vector(unlist(sapply(tracks,function(x){unlist(strsplit(x,'\\.'))[2]})))
		a_list <- list(dist_f=dist_f,a=a,dist=dist)
		save(a_list,file=file_f)
	} else {
		a_list <- get(load(file_f))
		a <- a_list[['a']]
		dist_f <- a_list[['dist_f']]
		dist <- a_list[['dist']]
	}

	a <- a[a$type %in% cellsToPlot,]
	hh <- hashmap (cells,colors)
	a$colors <- hh[[a$type]]
	if (is.null(disp_colors)){
		a$disp_colors <- as.vector(unlist(sapply(a$colors,function(x){makeTransparent(x,alpha=alpha)})))
	} else {
		hh <- hashmap (cells,disp_colors)
		a$disp_colors <- hh[[a$type]]
	}

	if (log_base==2){
		d_sums = ddply(a, .(type), function(x) {
			return(colSums(x[,2:(ncol(x)-3)]))
		})

		d_means = ddply(a, .(type), function(x) {
			means <- colMeans(x[,2:(ncol(x)-3)]/rowSums(x[,2:(ncol(x)-3)]))
			return(data.frame(t(means),x$colors[1],x$disp_colors[1]))
		})
		colnames(d_means)[(ncol(d_means)-1):ncol(d_means)] <- c('colors','disp_colors')
		is.na(d_means) <- sapply(d_means,is.infinite)

		d_sterr = ddply(a, .(type), function(x) {
			temp <- x[,2:(ncol(x)-3)]/rowSums(x[,2:(ncol(x)-3)])
			return(data.frame(t(apply(temp,2,function(x) sd(x))),x$colors[1],x$disp_colors[1]))
		})
		colnames(d_sterr)[(ncol(d_means)-1):ncol(d_means)] <- c('colors','disp_colors')
		is.na(d_sterr) <- sapply(d_sterr,is.infinite)

		max_mean <- max(d_means[,2:(ncol(d_means)-2)],na.rm=T)
		min_mean <- min(d_means[,2:(ncol(d_means)-2)],na.rm=T)
		max_sterr <- max(d_sterr[,2:(ncol(d_sterr)-2)], na.rm=T)
		min_sterr <- min(d_sterr[,2:(ncol(d_sterr)-2)], na.rm=T)

		axis_pos <- seq(min(dist),max(dist),by=1)
		axis_labels <- c(paste0(round_any(2^seq(min(dist),20,by=1),1000)/1000," kb"),paste0(round_any(2^seq(20,(max(dist)),by=1),10000)/1000000," mb"))

		pdf(out_f,width=10,height=8)
		plot(x=dist,y=d_means[1,2:(ncol(d_means)-2)],ylim=c((min_mean-min_sterr),(max_mean+max_sterr)),type='n',xlab="distance(log2)",ylab="obs/sum(obs)",xaxt='n',xlim=c(dist[findInterval(minDist,dist_f)+1],dist[findInterval(maxDist,dist_f)-1]))
		axis(1, at=axis_pos, labels = FALSE)
		text(axis_pos, par("usr")[3], labels = axis_labels, srt = 45,adj = c(1.1,1.1), xpd = TRUE)
		for (i in 1:nrow(d_means)){
			dispersion(x=t(dist),y=t(d_means[i,2:(ncol(d_means)-2)]),ulim=t(d_sterr[i,2:(ncol(d_means)-2)]),display.NA=FALSE,type='l',fill=as.vector(d_means$disp_colors)[i])
 			lines(dist,d_means[i,2:(ncol(d_means)-2)],col=as.character(d_means$colors[i]),lwd=1)
		}
		legend(x="topleft",d_means$type,cex=1.5,lwd=rep(2,length(cellsToPlot)),col=as.vector(d_means$colors))
		dev.off()

	} else if (log_base==10){

		a_log <- adply(a,1,function(x){
			if (norm.by.total.obs) {
				x[,2:(ncol(x)-3)] <- x[,2:(ncol(x)-3)]/sum(x[,2:(ncol(x)-3)],na.rm=T)
			}
			if (norm.by.binsize) {
				x[,2:(ncol(x)-3)] <- log10(x[,2:(ncol(x)-3)]/(10**dist * (sqrt(2)-1)))
			}

	 		coeff <- as.vector(as.numeric(x[,findInterval(coeff_range[1],dist_f):findInterval(coeff_range[2],dist_f)]))
			coeff <- coeff[!is.na(coeff)]
			is.na(coeff) <- sapply(coeff, is.infinite)
			if (sum(is.na(coeff)==0)){
	 			x$coeff <-   coef(lm ( coeff ~ log10(dist_f[findInterval(coeff_range[1],dist_f):findInterval(coeff_range[2],dist_f)])))[2]
			} else {
				x$coeff <- NA
			}
	#		x <- x[,c(1,idx,(ncol(x)-3):ncol(x))]
			return(x)
		})
    a_log$type <- factor(a_log$type,levels=cells)
		d_means = ddply(a_log, .(type), function(x) {
			means <- colMeans(x[,2:(ncol(x)-4)])
			return(data.frame(t(means),x$colors[1],x$disp_colors[1],mean(x$coeff,na.rm=T)))
		})
		colnames(d_means)[(ncol(d_means)-2):ncol(d_means)] <- c('colors','disp_colors','coeff')
		is.na(d_means) <- sapply(d_means,is.infinite)

		d_sterr = ddply(a_log, .(type), function(x) {
			temp <- x[,2:(ncol(x)-4)]
			return(data.frame(t(apply(temp,2,function(x) std.error(x))),x$colors[1],x$disp_colors[1],sd(x$coeff,na.rm=T)))
		})
		colnames(d_sterr)[(ncol(d_means)-2):ncol(d_means)] <- c('colors','disp_colors','coeff')
		is.na(d_sterr) <- sapply(d_sterr,is.infinite)

		max_mean <- max(d_means[,2:(ncol(d_means)-3)],na.rm=T)
		min_mean <- min(d_means[,2:(ncol(d_means)-3)],na.rm=T)
		max_sterr <- max(d_sterr[,2:(ncol(d_sterr)-3)], na.rm=T)
		min_sterr <- min(d_sterr[,2:(ncol(d_sterr)-3)], na.rm=T)

		axis_pos <- seq(min(dist),max(dist),by=1)
		axis_labels <- c(paste0(round_any(10^seq(min(dist)-0.01,5,by=1),1000)/1000," kb"),paste0(round_any(10^seq(6,(max(dist)),by=1),1000000)/1000000," mb"))

		pdf(out_f,width=8,height=8)
		plot(x=dist,y=d_means[1,2:(ncol(d_means)-3)],ylim=c((min_mean-min_sterr),(max_mean+max_sterr)),type='n',xlab="distance(log10)",ylab="Relative contact probability",xaxt='n',xlim=c(dist[findInterval(minDist,dist_f)+1],dist[findInterval(maxDist,dist_f)-1]))
		minor.ticks.axis(1,9,mn=3,mx=8)
		for (i in 1:nrow(d_means)){
	#		dispersion(x=t(dist),y=t(d_means[i,2:(ncol(d_means)-3)]),ulim=t(d_sterr[i,2:(ncol(d_means)-3)]),display.NA=FALSE,type='l',fill=as.vector(d_means$disp_colors)[i])
			dispersion(x=t(dist),y=t(d_means[i,2:(ncol(d_means)-3)]),ulim=t(d_sterr[i,2:(ncol(d_means)-3)]),display.NA=FALSE,type='a',col=as.vector(d_means$disp_colors)[i])
		#	lines(dist,d_means[i,2:(ncol(d_means)-3)],col=as.character(d_means$colors[i]),lwd=2)
		}
		for (i in 1:nrow(d_means)){
		  #		dispersion(x=t(dist),y=t(d_means[i,2:(ncol(d_means)-3)]),ulim=t(d_sterr[i,2:(ncol(d_means)-3)]),display.NA=FALSE,type='l',fill=as.vector(d_means$disp_colors)[i])
		 # dispersion(x=t(dist),y=t(d_means[i,2:(ncol(d_means)-3)]),ulim=t(d_sterr[i,2:(ncol(d_means)-3)]),display.NA=FALSE,type='a',col=as.vector(d_means$disp_colors)[i])
		  lines(dist,d_means[i,2:(ncol(d_means)-3)],col=as.character(d_means$colors[i]),lwd=2)
		}
	#	abline(v=log10(coeff_range[1]),lty=5,lwd=2,col='gray80')
	#	abline(v=log10(coeff_range[2]),lty=5,lwd=2,col='gray80')
	#	legend(x="topright",paste0(d_means$type,' coeff =',round(d_means$coeff,2),'±',round(d_sterr$coeff,2)),cex=1,lwd=rep(2,length(cellsToPlot)),col=as.vector(d_means$colors))
		legend(x="topright",gsub('E14_','',d_means$type),cex=1.5,lwd=rep(2,length(cellsToPlot)),col=as.vector(d_means$colors),seg.len = 1)
		
		dev.off()
	}
}

##############################################################################################################################################
#' Calculates the cis decay profiles per compartment based on A-B definition. Data can then be plotted using
#'
#' \code{comp_cisDecay}
#'
#' This function will calculate the distribution of contacts across range of distances separately for A-A, B-B and A-B pairs of regions. Compartments have to be defined and saved as gtrack where compartment identity is given in the cluster column - A domains are designated as 2,and B compartments - as 1.
#'
#' @param comp_tracks Which hic tracks to work on. Defaults to all tracks listed in the config file. Only tracks for the selected conditions will be further considered.
#' @param cells Which condition to work on. Vector of conditions is also accepted. Defaults to conditions listed in the config file.
#' @param path Where to store the data.
#' @param ins Insulation window which has been used to identify compartments. \code{analyzeCompartments} must have been run before using this insulation window parameter.
#' @param file_f Name of the file which stores the data. It will be saved in the directory path defined above.
#' @param minDist Minimum distance to consider when counting contacts. Default value - 1000.
#' @param maxDist Maximum distance to consider when counting contacts. Default value - 2e8. It is recommended to extract the data using the length of the biggest chromosome as a guide for this value. It can be later adjusted for plotting.
#'
#' @examples
#'
#' comp_cisDecay(cells=cells,path,ins=250,minDist=4e6,maxDist=2e8)    # Correlation between eigenvector value and chipseq signal enrichment in 100kb bins.
#'
#' @export
##########################################################################################################
comp_cisDecay <-function(comp_tracks,cells,path,ins=250,file_f,minDist=4e6,maxDist=2e8){
	log_base=10
	file_f <- paste0(path,file_f)
	comp_tracks <- as.vector(unlist(sapply(cells,function(x){comp_tracks[grep(x,comp_tracks)]})))
	start_n <- round(log(minDist) / log(log_base),2)
	end_n <- round(log(maxDist) / log(log_base),2)
	dist_f = log_base^seq(start_n, end_n, by=0.2/log_base)
	dist = seq(start_n, end_n, by=0.2/log_base)[-1]
	cis_decay <- list()
	for (cell in cells){
		name <- paste0(cell,'_compartments')
		domains <- paste0("hic.",cell,".ins",ins,"_k2_domains")
		if (!gintervals.exists(domains)) { next}
		domain <- gintervals.load(domains)
		intervals1 <- domain[domain$cluster==2,]
		intervals2 <- domain[domain$cluster==1,]
		AvsA <- construct.grid.comp(intervals1,intervals1,minDist,maxDist)
		AvsB <- construct.grid.comp(intervals1,intervals2,minDist,maxDist)
		BvsA <- construct.grid.comp(intervals2,intervals1,minDist,maxDist)
		AvsB <- rbind(AvsB,BvsA)
		BvsB <- construct.grid.comp(intervals2,intervals2,minDist,maxDist)
		tracks = comp_tracks[grep(cell,comp_tracks)]
		shuffled_tracks = paste0(tracks,'_shuffle')
		if (length(tracks)!=length(shuffled_tracks)){
			message('Check your tracks input - something went wrong with:', cell)
			next
		}
		temp_list <- list()
		for (i in 1:length(tracks)){
			obs_cis_decay1 <- as.data.frame(gcis_decay(tracks[i], log_base^dist, gintervals.all(), domain, AvsA))
			obs_cis_decay2 <- as.data.frame(gcis_decay(tracks[i], log_base^dist, gintervals.all(), domain, AvsB))
			obs_cis_decay3 <- as.data.frame(gcis_decay(tracks[i], log_base^dist, gintervals.all(), domain, BvsB))
			obs_cis_decay <- as.data.frame(cbind(obs_cis_decay1$inter,obs_cis_decay2$inter,obs_cis_decay3$inter))
			colnames(obs_cis_decay) <- c("AvsA",'AvsB','BvsB')
			obs_cis_decay$distance <- row.names(obs_cis_decay1)

			exp_cis_decay1 <- as.data.frame(gcis_decay(shuffled_tracks[i], log_base^dist, gintervals.all(), domain, AvsA))/2
			exp_cis_decay2 <- as.data.frame(gcis_decay(shuffled_tracks[i], log_base^dist, gintervals.all(), domain, AvsB))/2
			exp_cis_decay3 <- as.data.frame(gcis_decay(shuffled_tracks[i], log_base^dist, gintervals.all(), domain, BvsB))/2
			exp_cis_decay <- as.data.frame(cbind(exp_cis_decay1$inter,exp_cis_decay2$inter,exp_cis_decay3$inter))
			colnames(exp_cis_decay) <- c("AvsA",'AvsB','BvsB')
			exp_cis_decay$distance <- row.names(exp_cis_decay1)

			merged <- merge(obs_cis_decay,exp_cis_decay,by="distance",sort=FALSE)
			colnames(merged) <- c("distance","obs_AvsA","obs_AvsB","obs_BvsB","exp_AvsA","exp_AvsB","exp_BvsB")
			temp_list[[i]] <- merged
			names(temp_list)[i] <- unlist(strsplit(tracks[i],'\\.'))[2]
		}
		cis_decay[[cell]] <- temp_list
		cis_decay[['dist']] <- dist_f
		save(cis_decay,file=file_f)
	}
}

##############################################################################################################################################
#' Identifies conserved and differential domain boundaries across conditions based on insulation score.
#'
#' \code{analyzeInsulation}
#'
#' This function uses the insulation score per condition to identify domain boundaries considering replicates. In addition it determines conserved and differential boundaries given more than one conditions. Output can be insulation tracks as bedGraph files and borders as text files.
#'
#' @param tracks Which combined insulation tracks to work on. Defaults to all insulation tracks listed in the config file. Only tracks for the selected conditions will be further considered.
#' @param rep_tracks Insulation tracks per replicate to work on. Defaults to all insulation tracks per replicate listed in the config file. Only tracks for the selected conditions will be further considered.
#' @param cells Which conditions to work on. Vector of conditions is also accepted. Defaults to conditions listed in the config file.
#' @param border_max_distance Boundaries will be merged into one interval if separated by less than this distance.
#' @param path Where to save the resulting data and files.
#' @param ins_dom_thresh Insulation threshold to consider regions as boundaries. Use 0.1 for top 10 percent.
#' @param domain_size_range Vector of min and max size. All detected domains smaller or bigger than this will be discarded.
#' @param sig_thresh Quantile to consider the difference between max and min insulation as significant. For example for value of 0.8 only top 20 percent of the borders based on difference in insulation will be considered as significantly different.
#' @param no_sig_thresh Similar to sig_thresh but to determine a set of conserved boundaries. 0.1 will result in the returning the bottom 10 percent of the borders, where differences in insulation are minimal.
#' @param use.normFACS Should normalization based on genome average insulation be applied. Set to true if there are substantial difference in the average insulation across the whole genome (for example due to cell cycle differences).
#' @param use.sigThresh Boolean to use sig_thresh(respectively no_sig_thresh) as defined above. If FALSE all intervals are used.
#' @param write.insTracks Boolean whether to write the insulation tracks (normalized if use.normFACS=T) as bedGraph files.
#' @param write.borders Boolean whether to write the boundaries per condition as bed files.
#' @param anyDiff Boolean variable. If FALSE it will consider only regions where the difference between min and max insulation is across boundary score and is within the specified percentile given by sig_thresh (more stringent). Otherwise, the difference between min and max insulation has to be within the specified percentile given by sig_thresh but doesnt have to be across boundary score (less stringent).
#'
#' @return Returns a list with differential and conserved boundaries.
#'
#' @examples
#'
#' res <- analyzeInsulation(tracks=ins_tracks,rep_tracks=ins_rep_tracks,path=paste0(main_f,'analysis/insulation/'),cells=c('D0','D2','D4','D6','D10'),ins_dom_thresh=0.1,border_max_distance=1e4,sig_thresh=0.8,no_sig_thresh=0.2,use.normFACS=TRUE,use.sigThresh=TRUE,write.insTracks=T,write.borders=T,anyDiff=T)
#'
#' @export
##########################################################################################################
analyzeInsulation <- function(tracks,rep_tracks,cells,border_max_distance=1e4,path=path,ins_dom_thresh=0.1,domain_size_range=c(5e4, 4e6),sig_thresh=0.8,no_sig_thresh=0.1,use.normFACS=TRUE,use.sigThresh=TRUE,write.insTracks=F,write.borders=F,anyDiff=F){
	borders = c()
	tracks <- as.vector(unlist(sapply(cells,function(x){tracks[grep(x,tracks)]})))
	rep_tracks <- as.vector(unlist(sapply(cells,function(x){rep_tracks[grep(x,rep_tracks)]})))
	threshs <- as.vector(unlist(sapply(cells,function(x){ gquantiles(tracks[grep(x,tracks)],ins_dom_thresh) })))
	v_tracks <- paste("v_", tracks, sep="")
	for (i in 1:length(v_tracks)) {
		 gvtrack.create(v_tracks[i], tracks[i], 'min')
	}

	v_rep_tracks <- paste("v_", rep_tracks, sep="")
	for (i in 1:length(v_rep_tracks)) {
		 gvtrack.create(v_rep_tracks[i], rep_tracks[i], 'min')
	}

	rep_stats <- list()
	number_matrix <- matrix(nrow=length(tracks),ncol=nrow(sapply(cells,function(x){rep_tracks[grep(x,rep_tracks)]})))
	size_matrix <- matrix(nrow=length(tracks),ncol=nrow(sapply(cells,function(x){rep_tracks[grep(x,rep_tracks)]})))
	row.names(number_matrix) <- cells
	row.names(size_matrix) <- cells
	for (cell in cells) {
		i=grep(cell,cells)
		ins_track <- tracks[grep(cell,tracks)]
		ins_rep_tracks <- rep_tracks[grep(cell,rep_tracks)]
		for (j in 1:length(ins_rep_tracks)){
			rep_track <- ins_rep_tracks[j]
			rep_borders = get_borders(ins_track=rep_track,ins_dom_thresh=ins_dom_thresh,extend_border_region=border_max_distance)
			ins <- gextract(paste0('v_',rep_track),rep_borders,iterator=1000)
			ins <- ins[complete.cases(ins),]
			rep_borders <- ddply(ins,.(intervalID),function(x){
			return(x[which.min(x[,4]),])
			})
			is.na(rep_borders) <- sapply(rep_borders,is.infinite)
			rep_borders <- rep_borders[complete.cases(rep_borders),]
			domains <- gintervals.diff(gintervals.all(),rep_borders)
			domains$len = domains$end - domains$start
			domains = domains[domains$len >= domain_size_range[1], ]
			rep_stats[[rep_track]]$domains <- domains
			rep_stats[[rep_track]]$borders <- borders
			number_matrix[i,j] <- nrow(domains)
			size_matrix[i,j] <- round(mean(domains$len,na.rm=T),0)
		}
		track_borders = get_borders(ins_track=ins_track,ins_dom_thresh=ins_dom_thresh,extend_border_region=border_max_distance)
		message("Track:",ins_track,", initial N borders:",nrow(track_borders))
		ins <- gextract(paste0('v_',ins_track),track_borders,iterator=1000)
		ins <- ins[complete.cases(ins),]
		track_borders <- ddply(ins,.(intervalID),function(x){
			return(x[which.min(x[,4]),])
		})
		is.na(track_borders) <- sapply(track_borders,is.infinite)
		track_borders <- track_borders[complete.cases(track_borders),]
		track_borders$source = ins_track
#		rep_ins <- gextract(paste0('v_',ins_rep_tracks),track_borders,iterator=track_borders)
#		rep_ins <- rep_ins[complete.cases(rep_ins),]
#		quantiles <- as.vector(sapply(ins_rep_tracks,function(x) gquantiles(x, ins_dom_thresh)))
#		track_borders <- adply(rep_ins,1,function(x){
#			if(sum(x[,4:(ncol(x)-1)]<quantiles)==2) {return(x)}
#		})
#		message("Track:",ins_track,", N borders in at least 2 replicates:",nrow(track_borders))
		name_d <- paste0(ins_track,'_borders')
		if (write.borders) { write.table(track_borders[,1:3],paste0(path,cell,'_borders.bed'),col.names=F,row.names=F,quote=F,sep='\t')}
		if (!gintervals.exists(name_d)) {gintervals.save(name_d,track_borders)}
		name_d <- paste0(ins_track,'_domains_expanded')
		domains <- gintervals.diff(gintervals.all(),track_borders)
		domains$len = domains$end - domains$start
		domains = domains[domains$len >= domain_size_range[1], ]
		if (!gintervals.exists(name_d)) {gintervals.save(name_d,domains)}
		borders = rbind(borders[,1:3], track_borders[,1:3])
	}
	write.table(size_matrix,paste0(path,'Domain_size.txt'),quote=F,sep='\t',row.names=T,col.names=T)
	write.table(number_matrix,paste0(path,'Domain_number.txt'),quote=F,sep='\t',row.names=T,col.names=T)

	bb = gintervals.neighbors(borders, borders, mindist=1000, maxdist=border_max_distance, maxneighbors=100, na.if.notfound=TRUE)
	b = gintervals.canonic(rbind(bb[is.na(bb$dist),1:3],
	gintervals(bb[!is.na(bb$dist),1],
	apply(bb[!is.na(bb$dist),grep("start", colnames(bb))], 1, min),
	apply(bb[!is.na(bb$dist),grep("end", colnames(bb))], 1, max))))
	borders = intervals.normalize(b, 2000)

	allScores <- gextract(v_tracks, gintervals.all(), iterator=2000)
	allScores <-allScores[complete.cases(allScores),]
	normFACS <- mean(colMeans(allScores[,grep('ins',colnames(allScores))]))/colMeans(allScores[,grep('ins',colnames(allScores))])
	normFACS_df <- matrix(NA,ncol=length(cells),nrow=nrow(allScores))
	for (i in 1:length(normFACS)){normFACS_df[,i] <- normFACS[i]}
	if (use.normFACS) {allScores[,grep('ins',colnames(allScores))] <- allScores[,grep('ins',colnames(allScores))]*normFACS_df}
	if (write.insTracks){
		for (cell in cells){
			temp_df <- allScores
			temp_df[,grep('ins',colnames(temp_df))] <- temp_df[,grep('ins',colnames(temp_df))]*(-1)
			idx <- grep(cell,colnames(temp_df))
			write.table(temp_df[,c(1:3,idx)],paste0(path,cell,'_ins.bedGraph'),col.names=F,row.names=F,quote=F,sep='\t')
		}
	}
	borderScores <- gintervals.neighbors(borders,allScores)[,-c(4:6)]
	borderScores <- borderScores[,-c(ncol(borderScores)-1,ncol(borderScores))]

	variable_insulation <- adply(borderScores,1,function(x){
		scoreDiff <- max(x[,4:ncol(x)])-min(x[,4:ncol(x)])
		return(scoreDiff)
	})
	if (use.sigThresh){
		sig_thresh <- quantile(variable_insulation$V1,sig_thresh)
		no_sig_thresh <- quantile(variable_insulation$V1,no_sig_thresh)
		pdf(paste0(path,'diffInsulation.pdf'),width=5,height=5)
		plot(density(variable_insulation$V1,na.rm=T),lwd=2,main='Max(insulation) - min(insulation)')
		abline(v=sig_thresh,col='red',lty=5)
		abline(v=no_sig_thresh,col='blue',lty=5)
		dev.off()
	} else {
		sig_thresh <- 0
		no_sig_thresh <- 0
	}

	variable_ins = gscreen(paste0("max(",paste0('v_',tracks,collapse=','),')-min(',paste0('v_',tracks,collapse=','),') >= ',sig_thresh), borders, iterator=borders)
	variable_insulation_2k = gextract(paste0('v_',tracks), variable_ins, iterator=variable_ins)
	invariable_ins = gscreen(paste0("max(",paste0('v_',tracks,collapse=','),')-min(',paste0('v_',tracks,collapse=','),') < ',no_sig_thresh), borders, iterator=borders)
	invariable_insulation_2k = gextract(paste0('v_',tracks), invariable_ins, iterator=invariable_ins)

	temp <- variable_insulation_2k[,-c(1:3,ncol(variable_insulation_2k))]
	idx <- temp
	for (j in 1:ncol(temp)){
		if (anyDiff) {
			idx[,j] <- (temp[,j]<threshs[j]) & (adply(temp[,-j],1,function(x) { return(any(x>threshs[-j]))})$V1)
		} else {
			idx[,j] <- (temp[,j]<threshs[j]) & (rowMins(as.matrix(temp[,-j])) > min(threshs[-j]))
		}
	}
	variable_insulation_2k <- variable_insulation_2k[rowSums(idx)>0,]
	variable_insulation_2k[,grep('ins',colnames(variable_insulation_2k))] <- variable_insulation_2k[,grep('ins',colnames(variable_insulation_2k))]*normFACS_df
	res <- list(differential=variable_insulation_2k,constant=invariable_insulation_2k)
	message('N differential: ',nrow(variable_insulation_2k),' . N constant: ',nrow(invariable_insulation_2k))
	return(res)
}

##############################################################################################################################################
#' Analyzing differential or conserved boundaries based on gene expression.
#'
#' \code{analyzeDiffBorders}
#'
#' This function determines the relationship between gene expression and a set of boundaries as identified by \code{analyzeInsulation}. It returns the p.value for the difference in gene expression per cluster.
#'
#' @param input data frame generated by \code{analyzeInsulation}. Should be one of "differential" or "constant" elements of the list.
#' @param tss gtrack defining promoter coordinates. Usually defined in the configuration file.
#' @param minTss_dist Maximum distance between the boundary region and the gene promoter to be considered as overlapping.
#' @param expression_f Expression file containing gene expression information. Format should be a tab-separated file where each row is a separate gene and the columns are individual replicates. Column names should match the same names defined in the config file.
#' @param path Where to save the resulting data and files.
#' @param cells Conditions to work on. Should be same as in the call to \code{analyzeDiffBorders}.
#' @param scale Boolean whether to convert absolute insulation score values into Z-score.
#' @param clustering Clustering method to use. One of "kmeans" or "hierarchical".
#' @param method Clustering method to use. Same as pheatmap.
#' @param K Number of desired clusters for Kmeans or as input for cut tree.
#' @param write.borders Output files containing border coordinates and genes located in close proximity together with cluster identity.
#' @param cell_colors Vector of colors to use for conditions when plotting. Length should be same as the length of the conditions.
#'
#' @examples
#'
#' analyzeDiffBorders(input=res$differential,expression_f='/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt',path=paste0(main_f,'analysis/insulation/'),cells=c('D0','D2','D4','D6','D10'),scale=F,clustering='hierarchical',method='ward.D2',K=4,write.borders=T)
#'
#' @export
##########################################################################################################
analyzeDiffBorders <- function(input=res,tss=gintervals.load(tss_f),minTss_dist=1e4,expression_f='',path=paste0(main_f,'analysis/insulation/'),cells=c('D0','D2','D4','D6','D10'),scale=T,clustering='hierarchical',method='ward.D',K=2,write.borders=F,cell_colors=rainbow(length(cells))){
	#require(beanplot)
	peaks_df <- input[,1:3]
	if (clustering=='kmeans'){
		if (scale) {
			input_clusters = cluster_and_plot(input[,4:(ncol(input)-1)]-rowMeans(input[,4:(ncol(input)-1)]), paste0(path,"t"), paste0(path,'insulation_kmeans_k',K,'.png'), K, image_to_file=TRUE, zlim=c(-1, 1), fig.height=1600, fig.width=800, shades=colorRampPalette(c("red","white","blue")), shade_count=100)
		} else {
			input_clusters = cluster_and_plot(input[,4:(ncol(input)-1)]*(-1), paste0(path,"t"), paste0(path,'insulation_kmeans_k',K,'.png'), K, image_to_file=TRUE, zlim=c(1.5, 3.5), fig.height=1600, fig.width=800, shades=colorRampPalette(c("blue","white","red")), shade_count=100)
		}
		peaks_df$clust = input_clusters$clust
	} else {
		res <- pheatmap(-input[,4:(ncol(input)-1)],scale=ifelse(scale,'row','none'),cutree_rows=K,border_color=NA,color = colorRampPalette(c("blue", "white", "red"))(50),cluster_rows=T,cluster_cols=F,show_rownames=F,show_colnames=T,clustering_method=method,labels_col=cells,legend_labels='Insulation score',filename=paste0(path,"/insulation_hierarchical.pdf"),width=4,height=10)
		peaks_df$clust <- cutree(res$tree_row, k = K)
	}

	peaks_df_inters <- gintervals.neighbors(tss,peaks_df)
	tss_peaks <- peaks_df_inters[abs(peaks_df_inters$dist)<=minTss_dist,]
	peaks_df_inters <- gintervals.neighbors(peaks_df,tss)
	notTss_peaks <- peaks_df_inters[abs(peaks_df_inters$dist)>minTss_dist,]

	input$clust <- peaks_df$clust
	par(mfrow=c(K,1))
	for (i in c(rev(seq(1,K,by=1)))){
		pdf(paste0(path,"/cluster_boxIns_k",i,".pdf"),width=5,height=5)
		x <- input[input$clust==i,]
		x <- x[,grep('ins',colnames(x))]*(-1)
		boxplot(x,col=cell_colors,names=cells,notch=FALSE,outline=FALSE,ylab='Insulation Score',main=paste0('N: ',nrow(x)))
		points(1:length(cells),colMeans(log2(x+0.001),na.rm=T),pch='*',col='darkgrey',lwd=4)
		dev.off()
	}

#	fpkm_f <- read.table(expression_f,header=T)
#	clust_colors=rainbow(K)

#	fpkm <- merge(fpkm_f,tss_peaks[,c('geneName','clust')],by.x='gene_name',by.y='geneName')

	### Write borders
	if (write.borders){
		for (i in 1:K){
			write.table(tss_peaks[tss_peaks$clust==i,1:3],paste0(path,'cluster',i,'_tss.bed'),quote=F,row.names=F,col.names=F,sep='\t')
			write.table(peaks_df[peaks_df$clust==i,1:3],paste0(path,'cluster',i,'_all.bed'),quote=F,row.names=F,col.names=F,sep='\t')
			write.table(notTss_peaks[peaks_df$clust==i,1:3],paste0(path,'cluster',i,'_noTss.bed'),quote=F,row.names=F,col.names=F,sep='\t')
		}
		gene_inters <- gintervals.neighbors(tss,input)
		gene_inters <- gene_inters[abs(gene_inters$dist)<=minTss_dist,]
		gene_inters <- cbind(gene_inters[,c('geneName','clust')],gene_inters[,grep('ins',colnames(gene_inters))]*(-1))
		write.table(gene_inters,paste0(path,'diffBorder_genes.txt'),quote=F,row.names=F,col.names=F,sep='\t')
	}
	#####


	fpkm_df <- fpkm[,-1]
	fpkm_means <- as.data.frame(matrix(NA,ncol=length(cells),nrow=nrow(fpkm_df)))
	colnames(fpkm_means) <- cells
	for (i in 1:length(cells)){
		fpkm_means[,i] <- rowMeans(fpkm_df[,grep(cells[i],colnames(fpkm_df))],na.rm=T)
		}
	fpkm_means$clust=fpkm_df$clust
 	for (i in 1:K){
 		pdf(paste0(path,"/cluster_boxTss_k",i,".pdf"),width=5,height=5)
 		x <- fpkm_means[fpkm_means$clust==i,]
 		x <- x[,-ncol(x)]
 		boxplot(log2(x+0.001),col=cell_colors,names=cells,notch=FALSE,outline=FALSE,ylab='log2(FPKM+1)',main=paste0('N: ',nrow(x)))
 		points(1:length(cells),colMeans(log2(x+0.001),na.rm=T),pch='*',col='darkgrey',lwd=4)
 		dev.off()
 	}

	for (i in 1:K){
	pdf(paste0(path,"/cluster_beanplot_k",i,".pdf"),width=5,height=5)
	x <- fpkm_means[fpkm_means$clust==i,]
	x <- x[,-ncol(x)]
	beanplot(x+0.01,what=c(1,1,1,0),log='y', side='no', col=as.list(cell_colors),names=cells,ylab='FPKM',main=paste0('N: ',nrow(x)))
	dev.off()
	}

	stat_test <- as.data.frame(matrix(NA,nrow=K,ncol=2))
	colnames(stat_test) <- c(paste0(cells[1],'vs',cells[2]),paste0(cells[1],'vs',cells[length(cells)]))
	for (i in 1:K){
	x <- fpkm_means[fpkm_means$clust==i,]
		stat_test[i,1] <- wilcox.test(x[,1],x[,2],paired=T)$p.value
		stat_test[i,2] <- wilcox.test(x[,1],x[,length(cells)],paired=T)$p.value
	}
	return(stat_test)
}

##############################################################################################################################################
#' Calculates the average contact enrichment per TAD.
#'
#' \code{averageTAD}
#'
#' This function operated in two modes specified by the parameter stats_f. If true it will calculate the contact enrichment inside and outside TADs as shown in Bonev et al., 2017 - Figure S3A. If False it will extend each TAD by its length, split it in 100 bins and then calculate the average contact enrichment per bin to determine averageTAD profile as Bonev et al., 2017 - Figure 3A.
#'
#' @param tracks Which combined insulation tracks to work on. Defaults to all insulation tracks listed in the config file. Only tracks for the selected conditions will be further considered.
#' @param cells Which conditions to work on. Vector of conditions is also accepted. Defaults to conditions listed in the config file.
#' @param domains If not NULL, it needs to contain gintervals holding domain intervals with an additional column called cluster where domain type is specified, otherwise if NULL - domains that have been previously generated using \code{analyzeCompartments} will be used.
#' @param domain_size A vector of minimum and max size. Domains smaller or bigger than this will be discarded.
#' @param file_f Name of the output file.
#' @param path Folder where to store resulting data.
#' @param bins How many bins to split the extended domains. Use smaller number for sparse data and bigger number for high-resolution data. Recommended 50/100.
#' @param stats_f Controls whether to generate data to plot averageTAD (if FALSE) or calculate contact enrichment inside and outside domains (if TRUE)
#' @param min_dist only consider contacts separated by at least this distance. Set to at least 1kb to avoid digestion artefacts.
#'
#' @examples
#'
#' averageTAD(tracks=all_tracks,cells=cells,domains='sexton_TADs',domain_size=c(1e4,2e6),file_f='eGFP_SextonTADs',path=paste0(main_f,'analysis/averageTAD/'),bins=100,stats_f=F)    # Generate an average domain heatmap representation using an external bed file with domain coordinates and annotation
#' averageTAD(tracks=all_tracks,cells=cells,domains=NULL,domain_size=c(1e4,2e6),file_f='eGFP_SextonTADs',path=paste0(main_f,'analysis/averageTAD/'),bins=100,stats_f=F)    # Generate an average domain heatmap representation using an domain classification based on insulation score and A-B compartment classification generated by \code{analyzeCompartments}
#' averageTAD(tracks=all_tracks,cells=cells,domains='sexton_TADs',domain_size=c(1e4,2e6),file_f='eGFP_SextonTADs',path=paste0(main_f,'analysis/averageTAD/'),bins=100,stats_f=T)    # Quantifies the contact enrichment inside versus outside domains
#'
#' @export
##########################################################################################################
averageTAD <- function(tracks=all_tracks,cells,domains=NULL,domain_size=c(1e5,2e6),file_f,path,bins,stats_f=T,min_dist=1e4,score_tracks){
	max_dist <- 1e8
	tracks <- as.vector(unlist(sapply(cells,function(x){tracks[grep(x,tracks)]})))
	chrom_sizes <- gintervals.all()
	for (t in 1:length(tracks)){
	gvtrack.create(paste0('v_',tracks[t]), tracks[t], "area")
	}
	shuffled_tracks <- paste0(tracks,'_shuffle')
#	shuffled_tracks <- paste0(tracks,'_shuffle_500Small_1000High')
	for (t in 1:length(tracks)){
	gvtrack.create(paste0('v_',shuffled_tracks[t]), shuffled_tracks[t], "area")
	}
	#	shuffled_tracks <- paste0(tracks,'_shuffle_500Small_1000High')
	if(!is.null(score_tracks)){
	  for (t in 1:length(score_tracks)){
	    gvtrack.create(paste0('v_',score_tracks[t]), score_tracks[t], "avg")
	  }
	}
	file_f <- paste0(path,file_f)

	if(!stats_f){
		tad_list <- list()
		for (cell in cells){
			temp_list <- list()
	#		if (!gintervals.exists(domains[grep(cell,domains)])) { next}
			if (is.null(domains)) {
				domain <- paste0("hic.",cell,".ins250_k2_domains")
				domain <- gintervals.load(domain)
			} else {
				if (!grepl('bed',domains)) {
					domain <- gintervals.load(domains)
				} else {
					domain <- read.table(domains,sep='\t',header=F)
					domain <- domain[domain[,1]%in% gintervals.all()$chrom,]
					domain <- gintervals(domain[,1],domain[,2],domain[,3])
				}
			}
			domain$len <- domain[,3]-domain[,2]
			domain <- domain[domain$len>=domain_size[1]&domain$len<=domain_size[2],]
			domain <- domain[domain$start>1e6,]
			chr_ends <- ddply(domain,.(chrom),function(x){return(max(x$end))})
			domain <- domain[!(domain$end %in% chr_ends$V1),]
			for (i in 1:nrow(domain)){
				interval <- domain[i,]
				interval$start = round_any(interval$start - interval$len,1000)
				interval$end = round_any(interval$end + interval$len,1000)
				interval <- gintervals.force_range(interval)
				bin <- (interval$end-interval$start)/bins
				interval2d <- gintervals.2d(interval[,1],interval[,2],interval[,3],interval[,1],interval[,2],interval[,3])
				bins_2d <- giterator.intervals(intervals=interval2d, iterator=c(bin,bin), band=c(-max_dist,-0))
				df <- gextract(paste0('v_',tracks[grep(cell,tracks)],collapse='+'),intervals=bins_2d,iterator=bins_2d)
				df[,7]=df[,7]/sum(df[,7],na.rm=T)   #distribution
				df_exp <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)],collapse='+'),')/2'),intervals=bins_2d,iterator=bins_2d)
				df_exp[,7]=df_exp[,7]/sum(df_exp[,7],na.rm=T)
				df_exp[,7] <- log2(df[,7]/df_exp[,7])
				df <- df[,-ncol(df)]
				df$bin1 <- df[,2]/bin
				df$bin2 <- df[,5]/bin
				df <- df[,7:ncol(df)]
				colnames(df)[1] <- 'count'
				df_cast <- dcast(df,factor(bin1,levels=unique(bin1))~factor(bin2,levels=unique(bin2)),value.var='count')
				colnames(df_cast) <- NULL
				df_cast <- df_cast[,-1]
				df_cast <- data.matrix(df_cast)
				diag(df_cast)=0
				is.na(df_cast) <- sapply(df_cast, is.infinite)
				df_cast <- df_cast[-(bins+1),-(bins+1)]
				temp_list$obs[[i]] <- as.matrix(df_cast)
##### Score calculations #####
				if(!is.null(score_tracks)){
				  df <- gextract(paste0('v_',score_tracks[grep(cell,score_tracks)]),intervals=bins_2d,iterator=bins_2d)
				  df <- df[,-ncol(df)]
				  df$bin1 <- df[,2]/bin
				  df$bin2 <- df[,5]/bin
				  df <- df[,7:ncol(df)]
				  colnames(df)[1] <- 'count'
				  df_cast <- dcast(df,factor(bin1,levels=unique(bin1))~factor(bin2,levels=unique(bin2)),value.var='count')
				  colnames(df_cast) <- NULL
				  df_cast <- df_cast[,-1]
				  df_cast <- data.matrix(df_cast)
				  diag(df_cast)=0
				  is.na(df_cast) <- sapply(df_cast, is.infinite)
				  df_cast <- df_cast[-(bins+1),-(bins+1)]
				  temp_list$score[[i]] <- as.matrix(df_cast)
				}
#######################				
				df_exp <- df_exp[,-ncol(df_exp)]
				df_exp$bin1 <- df_exp[,2]/bin
				df_exp$bin2 <- df_exp[,5]/bin
				df_exp <- df_exp[,7:ncol(df_exp)]
				colnames(df_exp)[1] <- 'count'
				df_exp_cast <- dcast(df_exp,factor(bin1,levels=unique(bin1))~factor(bin2,levels=unique(bin2)),value.var='count')
				colnames(df_exp_cast) <- NULL
				df_exp_cast <- df_exp_cast[,-1]
				df_exp_cast <- data.matrix(df_exp_cast)
				diag(df_exp_cast)=0
				is.na(df_exp_cast) <- sapply(df_exp_cast, is.infinite)
				df_exp_cast <- df_exp_cast[-(bins+1),-(bins+1)]
				temp_list$norm[[i]] <- as.matrix(df_exp_cast)
			}
			tad_list[[cell]] <- temp_list
		}
		save(tad_list,file=file_f)
	} else {
	#### Combined calculation
		tad_list <- list()
		for (cell in cells){
#			if (!gintervals.exists(domains[grep(cell,domains)])) { next}
			if (is.null(domains)) {
				domain <- paste0("hic.",cell,".ins250_k2_domains")
				domain <- gintervals.load(domain)
			} else {
				if (!grepl('bed',domains)) {
					domain <- gintervals.load(domains)
					if(length(domain$cluster)==0){domain$cluster=1}
				} else {
					domain <- read.table(domains,sep='\t',header=F)
					domain <- domain[domain[,1]%in% gintervals.all()$chrom,]
					domain <- gintervals(domain[,1],domain[,2],domain[,3])
					if(length(domain$cluster)==0){domain$cluster=1}
				}
			}
			domain$len <- domain[,3]-domain[,2]
			domain <- domain[domain$len>=domain_size[1]&domain$len<=domain_size[2],]
			domain <- domain[domain$start>1e6,]
			domain <- domain[!(domain$end %in% gintervals.all()$end),]
			chr_ends <- ddply(domain,.(chrom),function(x){return(max(x$end))})
			x <- domain[!(domain$end %in% chr_ends$V1),]
			interval2d <- gintervals.2d(x[,1],x[,2],x[,3],x[,1],x[,2],x[,3])
			df_obs <- gextract(paste0('v_',tracks[grep(cell,tracks)]),intervals=interval2d,iterator=interval2d,band=-c(domain_size[2],10000))
			df_exp <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)]),')/2'),intervals=interval2d,iterator=interval2d,band=-c(domain_size[2],10000))
			df_score <- gextract(paste0('v_',score_tracks[grep(cell,score_tracks)]),intervals=interval2d,iterator=interval2d,band=-c(domain_size[2],10000))
			interval2d_1 <- gintervals.2d(x[,1],x[,2],x[,3],x[,1],x[,2],x[,3])
			interval2d_1$start2 <- x$start+round(x$len/10)
			interval2d_1$end2 <- round((x$start+x$end)/2)
			interval2d_1$start1 <- x$start-x$len
			interval2d_1$end1 <- x$start-round(x$len/10)
			interval2d_1 <- gintervals.force_range(interval2d_1)
			df_obs1 <- gextract(paste0('v_',tracks[grep(cell,tracks)]),intervals=interval2d_1,iterator=interval2d_1,band=-c(domain_size[2]*2,10000))
			df_score1 <- gextract(paste0('v_',score_tracks[grep(cell,score_tracks)]),intervals=interval2d_1,iterator=interval2d_1,band=-c(domain_size[2]*2,10000))
			df_exp1 <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)]),')/2'),intervals=interval2d_1,iterator=interval2d_1,band=-c(domain_size[2]*2,10000))

			interval2d_2 <- gintervals.2d(x[,1],x[,2],x[,3],x[,1],x[,2],x[,3])
			interval2d_2$start1 <- round((x$start+x$end)/2)
			interval2d_2$end1 <- x$end-round(x$len/10)
			interval2d_2$start2 <- x$end+round(x$len/10)
			interval2d_2$end2 <- x$end+x$len
			interval2d_2 <- gintervals.force_range(interval2d_2)
			df_obs2 <- gextract(paste0('v_',tracks[grep(cell,tracks)]),intervals=interval2d_2,iterator=interval2d_2,band=-c(domain_size[2]*2,10000))
			df_score2 <- gextract(paste0('v_',score_tracks[grep(cell,score_tracks)]),intervals=interval2d_2,iterator=interval2d_2,band=-c(domain_size[2]*2,10000))
			df_exp2 <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)]),')/2'),intervals=interval2d_2,iterator=interval2d_2,band=-c(domain_size[2]*2,10000))
			message(dim(df_obs1),' ',dim(df_obs2),' ',dim(df_exp1),' ',dim(df_exp2))
			df_obs1[,7:(ncol(df_obs1)-1)] <- df_obs1[,7:(ncol(df_obs1)-1)]+df_obs2[,7:(ncol(df_obs2)-1)]
			df_score1[,7:(ncol(df_score1)-1)] <- df_score1[,7:(ncol(df_score1)-1)]+df_score2[,7:(ncol(df_score2)-1)]
			df_exp1[,7:(ncol(df_exp1)-1)] <- df_exp1[,7:(ncol(df_exp1)-1)]+df_exp2[,7:(ncol(df_exp2)-1)]
			df_out <- cbind(x$cluster,df_obs[,7:(ncol(df_obs)-1)],df_exp[,7:(ncol(df_exp)-1)],df_obs1[,7:(ncol(df_obs1)-1)],df_exp1[,7:(ncol(df_exp1)-1)])
			df_score_out <- cbind(x$cluster,df_score[,7:(ncol(df_score)-1)],df_score1[,7:(ncol(df_score1)-1)])
			tad_list[[cell]]$res <- df_out
			tad_list[[cell]]$score_res <- df_score_out
		}
		file_f <- paste0(file_f,'_quantification')
		save(tad_list,file=file_f)
	}
}

















averageTrack <- function(tracks=all_tracks,regions=NULL,bins,anchor_middle=F,flank_l=1000,flank_r=1000){
  max_dist <- 1e8
  chrom_sizes <- gintervals.all()
  for (t in 1:length(tracks)){
    gvtrack.create(paste0('v_',tracks[t]), tracks[t], "avg")
  }
  if (is.null(regions)) {
    domain <- paste0("hic.",cell,".ins250_k2_domains")
    domain <- gintervals.load(domain)
  } else {
    if (!grepl('bed|narrowPeak',regions)) {
      domain <- gintervals.load(regions)
    } else {
      domain <- read.table(regions,sep='\t',header=F)
      domain <- domain[domain[,1]%in% gintervals.all()$chrom,]
      domain <- gintervals(domain[,1],domain[,2],domain[,3])
    }
  }
  domain$len <- abs(domain[,3]-domain[,2])
  domain <- domain[domain$start>1e6,]
  chr_ends <- ddply(domain,.(chrom),function(x){return(max(x$end))})
  domain <- domain[!(domain$end %in% chr_ends$V1),]
  if(!anchor_middle){
    res <- foreach(i=1:nrow(domain))%dopar%{
      interval <- domain[i,]
      interval$start = interval$start - interval$len
      interval$end = interval$end + interval$len
      interval <- gintervals.force_range(interval)
      bin <- (interval$end-interval$start)/bins
      if(bin<1){bin=1}
      bins_r <- giterator.intervals(intervals=interval, iterator=bin)
      df <- gextract(paste0('v_',tracks),intervals=bins_r,iterator=bins_r)
      temp <- rep(NA,(bins+1))
      if (!is.null(interval$strand)){
        if(interval$strand==1){
          temp[1:nrow(df)] <- df[,4]
        } else {
          temp[1:nrow(df)] <- rev(df[,4])
        }
      }else {
        temp[1:nrow(df)] <- df[,4]
      }
      return(temp[1:(bins+1)])
    }
  } else {
    res <- foreach(i=1:nrow(domain))%dopar%{
      interval <- domain[i,]
      interval$start = round((domain[i,'start']+domain[i,'end'])/2) - flank_l
      interval$end = round((domain[i,'start']+domain[i,'end'])/2) + flank_r
      interval <- gintervals.force_range(interval)
      bin <- (interval$end-interval$start)/bins
      bins_r <- giterator.intervals(intervals=interval, iterator=bin)
      df <- gextract(paste0('v_',tracks),intervals=bins_r,iterator=bins_r)
      temp <- rep(NA,(bins+1))
      if(bin<1){bin=1}
      if (!is.null(interval$strand)){
        if(interval$strand==1){
          temp[1:nrow(df)] <- df[,4]
        } else {
          temp[1:nrow(df)] <- rev(df[,4])
        }
      }else {
        temp[1:nrow(df)] <- df[,4]
      }
      return(temp[1:(bins+1)])
    }
  }
  names(res) <- 1:length(res)
  res <- as.data.frame(t(bind_rows(res)))
  res$id <- paste0(domain$chrom,':',domain$start,'_',domain$end)
  if (!is.null(domain$cluster)){
    res$cluster <- domain$cluster
  } else {
    res$cluster <- 1
  }
  return(res)
}













