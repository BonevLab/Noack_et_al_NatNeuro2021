intervals.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  return(inv)
}
intervals.2d.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  inv[,5:6] = floor((inv[,5]+inv[,6])/2)
  inv[,6] = inv[,5] + 1
  return(inv)
}

intervals.expand <- function(inv,expansion=100){
  inv[,2]<-inv[,2]-expansion
  inv[,3]<-inv[,3]+expansion
  return(gintervals.force_range(inv))
}

intervals.2d.expand <- function(inv,expansion1, expansion2){
  inv[,2]<-inv[,2]-expansion1
  inv[,3]<-inv[,3]+expansion1
  inv[,5]<-inv[,5]-expansion2
  inv[,6]<-inv[,6]+expansion2
  return(gintervals.force_range(inv))
}


intervals.size <- function(inv){
  return(sum(inv[,3]-inv[,2]))
}

intervals.normalize <- function(inv, size) {
  centers <- intervals.centers(inv)
  centers$end <- centers$end-1;
  return(intervals.expand(centers, floor(size/2)))
}

TGLKMeans_wrapper=function(data, fn, k)
{
  write.table(data, fn, sep="\t", quote=F, col.names=NA)
  system(sprintf("%sscripts/TGLKMeans_static %s %s euclid -allow_nas=1 &> %s.log",main_f,fn,k,fn))
  km = list()
  m.k = read.table(paste(fn, "kclust", sep="."), header=T)
  m.k$clust = m.k$clust + 1
  m.c = read.delim(paste(fn, "center", sep="."), header=F)
  km$size = tapply(m.k$id, m.k$clust, length)
  km$cluster = m.k$clust
  names(km$cluster) = m.k$id
  km$centers = m.c[, 2:ncol(m.c)]
  km
}

KRnorm = function(A) {
	A_orig <- A
  # remove any cols/rows of 0s
  zeros = unique(which(colSums(A) == 0), which(rowSums(A) == 0))
  if (length(zeros) > 0) {
    A <- A[-zeros, -zeros]
    message(paste0('Cols/Rows removed: '))
    message(paste(" ", zeros, sep = " "))
  }
  # initialize
  tol = 1e-6; delta = 0.1; Delta = 3; fl = 0;
  # change NAs in matrix to 0's
  NAlist = which(is.na(A), arr.ind = TRUE)
  A[is.na(A)] = 0
  n = nrow(A)
  e = matrix(1, nrow=n, ncol = 1)
  x0 = e
  res = matrix(nrow = n, ncol=1)
  # inner stopping criterior
  g=0.9; etamax = 0.1;
  eta = etamax; stop_tol = tol*.5;
  x = x0; rt = tol^2; v = x*(A %*% x); rk = 1 - v;
  rho_km1 = t(rk) %*% rk; rout = rho_km1; rold = rout;
  MVP = 0; # We'll count matrix vector products.
  i = 0; # Outer iteration count.
  while(rout > rt) { # Outer iteration
    i = i + 1; k = 0; y = e;
    innertol = max(c(eta^2 %*% rout, rt));
    while( rho_km1 > innertol ) { #Inner iteration by CG
      k = k +1
      if( k ==1) {
        z = rk / v; p=z; rho_km1 = t(rk)%*%z;
      }else {
        beta = rho_km1 %*% solve(rho_km2)
        p = z + beta*p
      }

      # update search direction efficiently
      w = x * (A%*%(x*p)) + v*p
      alpha = rho_km1 %*% solve(t(p) %*% w)
      ap = c(alpha) * p
      # test distance to boundary of cone
      ynew = y + ap;
      if(min(ynew) <= delta) {
        if(delta == 0) break()
        ind = which(ap < 0);
        gamma = min((delta - y[ind])/ap[ind]);
        y = y + gamma %*% ap;
        break()
      }
      if(max(ynew) >= Delta) {
        ind = which(ynew > Delta);
        gamma = min((Delta-y[ind])/ap[ind]);
        y = y + gamma %*% ap;
        break()
      }
      y = ynew;
      rk = rk - c(alpha) * w; rho_km2 = rho_km1;
      Z = rk/v; rho_km1 = t(rk) %*% z;
    }
    x = x*y; v = x*(A %*% x);
    rk = 1 - v; rho_km1 = t(rk) %*% rk; rout = rho_km1;
    MVP = MVP + k + 1;
    # Update inner iteration stopping criterion.
    rat = rout %*% solve(rold); rold = rout; res_norm = sqrt(rout);
    eta_o = eta; eta = g %*% rat;
    if(g %*% eta_o^2 > 0.1) {
      eta = max(c(eta, g %*% eta_o^2));
    }
    eta = max(c(min(c(eta, etamax)), stop_tol/res_norm));
  }
  result = diag(c(x)) %*% A %*% diag(c(x))
  # reintroduce NAs in final matrix
  if(nrow(NAlist) > 0) {
    idx <- as.matrix(NAlist[, 1:2])
    result[idx] <- NA
  }
	test <- A_orig
	normed <- result
	test[1:ncol(normed),1:nrow(normed)] <- normed
	if (length(zeros)>0){
		for (i in 1:length(zeros)){
			idx <- zeros[i]
			if (idx!=ncol(test)){
				test[,(idx+1):ncol(test)] <- test[,idx:(ncol(test)-1)]
				test[(idx+1):nrow(test),] <- test[idx:(nrow(test)-1),]
			} else {
				test[,idx] <- 0
				test[idx,] <- 0
			}
		}
		test[,zeros] <- 0
		test[zeros,] <- 0
	}
  return(test)
}

gtrack.2d.get_insu_doms = function(insu_track, thresh, iterator=500)
{
	doms = gscreen(sprintf("is.na(%s) | %s > %f", insu_track, insu_track, thresh), iterator=iterator)
	return(doms)
}

image.scale <- function(z, zlim, col = heat.colors(12),breaks, axis.pos=1, add.axis=TRUE,cex.axis=2.4,cex.lab=10,axis_lab=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
 if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
 plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
 for(i in seq(poly)){
  if(axis.pos %in% c(1,3)){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(axis.pos %in% c(2,4)){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
 box()
 if(add.axis) {
 zlim=c(zlim,0)
 zlim <- zlim[order(zlim)]
 if(is.null(axis_lab)){
   axis_lab=zlim
 }
 axis(axis.pos,lab=axis_lab,at=zlim, las=2,cex.axis=cex.axis,cex.lab=cex.lab,lwd.ticks = ifelse(axis_lab,0,1))
 #mtext(side=2,text='log2(obs/exp',line=1,cex=2)
 }
}

image.scale.aggregate <- function(z, zlim, col = heat.colors(12),
breaks, axis.pos=1, add.axis=TRUE, label=NULL){
  par(mgp=c(3,0.1,0))
   if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-5)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-5)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
 if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
 plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
 for(i in seq(poly)){
  if(axis.pos %in% c(1,3)){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(axis.pos %in% c(2,4)){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
 box()
 if(add.axis) {
 axis(axis.pos,las=2,tck = -.01,labels = c(min(zlim),paste0('  ',0),paste0('+',max(zlim))),at=c(min(zlim),0,max(zlim)),cex.axis=1.2)
 mtext(side=2,line=0.5, text=label)
 
 }
}

rotate <- function(x) t(apply(x, 2, rev))

makeTransparent = function(..., alpha=0.5) {

  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")

  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)

  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }

  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)

  return(newColor)

}

minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,base=10,nticks=5,...){

  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]

  major.ticks <- pretty(lims,n=nticks)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)

  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]

  if (base==10){
  labels <- sapply(major.ticks,function(i)
            as.expression(bquote(10^ .(i)))
          )
    } else {
  labels <- sapply(major.ticks,function(i)
            as.expression(bquote(2^ .(i)))
          )
    }
  axis(ax,at=major.ticks,labels=labels,...)

  n <- n+2
  if (base==10) {
  	minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  	} else {
  minors <- log2(pretty(2^major.ticks[1:2],n))-major.ticks[1]
  }
  minors <- minors[-c(1,n)]

  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]


  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}



construct.grid.comp = function(interv1,interv2,min_dist,max_dist){
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

get_borders <- function(iterator=1e3, ins_track=ins_track,ins_dom_thresh,extend_border_region)
{
  ins = gscreen(sprintf("(!is.na(%s)) & (%s >=(-5)) & (%s <= %f) ", ins_track, ins_track, ins_track, gquantiles(ins_track, ins_dom_thresh)), iterator=iterator)
  chroms = gintervals.all()
  rownames(chroms) = as.character(chroms$chrom)
  ins$max = chroms[ as.character(ins$chrom), 'end']
  ins$start = pmax(ins$start - extend_border_region, 0)
  ins$end   = pmin(ins$end   + extend_border_region, ins$max)
  insc = gintervals.canonic(ins, unify=T)
  gvtrack.create("ins_min", ins_track, "min")
  bords = gextract("ins_min", intervals=insc, iterator=insc)
  bords <- bords[(!is.na(bords$ins_min))&(bords$ins_min>=(-5))&(bords$ins_min<gquantiles(ins_track, ins_dom_thresh)),]
  gvtrack.rm("ins_min")
  return(bords[,1:3])
}

get_domains = function(ins_track=ins_track,min_domainSize=5e4,ins_dom_thresh=0.1)
{
  d = gtrack.2d.get_insu_doms(ins_track, gquantiles(ins_track, ins_dom_thresh),iterator=1e3)
  d$len = d$end - d$start
  d = d[d$len >= min_domainSize, ]
  return(d[,1:3])
}

split_and_merge_matrix <- function(mat, grid_min, grid_max, grid_res, new_res) {
  if (!(new_res / grid_res == floor(new_res/grid_res))) {
	message(paste0(new_res, " is not a multiple of ", grid_res))
	return(integer(0))
  }
  expand = seq(grid_min, grid_max,by=grid_res)
  multi = new_res/grid_res
  m = melt(mat)
  m$i = as.numeric(m[,1])
  m$j = as.numeric(m[,2])
  if (identical(which(expand == 0), integer(0))) {
	#center spans the 0 need to take only a single bin
	if (multi/2 == floor(multi/2)) {
		message(paste0("cannot maintain center at 0 - ", multi, " is even, center spans 0"))
		return(integer(0))
	}
	center_bin = length(expand)/2
	new_pos_bin = c(rep(0, (multi-1)/2), rep(1:(max(m$i)/2/multi), each=multi))[1:(center_bin-1)]
	new_bin = c(rev(-new_pos_bin), 0, new_pos_bin)
  } else {
	if (!(multi/2 == floor(multi/2))) {
		message(paste0("cannot maintain center at 0 - ", multi, " is odd, center does not span 0"))
		return(integer(0))
	}
	center_bin_pos = which(expand == 0)
	new_pos_bin = c(rep(0, multi/2), rep(1:(max(m$i)/2/multi), each=multi))[1:(center_bin_pos-1)]
	new_bin = c(rev(-new_pos_bin), new_pos_bin)
  }
  m$n_i = new_bin[m$i]
  m$n_j = new_bin[m$j]
  d = ddply(m, .(n_i, n_j), summarise, sum=sum(value, na.rm=TRUE))
  return(dcast(d, n_i ~ n_j)[,-1])
}

extractDiagonal <- function(inMatrix)
{
	tempMatrix <- inMatrix
	nr <- nrow(inMatrix)
	nc <- ncol(inMatrix)
	countPasses <- 1
	outMatrix <- c()

	for (i in 1:nr)
	{
		diagonal <- c()
		if (i == nr)
		{
			diagonal <- inMatrix[1,nr]
			diagonal <- c(rep(NA,countPasses/2), diagonal, rep(NA,countPasses/2))
			if (ncol(inMatrix) %% 2 == 0){diagonal <- diagonal[1:nc]}
			outMatrix <- rbind(outMatrix, diagonal)
		}
		else
		{
			for (j in 1:nrow(tempMatrix))
			{
				for (k in 1:ncol(tempMatrix))
				{
					if (j==k){diagonal <- append(diagonal, tempMatrix[j,k])}
				}
			}

			if (length(diagonal) == ncol(inMatrix))
			{
				outMatrix=diagonal
			}
			else if (length(diagonal) %in% seq(1,ncol(inMatrix),by=2))
			{
				diagonal <- c(rep(NA,countPasses/2), diagonal, rep(NA,countPasses/2))
				if (ncol(inMatrix) %% 2 == 0){diagonal <- diagonal[1:nc]}
				outMatrix <- rbind(outMatrix, diagonal)
			}

			if (nrow(tempMatrix) > 2)
			{
				tempMatrix <- tempMatrix[,-1]
				tempMatrix <- tempMatrix[-(nrow(tempMatrix)),]
			}
			countPasses <- countPasses + 1
		}
	}
	return(outMatrix)
}

calc_domains_a_score <- function(tn="scell.nextera.pool_good_hyb_2i_all_es",domains=NULL, vfunc="weighted.sum", d=NULL, use_trans=T, rebuild=F)
{
  ifn = sprintf("%s/%s_a_score_%s_%dDoms_%s.rdata", sch_rdata_dir, tn, vfunc, ifelse(is.null(d), 0, nrow(d)), ifelse(use_trans, "trans", "cis"))

  if (!file.exists(ifn) | rebuild) {
    message(sprintf("calc fraction of contacts with A domains using %s ...", tn))
    rd = calc_domains_ab_frac(tn=tn, d=d, cis=F, vfunc=vfunc, calc_frac_a=T, rebuild=rebuild)
    rd2 = merge(rd$d, rd$stats[, c('chrom1', 'start1', 'A', 'B', 'pA')], by.x=c('chrom', 'start'), by.y=c('chrom1', 'start1'))
    rownames(rd2) = paste(rd2$chrom, rd2$start, sep="_")

    message(sprintf("calc weigthed fraction of contacts with A domains using %s ...", tn))
    dd = calc_domains_ab_frac(tn=tn, d=rd2, cis=!use_trans, vfunc=vfunc, calc_frac_a=F, rebuild=rebuild)
    dd = dd$stats
    rd2[paste(dd$chrom1, dd$start1, sep="_"), ifelse(use_trans, 'trans_A_score', 'cis_A_score')] = dd$pA

    write.table(rd2, ifn, quote=F, sep="\t")
  }

  read.table(ifn, header=T)
}

calc_domains_ab_frac <- function(tn=pool_tn, d=NULL, min_cis=2e6, vfunc="weighted.sum", cis=F, calc_frac_a=T, rebuild=F)
{
  if (is.null(d)) {
    d = sch_get_ab_tagged_domains(raw_tn=tn, vfunc=vfunc, rebuild=rebuild)
  }

  ints = gintervals.2d.all()
  if (cis) {
    ints = ints[as.character(ints$chrom1) == as.character(ints$chrom2),]
  }
  else {
    ints = ints[as.character(ints$chrom1) != as.character(ints$chrom2),]
  }
  d_chroms = unique(as.character(d$chrom))

  ints = ints[is.element(as.character(ints$chrom1), d_chroms) & is.element(as.character(ints$chrom2), d_chroms), ]

  commands = paste0("{source(\"/work/project/Cavalli-mammals/boyan/HiC/lscripts/analyzeHiC.r\"); calc_domains_ab_frac_for_chrom_pairs(d=d, tn=\"", tn, "\", chr1=\"", ints$chrom1, "\", chr2=\"", ints$chrom2, "\", min_cis=min_cis, vfunc=vfunc, calc_frac_a=calc_frac_a)}", collapse=",")

  res = eval(parse(text=paste("gcluster.run(",commands,",opt.flags = \"-l mem=20G -l h_vmem=20G\")")))
  #stopifnot(sum(unlist(lapply(res, function(r) typeof(r$retv))) != "list") == 0)
  res_s = do.call("rbind", lapply(res, function(r) {if (ncol(r$retv)==4){return(r$retv)}}))

  if (calc_frac_a) {
    stats = as.data.frame(summarize(group_by(res_s, chrom1, start1), A=sum(A), B=sum(B)))
    stats$pA = stats$A / (stats$A + stats$B)
  #  stats$pB = stats$B / (stats$A + stats$B)
  }
  else {
    stats = as.data.frame(summarize(group_by(res_s, chrom1, start1), wA=sum(wA), n=sum(n)))
    stats$pA = stats$wA / stats$n
   # stats$pB = stats$B / (stats$A + stats$B)
  }

  list(d=d, stats=stats)
}

calc_domains_ab_frac_for_chrom_pairs <- function(d, tn=pool_tn, chr1="chr1", chr2="chr2", min_cis=2e6, vfunc="weighted.sum", calc_frac_a=T)
{
  stopifnot(calc_frac_a && is.element("ab", colnames(d)) || !calc_frac_a && is.element("pA", colnames(d)))

  fstr = ifelse(calc_frac_a, 'ab', 'pA')

  library(plyr)
  gvtrack.create("tnv", tn, vfunc)
  d1 = d[as.character(d$chrom) == chr1, ]
  d2 = d[as.character(d$chrom) == chr2, ]

  comb = expand.grid(1:nrow(d1), 1:nrow(d2))
  ints = gintervals.2d(d1$chrom[comb[,1]], d1$start[comb[,1]], d1$end[comb[,1]], d2$chrom[comb[,2]], d2$start[comb[,2]], d2$end[comb[,2]])

  if (chr1 == chr2) {
    x = gextract("ifelse(is.na(tnv), 0, tnv)", intervals=ints, iterator=ints, colnames='count', band=c(-1e9, -min_cis))
    x_non_trim1 = data.frame(chrom1=x$chrom1, start1=ints[x$intervalID, 'start1'], end1=ints[x$intervalID, 'end1'], chrom2=x$chrom1, start2=ints[x$intervalID, 'start2'], end2=ints[x$intervalID, 'end2'], count=x$count, intervalID=x$intervalID, d2=d2[comb[x$intervalID, 2], fstr])
    x_non_trim2 = data.frame(chrom1=x$chrom1, start1=ints[x$intervalID, 'start2'], end1=ints[x$intervalID, 'end2'], chrom2=x$chrom1, start2=ints[x$intervalID, 'start1'], end2=ints[x$intervalID, 'end1'], count=x$count, intervalID=x$intervalID, d2=d2[comb[x$intervalID, 1], fstr])

    x = rbind(x_non_trim1, x_non_trim2)

  }
  else {
    x = gextract("ifelse(is.na(tnv), 0, tnv)", intervals=ints, iterator=ints, colnames='count')
    x$d2 = d2[comb[,2], fstr]

  }

  if (calc_frac_a) {
    res = ddply(x, .(chrom1, start1), function(v) { tapply(v$count, v$d2, sum) } )
  }
  else {
    res = ddply(x, .(chrom1, start1), function(v) { data.frame(wA=sum(v$count * v$d2, na.rm=T), n=sum(v$count)) } )
  }

  res
}

###################################################################################################
#' Quantify contact enrichment between pairs of intervals.
#'
#' \code{submit_cisDecayIntervals}
#'
#' This function creates pairs of intervals given 1D gintervals saved as Rdata file or bed files and calculates the contact enrichment (log2(observed/expected)) with a certain distance of this points. Data is saved and then can be plotted using \code{plot_cis_decay}
#'
#' @param cells Conditions to work on.
#' @param domains Which domain tracks to use when constructing intra- vs inter-domain pairs. Two options - either gintervals format which will be used for all datasets and exist in the database or string which will be searched for in each condition folder. Usually generated by \code{analyzeCompartments} or \code{analyzeInsulation}. If not explicitly specified will use the domains generated by \code{analyzeInsulation}.
#' @param window_i If bigger than 0 create pairwise intervals as points and expand them in each direction using this value (this creates a square). If it is 0, it will just create pairwise intervals as is - without expanding. This mode is useful for big intervals, for example H3K9me3, H3K27me3 or full domains.
#' @param intervals1 Either gintervals or bed files with peak coordinates of the first regions.
#' @param intervals2 Either gintervals or bed files with peak coordinates of the second regions.
#' @param grid_mode One of '1D','2D','merged'. If 1D - it will create only pairs between first intervals as upstream and second intervals as downstream anchors. If 2D, it will also include pairs where interval2 is upstream. if merged - it will create pairs by matching each row in the intervals1 file with the corresponding row in the intervals2. This mode is useful for premade pairs - for example per gene, or using juicebox.
#' @param mem Memory available - increase if memory issues.
#'
#' @examples
#'
#' submit_cisDecayIntervals(cells=c('D0','D2','D4','D6','D10'),window_i=10000,intervals1='/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/D0_CTCF_For.bed',intervals2='/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/D0_CTCF_Rev.bed',grid_mode='1D',mem=80)   #Will generate pairs based on for-rev CTCF peaks in the D0 dataset and examine the contact enrichment in 10x10kb square in the listed conditions (D0-D10)
#' submit_cisDecayIntervals(cells='D4',window_i=0,intervals1='/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/D4_H3K9me3.bed',intervals2='/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/D4_H3K9me3.bed',grid_mode='1D',mem=80)   #Will generate pairs based on H3K9me3 peaks in the D4 dataset and examine the contact enrichment within the specified regions (no extension) in the D4 condition.
#' submit_cisDecayIntervals(cells='D0',window_i=5000,intervals1='/work/project/Cavalli-mammals/satish/HiC/data/rna/D0_activeTSS',intervals2='/work/project/Cavalli-mammals/satish/HiC/data/rna/D0_activeTTS_matched',grid_mode='merged',mem=80)   #For each entry in the first file it will create a pair by taking the corresponding row from the second interval, expand it by 5x5kb and calculate contact enrichment.
#'
#' @export
##########################################################################################################
submit_cisDecayIntervals <- function(cells,domains='ins_250_domains_expanded',window_i=5000,intervals1,intervals2,grid_mode=c('1D','2D','merged'),mem=100){
	f1=basename(intervals1)
	f2=basename(intervals2)
	command_f <- paste0('Rscript ',main_f,'/scripts/cis_decay_intervals.R ',intervals1,' ',intervals2,' ', paste0(f1,'_',f2,'.',grid_mode,'.',window_i),' ', grid_mode,' ',window_i,' ',paste0(cells,collapse=','), ' ',main_f,' ',domains,' &>',main_f,'/logs/cis_decay/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'.log'))
	cat(command_f,file=paste0(main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_cisDecay.sh')))
	if (sunGrid){
		command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/cis_decay/',paste0(f1,'-',f2,'.',grid_mode,'.',window_i,'.e'),' -o ',main_f,'logs/cis_decay/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'.o '),main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_cisDecay.sh'))
 	} else {
 		command_f <- paste0('nohup sh ',main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_cisDecay.sh &>/dev/null &'))
 	}
	system(command_f)
}

###################################################################################################
#' Generate aggregate Hi-C plot to quantify the contact enrichment between two intervals.
#'
#' \code{submit_aggregateHiC}
#'
#' This function creates pairs of intervals given 1D gintervals saved as Rdata file or bed files and calculates the aggregate Hi-C profile (log2(observed/expected)) with a certain distance of this points. Data is saved and then can be plotted using \code{plot_aggregateHiC}. It will check if data files have already been generated and skip them - delete those files if you want to rerun.
#'
#' @param cells Conditions to work on.
#' @param tracks Which hi-c tracks to work on. Defaults to those defined in the configuration file.
#' @param domains Which domain tracks to use when constructing intra- vs inter-domain pairs. Must be in gintervals format and exist in the database. Usually generated by \code{analyzeCompartments} or \code{analyzeInsulation}. If not explicitly specified will use the domains generated by \code{analyzeInsulation}.
#' @param range_f Distance of the square emanating from the center of the pair to be examined. Increase this and decrease resolution for sparse data. Not recommended to exceed 80kb.
#' @param res_f Size of the bin to estimate contact enrichment when extracting the data. Recommended value for sparse data - 1kb, for high-resolution data - 500bp. This can be increased but not decreased when plotting (i.e. combining bins).
#' @param filter_f If not zero, all the pairs will first be filtered based on the maximum hi-c score in them and only pairs with higher than this value will be retained. Not recommended.
#' @param intervals1 Either gintervals or bed files with peak coordinates of the first regions.
#' @param intervals2 Either gintervals or bed files with peak coordinates of the second regions.
#' @param grid_mode One of '1D','2D','merged'. If 1D - it will create only pairs between first intervals as upstream and second intervals as downstream anchors. If 2D, it will also include pairs where interval2 is upstream. if merged - it will create pairs by matching each row in the intervals1 file with the corresponding row in the intervals2. This mode is useful for premade pairs - for example per gene, or using juicebox.
#' @param mem Memory available - increase if memory issues.
#' @param k Value of k when the score track was generated if filter_f is not zero.
#'
#' @examples
#'
#' submit_aggregateHiC(cells=c('D0','D2'),tracks=all_tracks,range_f=20000,res_f=1000,filter_f=0,intervals1='/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/D0_CTCF_For.bed',intervals2='/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/D0_CTCF_Rev.bed',grid_mode='1D',mem=80)
#'
#' @export
##########################################################################################################
submit_aggregateHiC <- function(cells,tracks=all_tracks,domains=gintervals.ls('expanded'),range_f,res_f,filter_f=0,intervals1,intervals2,grid_mode=c('1D','2D','merged'),mem=100,k=100,use_score=FALSE){
	f1=basename(intervals1)
	f2=basename(intervals2)
	for (cell in cells){
		n_tracks <- length(tracks[grep(cell,tracks)])
		for (i in 1:n_tracks){
			track <- tracks[grep(cell,tracks)]
			file_f <- paste0(cell,'_',i,'.',f1,'.',f2,'.',range_f,'.',res_f,'.',grid_mode,'.',filter_f)
			if (!file.exists(paste0(main_f,'data/aggregateHiC/',file_f))){
				command_f <- paste0('Rscript ',main_f,'/scripts/aggregateHiC.R ',main_f,' ',track[i],' k',k,' ',filter_f,' ',intervals1,' ',intervals2,' ', paste0(main_f,'data/aggregateHiC/',file_f),' ',range_f,' ',res_f,' ',domains,' ',grid_mode,' ',f1,'_',f2,' ',use_score,' &>',main_f,'/logs/aggregateHiC/',file_f,'.log')
				cat(command_f,file=paste0(main_f,'/qsub/',file_f,'_aggregate.sh'))
				if (sunGrid){
					command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/aggregateHiC/',file_f,'.e',' -o ',main_f,'logs/aggregateHiC/',file_f,'.o ',paste0(main_f,'/qsub/',file_f,'_aggregate.sh'))
				} else {
					command_f <- paste0('nohup sh ',main_f,'/qsub/',file_f,'_aggregate.sh 2>/dev/null &')
				}
				system(command_f)
			} else {message ('File ',file_f,' exists. Skipping...')}
		}
	}
}

###################################################################################################
#' Generate aggregate Hi-C plot to quantify the insulation at certain intervals.
#'
#' \code{submit_aggregateDiagonal}
#'
#' This function calculates the aggregate Hi-C profile (log2(observed/expected)) within a certain distance of interval points. Data is saved and then can be plotted using \code{plot_aggregateDiagonal}. It will check if data files have already been generated and skip them - delete those files if you want to rerun.
#'
#' @param cells Conditions to work on.
#' @param tracks Which hi-c tracks to work on. Defaults to those defined in the configuration file.
#' @param minDist Minimum distance to be considered when calculating enrichment.
#' @param maxDist  Maximum distance to be considered when calculating enrichment.
#' @param res Size of the bin to estimate contact enrichment when extracting the data. Recommended value for sparse data - 1kb, for high-resolution data - 500bp. This can be increased but not decreased when plotting (i.e. combining bins).
#' @param filter_f If not zero, all the pairs will first be filtered based on the maximum hi-c score in them and only pairs with higher than this value will be retained. Not recommended.
#' @param intervals_f Either gintervals or bed files with peak coordinates of the regions.
#' @param consider_strand Boolean indicating whether to reorient intervals based on strand of the nearest tss.
#' @param mem Memory available - increase if memory issues.
#' @param k Value of k when the score track was generated if filter_f is not zero.
#'
#' @examples
#'
#' submit_aggregateDiagonal(cells=c('D0','D2'),tracks=all_tracks,maxDist=2e5,res=500,filter_f=0,intervals_f='/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/D0_CTCF.bed',consider_strand=F,mem=140)
#'
#' @export
##########################################################################################################
submit_aggregateDiagonal <- function(cells,tracks=all_tracks,minDist=1e3,maxDist=2e5,res=500,filter_f=0,intervals_f,consider_strand='FALSE',mem=100,k=200){
	f1=basename(intervals_f)
	for (cell in cells){
		n_tracks <- length(tracks[grep(cell,tracks)])
		for (i in 1:n_tracks){
			track <- tracks[grep(cell,tracks)]
			file_f <- paste0(cell,'_',i,'.',f1,'.',maxDist,'.',res,'.',consider_strand,'.',filter_f)
			if (!file.exists(paste0(main_f,'data/aggregateDiagonal/',file_f))){
				command_f <- paste0('Rscript ',main_f,'/scripts/diagonal_featureGrid.R ',main_f,' ',track[i],' k',k,' ',filter_f,' ',intervals_f,' ', paste0(main_f,'data/aggregateDiagonal/',file_f),' ',minDist,' ',maxDist,' ',res,' ',consider_strand,' ',f1, ' &>',main_f,'/logs/aggregateDiagonal/',file_f,'.log')
	# 			print(command_f)
				cat(command_f,file=paste0(main_f,'/qsub/',file_f,'_aggregate.sh'))
				if (sunGrid){
					command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/aggregateDiagonal/',file_f,'.e',' -o ',main_f,'logs/aggregateDiagonal/',file_f,'.o ',paste0(main_f,'/qsub/',file_f,'_aggregate.sh'))
				} else {
					command_f <- paste0('sh ',main_f,'/qsub/',file_f,'_aggregate.sh')
				}
				system(command_f)
			} else {message ('File ',file_f,' exists. Skipping...')}
		}
	}
}

pool_aggregateHiC <- function(cells,intervals1,intervals2,range_f,filter_f=0,res_f,grid_mode,plot_res){
	f1=basename(intervals1)
	f2=basename(intervals2)
	for (cell in cells){
		file_idx <- paste0(f1,'.',f2,'.',range_f,'.',res_f,'.',grid_mode,'.',filter_f)
		files <- list.files(paste0(main_f,'data/aggregateHiC/'),pattern=paste0('.',file_idx),full.names=T)
		files <- files[grep(cell,files)]
		output = paste0(main_f,'data/aggregateHiC/',file_idx,'.',cell)
		load(files[1])
		g = grid;
		enrich_v <- list()
		for (i in 1:length(g$obs)) {
			g$obs[[i]] = g$obs[[i]] * g$total_obs[[i]]
			g$exp[[i]] = g$exp[[i]] * g$total_exp[[i]]
			o = split_and_merge_matrix(grid$obs[[i]] * grid$total_obs[[i]], -range_f, range_f, res_f, plot_res)
			e = split_and_merge_matrix(grid$exp[[i]] * grid$total_exp[[i]], -range_f, range_f, res_f, plot_res)
		  center_v <- sum(o[(which(colnames(o)==0)-1):(which(colnames(o)==0)+1),(which(colnames(o)==0)-1):(which(colnames(o)==0)+1)],na.rm=T)/(sum(e[(which(colnames(e)==0)-1):(which(colnames(e)==0)+1),(which(colnames(e)==0)-1):(which(colnames(e)==0)+1)],na.rm=T)/2)
      left_top <- sum(o[2:4,2:4],na.rm=T)/(sum(e[2:4,2:4],na.rm=T)/2)
	    left_bottom <- sum(o[(nrow(o)-3):(nrow(o)-1),2:4],na.rm=T)/(sum(e[(nrow(o)-3):(nrow(o)-1),2:4],na.rm=T)/2)
	    right_top <- sum(o[2:4,(ncol(o)-3):(ncol(o)-1)],na.rm=T)/(sum(e[2:4,(ncol(o)-3):(ncol(o)-1)],na.rm=T)/2)
	    right_bottom <- sum(o[(nrow(o)-3):(nrow(o)-1),(ncol(o)-3):(ncol(o)-1)],na.rm=T)/(sum(e[(nrow(o)-3):(nrow(o)-1),(ncol(o)-3):(ncol(o)-1)],na.rm=T)/2)
#			center_v <- sum(o[(which(colnames(o)==0)),(which(colnames(o)==0))],na.rm=T)/(sum(e[(which(colnames(o)==0)),(which(colnames(o)==0))],na.rm=T)/2)
#			left_top <- sum(o[3,3],na.rm=T)/(sum(e[3,3],na.rm=T)/2)
#			left_bottom <- sum(o[(nrow(o)-2),3],na.rm=T)/(sum(e[(nrow(o)-2),3],na.rm=T)/2)
#			right_top <- sum(o[3,(ncol(o)-2)],na.rm=T)/(sum(e[3,(ncol(o)-2)],na.rm=T)/2)
#			right_bottom <- sum(o[(nrow(o)-2),(ncol(o)-2)],na.rm=T)/(sum(e[(nrow(o)-2),(ncol(o)-2)],na.rm=T)/2)
		#	enrich_v[[i]] <- c(enrich_v[[i]],center_v/mean(c(left_top,left_bottom,right_top,right_bottom)))
			enrich_v[[i]] <- center_v/left_bottom
	    }
		for (file_f in files[2:length(files)]) {
		 load(file=file_f)
		 for (i in 1:length(g$obs)) {
		   o = split_and_merge_matrix(grid$obs[[i]] * grid$total_obs[[i]], -range_f, range_f, res_f, plot_res)
		   e = split_and_merge_matrix(grid$exp[[i]] * grid$total_exp[[i]], -range_f, range_f, res_f, plot_res)
		   center_v <- sum(o[(which(colnames(o)==0)-1):(which(colnames(o)==0)+1),(which(colnames(o)==0)-1):(which(colnames(o)==0)+1)],na.rm=T)/(sum(e[(which(colnames(e)==0)-1):(which(colnames(e)==0)+1),(which(colnames(e)==0)-1):(which(colnames(e)==0)+1)],na.rm=T)/2)
		   left_top <- sum(o[2:4,2:4],na.rm=T)/(sum(e[2:4,2:4],na.rm=T)/2)
		   left_bottom <- sum(o[(nrow(o)-3):(nrow(o)-1),2:4],na.rm=T)/(sum(e[(nrow(o)-3):(nrow(o)-1),2:4],na.rm=T)/2)
		   right_top <- sum(o[2:4,(ncol(o)-3):(ncol(o)-1)],na.rm=T)/(sum(e[2:4,(ncol(o)-3):(ncol(o)-1)],na.rm=T)/2)
		   right_bottom <- sum(o[(nrow(o)-3):(nrow(o)-1),(ncol(o)-3):(ncol(o)-1)],na.rm=T)/(sum(e[(nrow(o)-3):(nrow(o)-1),(ncol(o)-3):(ncol(o)-1)],na.rm=T)/2)
#		   center_v <- sum(o[(which(colnames(o)==0)),(which(colnames(o)==0))],na.rm=T)/(sum(e[(which(colnames(o)==0)),(which(colnames(o)==0))],na.rm=T)/2)
#		   left_top <- sum(o[3,3],na.rm=T)/(sum(e[3,3],na.rm=T)/2)
#		   left_bottom <- sum(o[(nrow(o)-2),3],na.rm=T)/(sum(e[(nrow(o)-2),3],na.rm=T)/2)
#		   right_top <- sum(o[3,(ncol(o)-2)],na.rm=T)/(sum(e[3,(ncol(o)-2)],na.rm=T)/2)
#		   right_bottom <- sum(o[(nrow(o)-2),(ncol(o)-2)],na.rm=T)/(sum(e[(nrow(o)-2),(ncol(o)-2)],na.rm=T)/2)
		   enrich_v[[i]] <- c(enrich_v[[i]],center_v/mean(c(left_top,left_bottom,right_top,right_bottom)))
#		   enrich_v[[i]] <- c(enrich_v[[i]],center_v/left_bottom)
		   g$obs[[i]] = g$obs[[i]] + grid$obs[[i]] * grid$total_obs[[i]]
		   g$exp[[i]] = g$exp[[i]] + grid$exp[[i]] * grid$total_exp[[i]]
		   g$total_obs[[i]] = g$total_obs[[i]] + grid$total_obs[[i]]
		   g$total_exp[[i]] = g$total_exp[[i]] + grid$total_exp[[i]]
		 }
		}
		for (i in 1:length(g$obs)) {
		  g$obs[[i]] = g$obs[[i]] / g$total_obs[[i]]
		  g$exp[[i]] = g$exp[[i]] / g$total_exp[[i]]
		  g$enrich[[i]] <- enrich_v[[i]]
		}
		grid = g
		save(grid, file=output)
	}
}

pool_aggregateDiagonal <- function(cells,intervals1,maxDist=2e5,res,filter_f=0,consider_strand='FALSE'){
	f1=basename(intervals1)
	for (cell in cells){
		file_idx <- paste0(f1,'.',maxDist,'.',res,'.',consider_strand,'.',filter_f)
		files <- list.files(paste0(main_f,'data/aggregateDiagonal/'),pattern=paste0('.',file_idx),full.names=T)
		files <- files[grep(cell,files)]
		output = paste0(main_f,'data/aggregateDiagonal/',file_idx,'.',cell)
		load(files[1])
		g = grid;
		for (i in 1:length(g$obs)) {
			g$obs[[i]] = g$obs[[i]] * g$total_obs[[i]]
			g$exp[[i]] = g$exp[[i]] * g$total_exp[[i]]
		}
		for (file_f in files[2:length(files)]) {
		 load(file=file_f)
		 for (i in 1:length(g$obs)) {
			g$obs[[i]] = g$obs[[i]] + grid$obs[[i]] * grid$total_obs[[i]]
			g$exp[[i]] = g$exp[[i]] + grid$exp[[i]] * grid$total_exp[[i]]
			g$total_obs[[i]] = g$total_obs[[i]] + grid$total_obs[[i]]
			g$total_exp[[i]] = g$total_exp[[i]] + grid$total_exp[[i]]
		 }
		}
		for (i in 1:length(g$obs)) {
		  g$obs[[i]] = g$obs[[i]] / g$total_obs[[i]]
		  g$exp[[i]] = g$exp[[i]] / g$total_exp[[i]]
		}
		grid = g
		save(grid, file=output)
	}
}

get.ctcf <- function(ctcf_chip_tracks,peaks,max_peaks=NULL,motif_v_thresh=13, bin_size=20) {
	#finding chip peaks
	ctcf_chip_v_tracks = paste0("v_", ctcf_chip_tracks)
	for (i in seq_along(ctcf_chip_tracks)) {
	  gvtrack.create(ctcf_chip_v_tracks[i], ctcf_chip_tracks[i], "global.percentile.max")
	}
	if (!is.null(max_peaks)){
		ctcf_rank <- gextract(paste0('-log2(1-',ctcf_chip_v_tracks,')'),peaks,iterator=peaks,colnames='enrichment')
		ctcf_rank <- ctcf_rank[order(ctcf_rank$enrichment,decreasing=T),]
		chip <- ctcf_rank[1:min(max_peaks,nrow(ctcf_rank)),]
	} else {
		chip <- peaks
		}
	#finding motif
	gvtrack.create("ctcf_e_forward", "motifs.CTCF_forward", "global.percentile.max")
	gvtrack.create("ctcf_e_reverse", "motifs.CTCF_reverse", "global.percentile.max")
	ctcf = gextract("-log2(1 - ctcf_e_forward)", "-log2(1 - ctcf_e_reverse)", chip, iterator=bin_size, colnames=c("forward", "reverse"))
	ctcf$pos_m = ctcf$forward >= motif_v_thresh
	ctcf$neg_m = ctcf$reverse >= motif_v_thresh

  ctcf$type = "C"
  ctcf[!is.na(ctcf$pos_m) & ctcf$pos_m, 'type'] = 'F'
  ctcf[!is.na(ctcf$neg_m) & ctcf$neg_m, 'type'] = 'R'
  ctcf[!is.na(ctcf$pos_m) & !is.na(ctcf$neg_m) & ctcf$pos_m & ctcf$neg_m, 'type'] = "B"
  ctcf_df <- ddply(ctcf,.(intervalID),function(x){
  	if ('F' %in% x$type & !('R' %in% x$type)){
  			return(c(as.character(as.vector(x$chrom[1])),min(x$start),max(x$end),'F'))
  		} else if ('R' %in% x$type & !('F' %in% x$type)) {
  			return(c(as.character(as.vector(x$chrom[1])),min(x$start),max(x$end),'R'))
  		} else if ('F' %in% x$type & ('R' %in% x$type)) {
  			return(c(as.character(as.vector(x$chrom[1])),min(x$start),max(x$end),'B'))
  		} else {
  			return(c(as.character(as.vector(x$chrom[1])),min(x$start),max(x$end),'C'))
  		}
  	})[,-1]
  	colnames(ctcf_df) <- c(colnames(peaks),'type')
  	ctcf_df$chrom <- factor(ctcf_df$chrom,levels=levels(gintervals.all()$chrom))
  return(ctcf_df)
  gvtrack.rm("ctcf_e_forward")
  gvtrack.rm("ctcf_e_reverse")
}

get.peaks <- function(chips=chip,peaks_f=peaks,max_peaks=NULL,bin_size=peaks) {
	if (!is.null(max_peaks)){
		gvtrack.create("v_chip", chips, "global.percentile")
		rank <- gextract(paste0('-log2(1-v_chip)'),peaks_f,iterator=bin_size,colnames='enrichment')
		rank <- rank[order(rank$enrichment,decreasing=T),]
		chips <- rank[1:min(max_peaks,nrow(rank)),]
		gvtrack.rm("v_chip")
	} else {
		chips <- peaks_f
		}
	return(chips)
}

sampleMatchedExpression <- function(genes,allGenes=allGenes, expression, eMargin=0.2){
rownames(allGenes) <- allGenes$geneName
genes[is.na(genes)] <- 0
genes$expression <- expression[match(genes$geneName,expression$gene_name),'FPKM']
fullSet <- allGenes
fullSet <- fullSet[!(fullSet$geneName %in% genes$geneName),]
fullSet$expression <- expression[match(fullSet$geneName,expression$gene_name),'FPKM']
fullSet[is.na(fullSet)] <- 0
rGenes <- c()
for (i in 1:nrow(genes))
{
	#print(i)
	#print(as.character(genes[i,'geneName']))
	eLevel <- as.numeric(as.character(genes[i,'expression']))
	### Select set of suitable genes
	availableGenes <- fullSet[fullSet$expression >= eLevel-eMargin*eLevel & fullSet$expression <= eLevel+eMargin*eLevel, ]
	if(nrow(availableGenes) == 0)
	{
		availableGenes <- fullSet[fullSet$expression >= eLevel-eMargin*eLevel*3 & fullSet$expression <= eLevel+eMargin*eLevel*3, ]
		rGenes <- rbind(rGenes, availableGenes[sample(nrow(availableGenes), 1),])
	}else{
		rGenes <- rbind(rGenes, availableGenes[sample(nrow(availableGenes), 1),])
	}
	### Remove selected genes from pool
	fullSet <- fullSet[-(fullSet$geneName %in% rGenes$geneName),]
# 	print(dim(fullSet))
}
row.names(rGenes) <- seq(1:nrow(rGenes))
return(rGenes)
}

###################################################################################################
#' Calculate, process and import chipseq peak files (bed, narrowPeak) into the genomic database.
#'
#' \code{process_peaks}
#'
#' This function can be used to determine genomic regions with highest enrichment in a rudiment peak calling mode, process already generated peaks (including subsetting them based on signal enrichment) and filtering for CTCF peaks if available.
#'
#' @param peaks File name of the ChIPseq interval in the format "condition"_"mark" (simplify it as much as possible - for example ES_H3K4me1). Names of conditions should match those in the hic database. Rename narrowPeak to bed.
#' @param peaks_mode One of "sharp" or "broad". When in "sharp" mode it will process user-generated peaks (such as from MACS2 or other), when in "broad" mode it will attempt to generate highly enriched regions based on the cut_off parameter. If "broad" chiptrack with the same name must exist in the database.
#' @param binSize Size of the genomic bins to consider when calculating enrichment.
#' @param cut_off  Minimum value of global percentile enrichment to consider for enriched regions. For example, value of 0.98 will return the 2 percent genomic bins with highest enrichment. Should be manually optimised for the desired set.
#' @param path_f Folder where the chipseq peaks are stored.
#' @param chip_max_peaks Number of peaks to consider when working with chip tracks. Set NULL to take all. If not NULL chiptrack with the same name must exist in the database.
#' @param chip misha location of the chiptrack to calculate enrichment
#' @param blacklist_f blacklist locations in gintervals format
#' @param ctcf_peaks full path to CTCF bed file if you want to filter for CTCF positions
#'
#'
#' @examples
#'
#' cells=c('D2')
#' for (cell in cells){
#' 	process_peaks(peaks=paste0(cell,'_CTCF.bed'),chip_max_peaks=NULL,peaks_mode='sharp')     # Imports and processes all CTCF peaks
#' 	process_peaks(peaks=paste0(cell,'_H3K4me3.bed'),chip_max_peaks=10000,peaks_mode='sharp',binSize=500)    # Imports H3K4me3 peaks and identifies top 10000 peaks based on signal enrichment.
#' 	process_peaks(peaks=paste0(cell,'_H3K27me3.bed'),chip_max_peaks=NULL,peaks_mode='broad',binSize=1000,cut_off=0.98)   # Identifies H3K27me3 enriched regions based on signal enrichment (global percentile).
#' }
#'
#' @export
##########################################################################################################
process_peaks <- function(peaks,peaks_mode='sharp',binSize=1000,cut_off=0.98,path_f=paste0(main_f,'data/peaks/'),chip_max_peaks=NULL,chip=paste0('chipseq_RPM.',gsub('.bed','',peaks)),blacklist_f=blacklist,ctcf_peaks=NULL){
  tss <- gintervals.load(tss_full)
  cell=unlist(strsplit(peaks,'_'))[1]
  files <- paste0(path_f,peaks)
  path <- paste0(path_f,'processed/')
  dir.create(path, showWarnings = FALSE,recursive=T)
  name <- peaks
  name <- gsub('.bed','',name)
  print(name)
  if (peaks_mode!='broad'){
    peaks <- read.table(files)
    peaks <- peaks[peaks[,1] %in% gintervals.all()[,1],]
    peaks <- gintervals(peaks[,1],peaks[,2],peaks[,3])
  } else {
    vName <- paste("v_",name,sep="")
    gvtrack.create(vName, chip, "global.percentile")
    peaks = gscreen(paste(vName, '> ',cut_off),iterator=binSize)
  }
  if (!is.null(blacklist_f)){
    peaks <- gintervals.neighbors(peaks,blacklist)
    peaks <- peaks[is.na(peaks$dist) | abs(peaks$dist)>=5000,1:3]
  }
  if (!is.null(chip_max_peaks)){
    peaks <- intervals.expand(peaks,binSize)
    peaks <- gintervals.canonic(peaks)
  }
  #	print(chip)
  peaks <- get.peaks(chips=chip,peaks_f=peaks,max_peaks=chip_max_peaks,bin_size=peaks)
  write.table(peaks[,1:3],paste0(path,name,if(is.null(chip_max_peaks)){'.bed'} else {paste0('_top',chip_max_peaks/1000,'k.bed')}),quote=F,col.names=F,row.names=F,sep='\t')
  peaks_noTss <- gintervals.neighbors(peaks,tss)
  peaks_noTss <- peaks_noTss[is.na(peaks_noTss$dist) | abs(peaks_noTss$dist)>5000,1:3]
  write.table(peaks_noTss[,1:3],paste0(path,name,if(is.null(chip_max_peaks)){'_noTss.bed'} else {paste0('_top',chip_max_peaks/1000,'k_noTss.bed')}),quote=F,col.names=F,row.names=F,sep='\t')
  peaks_Tss <- gintervals.neighbors(peaks,tss)
  peaks_Tss <- peaks_Tss[abs(peaks_Tss$dist)<=5000,1:3]
  write.table(peaks_Tss[,1:3],paste0(path,name,if(is.null(chip_max_peaks)){'_Tss.bed'} else {paste0('_top',chip_max_peaks/1000,'k_Tss.bed')}),quote=F,col.names=F,row.names=F,sep='\t')
  message(name)
  message('peaks:',nrow(peaks))
  message('peaks noTSS:',nrow(peaks_noTss))
  message('peaks TSS:',nrow(peaks_Tss))
  if (!is.null(ctcf_peaks)) {
    ctcf <- read.table(ctcf_peaks)
    ctcf <- gintervals(ctcf[,1],ctcf[,2],ctcf[,3])
    peaks_noCTCF <- gintervals.neighbors(peaks,ctcf)
    peaks_noCTCF <- peaks_noCTCF[is.na(peaks_noCTCF$dist) | abs(peaks_noCTCF$dist)>5000,1:3]
    peaks_noCTCF_noTss <- gintervals.neighbors(peaks_noTss,ctcf)
    peaks_noCTCF_noTss <- peaks_noCTCF_noTss[is.na(peaks_noCTCF_noTss$dist) | abs(peaks_noCTCF_noTss$dist)>5000,1:3]
    peaks_noCTCF_Tss <- gintervals.neighbors(peaks_Tss,ctcf)
    peaks_noCTCF_Tss <- peaks_noCTCF_Tss[is.na(peaks_noCTCF_Tss$dist) | abs(peaks_noCTCF_Tss$dist)>5000,1:3]
    message('peaks noCTCF:',nrow(peaks_noCTCF))
    message('peaks noCTCF noTSS:',nrow(peaks_noCTCF_noTss))
    message('peaks noCTCF TSS:',nrow(peaks_noCTCF_Tss))
    write.table(peaks_noCTCF[,1:3],paste0(path,name,if(is.null(chip_max_peaks)){'_noCTCF.bed'} else {paste0('_top',chip_max_peaks/1000,'k_noCTCF.bed')}),quote=F,col.names=F,row.names=F,sep='\t')
    write.table(peaks_noCTCF_noTss[,1:3],paste0(path,name,if(is.null(chip_max_peaks)){'_noCTCF_noTss.bed'} else {paste0('_top',chip_max_peaks/1000,'k_noCTCF_noTss.bed')}),quote=F,col.names=F,row.names=F,sep='\t')
    write.table(peaks_noCTCF_Tss[,1:3],paste0(path,name,if(is.null(chip_max_peaks)){'_noCTCF_Tss.bed'} else {paste0('_top',chip_max_peaks/1000,'k_noCTCF_Tss.bed')}),quote=F,col.names=F,row.names=F,sep='\t')
  }
}

###################################################################################################
#' Ranks and outputs gene promoters based on expression data. Also filters for CTCF and neighboring Tss.
#'
#' \code{rank_Tss}
#'
#' This function uses gene expression matrix to classify and rank gene promoters based on their expression strength. Resulting files can then be used in \code{submit_aggregateHiC} and \code{submit_aggregateDiagonal} calls.
#'
#' @param cells Conditions to work on.
#' @param tss Gene promoter coordinates in gintervals format.
#' @param path Folder where the rna-based data is stored.
#' @param fpkm_f Expression file containing gene expression information. Format should be a tab-separated file where each row is a separate gene and the columns are individual replicates. Column names should match the same names defined in the config file.
#' @param fpkm_thresh Threshold to consider gene as expressed.
#' @param ranks Quantiles in which to separate active genes based on their expression
#'
#' @examples
#'
#' rank_Tss(cells=c('D0','D2'),fpkm_thresh=1,ranks=4)
#'
#' @export
##########################################################################################################
rank_Tss <- function(cells=cells,tss=gintervals.load(tss_f),path=paste0(main_f,'data/rna/'),fpkm_f='/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt',fpkm_thresh=1,ranks=4){
	fpkm_df <- read.table(fpkm_f,header=T)
	dir.create(path, showWarnings = FALSE,recursive=T)
	for (cell in cells){
		for (state in c('active','inactive')){
			message('Working on ',cell,'_',state)
			fpkm <- fpkm_df[,grep(cell,colnames(fpkm_df),fixed=TRUE)]
			if (state=='active'){
				active <- fpkm_df[(fpkm[,1]>=fpkm_thresh)&(fpkm[,2]>=fpkm_thresh), ]
				} else {
				active <- fpkm_df[(fpkm[,1]<fpkm_thresh)&(fpkm[,2]<fpkm_thresh), ]
				}
			active_tss_byFPKM <- unique(tss[tss$geneName %in% active$gene_name,c(1:5)])
			if (file.exists(paste0(main_f,'data/peaks/processed/',cell,'_CTCF.bed'))){
				ctcf <- read.table(paste0(main_f,'data/peaks/processed/',cell,'_CTCF.bed'))
			} else {
				ctcf <- read.table(paste0(main_f,'data/peaks/processed/',cells[1],'_CTCF.bed'))
			}
			ctcf <- gintervals(ctcf[,1],ctcf[,2],ctcf[,3])
			save(active_tss_byFPKM,file=paste0(path,cell,'_',state,'Tss'))
			if (state=='active'){
				active$rank <- ntile(rowMeans(fpkm[(fpkm[,1]>=fpkm_thresh)&(fpkm[,2]>=fpkm_thresh), ],na.rm=T), ranks)
				for (i in 1:4){
					active_Q <- subset(active,rank==i)
					active_tss_Q <- unique(tss[tss$geneName %in% active_Q$gene_name,c(1:5)])
					ctcfIntersect = gintervals.neighbors(active_tss_Q, ctcf)
					peaksNoCtcf = ctcfIntersect[is.na(ctcfIntersect$dist) | abs(ctcfIntersect$dist) > 5000,1:5]
					save(active_tss_Q,file=paste0(path,cell,'_activeTss_Q',i))
					save(peaksNoCtcf,file=paste0(path,cell,'_activeTss_Q',i,'_noCTCF'))
				}
			}
			tss_filtered <- unique(adply(active_tss_byFPKM,1,function(x){
				tssIntersect = gintervals.neighbors(x,tss,mindist=-20000,maxdist=20000,maxneighbors=2)
				if(nrow(tssIntersect)==1){return(x)}
			}))
			save(tss_filtered,file=paste0(path,cell,'_',state,'Tss_noTss20kb'))
			ctcfIntersect = gintervals.neighbors(active_tss_byFPKM, ctcf)
			peaksNoCtcf = ctcfIntersect[is.na(ctcfIntersect$dist) | abs(ctcfIntersect$dist) > 5000,1:5]
			save(peaksNoCtcf,file=paste0(path,cell,'_',state,'Tss_noCTCF'))
			tss_filtered <- unique(adply(peaksNoCtcf,1,function(x){
				tssIntersect = gintervals.neighbors(x,tss,mindist=-20000,maxdist=20000,maxneighbors=2)
				if(nrow(tssIntersect)==1){return(x)}
			}))
			save(tss_filtered,file=paste0(path,cell,'_',state,'Tss_noCTCF_noTss20kb'))
			peaksNoCtcf = ctcfIntersect[is.na(ctcfIntersect$dist) | abs(ctcfIntersect$dist) > 20000,1:5]
			save(peaksNoCtcf,file=paste0(path,cell,'_',state,'Tss_noCTCF_20kb'))
			tss_filtered <- unique(adply(peaksNoCtcf,1,function(x){
				tssIntersect = gintervals.neighbors(x,tss,mindist=-20000,maxdist=20000,maxneighbors=2)
				if(nrow(tssIntersect)==1){return(x)}
			}))
			save(tss_filtered,file=paste0(path,cell,'_',state,'Tss_noCTCF20kb_noTss20kb'))
		}
	}
}

###################################################################################################
#' Ranks and outputs gene coordinates based on expression data and exon count.
#'
#' \code{rank_Exons}
#'
#' This function uses gene expression matrix to classify and rank genes based on expression data and exon count. Output is the center of the gene. Resulting files can then be used in \code{submit_aggregateHiC} and \code{submit_aggregateDiagonal} calls.
#'
#' @param cells Conditions to work on.
#' @param path Folder where the rna-based data is stored.
#' @param exon_count_f KnownGenes table from UCSC where the 1st column is the uID and 7th column represents the number of exons.
#' @param fpkm_f Expression file containing gene expression information. Format should be a tab-separated file where each row is a separate gene and the columns are individual replicates. Column names should match the same names defined in the config file.
#' @param fpkm_thresh Threshold to consider gene as expressed.
#' @param fpkm_ranks Quantiles in which to separate active genes based on their expression.
#' @param exon_ranks Quantiles in which to separate active genes based on the number of exons.
#'
#' @examples
#'
#' rank_Exons(cells=cells,fpkm_thresh=1,fpkm_ranks=4,exon_ranks=6)
#'
#' @export
##########################################################################################################
rank_Exons <- function(cells,path=paste0(main_f,'data/exons/'),exon_count_f='/work/project/Cavalli-mammals/satish/HiC/data/known_genes.txt',fpkm_f='/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt',fpkm_thresh=1,fpkm_ranks=4,exon_ranks=6){
	dir.create(path, showWarnings = FALSE,recursive=T)
	tss <- gintervals.load(tss_f)
	genes <- gintervals.load(genes_f)
	exon_count <- read.table(exon_count_f,sep='\t')

	genes <- merge(genes,exon_count[,c(1,7)],by.x='uID',by.y='V1')
	colnames(genes)[8] <- 'exon_count'
	genes <- genes[,-1]
	tss <- genes
	fpkm_df <- read.table(fpkm_f,header=T)
	for (cell in cells){
		for (state in c('active','inactive')){
			message('Working on ',cell,'_',state)
			fpkm <- fpkm_df[,grep(cell,colnames(fpkm_df),fixed=TRUE)]
			if (state=='active'){
				active <- fpkm_df[(fpkm[,1]>=fpkm_thresh)&(fpkm[,2]>=fpkm_thresh), ]
				active$FPKMrank <- ntile(rowMeans(fpkm[(fpkm[,1]>=fpkm_thresh)&(fpkm[,2]>=fpkm_thresh), ],na.rm=T), fpkm_ranks)
				} else {
				active <- fpkm_df[(fpkm[,1]<fpkm_thresh)&(fpkm[,2]<fpkm_thresh), ]
				active$FPKMrank <- 'NE'
				}
			tss_byFPKM <- unique(tss[tss$geneName %in% active$gene_name,c(1:5,7)])
			tss_byFPKM$length <- abs(tss_byFPKM[,3]-tss_byFPKM[,2])
			tss_byFPKM$rank <- ntile(tss_byFPKM$exon_count,exon_ranks)
			stats <- ddply(tss_byFPKM,.(rank),function(x){return(c(mean(x$exon_count,na.rm=T),mean(x$length,na.rm=T)))})
			tss_byFPKM <- merge(tss_byFPKM,active[,c('gene_name','FPKMrank')],by.x='geneName',by.y='gene_name',sort=F)
			for (i in unique(tss_byFPKM$rank)){
				for (j in unique(tss_byFPKM$FPKMrank)){
					temp <- tss_byFPKM[tss_byFPKM$rank==i&tss_byFPKM$FPKMrank==j,c(2:5,1)]
					message(cell,'_',state,'_fpkmRank:',j,', exonRank:',i,' : ',nrow(temp))
					save(temp,file=paste0(path,cell,'_',state,'Tss_exonQ',i,'_fpkmQ',j))
				}
			}
		}
	}
}


###################################################################submit_aggregateHiC######################################
#' Outputs equal number of genes with similar expression levels given a subset of genes and expression matrix.
#'
#' \code{sampleMatchedExpression}
#'
#' This function uses gene expression matrix to output equal number of genes with similar expression levels given a subset of genes and expression matrix. Resulting files can then be used in \code{submit_aggregateHiC} and \code{submit_aggregateDiagonal} calls.
#'
#' @param genes Genes (subset of promoters in gintervals format) to use as starting point.
#' @param allGenes Pool of genes from which to choose for matched by expression. Default is all genes.
#' @param expression Expression file containing gene expression information. Format should be a tab-separated file where each row is a separate gene and the columns are individual replicates. Column names should match the same names defined in the config file. Column "FPKM" should contains values to be used when matching expression.
#' @param eMargin Fold difference margin to accept gene substitute.
#'
#' @examples
#'
#' path <- '/work/project/Cavalli-mammals/satish/HiC/data/rna/'
#' file_f <- '/work/project/Cavalli-mammals/satish/HiC/data/peaks/senescence_genes.txt'
#' which_column <- 1
#' cellToMatch <- 'D0'
#' int_name <- paste0('SASP_matched',cellToMatch)
#'
#' tss <- gintervals.load('glp.intervals.ucscCanTSS')
#' allGenes <- gintervals.load('glp.intervals.ucscCanTSS')
#' if (grepl('txt',file_f)|grepl('bed',file_f)) {
#' 	peaks <- read.table(file_f,header=F)
#' 	genes <- tss[tss$geneName %in% peaks[,which_column],]
#' } else {
#' 	peaks <- get(load(file_f))
#' 	genes <- tss[tss$geneName %in% peaks[,which_column],]
#' }
#' expression <- read.table("/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt",header=T)
#' expression$FPKM <- rowMeans(expression[,grep(cellToMatch,colnames(expression))],na.rm=T)
#' df <- sampleMatchedExpression(genes,allGenes,expression,0.2)
#' temp <- df[,1:5]
#' save(temp, file=paste0(path,int_name))
#'
#' @export
##########################################################################################################
sampleMatchedExpression <- function(genes,allGenes=gintervals.load(tss_f), expression='/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt', eMargin=0.2){
rownames(allGenes) <- allGenes$geneName
genes[is.na(genes)] <- 0
genes$expression <- expression[match(genes$geneName,expression$gene_name),'FPKM']
fullSet <- allGenes
fullSet <- fullSet[!(fullSet$geneName %in% genes$geneName),]
fullSet$expression <- expression[match(fullSet$geneName,expression$gene_name),'FPKM']
fullSet[is.na(fullSet)] <- 0
rGenes <- c()
for (i in 1:nrow(genes))
{
	#print(i)
	#print(as.character(genes[i,'geneName']))
	eLevel <- as.numeric(as.character(genes[i,'expression']))
	### Select set of suitable genes
	availableGenes <- fullSet[fullSet$expression >= eLevel-eMargin*eLevel & fullSet$expression <= eLevel+eMargin*eLevel, ]
	if(nrow(availableGenes) == 0)
	{
		availableGenes <- fullSet[fullSet$expression >= eLevel-eMargin*eLevel*3 & fullSet$expression <= eLevel+eMargin*eLevel*3, ]
		rGenes <- rbind(rGenes, availableGenes[sample(nrow(availableGenes), 1),])
	}else{
		rGenes <- rbind(rGenes, availableGenes[sample(nrow(availableGenes), 1),])
	}
	### Remove selected genes from pool
	fullSet <- fullSet[-(fullSet$geneName %in% rGenes$geneName),]
# 	print(dim(fullSet))
}
row.names(rGenes) <- seq(1:nrow(rGenes))
return(rGenes)
}

#########################################################################################################
#' Generates h5 format files suitable for import in diffHiC
#'
#' \code{generate_diffHiC}
#'
#' This function serves as a wrapper between misha database and diffHiC R package in order to determine differential contacts between two (or more) HiC samples.
#'
#' @param cell Which condition to work on.
#' @param chr Which chromosome to work on. Trans contacts are not supported.
#' @param mode One of "none","exp" and "norm". If "none" processes only observed contacts; if "exp" also generates h5 file from the shuffled contacts; if "norm" only considers observed contacts with score above a given value.
#' @param score Minimum score value to consider when running in norm mode.
#' @param tracks Which hic tracks to consider. Only tracks containing the names given in the cell paramater will be processed.
#' @param path Folder where to store generated h5 files.
#' @param mem Available memory. Increase in memory-related errors.
#'
#' @examples
#'
#' cells <- c('D0','D2','D4','D6','D10')
#' for (cell in cells){
#' 	for (chr in chrom_sizes[,1]){
#' 		generate_diffHiC(cell=cell,chr=chr,mode='none',mem=30)
#' 	}
#' }
#' @export
##########################################################################################################
generate_diffHiC <- function(cell,chr,mode='exp',score=30,tracks=all_tracks,path=paste0(main_f,'data/diffHiC/'),mem=20){
	n_tracks <- length(tracks[grep(paste0('.',cell),tracks,fixed=T)])
	for (n_track in 1:n_tracks){
		track <- tracks[grep(paste0('.',cell),tracks,fixed=T)]
		file_f <- paste0(path,cell,'/',cell,'_',n_track,'/',as.character(chr),'.h5')
		if (!file.exists(paste0(path,cell,'/',cell,'_',n_track,'/',as.character(chr),'_shuffle.h5'))){
			command_f <- paste0('Rscript ',main_f,'/scripts/generate_diffHiC.R ',cell,' ',track[n_track],' ',n_track,' ',chr,' ',path,' ',mode,' ',score,' ',main_f,' &>',main_f,'/logs/diffHiC/',cell,n_track,'_',chr,'.log')
# 			print(command_f)
			cat(command_f,file=paste0(main_f,'/qsub/',cell,n_track,'_',chr,'_genDiffHiC.sh'))
			if (sunGrid){
				command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/diffHiC/',cell,n_track,'_',chr,'.e',' -o ',main_f,'logs/diffHiC/',cell,n_track,'_',chr,'.o ',paste0(main_f,'/qsub/',cell,n_track,'_',chr,'_genDiffHiC.sh'))
			} else {
				command_f <- paste0('nohup sh ',main_f,'/qsub/',cell,n_track,'_',chr,'_genDiffHiC.sh 2>',main_f,'/logs/diffHiC/',cell,n_track,'_',chr,'.log &')
			}	
			system(command_f)
		} else {message ('File ',file_f,' exists. Skipping...')}
	}
}

#########################################################################################################
#' Determines differential contacts between two or more hic samples using diffHiC (Lun et al.,2015)
#'
#' \code{analyze_diffHiC}
#'
#' This function serves as a wrapper between misha database and diffHiC R package in order to determine differential contacts between two (or more) HiC samples. In addition to the standard normalization described in the diffhic package (loess), the function will also attempt to normalize using expected tracks (caution - experimental).
#'
#' @param cell Which condition to work on.
#' @param chr Which chromosome to work on. Trans contacts are not supported.
#' @param binSize Resolution at which to analyze differential contacts.
#' @param path Folder where to store generated data.
#' @param read_thresh Only consider bins with total reads more than this threshold.
#' @param shuffle_thresh Only consider bins where the ratio of obs/exp reads per bin is higher than this value.
#' @param name Desired name of run.
#' @param shuffle.keep Boolean to consider enrichment over expected (paramater shuffle_thresh) when outputting results.
#' @param mem Available memory. Increase in memory-related errors.
#'
#' @examples
#'
#' binSizes=c(5e4,1e5)
#' for (chr in chrom_sizes[,1]){
#' 	for (binSize in binSizes){
#' 		analyze_diffHiC(cells=c('D0','D2','D4','D6','D10'),binSize=binSize,chr=chr,name='Oncogene',mem=40)
#' 	}
#' }
#'
#' @export
##########################################################################################################
analyze_diffHiC <- function(cells,chr,binSize=1e5,path=paste0(main_f,'analysis/diffHiC/'),read_thresh=5,shuffle_thresh=1.3,name='test',shuffle.keep=T,mem=40){
#	command_f <- "export R_LIBS=/work/bbonev/software/R-3.4.1\n"
  command_f <- ''
# 	print(command_f)
	cat(command_f,file=paste0(main_f,'/qsub/',name,'_',binSize,'_',chr,'_DiffHiC.sh'))
	command_f <- paste0('Rscript ',main_f,'/scripts/diffHiC.R ',paste0(cells,collapse=','),' ',chr,' ',binSize,' ',path,' ',read_thresh,' ',shuffle_thresh,' ',shuffle.keep,' ',name,' ',main_f,' &>',main_f,'/logs/diffHiC/',name,'_',binSize,'_',chr,'.log')
	cat(command_f,file=paste0(main_f,'/qsub/',name,'_',binSize,'_',chr,'_DiffHiC.sh'),append=T)
	print(paste0(main_f,'/logs/diffHiC/',name,'_',binSize,'_',chr,'.log'))
	if (sunGrid){
		command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/diffHiC/',name,'_',binSize,'_',chr,'.e',' -o ',main_f,'logs/diffHiC/',name,'_',binSize,'_',chr,'.o ',paste0(main_f,'/qsub/',name,'_',binSize,'_',chr,'_DiffHiC.sh'))
	} else {
		command_f <- paste0('nohup sh ',main_f,'/qsub/',name,'_',binSize,'_',chr,'_DiffHiC.sh 2>',paste0(main_f,'/logs/diffHiC/',name,'_',binSize,'_',chr,'.log'),' &')
	}	
 	system(command_f)
}

#########################################################################################################
#' Function to convert differential contacts from diffhic to Juicer 2D format
#'
#' \code{txtToJuicer}
#'
#' Converting text files with Diffhic results to juicer format for inspection.
#'
#' @param path Folder where to store generated data.
#' @param type One of "direct","independent","best","boxed" according to the diffhic manual.
#' @param logFC Filter for desired minimum absolute logFC change.
#' @param FDR Filter for desired minimum FDR value.
#'
#' @examples
#'
#' txtToJuicer(path=paste0(main_f,'analysis/diffHiC/'),type='independent',logFC=2,FDR=0.001)
#'
#' @export
##########################################################################################################
txtToJuicer <- function(path,type='independent',logFC=2,FDR=0.001){
	files <- list.files(path,recursive=T,full.names=T)
	files <- files[grep('chr',files)]
	files <- files[grep('DI',files,invert=T)]
	files <- files[grep(type,files,fixed=T)]
	df <- read.table(files[1],header=T)
	for (i in 1:length(files)){
		df <- rbind(df,read.table(files[i],header=T))
	}
	if (grepl('direct',files[1])){
		juicer_df <- df[,1:6]
		juicer_df[df$logFC>0,7]='0,255,0'
		juicer_df[df$logFC<0,7]='0,0,255'
		colnames(juicer_df) <- c('chr1','x1','x2','chr2','y1','y2','color')
		juicer_df$logFC=df$logFC
		juicer_df$FDR=df$FDR
		juicer_df$ID=row.names(df)
	} else if (grepl('independent',files[1])) {
		juicer_df <- df[,1:6]
		juicer_df[df$direction=='up',7] <- '0,255,0'
		juicer_df[df$direction=='down',7] <- '0,0,255'
		juicer_df[df$direction=='mixed',7] <- '0,255,255'
		colnames(juicer_df) <- c('chr1','x1','x2','chr2','y1','y2','color')
		juicer_df$logFC <- df$logFC.up
		juicer_df$logFC[df$logFC.down>0] <- df$logFC.down[df$logFC.down>0]
		juicer_df$FDR <- df$FDR
		juicer_df$ID <- row.names(df)
	} else if (grepl('best',files[1])) {
		juicer_df <- df[,1:6]
		juicer_df[df$direction=='up',7] <- '0,255,0'
		juicer_df[df$direction=='down',7] <- '0,0,255'
		juicer_df[df$direction=='mixed',7] <- '0,255,255'
		colnames(juicer_df) <- c('chr1','x1','x2','chr2','y1','y2','color')
		juicer_df$logFC <- df$logFC.up
		juicer_df$logFC[df$logFC.down>0] <- df$logFC.down[df$logFC.down>0]
		juicer_df$FDR <- df$FDR
		juicer_df$ID <- row.names(df)
	} else if (grepl('boxed',files[1],fixed=T)) {
		juicer_df <- df[,1:6]
		juicer_df[df$direction=='up',7] <- '0,255,0'
		juicer_df[df$direction=='down',7] <- '0,0,255'
		juicer_df[df$direction=='mixed',7] <- '0,255,255'
		colnames(juicer_df) <- c('chr1','x1','x2','chr2','y1','y2','color')
		juicer_df$logFC <- df$logFC.up
		juicer_df$logFC[df$logFC.down>0] <- df$logFC.down[df$logFC.down>0]
		juicer_df$FDR <- df$FDR
		juicer_df$ID <- row.names(df)
		if (ncol(df)>12) {
			juicer_df[!is.na(df$best.start1),2:3] <- df[!is.na(df$best.start1),c('best.start1','best.end1')]
			juicer_df[!is.na(df$best.start1),5:6] <- df[!is.na(df$best.start1),c('best.start2','best.end2')]
		}
	}
	juicer_df <- juicer_df[abs(juicer_df$logFC)>logFC,]
	juicer_df <- juicer_df[juicer_df$FDR<FDR,]
	output=gsub('chr10','Juicer',files[1])
	write.table(juicer_df,output,col.names=T,row.names=F,sep='\t',quote=F)
	output=gsub('chr10','Combined',files[1])
	write.table(df,output,col.names=T,row.names=F,sep='\t',quote=F)
}

methylation_quantiles <- function(tracks=meth_tracks,cells=c('D0','D2'),gpath=path,bin=1e4,inters=gintervals.all(), quantile_cutoff=c(95,100),method='avg',file_f='test.pdf',fig.height=5,fig.width=6,cex.legend=1){
	tracks <- as.vector(unlist(sapply(cells,function(x){tracks[grep(x,tracks)]})))
	for (t in 1:length(tracks)){
		gvtrack.create(paste0('v_',tracks[t]), tracks[t], method)
	}
	df <- gextract(tracks,inters,iterator=bin)
	df <- df[complete.cases(df),]
	stats_df <- matrix(ncol=length(cells),nrow=nrow(df))
	colnames(stats_df) <- cells
	for (i in 1:length(cells)){
		stats_df[,i] <- rowMeans(df[,grep(cells[i],colnames(df))],na.rm=T)
		}
	colors=rainbow(length(cells))
	pdf(paste0(path,file_f),height=fig.height,width=fig.width)
	plot(density(stats_df[,1],na.rm=T),col=colors[1],main='',xlab='% methylation',lwd=2)
	for (i in 2:length(cells)){
		lines(density(stats_df[,i],na.rm=T),col=colors[i],lwd=2)
	}
	legend(x='topright',legend=cells,cex=cex.legend,lwd=2,col=colors,bg='white')
	dev.off()
	for (cell in cells){
		df1 <- cbind(df[,1:3],stats_df[,grep(cell,colnames(stats_df))])
		df1$rank <- ntile(df1[,4],100)
		chosen <- df1[df1$rank>quantile_cutoff[1]&df1$rank<=quantile_cutoff[2],]
		write.table(chosen[,1:3],paste0(path,cell,'_',paste0(quantile_cutoff,collapse='_'),'.bed'),quote=F,col.names=F,row.names=F,sep='\t')
		write.table(chosen[,1:4],paste0(path,cell,'_',paste0(quantile_cutoff,collapse='_'),'.bedGraph'),quote=F,col.names=F,row.names=F,sep='\t')
	}
}

###################################################################################################
#'  prepare hic tracks for downstream processing
#'
#' \code{submit_prepareTracks}
#'
#' This function serves as a wrapper for the shaman package and generates a background model and a score track for each hic dataset.
#' In addition it also calculates the insulation score and the eigenvector value per hic sample and per condition (combined replicates).
#' If the Sungrid option is true then SGE job submission is used, otherwise jobs are run locally. See prepare_Tracks.R for more information.
#'
#' @param ins_window Genomic window to calculate insulation score per 1kb. Defaults are: 250kb for mammalian genomes and 50kb for drosophila.
#' @param k Number of nearest neighbour points to calculate hic score as detailed in the shaman package
#' @param cell Which condition to work on. Should be listed in the config file.
#' @param near_cis Size of the matrix in grid as in the shaman package.
#' @param expand_cis How much to expand the matrix in order to calculate score for border points (recommended value 1e6).
#' @param eigenBin Size of the bin to calculate eigenvector values
#' @param mem Controlling memory of the the submitted cluster job (increase of memory problems)
#'
#' @examples
#'
#'	for (cell in cells){
#'		submit_prepareTracks(ins_window=50,k=100,cell=cell,eigenBin=100,mem=80)
#'	}
#' gdb.reload()
#' gtrack.ls("hic","shuffle") #new shuffled track that was created
#' gtrack.ls("hic","score") #new score track that was created
#' gtrack.ls("hic","ins") #new insulation track that was created
#' gtrack.ls("eigen") #new eigenvector track that was created
#'
#' @export
##########################################################################################################
submit_prepareTracks <- function(ins_window=250,k=100,cell=cells[1],near_cis=5e6,expand_cis=near_cis/5,eigenBin=100,mem=40){
  command_f <- paste0('Rscript ',main_f,'/scripts/prepare_Tracks.R ',ins_window,' ',k,' ', cell,' ', near_cis,' ',main_f,' ',eigenBin,' ',expand_cis,' &>',main_f,'/logs/prepareTracks/',paste0(cell,'_',ins_window,'_k',k,'_eigen',eigenBin,'.log'))
  cat(command_f,file=paste0(main_f,'/qsub/',paste0(cell,'_',ins_window,'_k',k,'_eigen',eigenBin,'_prepareTracks.sh')))
  if (sunGrid) {
    command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/prepareTracks/',paste0(cell,'_',ins_window,'_k',k,'_eigen',eigenBin,'_prepareTracks.e'),' -o ',main_f,'logs/prepareTracks/',paste0(cell,'_',ins_window,'_k',k,'_eigen',eigenBin,'_prepareTracks.o '),main_f,'/qsub/',paste0(cell,'_',ins_window,'_k',k,'_eigen',eigenBin,'_prepareTracks.sh'))
  } else {
    command_f <- paste0('nohup sh ',main_f,'/qsub/',paste0(cell,'_',ins_window,'_k',k,'_eigen',eigenBin,'_prepareTracks.sh &'))
  }
  dir.create(paste0(main_f,'/logs/prepareTracks/'),showWarnings = FALSE,recursive=T)
  #system(command_f)
  print(command_f)
}

mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


###################################################################################################
#'  prepare hic tracks for downstream processing
#'
#' \code{submit_map3c}
#'
#' This function can map hi-c data starting from either raw uncompressed fastq files (1), gz compressed fastq files or even directly based on SRA number.
#' It runs in two major modes which are controlled by the parameter "adj". If true, reads will be split based on the ligation junction and each part of the read will be mapped separately in end-to-end Bowtie2 mode, counting also multiple contacts (the process is called chaining). If adj=FALSE the reads will instead be mapped in a local mode counting only pairwise interactions. In general we recommend to use the first mode if your reads are longer than 50bp.
#' If the Sungrid option is true then SGE job submission is used, otherwise jobs are run locally. See mapHiC.R for more information.
#'
#' @param genome Which genome to use - specified in the config file. Currently accepts mm10, hg19 and dm3.
#' @param main_f Read main hic folder from config file.
#' @param map3c Folder where map3c scripts are located. Directory map3c should be inside. Normally taken from the config file.
#' @param bowtieIndx_dir Location where bowtie2 index files are located. Index names should start with the genome name.
#' @param REseq Sequence recognized by the restriction enzyme used in this experiment. For example for DpnII use "GATC"
#' @param adj Boolean variable indicating whether to split the reads based on the presence of the ligation junction and map independently or map them in a bowtie2 local mode with soft trimming. Recommended value:TRUE - takes longer and is only really necessary for reads > 50bp but increases the number of aligned reads.
#' @param mapq All reads with mapping quality below this value will not be considered.
#' @param fastq full path to the directory containing the fastq files. If SRA entries are to be mapped instead - location of the directory where the individual SRA subfolders will be created.
#' @param track_f track prefix where the tracks will be generated. Default is hic, and all the sample folders will be generated inside this main directory.
#' @param sample_name The name of the sample which will be the same across replicates. For example "ES", "GM12878" or "eGFP".
#' @param rep_nm Name of the current replicate to map - best practise is to use "rep1-n"
#' @param mem Controlling memory of the the submitted cluster job (increase in case of memory problems)
#' @param cpu Controlling number of processors of the the submitted cluster job (increase for faster mapping)
#' @param descr Description of the track for human use. No spaces allowed, use underscore instead.
#' @param fmode Indicates what is the expected format of the raw data to map. Takes the following values: "gz" - gzipped fastq files, "SRA" - SRA record, in which case the paramater SRA needs to be also indicated. Anything else will assume that the raw data is in fastq/fq format.
#' @param SRA Number of the SRA entry which will be processed. Required if fmode=SRA
#' @param import_adj boolean indicating that if an adjacency file has been previously generated it should be just imported to the database.

#' @examples
#'
#'
#' submit_map3c(fastq='/work/project/Cavalli-mammals/boyan/HiC/data/fastq_raw/',adj=T,mapq=30,fmode='SRA',track_f='hic',sample_name='GM12878',rep_nm='rep4',SRA='SRR1658573',mem=120,cpu=8,descr='Rao2014_GM12878_HIC_rep4',import_adj=F)    # Directly map SRR1658573 record as track: hic.GM12878.GM12878_rep4
#' submit_map3c(fastq='/work/project/Cavalli-mammals/boyan/HiC/data/fastq_raw/ES_rep1/',adj=T,mapq=30,fmode='fastq',track_f='hic',sample_name='ES',rep_nm='rep1',SRA='',mem=120,cpu=8,descr='ES_rep1_G0G1',import_adj=F)    # Map a hic track starting from fastq files
#' submit_map3c(fastq='/work/project/Cavalli-mammals/boyan/HiC/data/fastq_raw/ES_rep2/',adj=T,mapq=30,fmode='gz',track_f='hic',sample_name='ES',rep_nm='rep2',SRA='',mem=120,cpu=8,descr='ES_rep2_G0G1',import_adj=F)    # Map a hic track starting from fastq gz files. Recommended
#'
#' @export
##########################################################################################################
submit_map3c <- function(genome=genome,map3c=map3c_f,bowtieIndx_dir=bowtieInd,REseq='GATC',adj=F,mapq=30,fastq=fastq_dir,track_f='hic',sample_name='test',rep_nm='rep1',mem=80,cpu=8,descr='',fmode='SRA',SRA='',import_adj=F){
	append_f=F
	if (fmode=='gz'){
		command_f <- paste0('cd ',fastq,' ;gzip -d *.gz\n')
		cat(command_f,file=paste0(main_f,'/qsub/',sample_name,'_',rep_nm,'_map3c.sh'))
		append_f=T
	} else if (fmode=='SRA'){
		fastq <- paste0(main_f,'data/fastq/',SRA)
		dir.create(fastq, showWarnings = FALSE,recursive=T)
		if (length(list.files(fastq))==0){
			wget_cmd <- paste0('ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/',substr(SRA,1,6),'/',SRA,'/',SRA,'.sra')
			command_f <- paste0('cd ',fastq,';wget ',wget_cmd,' 2>/dev/null; fastq-dump --split-files ',SRA,'.sra\n')
			cat(command_f,file=paste0(main_f,'/qsub/',sample_name,'_',rep_nm,'_map3c.sh'))
			append_f=T
		} else if (gsub("^.*\\.", "", list.files(fastq))=='sra'){
			command_f <- paste0('cd ',fastq,';fastq-dump --split-files ',list.files(fastq,pattern='sra',full.names=T),'\n')
			cat(command_f,file=paste0(main_f,'/qsub/',sample_name,'_',rep_nm,'_map3c.sh'))
			append_f=T
		}

	}
	config <- readLines(paste0(main_f ,'/scripts/map3c/',ifelse(adj,'adj_shared.conf','shared.conf')))
	config <- mgsub(c("GATC",8,36,'/work/project/Cavalli-mammals/hg19/hg19','temp_fastq','fastq_reg','sample_name','/work/project/Cavalli-mammals/boyan/HiC/src/tlsrc/temp','/work/project/Cavalli-mammals/boyan/HiC/src/tlsrc/pipeline/map3c','/work/project/Cavalli-mammals/boyan/HiC/src/tlsrc/pipeline/map3c/TG3C/import3C.pl','/work/project/Cavalli-mammals/boyan/HiC/src/tlsrc/pipeline/map3c/TG3C/combine_adjs.pl','/work/project/Cavalli-mammals/boyan/HiC/trackdb/hg19/seq/redb/','rep1'),c(REseq,cpu*2,mapq,bowtieIndx_dir,fastq,ifelse(fmode=='SRA','fastq',gsub("^.*\\.", "", gsub('.gz','',list.files(fastq)[1]))),sample_name,paste0(main_f,'temp/'),paste0(map3c,'/map3c'),paste0(map3c,'/map3c/TG3C/import3C.pl'),paste0(map3c,'/map3c/TG3C/combine_adjs.pl'),paste0(trackdb,'seq/redb/'),rep_nm),config)
	writeLines(config, con=paste0(main_f,'logs/map3c/',sample_name,'_',rep_nm,'.config'))
	command_f <- paste0('Rscript ',main_f,'/scripts/mapHiC.R ',main_f,' ',paste0(track_f,'.',sample_name),' ', paste0(main_f,'logs/map3c/',sample_name,'_',rep_nm,'.config'),' ', descr,' ',import_adj,' &>',main_f,'/logs/map3c/',sample_name,'_',rep_nm,'.log')
	cat(command_f,file=paste0(main_f,'/qsub/',sample_name,'_',rep_nm,'_map3c.sh'),append=append_f)
	if (sunGrid) {
		command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -pe parallel_fill ',cpu,' -e ',main_f,'logs/map3c/',sample_name,'_',rep_nm,'_map3c.e',' -o ',main_f,'logs/map3c/',sample_name,'_',rep_nm,'_map3c.e ',main_f,'/qsub/',sample_name,'_',rep_nm,'_map3c.sh')
 	} else {
 		command_f <- paste0('sh ',main_f,'/qsub/',sample_name,'_',rep_nm,'_map3c.sh')
 	}
 	#system(command_f)
 	dir.create(paste0(trackdb,'/tracks/',track_f,'/',sample_name),showWarnings = FALSE,recursive=T,mode = "0777")
}

###################################################################################################
#' Imports wig, bigWig or bedGraph file into the genomic database.
#'
#' \code{import_track}
#'
#' This function is used to import wig, bigWig or bedGraph file into the genomic database. Chromosome names must match with the database (generally contain the "chr" string).
#'
#' @param file_f Full path to the file to import.
#' @param track_folder Folder inside the genomic database where to store the track.
#' @param track_name Preferred name of the track inside the database. Ideally it should be in the following format: "Condition"_"Mark"
#' @param descr Optional description of the track. Use this to give more detailed description of the track.
#' @param binSize If >0 tracks are stored as continious coverage with this binsize, otherwise they are stored as sparse tracks with exact coordinates (slower). Use binSize=>10 for ChIPseq/RNAseq tracks and binSize=0 for methylation data.
#'
#' @examples
#'
#' import_track(file_f='D0_CTCF.bw',track_folder='chipseq_RPM',track_name='D0_CTCF',binSize=10)
#'
#' @export
##########################################################################################################
import_track <-function(file_f,track_folder,track_name,descr='',binSize=10){
  gdir.create(track_folder, showWarnings=F)
  name <- paste0(track_folder,'.',track_name)
  gtrack.import(track=name,description=descr,file=file_f,binsize=binSize)
}

###################################################################################################
#' Quantify contact enrichment between pairs of intervals.
#'
#' \code{submit_APA}
#'
#' This function creates pairs of intervals given 1D gintervals saved as Rdata file or bed files and calculates the contact enrichment based on the juicer .hic files. Data is saved and plotted.
#'
#' @param hic_file Full path to relevant .hic file.
#' @param domains Which domain tracks to use when constructing intra- vs inter-domain pairs. Two options - either gintervals format which will be used for all datasets and exist in the database or string which will be searched for in each condition folder. Usually generated by \code{analyzeCompartments} or \code{analyzeInsulation}. If not explicitly specified will use the domains generated by \code{analyzeInsulation}.
#' @param window_i Creates pairwise intervals as points and expand them in each direction using this value (this creates a square).
#' @param intervals1 Either gintervals or bed files with peak coordinates of the first regions.
#' @param intervals2 Either gintervals or bed files with peak coordinates of the second regions.
#' @param grid_mode Either cis or trans depending on the pairs interrogated.
#' @param min_dist Min distance to consider when creating pairs. Only applies when grid_mode is "cis".
#' @param max_dist Max distance to consider when creating pairs. Only applies when grid_mode is "cis".
#' @param res HiC resolution when calculating APA (juicer based) 
#' @param window_f How far to extend the window to calculate APA.
#' @param domain_mode Construct pairs either inside or outside TAD when grid_mode is "cis". Accepted values are either "intra" or "inter".
#' @param save_folder Folder where the output data should be saved.
#' @param mem Memory available - increase if memory issues.
#'
#' @examples
#'
#' submit_APA(hic_file='/work/project/Cavalli-mammals/satish/juicer/projects/all_hic/all_hic/OIS_D0_merged_mapq30.hic',domains='hic.D0.ins_250_domains_expanded',window_i=500000,intervals1=paste0('/work/project/Cavalli-mammals/satish/misc/SAHDs_TSS/SAHDs/SAHD_v1.bed'),intervals2=paste0('/work/project/Cavalli-mammals/satish/misc/SAHDs_TSS/SAHDs/SAHD_v1.bed'),grid_mode='trans',domain_mode='inter',min_dist=5e5,max_dist=2e6,window_f=5000000,res=25000,save_folder='SAHD_trans_250kb_D0',mem=100)
#'
#' @export
##########################################################################################################
submit_APA <- function(hic_file='',domains='ins_250_domains_expanded',window_i=5000,intervals1,intervals2,grid_mode='cis',min_dist=5e4,max_dist=2e8,res=5e3,window_f=20*res,domain_mode='inter',save_folder='test',mem=100){
	f1=basename(intervals1)
	f2=basename(intervals2)
	command_f <- paste0('Rscript ',main_f,'/scripts/apa.R ',intervals1,' ',intervals2,' ', paste0(f1,'_',f2,'.',grid_mode,'.',window_i),' ', grid_mode,' ',window_i,' ',hic_file, ' ',main_f,' ',domains,' ',min_dist,' ',max_dist,' ',res,' ',window_f,' ',domain_mode,' ',save_folder,' &>',main_f,'/logs/APA/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'.log'))
	cat(command_f,file=paste0(main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_apa.sh')))
	if (sunGrid) {
		command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/APA/',paste0(f1,'-',f2,'.',grid_mode,'.',window_i,'.e'),' -o ',main_f,'logs/APA/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'.o '),main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_apa.sh'))
 	} else {
 		command_f <- paste0('nohup sh ',main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_apa.sh 2>log.txt &'))
 	}	
	system(command_f)
}


###################################################################################################
#' Plotting function.
#'
#' \code{plotMisha}
#'
#' This function creates pairs of intervals given 1D gintervals saved as Rdata file or bed files and calculates the contact enrichment based on the juicer .hic files. Data is saved and plotted.
#'
#' @param hic_file Full path to relevant .hic file.
#' @param domains Which domain tracks to use when constructing intra- vs inter-domain pairs. Two options - either gintervals format which will be used for all datasets and exist in the database or string which will be searched for in each condition folder. Usually generated by \code{analyzeCompartments} or \code{analyzeInsulation}. If not explicitly specified will use the domains generated by \code{analyzeInsulation}.
#' @param window_i Creates pairwise intervals as points and expand them in each direction using this value (this creates a square).
#' @param intervals1 Either gintervals or bed files with peak coordinates of the first regions.
#' @param intervals2 Either gintervals or bed files with peak coordinates of the second regions.
#' @param grid_mode Either cis or trans depending on the pairs interrogated.
#' @param min_dist Min distance to consider when creating pairs. Only applies when grid_mode is "cis".
#' @param max_dist Max distance to consider when creating pairs. Only applies when grid_mode is "cis".
#' @param res HiC resolution when calculating APA (juicer based) 
#' @param window_f How far to extend the window to calculate APA.
#' @param domain_mode Construct pairs either inside or outside TAD when grid_mode is "cis". Accepted values are either "intra" or "inter".
#' @param save_folder Folder where the output data should be saved.
#' @param mem Memory available - increase if memory issues.
#'
#' @examples
#'
#' submit_APA(hic_file='/work/project/Cavalli-mammals/satish/juicer/projects/all_hic/all_hic/OIS_D0_merged_mapq30.hic',domains='hic.D0.ins_250_domains_expanded',window_i=500000,intervals1=paste0('/work/project/Cavalli-mammals/satish/misc/SAHDs_TSS/SAHDs/SAHD_v1.bed'),intervals2=paste0('/work/project/Cavalli-mammals/satish/misc/SAHDs_TSS/SAHDs/SAHD_v1.bed'),grid_mode='trans',domain_mode='inter',min_dist=5e5,max_dist=2e6,window_f=5000000,res=25000,save_folder='SAHD_trans_250kb_D0',mem=100)
#'
#' @export
##########################################################################################################
#plotMisha <- function(targetCell,targetRegion,outDir,domainsToPlot,chipTracksToExtract,rnaTracksToExtract,maxDist,zoomRegion=c(0,100),ann2D,chipYlim='',rnaYlim=c(0),pointCEX = 0.05,chipRes=10,plotScale=TRUE,cex.axis=2,conditions,binSize=1e5,scoreTrackToExtract,plotOrder=list(scores=TRUE, domains=TRUE, loops=FALSE, VP=FALSE, arcs=FALSE, genes=TRUE, rna=TRUE, chip=TRUE,axis=TRUE,ideogram=TRUE),plotRatios=list(unitHeight=120, scores=5, VP=2, loops=2.2, rna=1, chip=1, domains=0.15, genes=0.8, arcs=3, axis=0.4,ideogram=0.3),arcsToDraw,vTypeChip='avg',vTypeScore='max',main_f=main_f){

 # command_f <- paste0('Rscript ',main_f,'/scripts/plotMisha.R ',main_f,' ',targetCell,' ',paste0(targetRegion,collapse=','),' ',outDir,' ',paste0(conditions,collapse=','),' ',scoreTrackToExtract,' ',domainsToPlot,' ',chipTracksToExtract,' ',rnaTracksToExtract,' ',arcsToDraw,' ',maxDist,' ',chipRes,' ',paste0(zoomRegion,collapse=','),' ',pointCEX,' ',plotScale,' ',ann2D,' ',max_dist,' ',res,' ',window_f,' ',domain_mode,' ',save_folder,' &>',main_f,'/logs/APA/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'.log'))
#  cat(command_f,file=paste0(main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_apa.sh')))
#  if (sunGrid) {
#    command_f <- paste0('qsub -l mem=',mem,'G -l h_vmem=',mem,'G -e ',main_f,'logs/APA/',paste0(f1,'-',f2,'.',grid_mode,'.',window_i,'.e'),' -o ',main_f,'logs/APA/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'.o '),main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_apa.sh'))
 # } else {
#    command_f <- paste0('sh ',main_f,'/qsub/',paste0(f1,'_',f2,'.',grid_mode,'.',window_i,'_apa.sh'))
#  }	
#  system(command_f)
#}







