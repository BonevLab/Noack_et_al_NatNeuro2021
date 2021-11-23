library(Hmisc)
palette.breaks = function(n, colors, breaks){
  colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
  
  for(i in 2:(length(colors)) ){
    colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
  }
  colspec = c( colspec,
               colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)])
  )
  colspec
}
# col.nbreaks <<- c(15,25,45,60,85)
# col.neg <<- rev(palette.breaks(100 , c("lightgrey", "lightblue", "cyan","blue", "black"), col.nbreaks ))
# col.pbreaks <<- c(15,25,45,60,85)
# col.pos <<- palette.breaks(100 , c("lightgrey", "pink","red", "orange", "yellow"), col.pbreaks)
# col.scores <<- c(col.neg, "lightgrey", col.pos)

col.pbreaks <<- c(20,35,50,65,75,85,95)        #Original
#col.pbreaks <<- c(15,25,35,45,55,65,75)
col.pos <<- palette.breaks(200 , c("lightgrey","lavenderblush2","#f8bfda","lightcoral","red","orange","yellow"), col.pbreaks)
#col.pos <<- palette.breaks(200 , c("#ffe6e6","#ffcccc","#f8bfda","lightcoral","red","orange","yellow"), col.pbreaks)
col.nbreaks <<- c(20,35,50,65,75,85,95)
col.neg <<- rev(palette.breaks(100 , c("lightgrey", "powderblue", "cornflowerblue", "blue","blueviolet", "#8A2BE2", "#4B0082"), col.nbreaks ))
col.scores <<- c(col.neg, "lightgrey", col.pos)

linear.scores <- colorRampPalette(c(rev(c("lightgrey", "lightblue", "cyan","blue", "black")), c("pink","red", "orange", "yellow")))(200)
linear.scores.2 <- colorRampPalette(c('black','blue','lightgrey','red','yellow'))(200)

countOverlaps <- function(inData,eDistCutOff){
  
  overlaps <- apply(inData, 1, function(x)
  {
    x1 <- as.numeric(as.character(x['start1']))
    y1 <- as.numeric(as.character(x['start2']))
    x2 <- inData$start1
    y2 <- inData$start2
    dist <- sqrt((x1-x2)^2+(y1-y2)^2)
    if (min(dist[dist>0]) <= eDistCutOff)
    {
      length(dist[dist <= eDistCutOff & dist > 0])
    }else{0}
  })
  return(overlaps)
}

overlappingMates <- function(inData,eDistCutOff){
  
  mates <- apply(inData, 1, function(x)
  {
    x1 <- as.numeric(as.character(x['start1']))
    y1 <- as.numeric(as.character(x['start2']))
    x2 <- inData$start1
    y2 <- inData$start2
    dist <- sqrt((x1-x2)^2+(y1-y2)^2)
    which(dist == min(dist[dist>0]))
  })
  return(mates)
}

extractOverlappingPairs <- function(inData,eDistCutOff){
  
  mates <- apply(inData, 1, function(x)
  {
    x1 <- as.numeric(as.character(x['start1']))
    y1 <- as.numeric(as.character(x['start2']))
    x2 <- inData$start1
    y2 <- inData$start2
    dist <- sqrt((x1-x2)^2+(y1-y2)^2)
    
    overlaps <- x2[which(dist <= eDistCutOff & dist > 0)]
    
    paste(paste(x1,overlaps,sep='_'),collapse=';')
    
  })
  return(mates)
}

minDistance2D <- function(inData){
  
  distances <- apply(inData, 1, function(x)
  {
    x1 <- as.numeric(as.character(x['start1']))
    y1 <- as.numeric(as.character(x['start2']))
    x2 <- inData$start1
    y2 <- inData$start2
    dist <- sqrt((x1-x2)^2+(y1-y2)^2)
    min(dist[dist>0])
  })
  return(distances)
}


removeOverlaps2d <- function(inData){
  if(nrow(inData) == 0){return()}
  removeLines <- apply(inData, 1, function(x){
    overlappingPairs <- which(x[2] > as.numeric(as.character(inData$start1)) & x[2] < as.numeric(as.character(inData$end1)) & x[5] >= as.numeric(as.character(inData$start2)) & x[5] <= as.numeric(as.character(inData$end2)))
    if (length(overlappingPairs) > 0){return(overlappingPairs)}
  })
  if(!(is.null(removeLines))){
    return(inData[-unique(unlist(removeLines)),])
  }else{return(inData)}
}

makeIter2D_TSS <- function(currentCoordinates,maxDist,binSize){
  center <- gintervals.force_range(gintervals(currentCoordinates$chrom, currentCoordinates$start-binSize/2, currentCoordinates$start+binSize/2))
  upDiag <- giterator.intervals(intervals=gintervals.force_range(data.frame(chrom=center$chrom,start=center$start-maxDist, end=center$start)), iterator=binSize)
  downDiag <- giterator.intervals(intervals=gintervals.force_range(data.frame(chrom=center$chrom,start=center$end, end=center$end+maxDist)), iterator=binSize)
  ###########################
  marginals <- gintervals.force_range(gintervals.rbind(upDiag,downDiag))
  ###########################
  iter2d <- gintervals.force_range(gintervals.2d(center$chrom, center$start, center$end, marginals$chrom, marginals$start, marginals$end))
  iter2d <- gintervals.canonic(iter2d)
  return(iter2d)
}

makeIter2D <- function(interv,binSize){
  tIntrv <- gintervals.2d(chr)
  iter2d <- giterator.intervals(intervals=tIntrv, iterator=c(binSize,binSize))
  return(iter2d)
}

plot2D_binned <- function(set1){
  plot((set1$start1 + set1$start2)/2, (set1$start2 - set1$start1), pch=19, col=col.scores[101+set1[,7]], xlab='', ylab='', yaxt='n', xaxt='n',frame=F, cex=10)
}

plotScoreLegend <- function(x){
  
  plot(1:200,y=rep(1,200),type='p', lwd=3, col=col.scores[101+-100:100], pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F)
  axis(1,labels=c(-100,0,100), at=c(1,100,200), tick=F, padj=-2)
}

plotSingle <- function(set1, header, plotGenes, plotAnn=F,pointCEX=0.05){
  if(is.null(set1)){return()}
  set1 <- subset(set1, set1$start1 < set1$start2)
  ########################################################
  ########################################################
  plotLim <- c(min(set1$start1),max(set1$end2))
  yLim <- c(0, max((set1$start2 - set1$start1)))
  ########################################################
  winSize <- plotLim[2]-plotLim[1];
  set1 <- set1[order(set1[,7]),]
  
  #if(nrow(set1) < 50000){pointCEX <- 0.1;}
  plot((set1$start1 + set1$start2)/2, (set1$start2 - set1$start1), xlim=c(plotLim[1],plotLim[2]), ylim=yLim, pch=19, col=col.scores[101+set1[,7]], xlab='', ylab='', yaxt='n', xaxt='n',frame=F, cex=pointCEX)
  lines(x=c(plotLim[1],plotLim[2]),y=c(0,0), col=1, lwd=2)
  text(x=(plotLim[1]+(plotLim[2]-plotLim[1])/5),y=yLim[2]-yLim[2]/4, labels=header, cex=2.4, col=1, pos=3, offset=1)
  axis(1, cex.axis=1.6)
  
  lineCheck <- which(plotGenes$geneName == header)
  if(length(lineCheck) > 0)
  {
    row <- plotGenes[lineCheck,]
    if(row$strand == 1){currentX=row$start;
    }else{currentX=row$end;}
    
    x1 <- currentX
    x2 <- x1+(plotLim[2]-x1)/2
    y2 <- (x2-x1)*2
    segments(x1,0,x2,y2,lwd=1.5)
    
    x2 <- plotLim[1]+ (x1-plotLim[1])/2
    y2 <- (x1-x2)*2
    segments(x1,0,x2,y2,lwd=1.5)
    
    print(row)
  }
  
  if(is.data.frame(plotAnn))
  {
    for(g in 1:nrow(plotAnn))
    {
      row <- plotAnn[g,]
      x1 <- row$start
      x2 <- x1+(plotLim[2]-x1)/2
      y2 <- (x2-x1)*2
      segments(x1,0,x2,y2,lwd=0.75)
      
      x2 <- plotLim[1]+ (x1-plotLim[1])/2
      y2 <- (x1-x2)*2
      segments(x1,0,x2,y2,lwd=0.75)
    }
  }
  
  if (nrow(plotGenes) > 0){genesPlot(plotGenes, plotLim, 'F', header)
  }else{genesFake1(plotLim,unique(set1$chrom1))}
  ########################################################
  ########################################################
}

plotSimpleBED <- function(set, plotMar=c(3.5,3), plotLim,cols=c('purple','#D4AF37')){
  par(mar = c(0, plotMar[1], 0, plotMar[2]))
  
  plot(0 ,ylim=c(0,50), xlim=plotLim, xlab='', ylab='', yaxt='n', xaxt='n',frame=F)
  
  for (i in 1:nrow(set))
  {
    rect(set[i,2], 0, set[i,3], 50, col=cols[i], border = NA)
    text(set[i,4], x=set[i,2]+(set[i,3]-set[i,2]),y=20, pos=4, cex=1.3, las=2)
  }
}


plotBedGraph_hm <- function(df,track, plotMar=c(3.5,3), plotLim,yLim,cex.axis,setName,cols=c('purple','#D4AF37'),idx,figure_mode=FALSE){
  set <- df[,c('chrom','start','end',track)]
  gen_color <- scales::col_numeric(palette = cols,domain = yLim,na.color = 'white')
  set$col <- gen_color(set[,4])
  
  par(mar = c(0.1,plotMar[1],0.1,plotMar[2]))
  
  #par(mar = c(0, plotMar[1], 0, plotMar[2]))
  
  plot(0 ,ylim=c(0,50), xlim=plotLim, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,yaxs='i',xaxs='i')
  
  for (i in 1:nrow(set))
  {
    rect(set[i,2], 0, set[i,3], 50, col=set[i,'col'], border = NA)
    #text(set[i,4], x=set[i,2]+(set[i,3]-set[i,2]),y=20, pos=4, cex=1.3, las=2)
  }
  mtext(setName, side = 4, outer = F,line=0,at=25,adj=1,las=2,cex=cex.axis)
}


plotSingleKNN <- function(set1, plotAnn=F, plotMar=c(3.5,3), plotLim, plot2Dann=F, plotThr=-101, highlightCircle=FALSE, highRad=25e3,pointCEX=0.05,plotScale=TRUE,cex.axis=2,window_scale=1){
  par(mar = c(0, plotMar[1], 0, plotMar[2]))
  if(is.null(set1)){return()}
  if(window_scale==1){
    set1 <- subset(set1, set1$start1 > plotLim[1] & set1$start2 < plotLim[2])
  } else {
    set1 <- subset(set1, (set1$start1 + set1$start2)/2 > plotLim[1] & (set1$start1 + set1$start2)/2 < plotLim[2])
  }
  set1 <- subset(set1, set1$start1 < set1$start2)
  set1 <- subset(set1, set1[,7] > plotThr)
  ########################################################
  ########################################################
  winSize <- plotLim[2]-plotLim[1]
  yLim <- c(0,winSize)
  # 	yLim <- c(0, max((set1$start2 - set1$start1)))
  ########################################################
  set1 <- set1[order(set1[,7]),]
  #	if(nrow(set1) < 50000){pointCEX <- 0.1;}
  
  # 	if(nrow(set1) > 50000)
  # 	{
  # 		set1 <- set1[sample(nrow(set1), 50000),]
  # 	}
  plot_label <- colnames(set1)[7]
  plot_label <- gsub('.score_k200','',plot_label)
  plot_label <- gsub('.score_k100','',plot_label)
  plot_label <- gsub('hic.','',plot_label)
  plot_label <- gsub('E14_','',plot_label)
  plot((set1$start1 + set1$start2)/2, (set1$start2 - set1$start1), xlim=c(plotLim[1],plotLim[2]), ylim=yLim/window_scale, pch=19, col=col.scores[101+set1[,7]], xlab='', ylab='', yaxt='n', xaxt='n',frame=F, cex=pointCEX,yaxs='i',xaxs='i')
  lines(x=c(plotLim[1],plotLim[2]),y=c(0,0), col=1, lwd=2)
  #axis(2, labels=plot_label, at=mean(yLim), las=2, tick=FALSE, cex.axis=1.6,hadj=0.5,line=1,font=2)
  text(plot_label,x=plotLim[1]+0.9*(plotLim[2]-plotLim[1]),y=0.75*winSize,font=2,cex=3)
  if(plotScale){	
    points(y=seq(yLim[2]/2,yLim[2], length.out=200),x=rep(plotLim[1],200),type='p', cex=3, col=col.scores[101+(-100:100)], pch=15)
    text('-100', x=plotLim[1],y=yLim[2]*0.49, pos=4, cex=cex.axis,offset=1)
#    text('-50', x=plotLim[1],y=yLim[2]*0.625, pos=4, cex=cex.axis,offset=1)
    text('0', x=plotLim[1],y=yLim[2]*0.75, pos=4, cex=cex.axis,offset=2)
#    text('+50', x=plotLim[1],y=yLim[2]*0.875, pos=4, cex=cex.axis,offset=1)
    text('+100', x=plotLim[1],y=yLim[2], pos=4, cex=cex.axis,offset=1)
  }
  # 	axis(1, cex.axis=1.6)
  
  ############ ANNOTATION
  if(is.data.frame(plotAnn))
  {
    plotAnn <- plotAnn[as.character(plotAnn$chrom) == as.character(set1[1,'chrom1']) & plotAnn$start > plotLim[1] & plotAnn$end < plotLim[2],]
    
    for(g in 1:nrow(plotAnn))
    {
      row <- plotAnn[g,]
      x1 <- row$start
      x2 <- x1+(plotLim[2]-x1)/2
      y2 <- (x2-x1)*2
      segments(x1,0,x2,y2,lwd=0.75)
      
      x2 <- plotLim[1]+ (x1-plotLim[1])/2
      y2 <- (x1-x2)*2
      segments(x1,0,x2,y2,lwd=0.75)
    }
  }
  
  if(is.data.frame(plot2Dann))
  {
    if(nrow(plot2Dann) > 0){annotate2D(plot2Dann,plotLim)} #function(track, plotLim, vertical=FALSE)
  }
  
  if(is.data.frame(highlightCircle))
  {
    if(nrow(highlightCircle) > 0)
    {
  #    print('trying to plot...')
      for(g in 1:nrow(highlightCircle))
      {
        row <- highlightCircle[g,]
        x <- row$start1+(row$start2-row$start1)/2
        y <- abs(row$start2-row$start1)
        draw.circle(x, y, highRad, nv = 1000, border = 'darkslategrey', col = NA, lty = 3, lwd = 2)
      }
    }
  }
  
}

plotBinnedKNN <- function(set1, plotAnn=F, plotMar=c(3.5,3), plotLim){
  par(mar = c(1, plotMar[1], 2, plotMar[2]))
  if(is.null(set1)){return()}
  set1 <- subset(set1, set1$start1 > plotLim[1] & set1$start2 < plotLim[2])
  set1 <- subset(set1, set1$start1 < set1$start2)
  
  # 	require(reshape2)
  
  try(library(reshape2), silent=TRUE)
  
  sqMat <- acast(set1, start1~start2, value.var=colnames(set1)[7])
  sqMat <- sqMat+101
  
  flatMat <- sqMat
  flatMat[is.na(flatMat)] <- 1
  flatMat[flatMat > 0] <- 0
  
  surf.colors <- function(x, col = colorRampPalette(c('black','blue','lightgrey','red','yellow'))(200)) {
    x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] + x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
    colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
    return(colors)
  }
  # 	require('plot3D')
  persp(as.matrix(flatMat), zlim=c(-0.1,0.1), r=1000, theta = 45, phi = 90, col = surf.colors(as.matrix(sqMat)), border=NA, box=T)
  # 	persp(as.matrix(flatMat), zlim=c(-0.1,0.1), r=1000, col = surf.colors(as.matrix(sqMat)), border=NA, box=T)
  
}

plotViewpoint <- function(hicNames, plotLim, targetGene, extReads, plotMar=c(3.5,3)){
  par(mar = c(0, plotMar[1], 0, plotMar[2]))
  yLim <- c(50, (length(hicNames)+1)*50)
  
  # 	plot(0,xlim=plotLim,ylim=yLim, main=targetGene, cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
  plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
  legend <- c()
  for(i in 1:length(hicNames))
  {
    hic <- hicNames[i]
    cell <- unlist(strsplit(gsub("ncx.","ncx_",hic),"\\."))[2]
    cell <- gsub("_score_k200",'',cell)
    cell <- gsub("ncx_Hes5",'ncxNPC',cell)
    cell <- gsub("ncx_Dcx",'ncxCN',cell)
  #  print(cell)
    legend <- append(legend, cell)
  #  print(paste('PLOTTING',hic,'....'))
    
    plotReads <- subset(extReads, extReads$start2 > plotLim[1] & extReads$end2 < plotLim[2])
    plotReads <- plotReads[order(abs(plotReads[,hic])),]
    rectCol <- col.scores[101+plotReads[,hic]]
    
    for(j in 1:nrow(plotReads))
    {
      row <- plotReads[j,]
      rect(row[5], i*50, row[6], i*50+50, col=rectCol[j], border = NA,yaxs='i',xaxs='i')
    }
    # 		abline(h=i*50, lwd=1.5, col=1)
    segments(plotLim[1],i*50,plotLim[2],i*50,lwd=1.75, col=1,yaxs='i',xaxs='i')
    # 		text(y=i*50+25, x=plotLim[1], label=cell, pos=2, offset=0, cex=1.6)
    
    # 		axis(3, cex.axis=1.6)
  }
  # 	axis(2)
  axis(2, labels=legend, at=seq(50,50*length(hicNames),length.out = length(hicNames))+25, las=2, tick=FALSE, cex.axis=1.6,hadj=0.5)
}

plotLoops <- function(hicNames, plotLim, targetGene, extReads, plotMar=c(3.5,3)){
  par(mar = c(0, plotMar[1], 0, plotMar[2]))
  plot(0,xlim=plotLim,ylim=c(0,1),cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
  
  for(i in 1:length(hicNames))
  {
    hic <- hicNames[i]
    cell <- unlist(strsplit(hic,"\\."))[2]
    plotReads <- subset(extReads, extReads$start2 > plotLim[1] & extReads$end2 < plotLim[2])
    plotReads[is.na(plotReads)] <- 0
    
    loops <- subset(plotReads, plotReads[,hic] > loopThr)
    if(nrow(loops) == 0){next}
    
    loopStart <- plotReads[1,'start1']+binSize/2
    loopEnds <- loops[,'start2']+binSize/2
    loopScores <- loops[,hic]/10
    lwdCex <- 3/(maxDist*1e-5)
    if(lwdCex < 0.1){lwdCex <- 0.5}
    
    for (loop in 1:length(loopEnds))
    {
      radius <- (loopEnds[loop]-loopStart)/2
      center <- loopStart + radius
      degree <- 180
      if(center < loopStart){degree <- -180}
      col <- chooseColor(cell)
      # 			draw.arc(x=center, y=0, radius=radius, deg2 = degree, lwd=lwdCex*loopScores[loop], col=col)
      draw.arc(x=center, y=0, radius=radius, deg2 = degree, lwd=lwdCex, col=col,yaxs='i',xaxs='i')
    }
  }
}

plotArcs <- function(loopSet, setOrder, plotLim, plotMar=c(3.5,3), loopThr=40, colData=FALSE,ltyData=FALSE, widthCorrection=20){
  par(mar = c(1, plotMar[1], 1, plotMar[2]))
  plot(0,xlim=plotLim,cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
  loopSet[is.na(loopSet)] <- 0
  
  plotSets <- unlist(strsplit(setOrder,'\\|'))
  
  binSize <- loopSet[1,3]-loopSet[1,2]
  
  print(plotSets)
  lSets <- c()
  for(hic in plotSets)
  {
    cell <- unlist(strsplit(hic,"\\."))[2]
    lSets <- append(lSets,cell)
    if(colData[1] != FALSE){col=colData[which(plotSets==hic)]
    }else{col <- chooseColor(hic)}
    
    if(ltyData[1] != FALSE){lty=ltyData[which(plotSets==hic)]
    }else{lty <- which(plotSets==hic)}
    
    setName <- paste0('v_',hic)
    loops <- subset(loopSet, loopSet[,setName] > loopThr)
    if(nrow(loops) == 0){next}
    
    for (i in 1:nrow(loops))
    {
      start <- loops[i,'start1']+binSize/2
      end <- loops[i,'start2']+binSize/2
      radius <- (end-start)/2
      center <- start+radius
      degree <- 180
      lwdCex <- 2
      rgbCol <- 'black'
      # 			lwdCex <- 2
      draw.arc(x=center, y=0, radius=radius, deg2 = degree, lwd=lwdCex, col=rgbCol, lty=lty,yaxs='i',xaxs='i',asp=1)
    }
  }
  legend("topleft", lty=ltyData, col=colData, lwd=4, cex=1.4, legend=lSets)
}

plotArcs2 <- function(arcIntervals,cols,plotLim){
  arcIntervals$center <- abs((arcIntervals$start1+arcIntervals$end1)/2 + (arcIntervals$start2+arcIntervals$end2)/2)/2
  arcIntervals$radius <- abs((arcIntervals$start2+arcIntervals$end2)/2-(arcIntervals$start1+arcIntervals$end1)/2)/2
  p <- ggplot(arcIntervals) + scale_x_continuous(limits=plotLim,expand = c(0,0)) + geom_arc(aes(x0 = center, y0 = 0, r = radius, start = 0, end = pi*2,color=Correlation)) 
  p <- p + scale_color_gradientn(colors=cols,name='',breaks=c(round(min(arcIntervals$Correlation),2),floor(max(arcIntervals$Correlation)*100)/100),labels=c(round(min(arcIntervals$Correlation),2),floor(max(arcIntervals$Correlation)*100)/100))
  p <- p + theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line=element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"))
  p <- p+ theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1, 1),legend.background =element_blank(),legend.key = element_blank())
  return(p + ylim(c(0,(plotLim[2]-plotLim[1]))))
}

plotArcs <- function(arcIntervals,plotMar=c(3.5,3),color_by,cols=colorRampPalette(c('blue','white','red')),plotLim,flip=F,plotScale=T,figure_mode){
  par(mar = c(ifelse(figure_mode,0.1,1), plotMar[1], ifelse(figure_mode,0.1,1), plotMar[2]))
  plotBedpe(bedpedata=arcIntervals, chrom=as.character(arcIntervals$chrom1[1]), chromstart=plotLim[1], chromend=plotLim[2], heights=1, colorby = arcIntervals[,color_by], colorbycol = cols, colorbyrange = c(min(arcIntervals[,color_by],na.rm=T),max(arcIntervals[,color_by],na.rm=T)), border = 'black', lwdby = arcIntervals[,color_by], lwdrange = c(1, 3), offset = 0, flip = flip, lwd = 1, xaxt = "n", yaxt = "n", bty = "n", plottype = "loops", maxrows = 10000, height = 1, ymax = 1.04,xaxs="i",yaxs="i")
  if(plotScale){
    pnts = cbind(x =c(plotLim[1]+(plotLim[2]-plotLim[1])*0.98,plotLim[1]+(plotLim[2]-plotLim[1])*0.995,plotLim[1]+(plotLim[2]-plotLim[1])*0.995,plotLim[1]+(plotLim[2]-plotLim[1])*0.98), y =c(0.85,0.85,0.2,0.2))
    SDMTools::legend.gradient(pnts=pnts, cols = cols(100),title ='',limits=c('',''))
    text(round(min(arcIntervals[,color_by],na.rm=T),2),x=plotLim[1]+(plotLim[2]-plotLim[1])*0.985,y=0.1,cex=0.8,offset=-2)
    text(round(max(arcIntervals[,color_by],na.rm=T),2), x=plotLim[1]+(plotLim[2]-plotLim[1])*0.985,y=0.95, cex=0.8,offset=-2)
  }
  }


loadDomains <- function(tracks, plotLim, intrv, plotMar=c(3.5,3),vertical=FALSE){
  for (set in unlist(strsplit(tracks, "\\|")))
  {
    domainCoordinates <- gintervals.load(set)
    domainCoordinates$cluster=1
    id <- unlist(strsplit(set, "\\."))[2]
    plotDomains <- subset(domainCoordinates, domainCoordinates$chrom == as.character(intrv$chrom) & domainCoordinates$end > intrv$start & domainCoordinates$start < intrv$end)
    domainsPlot(plotDomains, plotLim, id, plotMar, vertical)
  }
}

annotate2D <- function(track, plotLim, vertical=FALSE, points=FALSE){
  
  for(g in 1:nrow(track))
  {
    row <- track[g,]
    lwd=2
    col=rgb(0, 0, 0, alpha=0.6)
    
    if ('type' %in% names(row)){if(row$type == 'tss2dst'){col='3'}}
    # 		if(points==TRUE)
    # 		{
    # 			points(x=(row$start+(row$end-row$start)/2), y=(row$end-row$start), pch=1, col=1, cex=3)
    # 		}else{
    # 		next
    if(row$start < plotLim[1])
    {
      x1 <- row$end
      y1 <- 0
      x2 <- plotLim[1]+(row$end-plotLim[1])/2
      y2 <- (row$end-plotLim[1])
      segments(x1,y1,x2,y2,lwd=lwd, lty=5, col=col)
      segments(plotLim[1], 0, row$end, 0, lwd=lwd, lty=5, col=col)
    }else if(row$end > plotLim[2])
    {
      x1 <- row$start
      y1 <- 0
      x2 <- x1+(plotLim[2]-x1)/2
      y2 <- (plotLim[2]-row$start)
      segments(x1,y1,x2,y2,lwd=lwd, lty=5, col=col)
      segments(row$start,0,plotLim[2],0, lwd=lwd, lty=5, col=col)
      
    }else
    {
      segments(row$start,0,row$end,0,lwd=lwd, lty=5, col=col)
      x1 <- row$start
      y1 <- 0
      x2 <- x1+(row$end-x1)/2
      y2 <- (row$end-row$start)
      segments(x1,y1,x2,y2,lwd=lwd, lty=5, col=col)
      
      x1 <- row$end
      y1 <- 0
      x2 <- row$start+(row$end-row$start)/2
      y2 <- (row$end-row$start)
      segments(x1,y1,x2,y2,lwd=lwd, lty=5, col=col)
    }
    # 		}
  }
}


domainsPlot <- function(set1, plotLim, id, plotMar=c(3.5,3), vertical=FALSE){
  
  rgbCol1 <- col2rgb('orange')
  rgbCol2 <- col2rgb('purple')
  
 # print('domainsPlot called...')
  
  if(vertical==TRUE)
  {
    par(mar = c(plotMar[2], 0, plotMar[1],0))
    xLim=c(-50,50)
    plot(0, xlim=xLim, ylim=plotLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F,yaxs='i',xaxs='i')
    segments(0,plotLim[1],0,plotLim[2],lwd=1.75, col='darkgrey')
    for(g in 1:nrow(set1))
    {
      row <- set1[g,]
      # 			lwdCex <- 2
      if(row$cluster == 2){rect(-5, row$start, -50, row$end, col=rgb(rgbCol1[1,1]/255,rgbCol1[2,1]/255,rgbCol1[3,1]/255, alpha=0.9), border = 1); currentX=row$start; textY <- 50; textPos <- 3;
      }else{rect(5, row$start, 50, row$end, col=rgb(rgbCol2[1,1]/255,rgbCol2[2,1]/255,rgbCol2[3,1]/255, alpha=0.9), border = 1); currentX=row$end; textY <- -50; textPos <- 1;}
    }
    axis(3, labels=id, at=0, las=2, tick=FALSE)
  }else{
    par(mar = c(0, plotMar[1], 0, plotMar[2]))
    yLim=c(-50,50)
    plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F,yaxs='i',xaxs='i')
    segments(plotLim[1],0,plotLim[2],0,lwd=1.75, col='darkgrey')
    for(g in 1:nrow(set1))
    {
      row <- set1[g,]
      if(row$cluster == 2){rect(row$start, 5, row$end, 50, col=rgb(rgbCol1[1,1]/255,rgbCol1[2,1]/255,rgbCol1[3,1]/255, alpha=0.9), border = 1); currentX=row$start; textY <- 50; textPos <- 3;
      }else{rect(row$start, -5, row$end, -50, col=rgb(rgbCol2[1,1]/255,rgbCol2[2,1]/255,rgbCol2[3,1]/255, alpha=0.9), border = 1); currentX=row$end; textY <- -50; textPos <- 1;}
    }
    axis(2, labels=id, at=0, las=2, tick=FALSE)
  }
}

genesPlot <- function(set1, plotLim, vertical=FALSE, header='test', plotMar=c(3.5,3)){
  
  if(vertical == TRUE)
  {
    par(mar = c(plotMar[2], 0, plotMar[1], 0))
    xlim=c(-150,150)
    plot(0, xlab='', ylab='', yaxt='n', xaxt='n', frame=F, ylim=plotLim, xlim=xlim)
    print(plotLim)
    segments(0,plotLim[1],0,plotLim[2],lwd=1.75, col='darkgrey')
    if(nrow(set1) == 0){print(1);return()}
    for(g in 1:nrow(set1))
    {
      row <- set1[g,]
      #print(row)
      if(row$strand == -1){rect(5,row$start, 50, row$end, col=rgb(1,1,0, alpha=0.9), border = 1); textY=row$end; currentX <- 125; textPos <- 3;
      }else{rect(-5, row$start, -50, row$end, col=rgb(0,1,0, alpha=0.5), border = 1); textY =row$start; currentX <- -125; textPos <- 1;}
      ### text(x=currentX,y=textY, labels=row$geneName, pos=textPos, cex=1.2, srt = 90)
    }
    ### text(x=0,y=plotLim[2], labels=unique(set1$chrom), pos=3, cex=1.4, col='grey15')
  }else{
    par(mar = c(0, plotMar[1], 0, plotMar[2]))
    yLim=c(-150,150)
    plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F)
    segments(plotLim[1],0,plotLim[2],0,lwd=1.75, col='darkgrey')
    if(nrow(set1) == 0){print(2);return()}
    for(g in 1:nrow(set1))
    {
      row <- set1[g,]
      if(row$geneName == header){rect(row$start, yLim[1], row$end, yLim[2], col=rgb(0.66,0.66,0.66, alpha=0.4), border = 0)}
      if(row$strand == 1){rect(row$start, 5, row$end, 50, col=rgb(0,1,0, alpha=0.9), border = 1); currentX=row$start; textY <- 50; textPos <- 3;
      }else{rect(row$start, -5, row$end, -50, col=rgb(1,1,0, alpha=0.9), border = 1); currentX=row$end; textY <- -50; textPos <- 1;}
      ### text(x=currentX,y=textY, labels=row$geneName, pos=textPos, cex=1.2)
    }
    # 		text(x=plotLim[1],y=0, labels=unique(set1$chrom), pos=2, cex=1.4, col='grey15')
  }
}

plotAxis <- function(plotLim, plotMar, chr, orientation='DOWN', vertical=FALSE,cex.axis,figure_mode=FALSE){
  if(vertical == TRUE)
  {
    par(mar = c(plotMar[2], 0, plotMar[1], 0))
    xlim=c(0,0)
    plot(0, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F, ylim=plotLim, xlim=xlim)
    #		axis(2, labels=FALSE, cex.axis=1.6)
    
    xAxis <- round_any(seq(plotLim[1],plotLim[2],length.out=2),1e5)/1e6
    xPoints <- seq(plotLim[1],plotLim[2],length.out=length(xAxis))
    axis(2,labels=paste0(xAxis,'Mb'), at=xPoints, tick=T, las=1, cex.axis=cex.axis,padj=-0.3)
    
    # 		text(x=0,y=plotLim[2], labels=chr, pos=3, cex=1.4, col='grey15')
  }else{
    if(orientation == 'DOWN'){orientation=1; axMar <- c(1.3,0)
    }else{orientation=3; axMar <- c(0,1.3)}
    
    par(mar = c(axMar[1], plotMar[1], axMar[2], plotMar[2]))
    yLim=c(0,0)
    plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F)
    # 		axis(orientation, labels=FALSE, cex.axis=1.6)
    #  		axis(2, labels=chr, at=0, las=2, tick=FALSE, cex.axis=2.5)
    
    xAxis <- paste0(round_any(seq(plotLim[1],plotLim[2],by=1e5),1e5)/1e6,'Mb')
    xAxis[2:(length(xAxis)-1)] <- ''
    xPoints <- seq(plotLim[1],plotLim[2],length.out=length(xAxis))
    axis(1,labels=xAxis, at=xPoints, tck=0.4, las=1, cex.axis=cex.axis,mgp=c(3,0,0))
#    minor.tick(nx=10, ny=10, tick.ratio=3)
  #  axis(2,labels=xAxis, at=xPoints, tick=T, las=2)
  #  text(x=0,y=plotLim[1], labels=chr, pos=3, cex=1.4, col='grey15')
    axis(2, labels=chr, at=0, las=2, tick=FALSE, cex.axis=cex.axis,hadj=0.6,padj=1,col.axis='grey15')
    
  }
}

plotFake <- function(plotLim, chr, vertical=FALSE, plotMar, reps){
  if(reps == TRUE){reps=1}
  for(i in 1:reps)
  {
    print(reps)
    if(vertical == TRUE)
    {
    #  par(mar = c(plotMar[2], 0, plotMar[1], 0))
      xlim=c(-150,150)
      plot(0, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F, ylim=plotLim, xlim=xlim)
      segments(0,plotLim[1],0,plotLim[2],lwd=1.75, col='darkgrey')
      text(x=0,y=plotLim[2], labels=chr, pos=3, cex=1.4, col='grey15')
    }else{
    #  par(mar = c(0, plotMar[1], 0, plotMar[2]))
      yLim=c(-150,150)
      plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F,yaxs='i',xaxs='i')
      segments(plotLim[1],0,plotLim[2],0,lwd=1.75, col='darkgrey')
      text(x=plotLim[1],y=0, labels=chr, pos=2, cex=1.4, col='grey15')
    }
  }
}

plotSquare <- function(set1, header, plotGenes, interval, plotMar=c(3.5,3),pointCEX=0.05,cex.axis=2){
  par(mar = c(plotMar[2],plotMar[1], plotMar[1], plotMar[2]))
  if(is.null(set1)){return()}
  ########################################################
  ########################################################
  xLim <- c(interval$start1,interval$end1)
  yLim <- c(interval$start2,interval$end2)
  
  xDim <- xLim[2]-xLim[1]
  yDim <- yLim[2]-yLim[1]
  
  xMid <- xLim[1]+xDim/2
  yMid <- yLim[1]+yDim/2
  
  if(xDim > yDim){yLim <- c(yMid-xDim/2,yMid+xDim/2)
  }else if((xDim < yDim)){xLim <- c(xMid-yDim/2,xMid+yDim/2)}
  
#  print(paste(xDim, yDim))
  
  ########################################################
  set1 <- set1[order(abs(set1[,7])),]
  #if(nrow(set1) < 50000){pointCEX <- 0.5;}
  plot(set1$start1, set1$start2, ylim=yLim, xlim=xLim, pch=19, col=col.scores[101+set1[,7]], xlab='', ylab='', yaxt='n', xaxt='n',frame=F, cex=pointCEX)
  
  xAxis <- round_any(seq(xLim[1],xLim[2],length.out=2),1e5)/1e6
  xPoints <- seq(xLim[1],xLim[2],length.out=length(xAxis))
  yAxis <- round_any(seq(yLim[1],yLim[2],length.out=2),1e5)/1e6
  yPoints <- seq(yLim[1],yLim[2],length.out=length(xAxis))
  
  axis(1, labels=paste0(xAxis,'Mb'), at=xPoints, cex.axis=cex.axis)
  axis(4, labels=paste0(yAxis,'Mb'), at=yPoints, cex.axis=cex.axis)
  #	axis(4, labels=FALSE, cex.axis=1.6)
  ########################################################
  ########################################################
}

chooseColor <- function(name){
  
  if (length(grep('ES',name)) == 1){col=4
  }else if(length(grep('NPC',name)) == 1){col=3
  }else if(length(grep('CN',name)) == 1){col=2
  }else if(length(grep('Hes5',name)) == 1){col=3
  }else if(length(grep('Dcx',name)) == 1){col=2
  }else{col=1}
  return(col)
}

plotChIP <- function(chipData, plotOrder, plotMar=c(3.5,3), vertical=FALSE, chipRes=10, plotAnn=F, plotStart=F, yLim=0,cex.axis,setNames=NA,chipColors=NA,figure_mode=FALSE,chipIndex=1){
#  print(setNames)
#  print(plotOrder)
  if(vertical == TRUE)
  {
    par(mar = c(plotMar[2], ifelse(figure_mode,0.1,1), plotMar[1], ifelse(figure_mode,0.1,1)))
    cNames <- paste0('v_', plotOrder)
    for (chip in cNames)
    {
      i=1
      if(is.na(setNames)){
        setName=gsub('v_chipseq.|v_chipseq_RPM.|v_hic.','',chip)
      } else {
        setName=setNames
      }
      if (length(grep(chip,colnames(chipData))) == 0){next}
      if(is.na(chipColors)){
        col <- chooseColor(chip)
      } else {
        col <- chipColors
      }
      currentData <- chipData[,grep(chip, colnames(chipData))]
      currentData[is.na(currentData)] <- 0
      if(length(grep('.ins_',chip)) > 0){currentData <- -currentData}
      chipIndex <- which(cNames == chip)
      
      # 			if(!is.na(yLim[chipIndex,1]))
      # 			{
      # 				currentLim <- c(yLim[chipIndex,2],yLim[chipIndex,1])
      # 				barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, main=gsub('v_chipseq.|v_chipseq_RPM.|v_hic.','',chip), horiz=T, space=0)
      # 			}else{
      # 				currentLim <- c(max(currentData),0)
      # 				barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, main=gsub('v_chipseq.|v_chipseq_RPM.|v_hic.','',chip), horiz=T, space=0)
      # 			}
      # 			axis(1, cex.axis=1.6, las=1)
      if(!is.na(yLim[chipIndex,1]))
      {
        currentLim <- c(yLim[chipIndex,2],yLim[chipIndex,1])
        barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, horiz=T, space=0,main=ifelse(figure_mode,'',setName),yaxs='i',xaxs='i')
      }else{
        currentLim <- c(max(currentData),0)
        barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, horiz=T, space=0,main=ifelse(figure_mode,'',setName),yaxs='i',xaxs='i')
      }
      axis(1, labels=c("",""),at=currentLim, cex.axis=cex.axis, las=1,padj=c(0,1))
      
      if(is.data.frame(plotAnn))
      {
        for(g in 1:nrow(plotAnn))
        {
          row <- plotAnn[g,]
          rect(0,(row$start-plotStart)/chipRes, 1000, (row$end-plotStart)/chipRes, col=rgb(0.66,0.66,0.66, alpha=0.6), border = 0)
        }
      }
    i=i+1  
    }
  }else{
    par(mar = c(ifelse(figure_mode,0.1,1), plotMar[1], ifelse(figure_mode,0.1,1), plotMar[2]))
  #  par(mar = c(ifelse(figure_mode,0.1,1), 0, ifelse(figure_mode,0.1,1), 0))
    cNames <- paste0('v_', plotOrder)
    i=1
    for (chip in cNames)
    {
      if (length(grep(chip,colnames(chipData))) == 0){next}
      if(is.na(setNames)){
        setName=gsub('v_chipseq.|v_chipseq_RPM.|v_hic.|v_scATAC.','',chip)
      } else {
        setName=setNames[i]
      }     
      if(is.na(chipColors)){
        col <- chooseColor(chip)
      } else {
        col <- chipColors
      }
      
      currentData <- chipData[,grep(chip, colnames(chipData))]
      currentData[is.na(currentData)] <- 0
      if(length(grep('.ins_',chip)) > 0){currentData <- -currentData}
      if(is.na(setNames)){
        setName=gsub('v_chipseq.|v_chipseq_RPM.|v_hic.|v_scATAC.','',chip)
      } else {
        setName=setNames
      }      
      if(!is.na(yLim[chipIndex,1]))
      {
        currentLim <- c(yLim[chipIndex,1],yLim[chipIndex,2])
        barplot(currentData, cex.main=cex.axis, ylim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, space=0,main=ifelse(figure_mode,'',setName),yaxs='i',xaxs='i')
       # mtext(setName, side = 4, outer = F,line=0,adj=0.5,las=2,cex=cex.axis*1.2)
        #mtext(setName, side = 4, outer = F,line=0.5,adj=0,las=2,cex=cex.axis*1.2)
        mtext(setName, side = 4, outer = F,line=0,at=0.75*currentLim[2],adj=1,las=2,cex=cex.axis)
      }else{
        currentLim <- c(max(currentData),0)
        barplot(currentData, cex.main=cex.axis, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, space=0,main=ifelse(figure_mode,'',setName),yaxs='i',xaxs='i')
       # mtext(setName, side = 4, outer = F,line=3,adj=0.5,las=2,cex=cex.axis*1.2)
        #mtext(setName, side = 4, outer = F,line=-0.5,adj=0,las=2,cex=cex.axis*1.2)
        mtext(setName, side = 4, outer = F,line=0,at=0.75*currentLim[1],adj=1,las=2,cex=cex.axis)
      }
      if(cex.axis!=0) {axis_labels=round(currentLim,ifelse(min(currentLim)<0,2,1))} else {axis_labels=FALSE}
      # if (axis_labels==c(-0.06,0.06)){
      #   axis_labels=c('-0.06','+0.06')
      # }
      axis(2, labels=axis_labels,at=currentLim, cex.axis=ifelse(cex.axis!=0,cex.axis,1), las=2,padj=c(0,1))
      
      if(is.data.frame(plotAnn))
      {
        for(g in 1:nrow(plotAnn))
        {
          row <- plotAnn[g,]
          rect((row$start-plotStart)/chipRes, 0, (row$end-plotStart)/chipRes, 1000, col=rgb(0.66,0.66,0.66, alpha=0.6), border = 0)
        }
      }
      i=i+1  
    }
  }
}

plotMeth <- function(methData, plotOrder, plotMar=c(3.5,3),vertical=FALSE, plotAnn=F, plotStart=F, yLim=c(0,100),cex.axis,setNames=NA,methColors=NA,figure_mode=FALSE,methIndex=1,currentIntrv){
  #  print(setNames)
  #  print(plotOrder)
  if(vertical == TRUE)
  {
    par(mar = c(plotMar[2], ifelse(figure_mode,0.5,1), plotMar[1], ifelse(figure_mode,0.5,1)))
    cNames <- paste0('v_', plotOrder)
    for (meth in cNames)
    {
      i=1
      if(is.na(setNames)){
        setName=gsub('v_methseq.|v_methseq_RPM.|v_hic.','',meth)
      } else {
        setName=setNames
      }
      if (length(grep(meth,colnames(methData))) == 0){next}
      if(is.na(methColors)){
        col <- chooseColor(meth)
      } else {
        col <- methColors
      }
      currentData <- methData[,grep(meth, colnames(methData))]
      currentData[is.na(currentData)] <- 0
      if(length(grep('.ins_',meth)) > 0){currentData <- -currentData}
      methIndex <- which(cNames == meth)
      
      # 			if(!is.na(yLim[methIndex,1]))
      # 			{
      # 				currentLim <- c(yLim[methIndex,2],yLim[methIndex,1])
      # 				barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, main=gsub('v_methseq.|v_methseq_RPM.|v_hic.','',meth), horiz=T, space=0)
      # 			}else{
      # 				currentLim <- c(max(currentData),0)
      # 				barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, main=gsub('v_methseq.|v_methseq_RPM.|v_hic.','',meth), horiz=T, space=0)
      # 			}
      # 			axis(1, cex.axis=1.6, las=1)
      if(!is.na(yLim[methIndex,1]))
      {
        currentLim <- c(yLim[methIndex,2],yLim[methIndex,1])
        barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, horiz=T, space=0,main=ifelse(figure_mode,'',setName),yaxs='i',xaxs='i')
      }else{
        currentLim <- c(max(currentData),0)
        barplot(currentData, cex.main=1.6, xlim=currentLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, horiz=T, space=0,main=ifelse(figure_mode,'',setName),yaxs='i',xaxs='i')
      }
      axis(1, labels=c("",""),at=currentLim, cex.axis=cex.axis, las=1,padj=c(0,1))
      
      if(is.data.frame(plotAnn))
      {
        for(g in 1:nrow(plotAnn))
        {
          row <- plotAnn[g,]
          rect(0,(row$start-plotStart)/methRes, 1000, (row$end-plotStart)/methRes, col=rgb(0.66,0.66,0.66, alpha=0.6), border = 0)
        }
      }
      i=i+1  
    }
  }else{
    par(mar = c(ifelse(figure_mode,0,1), plotMar[1], ifelse(figure_mode,0.1,1), plotMar[2]))
    #  par(mar = c(ifelse(figure_mode,0.1,1), 0, ifelse(figure_mode,0.1,1), 0))
    cNames <- plotOrder
    i=1
    for (meth in cNames)
    {
      if (length(grep(meth,colnames(methData))) == 0){next}
      if(is.na(setNames)){
        setName=gsub('methylation.','',meth)
      } else {
        setName=setNames[i]
      }     
      if(is.na(methColors)){
        col <- chooseColor(meth)
      } else {
        col <- methColors
      }
      
      currentData <- methData[,c(1:3,grep(meth, colnames(methData)))]
      currentData[is.na(currentData[,4]),4] <- 0
      if(is.na(setNames)){
        setName=gsub('methylation.','',meth)
      } else {
        setName=setNames
      }      
      
      currentLim <- yLim
      plot(x = currentData$start+1,y=currentData[,4],xlim=c(currentIntrv$start,currentIntrv$end),ylim=c(-5,100),pch=19,cex=1.2, cex.main=cex.axis, xlab='', ylab='', yaxt='n', xaxt='n', col=col,main=ifelse(figure_mode,'',setName),xaxs='i',bty="n")
      mtext(setName, side = 4, outer = F,line=0.5,adj=0,las=2,cex=cex.axis*1.2)

      axis(2, labels=yLim,at=yLim, cex.axis=ifelse(cex.axis!=0,cex.axis,1), las=2,padj=c(0,1))
      
      if(is.data.frame(plotAnn))
      {
        for(g in 1:nrow(plotAnn))
        {
          row <- plotAnn[g,]
          rect((row$start-plotStart)/methRes, 0, (row$end-plotStart)/methRes, 1000, col=rgb(0.66,0.66,0.66, alpha=0.6), border = 0)
        }
      }
      i=i+1  
    }
  }
}

scorePair <- function(gN1,gN2,set,radius,tssCoordinates){
  window <- radius*25
  g1 <- tssCoordinates[gN1,]
  g2 <- tssCoordinates[gN2,]
  
  chrSize <- gintervals.all()
  chr <- as.character(tssCoordinates[g1$geneName,'chrom'])
  chrLim <- chrSize[chrSize$chrom == chr, 'end']
  center <- gintervals.2d(tssCoordinates[gN1,'chrom'], tssCoordinates[gN1,'start'], tssCoordinates[gN1,'end'],
                          tssCoordinates[gN2,'chrom'], tssCoordinates[gN2,'start'], tssCoordinates[gN2,'end'])
  
  rSet <- c()
  for (i in 1:100)
  {
    xShift <- sample(c((-window):window),1);
    yShift <- sample(c((-window):window),1);
    rC <- gintervals.2d(tssCoordinates[gN1,'chrom'], tssCoordinates[gN1,'start']+xShift, tssCoordinates[gN1,'end']+xShift,
                        tssCoordinates[gN2,'chrom'], tssCoordinates[gN2,'start']+yShift, tssCoordinates[gN2,'end']+yShift)
    rSet <- rbind(rSet, rC)
  }
  
  allPoints <- rbind(center,rSet)
  allPoints$eD <- sqrt((center$start1-allPoints$start1)^2+(center$start2-allPoints$start2)^2)
  allPoints <- allPoints[allPoints$eD==0  | (allPoints$eD > radius & allPoints$eD < window*0.8), ]
  allPoints <- subset(allPoints, allPoints$end1 < chrLim & allPoints$end2 < chrLim)
  
  ### GSCREEN INTERVALS
  ### Number of reads
  vType <- 'area'
  tName <- paste0("v_",set,'_',vType)
  gvtrack.create(tName, set, vType)
  gvtrack.iterator.2d(tName, sshift1 = -radius, eshift1 = radius, sshift2 = -radius, eshift2 = radius)
  areaIntervals <- gextract(set, tName, allPoints, iterator=allPoints, colnames=c('score','area'))
  gvtrack.rm(tName)
  ### Max Score
  vType='max'
  tName <- paste0("v_",set,'_',vType)
  gvtrack.create(tName, set, vType)
  gvtrack.iterator.2d(tName, sshift1 = -radius, eshift1 = radius, sshift2 = -radius, eshift2 = radius)
  scoreIntervals <- gextract(set, tName, allPoints, iterator=allPoints, colnames=c('score','maxScore'))
  gvtrack.rm(tName)
  ### Average Score
  vType='avg'
  tName <- paste0("v_",set,'_',vType)
  gvtrack.create(tName, set, vType)
  gvtrack.iterator.2d(tName, sshift1 = -radius, eshift1 = radius, sshift2 = -radius, eshift2 = radius)
  avgIntervals <- gextract(set, tName, allPoints, iterator=allPoints, colnames=c('score','meanScore'))
  gvtrack.rm(tName)
  ########################################################
  outTable <- scoreIntervals
  outTable$meanScore <- avgIntervals$meanScore
  outTable$area <- areaIntervals$area
  return(outTable)
}

combineRNA <- function(inData, tracks, plotMar=c(3.5,3), vertical=FALSE, yLim=NA,cex.axis=cex.axis,setNames=NA,rnaColors=NA,figure_mode=FALSE){
  for (set in unlist(strsplit(tracks, "\\|")))
  { 
    i=1
    r1_For <- paste0('v_',set,'_rep1_For')
    r2_For <- paste0('v_',set,'_rep2_For')
    
    r1_Rev <- paste0('v_',set,'_rep1_Rev')
    r2_Rev <- paste0('v_',set,'_rep2_Rev')
    
    plusStrand <- paste0('v_',set,'_For')
    minusStrand <- paste0('v_',set,'_Rev')
    
    inData[,plusStrand] <- rowMeans(inData[,c(r1_For, r2_For)])
    inData[,minusStrand] <- rowMeans(inData[,c(r1_Rev, r2_Rev)])
    if (is.na(setNames)){
      setName <- gsub('_str1|_For','',plusStrand)
      setName <- gsub('v_rnaseq_RPM.','rnaseq_',setName)
    } else {
      setName <- setNames
    }
    if(is.na(rnaColors)){
      col <- chooseColor(plusStrand)
    } else {
      col <- rnaColors
    }    
    plotRNA(inData, plusStrand, minusStrand,plotMar,vertical,yLim,cex.axis=cex.axis,setName=setName,col=col,figure_mode=figure_mode)
    i=i+1
  }
}

combineRNA_total <- function(inData, tracks, plotMar=c(3.5,3), vertical=FALSE, yLim=NA,cex.axis=cex.axis,setNames=NA,rnaColors=NA,figure_mode=FALSE){
  for (set in unlist(strsplit(tracks, "\\|")))
  { 
    i=1
    r1_For <- paste0('v_',set,'_rep1')
    r2_For <- paste0('v_',set,'_rep2')
    #   r3_For <- paste0('v_',set,'_rep3')
    
    plusStrand <- paste0('v_',set,'_For')
    minusStrand <- paste0('v_',set,'_Rev')
    
    inData[,plusStrand] <- rowMeans(inData[,c(r1_For, r2_For)])
    inData[,minusStrand] <- 0
    if (is.na(setNames)){
      setName <- gsub('_str1|_For','',plusStrand)
      setName <- gsub('v_rnaseq_RPM.','rnaseq_',setName)
    } else {
      setName <- setNames
    }
    if(is.na(rnaColors)){
      col <- chooseColor(plusStrand)
    } else {
      col <- rnaColors
    }    
    plotRNA(inData, plusStrand, minusStrand,plotMar,vertical,yLim,cex.axis=cex.axis,setName=setName,col=col,figure_mode=figure_mode)
    i=i+1
  }
}

plotRNA <- function(inData, plusStrand, minusStrand, plotMar=c(3.5,3), vertical=FALSE, yLim,cex.axis,setName='',col,figure_mode=figure_mode){
  if(vertical==TRUE)
  {
    par(mar = c(plotMar[2], ifelse(figure_mode,0.1,1), plotMar[1],ifelse(figure_mode,0.1,1)))
    plusData <- inData[,plusStrand]
    minusData <- -inData[,minusStrand]
    if(length(yLim) == 1){
      xLim <- c(max(plusData),min(minusData))
    }else{xLim <- c(yLim[2],yLim[1])}
   # print(xLim)
 #   col <- chooseColor(plusStrand)
    barplot(plusData, horiz=T, xlim=xLim, space=0, cex.main=1.6, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col,main=setName,yaxs='i',xaxs='i')
    barplot(minusData, horiz=T, space=0, col=col, border=col, add=TRUE, xlab='', ylab='', yaxt='n', xaxt='n',yaxs='i',xaxs='i')
    segments(0,0,nrow(rnaData),0,lwd=1.75, col='darkgrey')		
    axis(1,labels=round(xLim,1),at=xLim, cex.axis=cex.axis, las=1)
  }else{
    par(mar = c(ifelse(figure_mode,0.1,1), plotMar[1], ifelse(figure_mode,0.1,1), plotMar[2]))
    plusData <- inData[,plusStrand]
    minusData <- -inData[,minusStrand]
    if (is.na(setName)){
      setName <- gsub('_str1|_For','',plusStrand)
      setName <- gsub('v_rnaseq_RPM.','rnaseq_',setName)
    }
    if(length(yLim) == 1){
      yLim <- c(min(minusData),max(plusData))
    }
 #   print(yLim)
#    col <- chooseColor(plusStrand)
    barplot(plusData, horiz=F, ylim=yLim, space=0, cex.main=cex.axis, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col,main=ifelse(figure_mode,'',setName),yaxs='i',xaxs='i')
    barplot(minusData, horiz=F, space=0, col=col, border=col, add=TRUE, xlab='', ylab='', yaxt='n', xaxt='n',yaxs='i',xaxs='i')
    segments(0,0,nrow(rnaData),0,lwd=1.75, col='darkgrey')
    if(cex.axis!=0){axis_labels=round(yLim,1)} else {axis_labels=FALSE}
    axis(2,labels=axis_labels,at=yLim, cex.axis=ifelse(cex.axis!=0,cex.axis,1), las=2,padj=c(0,1))
  }
}

meanNorm <- function(in_data)
{
  norm_dataset <- as.matrix(in_data)
  for(i in 1:ncol(norm_dataset))
  {
    current_set <- subset(norm_dataset, norm_dataset[,i] != 0)
    mean_value <- mean(current_set[,i])
    for (j in 1:nrow(norm_dataset))
    {
      if (norm_dataset[j,i] != 0)
      {
        ######## MEAN NORM
        norm_dataset[j,i] = norm_dataset[j,i]/mean_value
      }
    }
  }
  return(norm_dataset)
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

get_ucsc <- function(genome, interv, stacking="dense") {
  track="RefSeq Genes"
  table = "refGene"
  if (genome == "mm10" | genome == "hg19") {
    track = "NCBI RefSeq"
  }
  genetrack=Gviz::UcscTrack(genome=genome, track=track, table="refGene",
                            trackType="GeneRegionTrack",chromosome=as.character(interv$chrom), 
                            rstart="exonStarts", rends="exonEnds", gene="name2", symbol="name2",
                            transcript="name2", strand="strand", name="RefSeq Genes",
                            feature="name2", stacking=stacking,
                            showId=T, from=interv$start, to=interv$end)
  return(genetrack)
}


plotIdeogram <- function(genome,plotLim,chrom){
  ideo_track <- IdeogramTrack(genome=genome, chromosome=as.character(chrom),name=as.character(chrom))
 # options(ucscChromosomeNames=FALSE)
  plotTracks(ideo_track, from=plotLim[1], to=plotLim[2], panel.only=T, labelPos="below",add=F,fontsize=16,fontcolor='black',showId=F)
 # grid.text(chrom)
}


plot_genes <- function(genome,gene_track=NULL, interval_range, gene_stacking,gene_color,gene_size,fontsize=20,targetGene,collapseTranscripts,geneNames=TRUE)
{
  if(is.null(gene_track)){
    gene_track <- get_ucsc(genome, interval_range, stacking=gene_stacking)
  }
  displayPars(gene_track) <- list(size=gene_size,col='black',fill=gene_color,fontsize.group=fontsize,collapseTranscripts =collapseTranscripts ,fontcolor.group='black',showId=geneNames)
  if (nrow(targetGene)>0){
    httrack <- HighlightTrack(trackList = list(gene_track,GenomeAxisTrack(col="black",fontcolor="black",labelPos='below', lwd=2,ticksAt=seq(ceiling(interval_range[,2]/10000+1)*10000,floor(interval_range[,3]/10000-1)*10000,length.out=3))), start = targetGene$start,end = targetGene$end, chromosome = targetGene$chrom)
    displayPars(httrack) <- list(col='lightgrey',fill='lightgrey',alpha=0.4,inBackground=T,fontsize=fontsize)
  } else {httrack <- gene_track}
    plotTracks(httrack, from=interval_range$start, to=interval_range$end, panel.only=T,add=F,frame=F, labelPos="below",margin=0,innerMargin=0)
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

scCoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 1,
  cluster_cols=hue_pal(length(idents)),
  cells = NULL,
  idents = NULL,
  singleCells=TRUE,
  average=FALSE,
  plotMar=c(5,3)
) {
  if (!is.null(cells)){
    cells <- cells[cells %in% colnames(object)]
  } else {
    cells <- colnames(x = object)
  }
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  if (!is(object = region, class2 = 'GRanges')) {
    region <- makeGRangesFromDataFrame(region)
  }
  reads <- GetReadsInRegion(
    object = object,
    region = region,
    cells = cells,
    group.by = group.by,
    fragment.path = fragment.path,
    verbose = FALSE
  )
  cells.per.group <- CellsPerGroup(
    object = object,
    group.by = group.by
  )
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  
  if (average){
    reads.per.group <- AverageCounts(
      object = object,
      group.by = group.by,
      verbose = FALSE
    )
    coverages <- suppressWarnings(CalculateCoverages(
      reads = reads,
      cells.per.group = cells.per.group,
      reads.per.group = reads.per.group,
      window = window/2,
      verbose = FALSE
    ))
    if (downsample > 1) {
      warning('Requested downsampling <0%, retaining all positions')
      downsample <- 1
    }
    
    stepsize <- 1 / downsample
    total_range <- end.pos - start.pos
    steps <- ceiling(x = (total_range / stepsize))
    retain_positions <- seq(from = start.pos, to = end.pos, by = stepsize)
    downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
    ymax <- signif(x = max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)
    downsampled_coverage <- downsampled_coverage[!is.na(x = downsampled_coverage$coverage), ]
    
    p <- ggplot(data = downsampled_coverage, mapping = aes(x = position, y = coverage, fill = group)) +
      geom_bar(stat = 'identity', width=1) +
      facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
      xlab(label = paste0(chromosome, ' position (bp)')) +
      ylab(label = paste0('0 - ', as.character(x = ymax))) +
      ylim(c(0, ymax)) +
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        legend.position = 'none',
        strip.text.y = element_text(angle = 0)
      )
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
  }
  
  if(singleCells){
    reads <- makeGRangesFromDataFrame(reads,keep.extra.columns = T)
    sliding_region <- slidingWindows(region,window,window)[[1]]
    mat <- t(countInsertions2(sliding_region,reads,cells=cells,by = 'cell',window=window))
    mat@x[mat@x > 0] <- 1
    df.long <- melt(as.matrix(mat))
    df.long$group <- Idents(object)[match(df.long$Var1,names(Idents(object)))]
    cells_group <- paste0(names(cells.per.group))
    names(cells_group) <- names(cells.per.group)
    cells_group <- ''
    plotLabel <- as.list(cells_group)
    plot_labeller <- function(variable,value){
      return(plotLabel[value])
    }
    df_sub <- subset(df.long, value == 1)
    df_sub$color <- df_sub$group
    if (sum(!(df.long$Var1%in%df_sub$Var1))>0){
      df_neg <- as.data.frame(df.long[!(df.long$Var1%in%df_sub$Var1),]%>% group_by(Var1) %>% dplyr::slice(1))
      df_neg$color <- NA
      df_comb <- rbind(df_sub,df_neg)
    } else {df_comb <- df_sub}
    p1 <- ggplot(data = df_comb, mapping = aes(x = Var2, y = reorder(Var1), fill = color,size=4)) +
      geom_tile() + scale_x_continuous(expand = c(0,0),limits=c(start(region),end(region))) +
      scale_fill_manual(name='',values=as.character(cluster_cols),na.value = 'white')+
      facet_grid(rows=vars(group), scales="free_y",space='free_y',labeller =plot_labeller,margins=F) +
      xlab(label = paste0(chromosome, ' position (bp)')) +
      theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line=element_blank(),
        legend.position = 'none',
        strip.text.y = element_text(hjust=-10),
        strip.background = element_rect(
          color="black", fill=cluster_cols, size=1, linetype="solid"
        ),
        plot.margin = unit(c(0,0,0,0), "cm")
      )
    g <- ggplot_gtable(ggplot_build(p1))
    stripr <- which(grepl('strip-r', g$layout$name))
    fills <- cluster_cols
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    p1 <- g
  } 
  
  if (!is.null(x = annotation)) {
    gr <- GRanges(
      seqnames = chromosome,
      IRanges(start = start.pos, end = end.pos)
    )
    hits <- DataFrame(findOverlaps(gr,transcripts(annotation),type='any'))
    hits <- transcripts(annotation)[hits$subjectHits]
    hits <- AnnotationDbi::select(org.Mm.eg.db,keys=gsub('\\..*','',hits$tx_name),keytype = "REFSEQ",columns = "SYMBOL")
    genes <- suppressMessages(expr = ggbio::autoplot(object = annotation, which=gr,names.expr = hits$SYMBOL))
    gene.plot <- genes@ggplot +
      xlim(start.pos, end.pos) +
      xlab(label = paste0(chromosome, ' position (bp)')) +
      theme_classic()
  }
  if (singleCells) {
    if (average){
      if(!is.null(x = annotation)){
        plotlist <- list(p1,p,gene.plot)
        rel.heights=c(2,1,0.6)
      }
      plotlist <- list(p1,p)
      rel.heights=c(2,1)
    }else if(!is.null(x = annotation)){
      plotlist <- list(p1,gene.plot)
      rel.heights=c(1,0.6)
    } else {
      plotlist <- list(p1)
      rel.heights=c(1)     
    }
  } else if(!is.null(x = annotation)){
    plotlist <- list(p,gene.plot)
    rel.heights=c(1,0.6)
  } else {
    plotlist <- list(p)
    rel.heights=c(1,0.6)    
  }
  
  p <- suppressWarnings(cowplot::plot_grid(
    plotlist = plotlist,
    ncol = 1,
    axis = 'btlr',
    rel_heights = rel.heights,
    align = 'v',
    greedy = FALSE
  ))
 # p <- p + theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  message('Succesfully extracted scATAC')
  return(p)
}



















extractMarginal <- function(intervals,obsTracks,binsize=0,orig_intervals=F){
  v_obsTracks <- paste("v_", obsTracks, sep="")
  for (i in 1:length(v_obsTracks)) {
    gvtrack.create(v_obsTracks[i], obsTracks[i], 'weighted.sum')
    gvtrack.iterator.2d(v_obsTracks[i], sshift1=-binsize/2, eshift1=binsize/2, sshift2=-binsize/2, eshift2=binsize/2)
    
  }
  if(!orig_intervals){
    intervals <- intervals[intervals$chrom %in% gintervals.all()$chrom,]
    intervals <- intervals[(as.character(intervals$chrom)!='chrM'|as.character(intervals$chrom)!='chrY'),]
    intervals <- intervals[complete.cases(intervals),]
  #  if (binsize!=0) {intervals <- intervals.expand(intervals.centers(intervals), binsize/2)}
    intervals <- gintervals.canonic(intervals)
  }
  chroms <- gintervals.all()
  interv_marg <- adply(intervals,1,function(x){
    return(chroms[chroms$chrom%in%x$chrom,])
  })
  interv_marg <- gintervals.2d(intervals[,'chrom'],intervals[,'start'],intervals[,'end'],interv_marg[,'chrom'],interv_marg[,'start'],interv_marg[,'end'])
  marginal = gextract(v_obsTracks, intervals=interv_marg, iterator = interv_marg, band=c(-max(gintervals.all()$end), -1024))
  marginal2 = gextract(v_obsTracks, intervals=interv_marg, iterator = interv_marg, band=c(1024,max(gintervals.all()$end)))
  for (i in 7:(6+length(v_obsTracks))){
    marginal[is.na(marginal)] <- 0
    marginal2[is.na(marginal2)] <- 0
    intervals[,i-3] <- marginal[,i]+marginal2[,i]
    colnames(intervals)[i-3] <- colnames(marginal)[i]
  }
  return(intervals)
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

extractPairs <- function(intervals1,intervals2,obsTracks,binsize=0,orig_intervals=F,min_dist=10000,max_dist=2e6){
  v_obsTracks <- paste("v_", obsTracks, sep="")
  for (i in 1:length(v_obsTracks)) {
    gvtrack.create(v_obsTracks[i], obsTracks[i], 'weighted.sum')
    gvtrack.iterator.2d(v_obsTracks[i], sshift1=-binsize/2, eshift1=binsize/2, sshift2=-binsize/2, eshift2=binsize/2)
    
  }
  if(!orig_intervals){
    intervals1 <- intervals1[intervals1$chrom %in% gintervals.all()$chrom,]
    intervals1 <- intervals1[(as.character(intervals1$chrom)!='chrM'|as.character(intervals1$chrom)!='chrY'),]
    intervals1 <- intervals1[complete.cases(intervals1),]
    if(binsize!=0){intervals1 <- intervals.centers(intervals1)}
    intervals1 <- gintervals.canonic(intervals1)
    
    intervals2 <- intervals2[intervals2$chrom %in% gintervals.all()$chrom,]
    intervals2 <- intervals2[(as.character(intervals2$chrom)!='chrM'|as.character(intervals2$chrom)!='chrY'),]
    intervals2 <- intervals2[complete.cases(intervals2),]
    if(binsize!=0){intervals2 <- intervals.centers(intervals2)}
    intervals2 <- gintervals.canonic(intervals2)
  }
  
  grid1 <- construct.grid(intervals1,intervals2,min_dist,max_dist)
  grid2 <- construct.grid(intervals2,intervals1,min_dist,max_dist)
  grid <- unique(rbind(grid1,grid2))
  
  marginal = gextract(v_obsTracks, intervals=grid, iterator = grid, band=c(-max(gintervals.all()$end), -1024))
 # marginal2 = gextract(v_obsTracks, intervals=grid, iterator = grid, band=c(1024,max(gintervals.all()$end)))
  return(marginal)
  for (i in 1:length(v_obsTracks)) {gvtrack.rm(v_obsTracks[i])}
}

averageTssTts <- function(tracks=all_tracks,cells,domains=NULL,domain_size=c(1e5,2e6),file_f,path,bins,min_dist=1e4){
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
  file_f <- paste0(path,file_f)
  tad_list <- list()
  for (cell in cells){
    #			if (!gintervals.exists(domains[grep(cell,domains)])) { next}
    if (is.null(domains)) {
      domain <- paste0("hic.",cell,".ins250_k2_domains")
      domain <- gintervals.load(domain)
    } else {
      if (!grepl('bed',domains)) {
        #		domain <- gintervals.load(domains)
        domain <- get(load(domains))[,1:4]
        domain <- gintervals.canonic(domain)
        domain$cluster <- gintervals.neighbors(domain,get(load(domains))[,1:4],na.if.notfound = T)[,'strand']
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
    #### 
    x <- adply(domain,1,function(x){
      if(x$cluster==1){
        x$end2 <- x$end
        x$end <- x$start + 10000
        x$start2 <- x$start + round(x$len*0.5)
      } else {
        x$start2 <- x$end - 10000 ; x$end2 <- x$end 
        x$end <- x$end - round(x$len*0.5)
      }
      return(x)
    })
    ####
    interval2d <- gintervals.2d(x$chrom,x$start,x$end,x$chrom,x$start2,x$end2)
    df_obs <- gextract(paste0('v_',tracks[grep(cell,tracks)]),intervals=interval2d,iterator=interval2d,band=-c(domain_size[2],10000))
    df_exp <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)]),')/2'),intervals=interval2d,iterator=interval2d,band=-c(domain_size[2],10000))
    x <- domain
    interval2d_1 <- gintervals.2d(x[,1],x[,2],x[,3],x[,1],x[,2],x[,3])
    interval2d_1$start2 <- x$start+round(x$len/10)
    interval2d_1$end2 <- round((x$start+x$end)/2)
    interval2d_1$start1 <- x$start-x$len
    interval2d_1$end1 <- x$start-round(x$len/10)
    interval2d_1 <- gintervals.force_range(interval2d_1)
    df_obs1 <- gextract(paste0('v_',tracks[grep(cell,tracks)]),intervals=interval2d_1,iterator=interval2d_1,band=-c(domain_size[2]*2,10000))
    df_exp1 <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)]),')/2'),intervals=interval2d_1,iterator=interval2d_1,band=-c(domain_size[2]*2,10000))
    
    interval2d_2 <- gintervals.2d(x[,1],x[,2],x[,3],x[,1],x[,2],x[,3])
    interval2d_2$start1 <- round((x$start+x$end)/2)
    interval2d_2$end1 <- x$end-round(x$len/10)
    interval2d_2$start2 <- x$end+round(x$len/10)
    interval2d_2$end2 <- x$end+x$len
    interval2d_2 <- gintervals.force_range(interval2d_2)
    df_obs2 <- gextract(paste0('v_',tracks[grep(cell,tracks)]),intervals=interval2d_2,iterator=interval2d_2,band=-c(domain_size[2]*2,10000))
    df_exp2 <- gextract(paste0('(',paste0('v_',shuffled_tracks[grep(cell,tracks)]),')/2'),intervals=interval2d_2,iterator=interval2d_2,band=-c(domain_size[2]*2,10000))
    message(dim(df_obs1),' ',dim(df_obs2),' ',dim(df_exp1),' ',dim(df_exp2))
    df_obs1[,7:(ncol(df_obs1)-1)] <- df_obs1[,7:(ncol(df_obs1)-1)]+df_obs2[,7:(ncol(df_obs2)-1)]
    df_exp1[,7:(ncol(df_exp1)-1)] <- df_exp1[,7:(ncol(df_exp1)-1)]+df_exp2[,7:(ncol(df_exp2)-1)]
    df_out <- cbind(x$cluster,df_obs[,7:(ncol(df_obs)-1)],df_exp[,7:(ncol(df_exp)-1)],df_obs1[,7:(ncol(df_obs1)-1)],df_exp1[,7:(ncol(df_exp1)-1)])
    tad_list[[cell]]$res <- df_out
  }
  file_f <- paste0(file_f,'tssTts_quantification')
  save(tad_list,file=file_f)
}