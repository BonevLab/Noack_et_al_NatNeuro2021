#### If Please cite as Noack et al., Nature Neuroscience 2022 """
#### Created by Boyan Bonev ####


main_f <- '/home/hpc/bonev/projects/SC/analysis/'        #Please set this to the main path where you wish to perform the analysis
library(pals)
### Global Directories ###
rna_sample_dir <- '/home/hpc/bonev/projects/SC/analysis/scRNA/'             # the files from the 10x scRNA run should be located here
atac_sample_dir <- '/home/hpc/bonev/projects/SC/analysis/scATAC/'           # the files from the 10x scATAC run should be located here
sample_names <- c('E14_rep1','E14_rep2')

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

### Colors
archr_colors <- ArchR::ArchRPalettes
sample_colors <-  function(x){return(brewer.pal(n = x, name = "Reds"))}
rep_colors <-  c('#F0F0F0','#5A5A5A')
#cluster_colors <- function(x){glasbey(x)}
RNA_cluster_colors <- c("#F51111","#FFCAC7","#16A810","#D6FFD1","#EBEB17","#4861F0","#F20DFA","#16EAF5","#7400D9","#9E560E","#FFA200")
ATAC_cluster_colors <- c("#F51111","#16A810","#EBEB17","#4861F0","#16EAF5","#7400D9","#F78F4ADF")
gene_colors <- function(x){colorRampPalette(rev(colorpalette('reds',10)))(x)}
#cluster_colors <- function(x){cols25(x)}
heatmap_colors <- colorpalette('matlablike',10)
col.pbreaks <<- c(20,35,50,65,75,85,95)        #Original
col.pos <<- palette.breaks(100 , c("lavenderblush2","#f8bfda","lightcoral","red","orange","yellow"), col.pbreaks)
col.nbreaks <<- c(20,35,50,65,75,85,95)
col.neg <<- rev(palette.breaks(100 , c("powderblue", "cornflowerblue", "blue","blueviolet", "#8A2BE2", "#4B0082"), col.nbreaks ))
hic.scores <<- c(col.neg, "lightgrey", col.pos)


### Marker Genes
diff_markers <- c('Hes1','Neurog2','Eomes','Dcx','Mapt')
layer_markers <- c('Tle4','Bcl11b','Rorb','Satb2','Cux1','Cux2')
other_markers <- c('C1qc','Vtn','Reln','Gad2','Gad1','Vip','Sst')
temp_markers <- c('Nusap1','Prc1','Hmga2','Slc1a3','Apoe','Mfge8','Fabp7')
regional_markers <- c('Foxg1','Emx2','Gsx2','Pax6','Sp8','Nr2f1')
mitotic_markers <- c('Cdk1','Ube2c','Top2a','Hist1h4e','Mki67','Pcna')
tf_markers <- c('Id4','Id1','Tcf4','Hey1','Hey2','Zeb1','Meis2','Smarca2','Tcf12','Ldb1','Ldb2','Lhx2','Lmo3','Lmo4','Isl1','Nr2f1','Smarca5','Chd4')

### Genomic databases ###

#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
