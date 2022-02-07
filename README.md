# Multimodal profiling of the transcriptional regulatory landscape of developing mouse cortex identifies Neurog2 as a key epigenome remodeler 

*<h3>Noack F., Vangelisti S., Raffl G., Carido M., Diwakar J., Chong F., Bonev B. Nature Neuroscience 2022</h3>*

*<h3>[Please find the publication here](https://www.nature.com/articles/s41593-021-01002-4)</h3>*


![](Gif1.gif)

# Description of Scripts and Analysis 

## scRNA-seq Analyses 

**scRNA1-QC** - Quality control and standard preprocessing of scRNA data using Seurat v3.1.5

**scRNA2-UMAP** - Louvain clustering + UMAP visualization of scRNA data following transformation.

**scRNA3-DE** - Assignment of cluster identities based DE genes (using MAST - Finak et al., 2015). 

**scRNA4-monocle3** - Pseudotime analysis using Monocle3. Analyses include model fitting to identify gene expression changes as a function of pseudotime.


## scATAC-seq Analyses 

**scATAC1-QC** - Quality control of scATAC data based on number of fragments and TSS enrichment. 

**scATAC2-aggregateBin** - Final QC and code for generating count matrices based on genomic bins (for initial clustering), gene bodies and promoters. 

**scATAC3-normalization** - Code for inital clustering based on fragments from fixed-size genome wide bins.

**scATAC4-peakNormalization** - Final peak calling based on intial clusters to generate high quality peak set, used for final clustering and visualization.

**scATAC5-chromVar_motifs** - Computing motif accessibility deviations using chromVAR (Schep et al., 2017) implemented in [Signac](https://github.com/timoast/signac).

**scATAC6-Compute_Gene_Scores** - Computing gene activity scores using Cicero, used for subsequent integration analyses (seen below).

**scATAC7-cluster_unique_peaks** - Identification of cluster specific peaks.


## Integration Analyses of scRNA & scATAC

*Credit for many of the functions in this part of the analysis goes to [Satpathy\*, Granja\* et al. 2019](https://github.com/GreenleafLab/MPAL-Single-Cell-2019) and [Granja et al. 2019](https://github.com/GreenleafLab/10x-scATAC-2019)*

**scRNA_scATAC_Integration_01_Align_scATAC_scRNA** - Integration of scRNA and scATAC data using Seurat CCA and identification of nearest neighbors (kNN).

**scRNA_scATAC_Integration_02_Create_Aggregate_scATAC_scRNA** - Aggregate scRNA and scATAC data using nearest neighbor information.

**scRNA_scATAC_Integration_03_Compute_Peak_to_Gene_links** - Identification of enhancer-gene pairs by correlating each pair of distal scATAC peak and gene promoter.

**scRNA_scATAC_Integration_04_P2G_analysis** - Further characterization of identified enhancer-gene pairs.

**scRNA_scATAC_Integration_05_P2G_monocle** - Pseudotime analysis on integrated scRNA-scATAC object using Monocle3. Analyses also include model fitting to identify changes of accessibility and motif deviations as a function of pseudotime.

**scRNA_scATAC_Integration_06_chromVar** - Motif analysis for integrated object using chromVar.

## List of Figures 

**Figure1** - scRNA-seq analysis of mouse E14.5 cortical development. Associated Extended Data Fig. 1-2. 

**Figure2** - scATAC-seq identifies dynamic TF motis and variable distal regulatory elements. Associated Extended Data Fig. 3. 

**Figure3** - Lineage dynamics of enhancer-gene pairs and transcription factor motifs. Associated Extended Data Fig. 4-5.

**Figure4** - _In vivo_ immunoMPRA validates cell type specific activity of identified CREs and their regulation by transcription factors. Associated Extended Data Fig. 6.

**Figure5** - ImmunoMethyl-HiC identifies DNA methylation-independent global changes in 3D genome architecture during cortical development. Associated Extended Data Fig. 7.

**Figure6** - Dynamic enhancer-promoter loops and DNA methylation levels at regulatory regions. Associated Extended Data Fig. 8. 

**Figure7** - Transcription factors are associated with changes in both chromatin looping and DNA methylation levels. Associated Extended Data Fig. 9. 

**Figure8** - Neurog2 is sufficient to induce multilayered epigenome changes in vivo. Associated Extended Data Fig. 10. 

**Extended Data Fig. 1** - Quality control and further validation of the scRNA-seq data 

**Extended Data Fig. 2** - Pseudo-temporal analysis and comparison to human data.

**Extended Data Fig. 3** - Quality control and validation of the scATAC-seq data.

**Extended Data Fig. 4** - scATAC-scRNA integration metrics and identification of enhancer-gene pairs. 

**Extended Data Fig. 5** - Properties of the identified enhancer-gene links.

**Extended Data Fig. 6** - In vivo (immune-)MPRA validates cell type specific activity of identified CREs and their regulation by transription factors.

**Extended Data Fig. 7** - Global reorganization of the 3D genome in vivo during differentiation.

**Extended Data Fig. 8** - Non-correlated enhancer-gene pairs are not associated with dynamic chromatin looping or changes at DNA methylation levels. 

**Extended Data Fig. 9** - Transcription factors associated with cell-type-specific looping based on positively correlated enhancer-gene pairs. 

**Extended Data Fig. 10** - No global reorganization of chromatin accessibility, DNA methylation and 3D genome architecture upon Neurog2 overexpression. 


# Raw Data Download 

### GEO Accession : GSE155677 [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155677]
