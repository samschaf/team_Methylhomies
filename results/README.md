
Results
========
in this folder we will discuss the imput and outputs of the code as well as interpretations. the breakdown of our codes can be found in the [src folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src/final_codes)

Overarcing aims:
1. To determine whether cell type correction between brain regions is necessary in the analysis of Illumina HumanMethylation450 BeadChip array data
2. To investigate probes differentially methylated between cerebellum and cortex regions

Processing
-------------

Overview of processing steps:

![](/Images/Pipeline_of_Methods.png)

**Aim**

**Imput: **

Dasen data = file: 

QN data = file:

meta data = file

**Output:**
  - files
  
--Dasen filtered data = file: GSE59685_batch_cor.RData

--Dasen filtered + cell corrected data = file:GSE59685_cell_cor.RData

--meta data = file: Meta_brain_cell_batch.RData

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

  - images
  
**conclusion**

Analysis
--------

*Sample to sample correlation:*

**Aim**
to determine whether the samples correlate well with eachother. If processing normalized the data we should see that samples will correlate with samples within the same cell type and outside of the same cell type.

**inputs**

--Dasen filtered data = file: GSE59685_batch_cor.RData

--Dasen filtered + cell corrected data = file:GSE59685_cell_cor.RData

--meta data = file: Meta_brain_cell_batch.RData

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

**outputs**
- images

Dasen

![](/results/Heatmaps/Non-cell_type_corrected_sample_cor_dasen.png)

![](/results/Heatmaps/Cell_type_corrected_sample_cor_dasen.png)

BMIQ

![](/results/Heatmaps/Non_cell_type_corrected_sample_cor_BMIQ.png)

![](/results/Heatmaps/Cell_type_corrected_sample_cor_BMIQ.png)

![](/results/Heatmaps/Batch_BMIQ_heatmap.png)

![](/results/Heatmaps/Cell_type_corrected_BMIQ.png)

**conclusions**
QN+BMIQ normalization increases sample to sample correlations compared to dasen normalization. Cell type correction increases inter tissue sample correlations in the dasen model, in line with the purpose of the normalization method to increase the comparability between tissues. This is unable to be seen in the BMIQ data as the samples are already highly correlated.


*DMP and DMR analysis*

**Aim**
to identify and compare differentially methylated probes (at an FDR 0.05 and absolute delta beta >=0.05) between tissue type, and differentially methylated regions in both the batch-corrected only data, and the both batch-corrected and cell-type corrected data.

**inputs**

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

**outputs**
Top Probes:

List of top probes significant in both datasets
chr [1:11095] "cg04197371" "cg18709737" "cg25960893" "cg00267142" ...

List of top probes significant in cell-type and batch corrected only 
chr [1:1768] "cg23112745" "cg12645236" "cg14165909" "cg07133729" ...

List of top probes significant in batch-corrected only 
chr [1:7419] "cg06826710" "cg25057705" "cg01219907" "cg04462567" ...

![Number of significant probes found in each type of analysis conducted](/results/Venndiagram/DMPs_venn.png)

**conclusions**
correcting for cell type when examining cortex and cerebellum methylation differences may alter the number of significant probes and regions identified.

*single probe comparisons and Wilcoxon tests*

**Aim**
to show statistical differences between tissues regions in probes found only in the cell type corrected analysis, only in the non cell type corrected analysis and in both.

**inputs**

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

**outputs**

- images

![Most significant DMP in cell-type corrected data](/results/boxplots/Rplotcell.png)

![Most significant DMP in non-cell-type corrected data](/results/boxplots/Rplotbatch.png)

![Most significant DMP found in both data sets](/results/boxplots/Rplotbatch.png)

- statistics: wilcoxon test was performed for the 3 probes (shown in boxplots above) on the cell-type corrected dataset and non cell-type corrected data

**statistics for cell-type corrected data**

probe IDs     | p value       | adjusted p value  |
------------- | ------------- | ----------------- |   
cg04197371    | 3.161739e-13  | 4.742609e-13      |
cg23112745    | 7.113933e-11  | 7.113933e-11      |  
cg06826710    | 4.018419e-14  | 1.205526e-13      |


**statistics for non cell-type corrected data**

probe IDs     | p value       | adjusted p value  |
------------- | ------------- | ----------------- |   
cg04197371    | 9.664758e-15  | 1.449714e-14      |
cg23112745    | 9.327963e-06  | 9.327963e-06      |  
cg06826710    | 3.652858e-15  | 1.095857e-14      |


**conclusions**
correcting for cell type can effect the relationship between cortex and cerebellum methylation differences. 

Overall conclusions:
======================
• BMIQ normalization is more suitable for cross-sample comparisons
• Technical and biological variability to account for are batch effects and cell type variation across
brain regions to allow for comparison between tissues
• Differential methylation analysis reveals DMPs that are attributed to cell type variation and those
that are independent of cell type variation
• Failing to correct for cell type may result in identification of false positive DMPs/DMRs, with
implications for biological conclusions concerning clinically relevant loci
