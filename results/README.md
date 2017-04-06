
Results
========
In this folder we will discuss the input and outputs of the code as well as interpretations. The breakdown of our codes can be found in the [src folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src/final_codes)

Overarching aims:
1. To determine whether cell type correction between brain regions is necessary in the analysis of Illumina HumanMethylation450 BeadChip array data
2. To investigate probes differentially methylated between cerebellum and cortex regions

Processing
-------------

Overview of processing steps:

![](/Images/Pipeline_of_Methods.png)
*Probe filtering:*

**Aim**   
To remove probes from the data which would confound our analysis, including probes with SNPs, probes with NAs, and cross-hybridizing probes. Data was also filtered to contain only brain samples (originally had brain + blood).

**Input:**

--Dasen data = file: GSE59685.RData

--meta data = file: Meta_matched_GSE59685.RData

--QN data = file: GSE43414_betaqn_geo_all_cohorts.csv

--meta data = file: Brain_meta_matched_GSE43414.RData

**Output:**
  
--Dasen filtered data = file: GSE59685_filtered.RData

--Dasen filtered meta data = file: Meta_matched_brain.txt

--QN filtered data = file: GSE43414_filtered.RData

--QN filtered meta data = file: Meta_uncor.RData
  
**Conclusion**
Prior to filtering, there were 485577 probes. Filtering removed 65 SNP probes, 20869 probes that had a SNP within, 11475 probes hybridzing to sex chromosomes, 10673 cross hybrizing X and Y probes, 27415 cross hybridizing autosomal probes, and 176 NA probes. A total of 414904 probes remained.

*Beta-mixture quantile (BMIQ) normalization:*

**Aim**   
To correct for intra-array differences between probe type I and II chemistries that cause differing beta value distributions. This was performed on the QN-normalized data we had decided to move forward with, in order to make direct comparisons between any CpG site across any two brain regions.

**Input:**

--QN data = file: GSE43414_filtered.RData

--meta data = file: Meta_uncor.RData

**Output:**

--BMIQ normalized data = file: GSE43414_BMIQ.RData

**Conclusion**
![](/Images/raw_normalized_ggplot.png)   
Applying BMIQ normalization made the type I and II probe distributions highly similar, and now much more appropriate for direct comparisons and further analysis.

*ComBat batch correction, cell type prediction, and PCA:*

**Aim**   
Multiple BeadChips were used to generate this data; we aimed to correct for variation due to the individual chip and position on the chip. Since we would be analyzing the effect of cell type, we also used to CETS package to predict proportions of neurons and glia based on the non-corrected and the corrected beta values. These neuronal propotions were attached to the meta data and used as continous variables in PCA analysis.

**Input:**

--BMIQ normalized data = file: GSE43414_BMIQ.RData

--meta data = file: Meta_uncor.RData

**Output:**

--ComBat corrected data = file: GSE43414_batch_cor.RData

--meta data = file: Meta_batch_cor.RData


**Conclusion**
![](/Images/PCA_uncor_zoom.png)   
Before batch correction, there was a large effect of chip and tissue; it appeared the samples had been batched by tissue type.   

![](/Images/PCA_batch_zoom.png)   
After batch correction, the chip and tissue effects plus much of the other variation was removed, leaving only a Neuron (cell type) effect.

*CETS cell type correction and PCA:*

**Aim**   
Linear regression and the predicted neuronal proportion was used to fit our beta values to reference profiles of neuronal and glial methylation, removing variation due to cell type.

**Input:**

--ComBat corrected data = file: GSE43414_batch_cor.RData

--meta data = file: Meta_batch_cor.RData

**Output:**

--CETS corrected data = file: GSE43414_cell_cor.RData


**Conclusion**
![](/Images/PCA_cell_zoom.png)   
After cell type correction, Neuron comprised a smaller propotion of total variance, and we unmasked some effects of braak stage, age, and AD disease status. 


Analysis
--------

*Sample to sample correlation:*

**Aim**   
To determine whether the samples correlate well with each other. If processing normalized the data we should see that samples will correlate with samples within the same cell type and outside of the same cell type.

**Inputs**

--Dasen filtered data = file: GSE59685_batch_cor.RData

--Dasen filtered + cell corrected data = file:GSE59685_cell_cor.RData

--meta data = file: Meta_brain_cell_batch.RData

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file: Meta_batch_cor.RData

**Outputs**

Dasen

![](/results/Heatmaps/Non-cell_type_corrected_sample_cor_dasen.png)

![](/results/Heatmaps/Cell_type_corrected_sample_cor_dasen.png)

BMIQ

![](/results/Heatmaps/Non_cell_type_corrected_sample_cor_BMIQ.png)

![](/results/Heatmaps/Cell_type_corrected_sample_cor_BMIQ.png)

![](/results/Heatmaps/Batch_BMIQ_heatmap.png)

![](/results/Heatmaps/Cell_type_corrected_BMIQ.png)

**Conclusions**   
QN+BMIQ normalization increases sample to sample correlations compared to dasen normalization. Cell type correction increases inter tissue sample correlations in the dasen model, in line with the purpose of the normalization method to increase the comparability between tissues. This is unable to be seen in the BMIQ data as the samples are already highly correlated.


*DMP and DMR analysis*

**Aim**   
to identify and compare differentially methylated probes (at an FDR 0.05 and absolute delta beta >=0.05) between tissue type, and differentially methylated regions in both the batch-corrected only data, and the both batch-corrected and cell-type corrected data.

**Inputs**

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

**Outputs**

- Top Probes (The full txt files are too large for upload, so a subset is presented here):

List of top probes significant in both datasets:

chr [1:11095] "cg04197371" "cg18709737" "cg25960893" "cg00267142" ...

List of top probes significant in cell-type and batch corrected only:

chr [1:1768] "cg23112745" "cg12645236" "cg14165909" "cg07133729" ...

List of top probes significant in batch-corrected only:

chr [1:7419] "cg06826710" "cg25057705" "cg01219907" "cg04462567" ...

![Batch-corrected DMPs](/results/Volcano plots/volcano_DB_batch.png)

![Cell-corrected DMPs](/results/Volcano plots/volcano_DB_cell.png)

![Number of significant probes found in each type of analysis conducted](/results/Venn diagram/DMPs_venn.png)

**Conclusions**   
Correcting for cell type when examining cortex and cerebellum methylation differences may alter the number of significant probes and regions identified.

*Single probe comparisons and Wilcoxon tests*

**Aim**   
To show statistical differences between tissues regions in probes found only in the cell type corrected analysis, only in the non cell type corrected analysis and in both.

**Inputs**

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

**Outputs**

- images

![Most significant DMP in cell-type corrected data](/results/boxplots/Rplotcell.png)

![Most significant DMP in non-cell-type corrected data](/results/boxplots/Rplotbatch.png)

![Most significant DMP found in both data sets](/results/boxplots/Rplotbatch.png)

- statistics: wilcoxon test was performed for the 3 probes (shown in boxplots above) on the cell-type corrected dataset and non cell-type corrected data

**Statistics for cell-type corrected data**

probe IDs     | p value       | adjusted p value  |
------------- | ------------- | ----------------- |   
cg04197371    | 3.161739e-13  | 4.742609e-13      |
cg23112745    | 7.113933e-11  | 7.113933e-11      |  
cg06826710    | 4.018419e-14  | 1.205526e-13      |


**Statistics for non cell-type corrected data**

probe IDs     | p value       | adjusted p value  |
------------- | ------------- | ----------------- |   
cg04197371    | 9.664758e-15  | 1.449714e-14      |
cg23112745    | 9.327963e-06  | 9.327963e-06      |  
cg06826710    | 3.652858e-15  | 1.095857e-14      |


**Conclusions**   
correcting for cell type can affect the relationship between cortex and cerebellum methylation differences. 

Overall conclusions:
======================
• BMIQ normalization is more suitable for cross-sample comparisons
• Technical and biological variability to account for are batch effects and cell type variation across
brain regions to allow for comparison between tissues
• Differential methylation analysis reveals DMPs that are attributed to cell type variation and those
that are independent of cell type variation
• Failing to correct for cell type may result in identification of false positive DMPs/DMRs, with
implications for biological conclusions concerning clinically relevant loci
