
This folder contains the steps and scripts used for processing the metadata and beta values 

Steps after [downloading the Data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/Acquiring%20GEO%20meta%20data.Rmd)



Step 1
- filtering probes - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/scr/Probe%20Filtering.Rmd)
- Remove blood samples from data - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/scr/Create.brain.only.Rmd)
- Normalizing beta values (BMIQ) - [here]NEED this one still
- Normalizing beta values (Dasen) - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/Dasen%20Normalization%20of%20Beta%20Data.Rmd)


Step 2
- Cell type Prediction, PCA, and Combat - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/PCA%20%26%20ComBat.Rmd) 


Step 4
- Plot and determine if the differences between corrected and uncorrected matter 

[Generating data with samples sorted by brain region](analysis_script_v2.Rmd) - this code was incorperated into the plotting codes.

[broad expression and correlation heatmaps, methylation levels of probes generated for corrected and uncorrected data with DASEN normalized data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/Correction_matters.md) 

[broad expression and correlation heatmaps, methylation levels of probes generated for corrected and uncorrected data with BMIQ normalized data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/Correction_matters.2.md)

[Wilcoxon Test script to calculate p-value for significantly differentially methylated probes between cerebellum and cortex](WilcoxTestProbes.R)
