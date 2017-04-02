
This folder contains the steps and scripts used for processing the metadata and beta values 

Steps after [downloading the Data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/processed_data/Acquiring%20GEO%20meta%20data.Rmd)



Step 1
- Normalizing beta values - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/processed_data/Normalization%20of%20Beta%20Data.Rmd)
- filtering probes - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/processed_data/Probe%20Filtering.Rmd)
- Remove blood samples from data - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/processed_data/Create.brain.only.Rmd)

Step 2
- Cell type Prediction, PCA, and Combat - [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/processed_data/PCA.Rmd) 
- we used a [code](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/PCA%20%26%20ComBat.Rmd) previously made as a reference 

Step 3 
- Rearrange data so all data sets have same order - here() #add this

Step 4
- Plot and determine if the differences between corrected and uncorrected matter 


[Generating data with samples sorted by brain region](analysis_script_v2.Rmd) 


[broad expression and correlation heatmaps, methylation levels of probes generated for corrected and uncorrected data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/processed_data/Correction_matters.md) 

**NEW** [Generating data with samples sorted by brain region on the new corrected data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/processed_data/Correction_matters.2.md)

[Wilcoxon Test script to calculate p-value for significantly differentially methylated probes between cerebellum and cortex](WilcoxTestProbes.R)
