This folder contains the steps and scripts used for our analysis

## Data acquisition

A description of the data and how it was obtained can be found [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/README.md)

We downloaded both the Dasen normalized data and the QN normalized data. Meta data was obtained from GEO using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Acquiring%20GEO%20meta%20data.Rmd) code

## Processing
- Dasen data was filtered using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering%20Dasen.Rmd) code and QN data was filtered using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering%20QN.Rmd) code   
- [Blood samples were then removed from the data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Create.brain.only.Rmd)
- The QN-normalized beta values were then [BMIQ normalized](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/BMIQ_final.md)
- The BMIQ data set went through [cell type prediction](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/Cell%20Type%20Prediction.Rmd) and [combat/PCA and cell type correction](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/PCA%20%26%20ComBat.md)

A overview of the processing steps are visualized in this pipeline:

![](/Images/Pipeline_of_Methods.png)

## Analysis
- Sample to sample correlation plots and heatmaps were created from various codes that can be found with explanations in the working_code [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md). This is the final code for the [BMIQ data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(BMIQ).Rmd) and [Dasen data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd)
- Differentially methylated probe (DMP) and differentially methylated region (DMR) analysis was conduced using various codes
  -- Linear modelling for delta betas and generating volcano plots was performed on batch-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/DMR%20batch%20cor.md) and cell-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/DMR%20cell%20cor.md)
  -- DMP analysis was conducted on the batch-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially%20Methylated%20Probe%20Analysis%20-%20Batch%20Corrected%20Only%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially_Methylated_Probe_Analysis_-_Batch_Corrected_Only__DMR_Setup__Final.md). The same analysis was also conducted on the batch and cell-type corrected dataset [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially%20Methylated%20Probe%20Analysis%20-%20Cell-Type%20Corrected%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially_Methylated_Probe_Analysis_-_Cell-Type_Corrected__DMR_Setup__Final.md), where comparisons of significant probes (at an FDR <= 0.05 and with an absolute value delta beta >= 0.05) found in these two datasets were also made. 
- Boxplots for the most differentially methylated probes were generated [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd). Preliminary versions of the code can be found in the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md)
- Wilcoxon test was performed to compare probe values between cerebellum and cortex for the most significantly differentially methylated probe. The most significantly differentially methylated probe that was found only in the list from cell type corrected, only in non cell-type corrected, and in the overlap was chosen for visualization with the boxplots. The code for the wilcoxon test can be found [here](WilcoxTestProbes_analysis.R)




