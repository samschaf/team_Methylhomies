This folder contains the steps and scripts used for our analysis

## Data acquisition

A discription of the data and how it was obtained can be found [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/README.md)

Steps after [downloading the Data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/Acquiring%20GEO%20meta%20data.Rmd)

We downloaded both the Dasen normalized data and the QN normalized data

## Processing
- Both data sets were filtered using this [code](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering_Revised.Rmd)
- [Blood samples were then removed from the data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Create.brain.only.Rmd)
- The beta values were then normaled for both the [QN/BMIQ]() **ADD** and [Dasen](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Dasen%20Normalization%20of%20Beta%20Data.Rmd) data sets that we downloaded 
- Both the BMIQ data sets and Dasen data sets went through [combat/PCA](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/PCA%20%26%20ComBat.Rmd) and [cell type correction]()**Add this**

A overview of the processing steps are visualized in this pipeline:

![](/Images/Pipeline_of_Methods.png)

## Analysis
- Sample to sample correlation plots and heatmaps were created from various codes that can be found with explanations in the working_code [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md). This is the final code for the [BMIQ data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(BMIQ).Rmd) and [Dasen data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd)
- Differentially methylated probe (DMP) and differentially methylated region (DMR) analysis was conduced using various codes
  -- DMP analysis was conducted on the batch-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially%20Methylated%20Probe%20Analysis%20-%20Batch%20Corrected%20Only%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially_Methylated_Probe_Analysis_-_Batch_Corrected_Only__DMR_Setup__Final.md). The same analysis was also conducted on the batch and cell-type corrected dataset [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially%20Methylated%20Probe%20Analysis%20-%20Cell-Type%20Corrected%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially_Methylated_Probe_Analysis_-_Cell-Type_Corrected__DMR_Setup__Final.md), where comparisons of significant probes (at an FDR <= 0.05 and with an absolute value delta beta >= 0.05) found in these two datasets were also made.
- Boxplots for the most differentially methylated probes were generated [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd). Preliminary versions of the code can be found in the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md)
- Wilcoxon test was performed to compare probe values between cerebellum and cortex. The most significantly differentially methylated probe that was found only in the list from cell type corrected, only in non cell-type corrected, and in the overlap was used. The code for the wilcoxon test can be found [here](WilcoxTestProbes_analysis.R)




