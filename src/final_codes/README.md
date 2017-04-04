This folder contains the steps and scripts used for our analysis

Steps after [downloading the Data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/Acquiring%20GEO%20meta%20data.Rmd)

We downloaded both the Dasen normalized data and the QN normalized data

Processing:
- Both data sets were filtered using this [code](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering_Revised.Rmd)
- [Blood samlples were then removed from the data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Create.brain.only.Rmd)
- The beta values were then normaled for both the [QN/BMIQ]() **ADD** and [Dasen](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Dasen%20Normalization%20of%20Beta%20Data.Rmd) data sets that we downloaded 
- Both the BMIQ data sets and Dasen data sets went through [combat/PCA](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/PCA%20%26%20ComBat.Rmd) and [cell type corection]()**Add this**

Analysis
- sample to sample correlation plots and heatmaps were created from various codes that can be found with explanations in the working_code [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md). This is thefinal code for the [BMIQ data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(BMIQ).Rmd) and [Dasen data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd)
- DMR analysis was conduced using various codes
  -- **WHat did each of these do** [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially%20Methylated%20Probe%20Analysis%20-%20Batch%20Corrected%20Only%20(DMR%20Setup)%20Final.Rmd)
  -- write [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially%20Methylated%20Probe%20Analysis%20-%20Cell-Type%20Corrected%20(DMR%20Setup)%20Final.Rmd)
  --more [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially_Methylated_Probe_Analysis_-_Batch_Corrected_Only__DMR_Setup__Final.md)
  -- here  [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Differentially_Methylated_Probe_Analysis_-_Cell-Type_Corrected__DMR_Setup__Final.md)
- Boxplots for the most differentially methylated probes were generated [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd). Pleminary versions of the code can be found in the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md)
Wilcoxon test was performed to compare probe values between cerebellum and cortex. The most significantly differentially methylated probe that was found only in the list from cell type corrected, only in non cell-type corrected, and in the overlap was used. The code for the wilcoxon test can be found [here] (WilcoxTestProbes_analysis.R)




