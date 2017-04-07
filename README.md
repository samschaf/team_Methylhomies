Normalizing for Cell Type: Does it matter?
============================================

Team
------------

**Table of group members with their background, degree, affiliations and job assignments.**

| Name | Degree/Affiliation | Speacilizations |Job Assignments | 
| ------------- | ------------- | ------------- | ------------- |
| Samantha Schaffner | Medical Genetics  | Neuroepigenetics of human development and Parkinson’s Disease | DMR analysis and help with data processing and normalization |
| Cassia Warren | Interdisciplinary Oncology   | Biomarker development for Vitamin C treatment for Pancreatic cancer | Maintain/organize Github and Cell type correction analysis|
| Hilary Brewis  | Medical Genetics  | Histone variants and chromatin structure | Normalization and data processing |
| Randip Gill  | Educational Psychology  | Social epigenetics | DMR analysis |
| Lisa Wei | Bioinformatics  | Cancer genomics, sequence analysis | Help maintain/organize Github and Overall methylation level differences between brain subregions |

Our project uses data that was generated for the paper [Interindividual methylomic variation across blood, cortex, and cerebellum: implications for epigenetic studies of neurological and neuropsychiatric phenotypes. Hannon *et al* 2015](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/Interindividual%20methylomic%20variation%20across%20blood%20cortex%20and%20cerebellum%20implications%20for%20epigenetic%20studies%20of%20neurological%20and%20neuropsychiatric.pdf)

An outlne of our project can be found in our [Project proposal](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md)

Introduction
------------

DNA methylation (DNAm), or covalent attachment of a methyl group to the 5’ position of cytosine bases located in CpG dinucleotides, plays a crucial role in maintaining patterns of gene expression during human development and aging. Along with other epigenetic modifications, DNAm is sensitive to environmental influences and can change over an individual’s lifespan. Thus, understanding the human methylome is important for determining both biomarkers for - and direct pathways implicating - health and disease [1](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). Aberrant DNAm patterns have been correlated with common neurodegenerative disorders, including Alzheimer’s Disease and Parkinson’s Disease [2](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md), as well as with mental disorders such as schizophrenia [3,4](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). However, much of the literature on either the diseased or healthy brain methylome fails to separate DNAm data by cell type composition - a major driver of DNAm variability - or by brain region [5-8](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). 

Motivation
---------
The [Hannon *et al* paper](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/Interindividual%20methylomic%20variation%20across%20blood%20cortex%20and%20cerebellum%20implications%20for%20epigenetic%20studies%20of%20neurological%20and%20neuropsychiatric.pdf) compared methylation levels across brain regions and the blood to determine differential methylation between the regions. However their analysis did not include normalizing for cell type which is a major driver of DNAm variability [5-8](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). Correcting for cell type would enable a comparison of differential methylation independant of cell type. Therefore, we wanted to further evaluate whether cell type correction is nessisary for this analysis and determine differentially methylated regions that would be missed if cell type correction was not included in the analysis.

Objectives
----------

1.To determine whether cell type correction between brain regions is necessary in the analysis of Illumina HumanMethylation450 BeadChip array data
2.To investigate probes differentially methylated between cerebellum and cortex regions

Data
-----
Samples of entorhinal cortex (EC), prefrontal cortex (PFC), superior temporal gyrus (STG), and cerebellum (CER) tissue, as well as blood samples, were collected from patients in archives of the MRC London Neurodegenerative Disease Brain Bank. Degradation and purity analysis were done on the samples using phenol-chloroform extraction and bisulfite conversion was performed using the Zymo EZ 96 DNA methylation kit (Zymo Research). A total of 531 samples from 122 individuals were included for further analysis

DNA methylation was determined using Illumina Infinium HumanMethylation450 BeadChip (Illumina) and an Illumina HiScan System (Illumina). Raw signal intensities of the probe were quantified using the Illumina Genome Studio software and converted into beta values using the methylumi package in R. 

A description of the data and how it was obtained can be found [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/README.md)

We downloaded both the Dasen normalized data and the QN normalized data. Meta data was obtained from GEO using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Acquiring%20GEO%20meta%20data.Rmd) code

[Analysis and tasks directory](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/README.md)
----------------------------
Methodology with division of labour and links to analysis. 

1. Processing
- Dasen data was filtered using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering%20Dasen.Rmd) code and QN data was filtered using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering%20QN.Rmd) code (**Hilary,Sam**)
- [Blood samples were then removed from the data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Create.brain.only.Rmd) (**Cassia**)
- The QN-normalized beta values were then [BMIQ normalized](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/BMIQ_final.md) (**Sam**)
- The BMIQ data set went through [cell type prediction](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/Cell%20Type%20Prediction.Rmd) and [combat/PCA and cell type correction](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/PCA%20%26%20ComBat.md) (**Sam**)

A overview of the processing steps are visualized in this pipeline:

![](/Images/Pipeline_of_Methods.png)

2. Analysis
- Sample to sample correlation plots and heatmaps were created from various codes that can be found with explanations in the working_code [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md). This is the final code for the [BMIQ data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(BMIQ).Rmd) and [Dasen data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd) (**Cassia, Lisa, Hilary**)
- Differentially methylated probe (DMP) and differentially methylated region (DMR) analysis was conduced using various codes
  +Linear modelling for delta betas and generating volcano plots was performed on batch-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/DMR%20batch%20cor.md) and cell-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/DMR%20cell%20cor.md)(**Randip**)
 + DMP analysis was conducted on the batch-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially%20Methylated%20Probe%20Analysis%20-%20Batch%20Corrected%20Only%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially_Methylated_Probe_Analysis_-_Batch_Corrected_Only__DMR_Setup__Final.md). The same analysis was also conducted on the batch and cell-type corrected dataset [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially%20Methylated%20Probe%20Analysis%20-%20Cell-Type%20Corrected%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially_Methylated_Probe_Analysis_-_Cell-Type_Corrected__DMR_Setup__Final.md), where comparisons of significant probes (at an FDR <= 0.05 and with an absolute value delta beta >= 0.05) found in these two datasets were also made. (**Randip**)
- Boxplots for the most differentially methylated probes were generated [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd). Preliminary versions of the code can be found in the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md) (**Cassia**)
- Wilcoxon test was performed to compare probe values between cerebellum and cortex for the most significantly differentially methylated probe. The most significantly differentially methylated probe that was found only in the list from cell type corrected, only in non cell-type corrected, and in the overlap was chosen for visualization with the boxplots. The code for the wilcoxon test can be found [here](WilcoxTestProbes_analysis.R) (**Lisa))

3. Administrative tasks and preparation of deliverables:
- [Project Proposal](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md)(**Lead:Sam, Everyone**)
- Repo Readme directories(**Lead: Cassia, everyone**)
- [Progress report](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/progress_report.md)(**Lead: Sam + Hilary, everyone**)
- [Poster](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/Poster.pdf)(**Lead: Cassia + Hilary + Sam, everyone**)
- [Bibliography](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md)(**Everyone**)

[A summary of major analysis, their motivations, and their main results](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/results)
----------------------------------------------------------------------

Discussion
------------

Deliverables
--------------

- [Project proposal](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md)
- [Progress report](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/progress_report.md)
- [Poster](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/Poster.pdf)

Repo Directory
----------

We aim to study the differences in methylation between cerebellum and cortex. For more information on our project please visit our project proposal [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md)

We will be using methylation data from Hannon *et al.* paper that is publically available. More information about our data set please check it out [here](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data)

Our [data folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data) consists of two subfolders
  - [raw data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/raw_data)
  - [processed data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/processed_data)
  our data set is too large to upload into the repo but we have included the raw and processed meta data.
  
Our [background information folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/background_information) is there to provide your with some relevant papers to our project

Our [src folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src) is where the codes used for our analysis are stored. 
  - in [working codes](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src/working_codes) you will find parts of codes that were created and then shared to help create code for our final analysis codes. Explain through the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md)
  - in [final codes](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src/final_codes) We have the final copies of our code used to optain our results. Detailed explanation of which step they were used for can be found in the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/README.md)
  
Our [results folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/results) contains the decriptions of the aims and conclusions of our analysis steps and the important final images derived from the analysis

**Other folders**

Our [resources folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/resources) includes codes that we modified for our analysis. These codes have been referenced in our src where appropriate

Our [images folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/Images) includes images we shared with eachother to keep updates on the project and help with trouble shooting

Our final poster can be found in the root directory or [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/Poster.pdf)





