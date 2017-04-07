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
| Lisa Wei | Bioinformatics  | Cancer genomics, sequence analysis | Help maintain/organize Github and overall methylation level differences between brain subregions |

Our project uses data that was generated for the paper [Interindividual methylomic variation across blood, cortex, and cerebellum: implications for epigenetic studies of neurological and neuropsychiatric phenotypes. Hannon *et al* 2015](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/Interindividual%20methylomic%20variation%20across%20blood%20cortex%20and%20cerebellum%20implications%20for%20epigenetic%20studies%20of%20neurological%20and%20neuropsychiatric.pdf)

An outline of our project can be found in our [Project proposal](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md).

Introduction
------------

DNA methylation (DNAm), or covalent attachment of a methyl group to the 5’ position of cytosine bases located in CpG dinucleotides, plays a crucial role in maintaining patterns of gene expression during human development and aging. Along with other epigenetic modifications, DNAm is sensitive to environmental influences and can change over an individual’s lifespan. Thus, understanding the human methylome is important for determining both biomarkers for - and direct pathways implicating - health and disease [1](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). Aberrant DNAm patterns have been correlated with common neurodegenerative disorders, including Alzheimer’s Disease and Parkinson’s Disease [2](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md), as well as with mental disorders such as schizophrenia [3,4](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). However, much of the literature on either the diseased or healthy brain methylome fails to separate DNAm data by cell type composition - a major driver of DNAm variability - or by brain region [5-8](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). 

Motivation
---------
The [Hannon *et al* paper](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/Interindividual%20methylomic%20variation%20across%20blood%20cortex%20and%20cerebellum%20implications%20for%20epigenetic%20studies%20of%20neurological%20and%20neuropsychiatric.pdf) compared methylation levels across brain regions and the blood to determine differential methylation between the regions. However their analysis did not include normalizing for cell type which is a major driver of DNAm variability [5-8](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). Correcting for cell type would enable a comparison of differential methylation independent of cell type. Therefore, we wanted to further evaluate whether cell type correction is necessary for this analysis and determine differentially methylated regions that would be missed if cell type correction was not included in the analysis.

Objectives
----------

1. To determine whether cell type correction between brain regions is necessary in the analysis of Illumina HumanMethylation450 BeadChip array data
2. To investigate probes differentially methylated between cerebellum and cortex regions

Data
-----
Samples of entorhinal cortex (EC), prefrontal cortex (PFC), superior temporal gyrus (STG), and cerebellum (CER) tissue, as well as blood samples, were collected from patients in archives of the MRC London Neurodegenerative Disease Brain Bank. DNA was extracted by phenol-chloroform, followed by degradation and purity analysis; bisulfite conversion was then performed using the Zymo EZ 96 DNA methylation kit (Zymo Research). A total of 531 samples from 122 individuals were included for further analysis.

DNA methylation was determined using Illumina Infinium HumanMethylation450 BeadChip (Illumina) and an Illumina HiScan System (Illumina). Raw signal intensities of the probe were quantified using the Illumina Genome Studio software and converted into beta values using the methylumi package in R. 

A description of the data and how it was obtained can be found [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/data/README.md).

We downloaded both the Dasen normalized data and the QN normalized data. Meta data was obtained from GEO using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Acquiring%20GEO%20meta%20data.Rmd) code.

[Analysis and tasks directory](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/README.md)
----------------------------
Methodology with division of labour and links to analysis. 

1. Processing
- Dasen data was filtered using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering%20Dasen.Rmd) code and QN data was filtered using [this](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Probe%20Filtering%20QN.Rmd) code (**Hilary,Sam**)
- [Blood samples were then removed from the data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Create.brain.only.Rmd) (**Cassia**)
- The QN-normalized beta values were then [BMIQ normalized](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/BMIQ_final.md) (**Sam**)
- The BMIQ data set went through [cell type prediction](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/Cell%20Type%20Prediction.Rmd) and [Principal Component Analsyis, ComBat, and CETS cell type correction](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/PCA%20%26%20ComBat.md) (**Sam**)

A overview of the processing steps are visualized in this pipeline:

![](/Images/Pipeline_of_Methods.png)

2. Analysis
- Sample to sample correlation plots and heatmaps were created from various codes that can be found with explanations in the working_code [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md). This is the final code for the [BMIQ data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(BMIQ).Rmd) and [Dasen data](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd) (**Cassia, Lisa, Hilary**)
- Differentially methylated probe (DMP) and differentially methylated region (DMR) analysis was conduced using various codes
  and linear modelling for delta betas and generating volcano plots was performed on batch-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/DMR%20batch%20cor.md) and cell-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/DMR%20cell%20cor.md) (**Randip, Sam**)
 + DMP analysis was conducted on the batch-corrected data [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially%20Methylated%20Probe%20Analysis%20-%20Batch%20Corrected%20Only%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially_Methylated_Probe_Analysis_-_Batch_Corrected_Only__DMR_Setup__Final.md). The same analysis was also conducted on the batch and cell-type corrected dataset [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially%20Methylated%20Probe%20Analysis%20-%20Cell-Type%20Corrected%20(DMR%20Setup)%20Final.Rmd) and [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/DMR%20Analysis/Differentially_Methylated_Probe_Analysis_-_Cell-Type_Corrected__DMR_Setup__Final.md), where comparisons of significant probes (at an FDR <= 0.05 and with an absolute value delta beta >= 0.05) found in these two datasets were also made (**Randip**)
- Boxplots for the most differentially methylated probes were generated [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/Heatmaps%20(dasen).Rmd). Preliminary versions of the code can be found in the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md) (**Cassia**)
- Wilcoxon test was performed to compare probe values between cerebellum and cortex for the most significantly differentially methylated probe. The most significantly differentially methylated probe that was found only in the list from cell type corrected, only in non cell-type corrected, and in the overlap was chosen for visualization with the boxplots. The code for the wilcoxon test can be found [here](WilcoxTestProbes_analysis.R) (**Lisa**)

3. Administrative tasks and preparation of deliverables:
- [Project Proposal](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md) (**Lead:Sam +Lisa, Everyone**)
- Repo README directories (**Lead: Cassia, everyone**)
- [Progress report](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/progress_report.md) (**Lead: Sam + Hilary, everyone**)
- [Poster](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/Poster.pdf) (**Lead: Cassia + Hilary + Sam, everyone**)
- [Bibliography](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md) (**Everyone**)

[A summary of major analysis, their motivations, and their main results](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/results)
----------------------------------------------------------------------

**Processing**

Probe filtering was performed to remove probes from the data which would confound our analysis, including probes with SNPs, probes with NAs, and cross-hybridizing probes. Data was also filtered to contain only brain samples (originally had brain + blood). Prior to filtering, there were 485577 probes. Filtering removed 65 SNP probes, 20869 probes that had a SNP within, 11475 probes hybridzing to sex chromosomes, 10673 cross hybrizing X and Y probes, 27415 cross hybridizing autosomal probes, and 176 NA probes. A total of 414904 probes remained. 

We had originally downloaded the dasen normalized data that was used in the Hannon paper. After probe filtering and normalization the PCA analysis was unusual. After further investigation into the literature we determined that QN and BMIQ normalization would be more appropriate for our analysis. BMIQ normalization was done on the QN data to correct for intra-array differences between probe type I and II chemistries that cause differing beta value distributions. This was performed on the QN-normalized data we had decided to move forward with, in order to make direct comparisons between any CpG site across any two brain regions. 

Multiple BeadChips were used to generate this data; we aimed to correct for variation due to the individual chip and position on the chip with Combat. Since we would be analyzing the effect of cell type, we also used the CETS package to predict proportions of neurons and glia based on the non-corrected and the corrected beta values. These neuronal propotions were attached to the meta data and used as continous variables in PCA analysis. Before batch correction, there was a large effect of chip and tissue; it appeared the samples had been batched by tissue type. After batch correction, the chip and tissue effects plus much of the other variation was removed, leaving only a Neuron (cell type) effect. Linear regression and the predicted neuronal proportion was used to fit our beta values to reference profiles of neuronal and glial methylation, removing variation due to cell type. After cell type correction, Neuron comprised a smaller propotion of total variance, and we unmasked some effects of braak stage, age, and AD disease status.

**Analysis**

Sample to sample correlations were done to determine whether the samples correlate well with each other. If processing normalized the data we should see that samples will correlate with samples both within the same cell type and outside of the same cell type. QN+BMIQ normalization increased sample to sample correlations compared to dasen normalization. Cell type correction increases inter-tissue sample correlations in the dasen model, in line with the purpose of the normalization method to increase the comparability between tissues. This is unable to be seen in the BMIQ data as the samples are already highly correlated.

DMP and DMR analyses were done to identify and compare differentially methylated probes (at an FDR <= 0.05 and absolute delta beta >= 0.05) between tissue type, and differentially methylated regions in both the batch-corrected only data, and the both batch-corrected and cell-type corrected data. We found that correcting for cell type when examining cortex and cerebellum methylation differences may alter the number of significant probes and regions identified. We then plotted individual probes found only within the non cell corrected, cell corrected, and both analyses and performed wilcoxon tests to show statistical differences between tissues regions in probes found only in the cell type corrected analysis, only in the non cell type corrected analysis, and in both. 

![](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/results/Venndiagram/Screen%20Shot%202017-04-07%20at%2012.29.33%20PM.png)

**Summary**   
• BMIQ normalization is more suitable for cross-sample comparisons   
• Technical and biological variability to account for are batch effects and cell type variation across brain regions to allow for comparison between tissues   
• Differential methylation analysis reveals DMPs that are attributed to cell type variation and those that are independent of cell type variation   
• Failing to correct for cell type may result in identification of false positive DMPs/DMRs, with implications for biological conclusions concerning clinically relevant loci

Discussion
------------

Our project analyzed different normalization methods to show that: 1) cell type has an effect on methylation; and 2) cell type correction is necessary to determine regions of differential methylation that are independent of neuronal composition. Differentially methylated regions independent of cell type that our analysis found include: TRAK1, ERCC5, METTL21EP, RGS6, and WSCD1, among others. Previous work has shown that methylation in certain brain regions may play a role in the development of Alzheimer's disease. Correcting for cell type may be a critical step in determining genes that are differentially methylated between brain regions and could be further used to compare between healthy and diseased individuals [13](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/background_information/README.md). For example, [RGS6](http://www.genecards.org/cgi-bin/carddisp.pl?gene=RGS6), which was only found in our CETS analysis, has been linked to Alzheimer's disease.


Deliverables
--------------

- [Project proposal](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md)
- [Progress report](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/progress_report.md)
- [Poster](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/Poster.pdf)

Repo Directory
----------

We aimed to study the differences in methylation between cerebellum and cortex. For more information on our project please visit our project proposal [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/project_proposal.md).

We used methylation data from [Hannon *et al.* paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844197/). More information about our data set please check it out [here](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data).

Our [data folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data) consists of two subfolders:
  - [raw data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/raw_data)
  - [processed data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/processed_data)   
Our data set is too large to upload into the repo, but we have included the raw and processed meta data.
  
Our [background information folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/background_information) is there to provide you with some relevant papers to our project.

Our [src folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src) is where the codes used for our analysis are stored: 
  - In [working codes](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src/working_codes) you will find parts of codes that were created and then shared to help create code for our final analysis codes. See the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/working_codes/README.md) for more information   
  - In [final codes](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src/final_codes) We have the final copies of our code used to obtain our results. Detailed explanations of which step each was used for can be found in the [README](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/src/final_codes/README.md)
  
Our [results folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/results) contains decriptions of the aims and conclusions of our analysis steps and the important final images derived from the analysis.

**Other folders**

Our [resources folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/resources) includes codes that we modified for our analysis. These codes have been referenced in our src where appropriate.

Our [images folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/Images) includes images we shared with each other to keep updates on the project and help with trouble shooting.

Our final poster can be found in the root directory or [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/Poster.pdf).





