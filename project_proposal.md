Final Project Proposal
======================

Background and rationale for the study. 
----------------------------------------
**Why is it important to answer this question? Identify knowledge gaps.**
DNA methylation (DNAm), or covalent attachment of a methyl group to the 5’ position of cytosine bases located in CpG dinucleotides, plays a crucial role in maintaining patterns of gene expression during human development and aging. Along with other epigenetic modifications, DNAm is sensitive to environmental influences and can change over an individual’s lifespan. Thus, understanding the human methylome is important for determining both biomarkers for - and direct pathways implicating - health and disease [1]. Aberrant DNAm patterns have been correlated with common neurodegenerative disorders, including Alzheimer’s Disease and Parkinson’s Disease [2], as well as with mental disorders such as schizophrenia [3,4]. However, much of the literature on either the diseased or healthy brain methylome fails to separate DNAm data by cell type composition - a major driver of DNAm variability - or by brain region [5-8]. A comprehensive characterization of differentially methylated regions (DMRs) in the brain both in normal individuals [9] and neurodegenerative disorders such as autism has been performed [10]. Using a publicly available methylation dataset for the cerebellum and cortex of 122 healthy individuals [11], we aim to confirm established DMRs between brain regions and discover new region-specific DMRs to add to the list.

What is the overall hypothesis and objectives for the project?
---------------------------------------------------------------
Our overall hypothesis is that there will be differences in methylation levels between the cerebellum and cortex, and that correcting for cell type is important for the accurate comparison of DNAm across brain regions. In addition, the presence of DMRs between regions has important biological consequences such as local activation or repression of promoters and regulatory elements, influencing gene expression patterns.
We aim to provide a baseline of variation in the normal human methylome between the cerebellum and cortex, and to determine DMRs in these regions. Ultimately, the project will reveal the magnitude and significance of regional-specific methylation differences, provide functional characterization of differentially methylated regions in the brain, and produce guidelines to analyze existing brain epigenetic literature and to plan robust experimental design moving forward.

Division of Labour
-------------------------------------------------------------------------------------------------------
**Who is going to do what? State assignment of tasks and projected contributions for each group members**.
Everyone will help with writing and data interpretation. Final poster organization will be done by 2 or 3 designated group members.

**Hilary (with backup support from Sam):**
- Data retrieval - GEOquery 
- Preprocessing data - merging metadata and raw beta values (methylation scores)
- Clean up data - filter probes (methylumi)
- Normalization - quantro (quantile) and BMIQ (probe type)
- Cell type correction - CETS
- PCA - identify factors most associated with each principal component and maybe remove unwanted effects/influences 
- Remove batch effects - ComBat  
- PCA analysis again after batch effect removal to check data
- Remove outliers - mvoutliers package, heatmap visualization for correlation

**Cassia and Lisa:**
- Answer question on importance of cell type correction:
    - Compare cell type corrected data set with the same data set where the cell type correction was not performed. 
    -  When cell type correction is applied, does it change the correlation between the cortex and cerebellum observed in Hannon et al. [11]?
- Maintain group GitHub Readme files and repositories:
    - Unify formatting of the different analyses
- Coordinate and design poster

**Lisa:**
- Compare the overall methylation level in the different brain regions (cortex vs cerebellum)
- Compare these results with Hannon et al. [11].

**Randip and Sam:**
- Look at the DMRs between brain regions. Where are they located? Is there enrichment at promoter/enhancer regions? Could this result in functional differences?
    - Filter probes for CpG island using the IlluminaHumanMethylation450k.db Bioconductor package
    - Use limma package to analyze differentially methylated regions
    - Map differentially methylated regions onto specific genomic locations using probe information
    - Perform interpretation and functional enrichment analysis using the IlluminaHumanMethylation450k.db package

**Table of group members with their background, degree, affiliations and job assignments.**

| Name | Degree/Affiliation | Speacilizations |Job Assignments | 
| ------------- | ------------- | ------------- | ------------- |
| Samantha Schaffner | Medical Genetics  | Neuroepigenetics of human development and Parkinson’s Disease | DMR analysis and help with data processing and normalization |
| Cassia Warren | Interdisciplinary Oncology   | Biomarker development for Vitamin C treatment for Pancreatic cancer | Maintain/organize Github and Cell type correction analysis|
| Hilary Brewis  | Medical Genetics  | Histone variants and chromatin structure | Normalization and data processing |
| Randip Gill  | Educational Psychology  | Social epigenetics | DMR analysis |
| Lisa Wei | Bioinformatics  | Cancer genomics, sequence analysis | Help maintain/organize Github and Overall methylation level differences between brain subregions |


Dataset
---------
**What kind of data are you working with? What is the general description and characteristics of the data? State the technology used to generate the data.**
We will use a publicly available dataset of n=122 healthy individuals aged 40-105 [11], consisting of raw intensity values obtained using the Illumina HumanMethylation450 BeadChip array platform. 

Samples of entorhinal cortex (EC), prefrontal cortex (PFC), superior temporal gyrus (STG), and cerebellum (CER) tissue were collected from patients in archives of the MRC London Neurodegenerative Disease Brain Bank. Blood samples were also used when available. Degradation and purity analyses were done on the samples following phenol-chloroform DNA extraction, and bisulfite conversion was performed using the Zymo EZ 96 DNA methylation kit (Zymo Research). A total of 531 samples were included for further analysis.
DNA methylation was determined using Illumina Infinium HumanMethylation450 BeadChip and an Illumina HiScan System. Raw signal intensities of the probe were quantified using the Illumina Genome Studio software. We attained this raw data and beta values through GEO (accession identifier GSE59685). Raw intensity values were converted into beta values using the methylumi package in R. Beta values range between 0 and 1 and are determined as a ratio of probe intensity over overall intensity [12].
The accompanying data file of beta values consists of 485579 total rows: 485577 probes (top 2 rows are sample identifiers) and 532 columns - 531 samples, 1st column for probe identifier (each column is a sample from a brain region from a specific patient, 122 patients overall). We are planning to start with this file of beta values for our subsequent analysis.

Aims and methodology:
----------------------

1) Determine if cell type correction will affect correlation between the cortex and cerebellum observed in Hannon et al. [11]. Compare cell type corrected data set with the same data set where the cell type correction was not performed. 
Statistical analysis:
- use CETS cell type correction package

2) Compare the overall methylation level in the different brain regions (cortex vs cerebellum) to establish a baseline of variation in normal individuals. 
Statistical analysis: 
- Boxplots for overall visualization (use data after cell type correction)
- Quantify correlation of methylation between regions, distribution, median/means, percent methylation of each region
    
3) Look at the DMRs between brain regions. Where are they located? Is there enrichment at promoter/enhancer regions? Could this result in functional differences? 
Statistical analysis: 
- Limma to identify DMRs, make age a covariate in the linear regression model
- IlluminaHumanMethylation450k.db Bioconductor package for grouping probes by CpG islands, mapping onto specific location on chromosome, DMR functional analysis  

References
------------

[1] Bernstein, B. E., Meissner, A., & Ladner, E. S. (2007). The Mammalian Epigenome. Cell. 128(4): 669-681. 

[2] Sanchez-Mut, J. V., Heyn, H., Vidal, E., Moran, S., Sayols, S., Delgado-Morales, R., … Esteller, M. (2016). Human DNA methylomes of neurodegenerative diseases show common epigenomic patterns. Transl Pysch. 6: e817.

[3] Huang, H. S., & Akbarian, S. (2007). GAD1 mRNA expression and DNA methylation in prefrontal cortex of subjects with schizophrenia. PloS one, 2(8), e809.

[4] Huang, H. S., Matevossian, A., Whittle, C., Kim, S. Y., Schumacher, A., Baker, S. P., & Akbarian, S. (2007). Prefrontal dysfunction in schizophrenia involves mixed-lineage leukemia 1-regulated histone methylation at GABAergic gene promoters. Journal of Neuroscience, 27(42), 11254-11262.

[5] Shin, J., Ming, G., & Song, H. (2014). Decoding neural transcriptomes and epigenomes via high-throughput sequencing. Nat Neurosci 17(11): 1463-1475.

[6] Jaffe, A. E. & Irizarry, R. A. (2014). Accounting for cellular heterogeneity is critical in epigenome-wide association studies. Genome Biology. 15: R31.

[7] Guintivano, J., Aryee, M. J., & Kaminsky, Z. A. (2013). A cell epigenotype specific model for the correction of brain cellular heterogeneity bias and its application to age, brain region and major depression. Epigenetics. 8(3): 290-302.

[8] Montano, C. M., Irizarry, R. A., Kaufmann, W. E., Talbot, K., Gur, R. E., Feinberg, A. P., & Taub, M. A. (2013). Measuring cell-type specific differential methylation in human brain tissue. Genome Biol. 14: R94.

[9] Ladd-Acosta, C. et al. (2007). DNA methylation signatures within the human brain. The 
American Journal of Human Genetics. 81(6): 1304-1315.

[10] Ladd-Acosta, C., Hansen, K.D., Briem, E., Fallin, M.D., Kaufmann, W.E., Feinberg, A.P. 
(2014). Common DNA methylation alterations in multiple brain regions in autism. Molecular Psychiatry. 19: 862-871.   

[11] Hannon, E., Lunnon, K., Schalkwyk, L., & Mill, J. (2015). Interindividual methylomic variation across blood, cortex, and cerebellum: implications for epigenetic studies of neurological and neuropsychiatric phenotypes. Epigenetics. 10(11): 1024-1032.

[12]  Du P, Zhang X, Huang C-C, Jafari N, Kibbe WA, Hou L, et al. (2010). Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. BMC Bioinformatics. 11(1): 587. 


