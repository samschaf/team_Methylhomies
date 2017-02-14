
This folder contains data used in the project and detailed description of thee data 

Subdivisions:
- Raw Data
- Processed Data

We are using publically available data that is used and described in the [Hannon *et al.* paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844197) 

Details
-------
Samples of entorhinal cortex (EC), prefrontal cortex (PFC), superior temporal gyrus (STG), and cerebellum (CER) tissue were collected from 117 patients in archives of the MRC London Neurodegenerative Disease Brain Bank. Degredation and purity analysies were done on the samples using phenol-chloroform extraction. 

DNA methylation was determined using Illumina Infinium HumanMethylation450 BeadChip (Illumina) and an Illumina HiScan System (Illumina). Raw signal intensities of the probe were quantified using the Illumina Genome Studio software. We attained this raw data through GEO accession identifier [GSE59685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59685)).

The following code was used to retreive the data:

```{rret, eval=FALSE}
source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
filePaths = getGEOSuppFiles("GSE59685")
filePaths
``` 
The raw data set has X rows and X columns
