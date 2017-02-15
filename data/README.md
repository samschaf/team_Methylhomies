
This folder contains data used in the project and detailed description of thee data 

Subdivisions:
- [Raw Data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/raw_data)
- [Processed Data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/processed_data)

We are using publically available data that is used and described in the [Hannon *et al.* paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844197) 

Details
-------
Samples of entorhinal cortex (EC), prefrontal cortex (PFC), superior temporal gyrus (STG), and cerebellum (CER) tissue, as well as blood samples, were collected from patients in archives of the MRC London Neurodegenerative Disease Brain Bank. Degredation and purity analysies were done on the samples using phenol-chloroform extraction. A total of 531 samples were included for further analysis.

DNA methylation was determined using Illumina Infinium HumanMethylation450 BeadChip (Illumina) and an Illumina HiScan System (Illumina). Raw signal intensities of the probe were quantified using the Illumina Genome Studio software and converted into beta values using the [*methylumi*](https://bioconductor.org/packages/release/bioc/html/methylumi.html) package in R. We attained this beta value data through GEO accession identifier [GSE59685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59685)). Explanations and example caclulations of beta values can be found in [this paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587).

The following code was used to retreive the data:

```{rret, eval=FALSE}
source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
filePaths = getGEOSuppFiles("GSE59685")
filePaths
``` 
To read file containg the beta values into R the following code was used:

```{r beta, eval=FALSE}
GSE59685_betas <- read.csv("~/path/to/file/GSE59685/GSE59685_betas.csv.gz", 
comment.char="#", stringsAsFactors=FALSE)
```
The raw data set of beta values has 485579 total rows, 2 are sample identifiers and the rest are the 485577 probes. The 532 columns reflect the 521 samples and 1 column for the probe identifiers.
