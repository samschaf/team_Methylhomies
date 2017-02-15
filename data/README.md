
This folder contains data used in the project and detailed description of the data 

Subdivisions:
- [Raw Data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/raw_data)
- [Processed Data](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/data/processed_data)

We are using publically available data that is used and described in the [Hannon *et al.* paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844197) 

Details
-------
Samples of entorhinal cortex (EC), prefrontal cortex (PFC), superior temporal gyrus (STG), and cerebellum (CER) tissue, as well as blood samples, were collected from patients in archives of the MRC London Neurodegenerative Disease Brain Bank. Degradation and purity analysis were done on the samples using phenol-chloroform extraction and bisulfite conversion was performed using the Zymo EZ 96 DNA methylation kit (Zymo Research). A total of 531 samples were included for further analysis

DNA methylation was determined using Illumina Infinium HumanMethylation450 BeadChip (Illumina) and an Illumina HiScan System (Illumina). Raw signal intensities of the probe were quantified using the Illumina Genome Studio software and converted into beta values using the [*methylumi*](https://bioconductor.org/packages/release/bioc/html/methylumi.html) package in R. We attained this beta value data through GEO accession identifier [GSE59685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59685)). Explanations and example calculations of beta values can be found in [this paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587).

The following code was used to retreive the data:

```{rret, eval=FALSE}
source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
filePaths = getGEOSuppFiles("GSE59685")
filePaths
``` 
To read file containingg the beta values into R the following code was used:

```{r beta, eval=FALSE}
GSE59685_betas <- read.csv("~/path/to/file/GSE59685/GSE59685_betas.csv.gz", 
comment.char="#", stringsAsFactors=FALSE)
```
The raw data set of beta values has 485579 total rows, 2 are sample identifiers and the rest are the 485577 probes. These probes relate to the different areas of the DNA that is being tested for DNA methylation. The 532 columns reflect the 521 samples and 1 column for the probe identifiers.

Descriptions
-------------

**Bisulfite conversion** is a technique used to determine methylation on DNA. Treatment of DNA with bisulfite causes cytosine residues to convert to uracil. However, when the cytosine is methylated it is protected from this conversion. Therefore, methylation can be detected through sequencing as reading either a C when it is methylated or a T when it was non-methylated.

**Illumina Infinium HumanMethylation450 BeadChip (Illumina)** is a probe based way to determine methylation rather than sequencing. This chip can detect over 480,000 methylation sites across the human genome. On the chip there are beads for both the methylated and unmethylated loci. The DNA hybridizes to the bead which represents its methylated status at that region. The ddCTP (methylated) loci are labeled with biotin while the other ddNTPs (unmethylated) are labeled with 2,4-dinitrophenol. The chip is then stained with antibodies that differentiate the labels and scanned to determine intensity values. A value of 0 indicates non-methylation, 0.5 indicates one copy of the locus is methylated, and 1 indicates both copies are methylated.



