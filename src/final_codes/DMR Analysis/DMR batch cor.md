---
title: "MAR 30 CELL-TYPE AND BATCH CORRECTED DMR"
output: github_document
---

```r
setwd("/home/sschaffner/team_Methylhomies")
#Loading necessary packages: Functions of packages are noted here also for HW Assignment Q1-5
#Supressing warnings that were appearing to make my markdown less cluttered
options(warn=-1) 
suppressWarnings(suppressMessages(library(ggfortify)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(plotly)))
suppressWarnings(suppressMessages(library(GEOquery)))
suppressWarnings(suppressMessages(library(wateRmelon)))
suppressWarnings(suppressMessages(library(minfi)))
suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19)))
suppressWarnings(suppressMessages(library(IlluminaHumanMethylation450kmanifest)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(missMethyl)))
suppressWarnings(suppressMessages(library(matrixStats)))
suppressWarnings(suppressMessages(library(minfiData)))
suppressWarnings(suppressMessages(library(Gviz)))
suppressWarnings(suppressMessages(library(DMRcate)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(colorspace)))
suppressWarnings(suppressMessages(library(VennDiagram)))
suppressWarnings(suppressMessages(library(qpcR)))
```

```
## Error: package 'rgl' could not be loaded
```

```r
library(wateRmelon)
library(Sushi)
```

```
## Loading required package: zoo
```

```
## 
## Attaching package: 'zoo'
```

```
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```
## Loading required package: biomaRt
```

```r
library(IlluminaHumanMethylation450k.db)
library(lumi)
library(ggfortify) #For plotting
library(gridExtra) #For plotting
library(plotly) #For plotting
library(GEOquery) #For downloading genetic data
library(wateRmelon) #For model/analyses
library(limma) #For model/analyses
library(minfi) #For model/analyses
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #Gene annotation for analyses
library(IlluminaHumanMethylation450kmanifest) #Gene annotation for analyses
library(RColorBrewer) #Colors for plotting
library(missMethyl)  #For model/analyses
library(matrixStats) #For data manipulation
library(minfiData) #For model/analyses
library(Gviz) #For DMR plotting
library(DMRcate) #For DMR model/analyses
library(gplots) #For plotting
library(ggplot2) #For plotting
library(stringr) #For data manipulation
library(tidyverse) #For data manipulation
library(data.table) #For data manipulation
library(colorspace) #Colors for plotting
library(VennDiagram) #Necessary for the overlap counts/venn diagram
library(qpcR) #Necessary for the overlap counts/venn diagram
```

```
## Loading required package: rgl
```

```
## Error: package 'rgl' could not be loaded
```

```r
library(methylumi)
library(parallel)
library(dplyr)
```


```r
#Assigning annotation information for CpG probes on the 450k array
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
```

```
## DataFrame with 6 rows and 33 columns
##                    chr       pos      strand        Name    AddressA
##            <character> <integer> <character> <character> <character>
## cg00050873        chrY   9363356           -  cg00050873    32735311
## cg00212031        chrY  21239348           -  cg00212031    29674443
## cg00213748        chrY   8148233           -  cg00213748    30703409
## cg00214611        chrY  15815688           -  cg00214611    69792329
## cg00455876        chrY   9385539           -  cg00455876    27653438
## cg01707559        chrY   6778695           +  cg01707559    45652402
##               AddressB                                          ProbeSeqA
##            <character>                                        <character>
## cg00050873    31717405 ACAAAAAAACAACACACAACTATAATAATTTTTAAAATAAATAAACCCCA
## cg00212031    38703326 CCCAATTAACCACAAAAACTAAACAAATTATACAATCAAAAAAACATACA
## cg00213748    36767301 TTTTAACACCTAACACCATTTTAACAATAAAAATTCTACAAAAAAAAACA
## cg00214611    46723459 CTAACTTCCAAACCACACTTTATATACTAAACTACAATATAACACAAACA
## cg00455876    69732350 AACTCTAAACTACCCAACACAAACTCCAAAAACTTCTCAAAAAAAACTCA
## cg01707559    64689504 ACAAATTAAAAACACTAAAACAAACACAACAACTACAACAACAAAAAACA
##                                                     ProbeSeqB        Type
##                                                   <character> <character>
## cg00050873 ACGAAAAAACAACGCACAACTATAATAATTTTTAAAATAAATAAACCCCG           I
## cg00212031 CCCAATTAACCGCAAAAACTAAACAAATTATACGATCGAAAAAACGTACG           I
## cg00213748 TTTTAACGCCTAACACCGTTTTAACGATAAAAATTCTACAAAAAAAAACG           I
## cg00214611 CTAACTTCCGAACCGCGCTTTATATACTAAACTACAATATAACGCGAACG           I
## cg00455876 AACTCTAAACTACCCGACACAAACTCCAAAAACTTCTCGAAAAAAACTCG           I
## cg01707559 GCGAATTAAAAACACTAAAACGAACGCGACGACTACAACGACAAAAAACG           I
##               NextBase       Color    Probe_rs Probe_maf      CpG_rs
##            <character> <character> <character> <numeric> <character>
## cg00050873           A         Red          NA        NA          NA
## cg00212031           T         Red          NA        NA          NA
## cg00213748           A         Red          NA        NA          NA
## cg00214611           A         Red          NA        NA          NA
## cg00455876           A         Red          NA        NA          NA
## cg01707559           A         Red          NA        NA          NA
##              CpG_maf      SBE_rs   SBE_maf           Islands_Name
##            <numeric> <character> <numeric>            <character>
## cg00050873        NA          NA        NA   chrY:9363680-9363943
## cg00212031        NA          NA        NA chrY:21238448-21240005
## cg00213748        NA          NA        NA   chrY:8147877-8148210
## cg00214611        NA          NA        NA chrY:15815488-15815779
## cg00455876        NA          NA        NA   chrY:9385471-9385777
## cg01707559        NA          NA        NA   chrY:6778574-6780028
##            Relation_to_Island
##                   <character>
## cg00050873            N_Shore
## cg00212031             Island
## cg00213748            S_Shore
## cg00214611             Island
## cg00455876             Island
## cg01707559             Island
##                                                                                                                        Forward_Sequence
##                                                                                                                             <character>
## cg00050873 TATCTCTGTCTGGCGAGGAGGCAACGCACAACTGTGGTGGTTTTTGGAGTGGGTGGACCC[CG]GCCAAGACGGCCTGGGCTGACCAGAGACGGGAGGCAGAAAAAGTGGGCAGGTGGTTGCAG
## cg00212031 CCATTGGCCCGCCCCAGTTGGCCGCAGGGACTGAGCAAGTTATGCGGTCGGGAAGACGTG[CG]TTAAAGGGCTGAAGGGGAGGGACGGAACTGACAGTCTCTGTGACAGCTCTGAGGTGGGAG
## cg00213748 TCTGTGGGACCATTTTAACGCCTGGCACCGTTTTAACGATGGAGGTTCTGCAGGAGGGGG[CG]ACCTGGGGTAGGAGGCGTGCTAGTGGTGGATGACATTGTGGCAGAGATGGAGGTGGTGGC
## cg00214611 GCGCCGGCAGGACTAGCTTCCGGGCCGCGCTTTGTGTGCTGGGCTGCAGTGTGGCGCGGG[CG]AGGAAGCTGGTAGGGCGGTTGTCGCAAGCTCCAGCTGCAGCCTCCGCCTACGTGAGAAGA
## cg00455876 CGCGTGTGCCTGGACTCTGAGCTACCCGGCACAAGCTCCAAGGGCTTCTCGGAGGAGGCT[CG]GGGACGGAAGGCGTGGGGTGAGTGGGCTGGAGATGCAGGCGCGCCCGTGGCTGTGCAGCC
## cg01707559 AGCGGCCGCTCCCAGTGGTGGTCACCGCCAGTGCCAATCCCTTGCGCCGCCGTGCAGTCC[CG]CCCTCTGTCGCTGCAGCCGCCGCGCCCGCTCCAGTGCCCCCAATTCGCGCTCGGGAGTGA
##                                                     SourceSeq Random_Loci
##                                                   <character> <character>
## cg00050873 CGGGGTCCACCCACTCCAAAAACCACCACAGTTGTGCGTTGCCTCCTCGC            
## cg00212031 CGCACGTCTTCCCGACCGCATAACTTGCTCAGTCCCTGCGGCCAACTGGG            
## cg00213748 CGCCCCCTCCTGCAGAACCTCCATCGTTAAAACGGTGCCAGGCGTTAAAA            
## cg00214611 CGCCCGCGCCACACTGCAGCCCAGCACACAAAGCGCGGCCCGGAAGCTAG            
## cg00455876 GACTCTGAGCTACCCGGCACAAGCTCCAAGGGCTTCTCGGAGGAGGCTCG            
## cg01707559 CGCCCTCTGTCGCTGCAGCCGCCGCGCCCGCTCCAGTGCCCCCAATTCGC            
##            Methyl27_Loci UCSC_RefGene_Name        UCSC_RefGene_Accession
##              <character>       <character>                   <character>
## cg00050873                  TSPY4;FAM197Y2        NM_001164471;NR_001553
## cg00212031                          TTTY14                     NR_001543
## cg00213748                                                              
## cg00214611                   TMSB4Y;TMSB4Y           NM_004202;NM_004202
## cg00455876                                                              
## cg01707559               TBL1Y;TBL1Y;TBL1Y NM_134259;NM_033284;NM_134258
##              UCSC_RefGene_Group     Phantom         DMR    Enhancer
##                     <character> <character> <character> <character>
## cg00050873         Body;TSS1500                                    
## cg00212031               TSS200                                    
## cg00213748                                                         
## cg00214611        1stExon;5'UTR                                    
## cg00455876                                                         
## cg01707559 TSS200;TSS200;TSS200                                    
##                     HMM_Island Regulatory_Feature_Name
##                    <character>             <character>
## cg00050873   Y:9973136-9976273                        
## cg00212031 Y:19697854-19699393                        
## cg00213748   Y:8207555-8208234                        
## cg00214611 Y:14324883-14325218     Y:15815422-15815706
## cg00455876   Y:9993394-9995882                        
## cg01707559   Y:6838022-6839951                        
##                          Regulatory_Feature_Group         DHS
##                                       <character> <character>
## cg00050873                                                   
## cg00212031                                                   
## cg00213748                                                   
## cg00214611 Promoter_Associated_Cell_type_specific            
## cg00455876                                                   
## cg01707559
```

```r
#Load data and meta files
load("GSE43414_batch_cor.RData", verbose=TRUE) 
```

```
## Loading objects:
##   GSE43414_batch_cor
```

```r
load("Meta_batch_cor.RData", verbose=TRUE)
```

```
## Loading objects:
##   meta
```

```r
## transpose data such that probe names are colnames, and rows are patient samples
transpose_GSE43414_batch_cor <- t(GSE43414_batch_cor)
## order metadata by brain region to remove lunnon NA cases
meta_order_by_braak <- meta %>% arrange(braak.stage)
meta_order_by_braakdf <- as.data.frame(meta_order_by_braak)

#Remordering samples in beta data based on lunnon et al.
matches_GSE43414_batch_cor <- match(meta_order_by_braakdf$barcode, rownames(transpose_GSE43414_batch_cor))
GSE43414_batch_cor_sorted_by_braak <- t(transpose_GSE43414_batch_cor[matches_GSE43414_batch_cor,])
GSE43414_batch_cor_sorted_by_braakdf <- as.data.frame(GSE43414_batch_cor_sorted_by_braak)

#Removing NA and braak exclude samples in meta 
meta_order_by_braakdf <- meta_order_by_braakdf[-c(280:432), ]

#Removing the NA and braak exclude samples in data
GSE43414_batch_cor_sorted_by_braakdf <- subset(GSE43414_batch_cor_sorted_by_braakdf, select = -c(280:432))
GSE43414_batch_cor_sorted_by_braakmat <- as.matrix(GSE43414_batch_cor_sorted_by_braakdf)

#Create broad regions column
meta_order_by_braakdf$broad_regions <- ifelse(meta_order_by_braakdf$Tissue == "cerebellum", "cerebellum","cortex")

#Making tissue names and braak stage syntactically valid for analysis later (removes space, adds ".", numeric stage)
meta_order_by_braakdf$braak.stage <- as.numeric(meta_order_by_braakdf$braak.stage)

str(meta_order_by_braakdf)
```

```
## 'data.frame':	279 obs. of  17 variables:
##  $ series_id        : chr  "GSE43414" "GSE43414" "GSE43414" "GSE43414" ...
##  $ gsm              : chr  "GSM1068965" "GSM1069176" "GSM1069080" "GSM1069136" ...
##  $ Subject          : chr  "1" "1" "1" "6" ...
##  $ barcode          : chr  "6042316048_R05C01" "6042316103_R06C02" "6969568118_R03C02" "6042316054_R04C01" ...
##  $ lunnon.et.al     : chr  "TRUE" "TRUE" "TRUE" "TRUE" ...
##  $ tissue.code      : chr  "A" "F" "E" "F" ...
##  $ braak.stage      : num  0 0 0 1 1 1 1 1 1 1 ...
##  $ Sex              : chr  "FEMALE" "FEMALE" "FEMALE" "MALE" ...
##  $ ad.disease.status: chr  "C" "C" "C" "C" ...
##  $ age.brain        : num  82 82 82 78 85 85 92 78 78 85 ...
##  $ age.blood        : chr  "79" "79" "79" "78" ...
##  $ Tissue           : chr  "frontal cortex" "superior temporal gyrus" "entorhinal cortex" "superior temporal gyrus" ...
##  $ Neuron           : num  0.442 0.486 0.332 0.408 0.306 ...
##  $ Glia             : num  0.558 0.514 0.668 0.592 0.694 ...
##  $ chip             : chr  "6042316048" "6042316103" "6969568118" "6042316054" ...
##  $ row              : chr  "05" "06" "03" "04" ...
##  $ broad_regions    : chr  "cortex" "cortex" "cortex" "cortex" ...
```

```r
str(GSE43414_batch_cor_sorted_by_braakdf)
```

```
## 'data.frame':	338359 obs. of  279 variables:
##  $ 6042316048_R05C01: num  0.565 0.75 0.345 0.798 0.501 ...
##  $ 6042316103_R06C02: num  0.571 0.814 0.155 0.771 0.469 ...
##  $ 6969568118_R03C02: num  0.646 0.828 0.267 0.827 0.525 ...
##  $ 6042316054_R04C01: num  0.603 0.806 0.226 0.801 0.429 ...
##  $ 6042316063_R05C01: num  0.556 0.777 0.356 0.783 0.412 ...
##  $ 6042316069_R03C01: num  0.626 0.806 0.266 0.827 0.461 ...
##  $ 6042316094_R05C01: num  0.583 0.859 0.339 0.837 0.554 ...
##  $ 6042316099_R01C01: num  0.567 0.782 0.288 0.808 0.395 ...
##  $ 6042316127_R03C01: num  0.681 0.762 0.354 0.831 0.488 ...
##  $ 6057825014_R02C02: num  0.621 0.954 0.387 0.815 0.504 ...
##  $ 6057825017_R06C01: num  0.694 0.945 0.31 0.832 0.509 ...
##  $ 6969568084_R02C01: num  0.542 0.793 0.219 0.814 0.425 ...
##  $ 6969568118_R01C02: num  0.583 0.795 0.473 0.762 0.593 ...
##  $ 7786923046_R03C02: num  0.532 0.864 0.183 0.868 0.487 ...
##  $ 7786923107_R04C01: num  0.598 0.865 0.258 0.862 0.513 ...
##  $ 6042316066_R05C01: num  0.694 0.801 0.291 0.84 0.472 ...
##  $ 6042316103_R02C01: num  0.614 0.81 0.22 0.812 0.496 ...
##  $ 6042316103_R03C01: num  0.491 0.78 0.182 0.792 0.41 ...
##  $ 6042316121_R04C02: num  0.542 0.799 0.257 0.837 0.42 ...
##  $ 6929718123_R02C02: num  0.665 0.54 0.497 0.746 0.426 ...
##  $ 6969568084_R03C02: num  0.581 0.78 0.292 0.818 0.424 ...
##  $ 6969568087_R03C01: num  0.542 0.797 0.244 0.802 0.417 ...
##  $ 6969568087_R03C02: num  0.693 0.815 0.26 0.841 0.588 ...
##  $ 7796806002_R04C02: num  0.557 0.782 0.225 0.861 0.522 ...
##  $ 7796806022_R04C02: num  0.665 0.671 0.399 0.788 0.54 ...
##  $ 7796806038_R03C02: num  0.546 0.869 0.197 0.858 0.463 ...
##  $ 6042316035_R01C01: num  0.555 0.879 0.277 0.855 0.468 ...
##  $ 6042316035_R05C02: num  0.618 0.81 0.194 0.786 0.58 ...
##  $ 6042316036_R04C02: num  0.592 0.817 0.275 0.865 0.479 ...
##  $ 6042316036_R05C02: num  0.634 0.774 0.249 0.763 0.512 ...
##  $ 6042316048_R03C01: num  0.677 0.788 0.207 0.882 0.47 ...
##  $ 6042316050_R01C02: num  0.627 0.802 0.199 0.835 0.439 ...
##  $ 6042316050_R04C02: num  0.518 0.768 0.261 0.812 0.542 ...
##  $ 6042316050_R05C02: num  0.534 0.831 0.3 0.752 0.411 ...
##  $ 6042316050_R06C01: num  0.658 0.816 0.425 0.797 0.489 ...
##  $ 6042316053_R01C01: num  0.56 0.799 0.251 0.825 0.499 ...
##  $ 6042316054_R03C01: num  0.519 0.827 0.235 0.825 0.409 ...
##  $ 6042316061_R03C02: num  0.636 0.797 0.207 0.83 0.531 ...
##  $ 6042316066_R01C01: num  0.514 0.782 0.207 0.82 0.442 ...
##  $ 6042316066_R06C02: num  0.575 0.812 0.385 0.825 0.583 ...
##  $ 6042316085_R01C01: num  0.502 0.565 0.49 0.878 0.537 ...
##  $ 6042316110_R06C01: num  0.69 0.79 0.223 0.845 0.532 ...
##  $ 6042316121_R03C01: num  0.606 0.787 0.282 0.878 0.514 ...
##  $ 6042316121_R05C01: num  0.583 0.817 0.242 0.792 0.454 ...
##  $ 6042316127_R01C02: num  0.636 0.663 0.463 0.806 0.571 ...
##  $ 6042316127_R04C01: num  0.64 0.801 0.251 0.781 0.461 ...
##  $ 6057825008_R02C02: num  0.604 0.663 0.309 0.863 0.532 ...
##  $ 6057825008_R04C01: num  0.586 0.811 0.26 0.859 0.48 ...
##  $ 6057825017_R04C02: num  0.6 0.937 0.115 0.867 0.479 ...
##  $ 6057825017_R05C02: num  0.662 0.758 0.176 0.807 0.547 ...
##  $ 6057825018_R03C02: num  0.519 0.771 0.408 0.82 0.48 ...
##  $ 6057825018_R04C02: num  0.557 0.682 0.253 0.858 0.544 ...
##  $ 6929718123_R03C02: num  0.631 0.83 0.228 0.846 0.514 ...
##  $ 6929718136_R02C02: num  0.61 0.832 0.25 0.819 0.455 ...
##  $ 6929718136_R03C01: num  0.621 0.827 0.306 0.852 0.526 ...
##  $ 6929718138_R04C01: num  0.597 0.827 0.29 0.864 0.532 ...
##  $ 6929718138_R06C02: num  0.603 0.805 0.322 0.863 0.529 ...
##  $ 6969568082_R06C01: num  0.542 0.805 0.288 0.827 0.393 ...
##  $ 6969568084_R04C01: num  0.586 0.779 0.249 0.825 0.543 ...
##  $ 6969568087_R06C02: num  0.553 0.77 0.272 0.799 0.432 ...
##  $ 6969568118_R02C01: num  0.556 0.806 0.183 0.834 0.554 ...
##  $ 6969568118_R03C01: num  0.555 0.757 0.22 0.81 0.466 ...
##  $ 7786923063_R03C01: num  0.693 0.699 0.215 0.78 0.483 ...
##  $ 7786923107_R05C01: num  0.505 0.765 0.161 0.818 0.469 ...
##  $ 7786923107_R06C01: num  0.628 0.787 0.262 0.79 0.54 ...
##  $ 7796806022_R05C02: num  0.5 0.871 0.216 0.852 0.481 ...
##  $ 6042316035_R03C02: num  0.543 0.895 0.293 0.846 0.547 ...
##  $ 6042316035_R04C01: num  0.6 0.78 0.222 0.82 0.488 ...
##  $ 6042316050_R02C02: num  0.772 0.757 0.289 0.864 0.527 ...
##  $ 6042316053_R02C02: num  0.607 0.845 0.244 0.855 0.524 ...
##  $ 6042316053_R03C02: num  0.578 0.823 0.274 0.845 0.412 ...
##  $ 6042316063_R03C01: num  0.591 0.797 0.27 0.862 0.527 ...
##  $ 6042316063_R03C02: num  0.686 0.798 0.201 0.842 0.542 ...
##  $ 6042316069_R02C01: num  0.616 0.772 0.255 0.796 0.417 ...
##  $ 6042316085_R05C01: num  0.543 0.754 0.27 0.771 0.474 ...
##  $ 6042316094_R06C01: num  0.644 0.748 0.233 0.816 0.543 ...
##  $ 6042316099_R02C01: num  0.538 0.811 0.25 0.852 0.488 ...
##  $ 6042316107_R03C01: num  0.62 0.791 0.404 0.805 0.539 ...
##  $ 6042316113_R06C01: num  0.613 0.739 0.252 0.789 0.52 ...
##  $ 6057825014_R03C01: num  0.611 0.742 0.376 0.861 0.473 ...
##  $ 6057825018_R05C02: num  0.64 0.749 0.183 0.83 0.546 ...
##  $ 6929718136_R05C01: num  0.586 0.841 0.17 0.83 0.512 ...
##  $ 6929718138_R01C01: num  0.495 0.832 0.266 0.833 0.467 ...
##  $ 6969568084_R05C02: num  0.61 0.82 0.227 0.805 0.507 ...
##  $ 6969568087_R02C01: num  0.731 0.815 0.239 0.826 0.519 ...
##  $ 6969568118_R06C02: num  0.68 0.805 0.254 0.842 0.531 ...
##  $ 6042316035_R02C02: num  0.603 0.894 0.29 0.84 0.513 ...
##  $ 6042316048_R02C01: num  0.576 0.83 0.23 0.798 0.561 ...
##  $ 6042316050_R06C02: num  0.605 0.805 0.214 0.818 0.503 ...
##  $ 6042316053_R06C02: num  0.628 0.794 0.236 0.81 0.497 ...
##  $ 6042316054_R02C02: num  0.619 0.827 0.403 0.832 0.502 ...
##  $ 6042316054_R04C02: num  0.478 0.82 0.19 0.859 0.419 ...
##  $ 6042316061_R03C01: num  0.481 0.811 0.213 0.821 0.488 ...
##  $ 6042316061_R04C01: num  0.61 0.817 0.248 0.837 0.572 ...
##  $ 6042316061_R05C01: num  0.678 0.82 0.269 0.825 0.526 ...
##  $ 6042316063_R02C01: num  0.577 0.797 0.308 0.805 0.451 ...
##  $ 6042316063_R04C02: num  0.596 0.808 0.224 0.847 0.496 ...
##  $ 6042316065_R06C02: num  0.667 0.741 0.311 0.85 0.456 ...
##  $ 6042316066_R02C02: num  0.646 0.811 0.202 0.829 0.516 ...
##   [list output truncated]
```

```r
#Creating beta and M value dataframes
B.norm <- GSE43414_batch_cor_sorted_by_braakmat
B.normdf <- as.data.frame(B.norm)

M.norm <- logit2(GSE43414_batch_cor_sorted_by_braakmat)
M.normdf <- as.data.frame(M.norm)

#Memory Cleaning - Necessary when data is large
rm(GSE43414_batch_cor)
rm(transpose_GSE43414_batch_cor)
rm(GSE43414_batch_cor_sorted_by_braak)
rm(GSE43414_batch_cor_sorted_by_braakmat)
rm(matches_GSE43414_batch_cor)
gc()
```

```
##             used   (Mb) gc trigger    (Mb)   max used   (Mb)
## Ncells   9035506  482.6   14442815   771.4    9431000  503.7
## Vcells 525711895 4010.9 1324813031 10107.6 1262370963 9631.2
```


```r
# Pre-DMRcate Setup
# create the design for the model
design <- model.matrix(~meta_order_by_braakdf$broad_regions+meta_order_by_braakdf$braak.stage+meta_order_by_braakdf$Sex+meta_order_by_braakdf$age.brain, row.names=T)
View(design)
```

```
## Error in .External2(C_dataviewer, x, title): unable to start data viewer
```

```r
# fit the linear model
fit <- lmFit(M.norm, design)

# fit the contrasts
fit2 <- eBayes(fit)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))
```

```
##    (Intercept) meta_order_by_braakdf$broad_regionscortex
## -1      150888                                     53038
## 0        27765                                    225205
## 1       159706                                     60116
##    meta_order_by_braakdf$braak.stage meta_order_by_braakdf$SexMALE
## -1                               504                          4437
## 0                             337330                        290621
## 1                                525                         43301
##    meta_order_by_braakdf$age.brain
## -1                             777
## 0                           334903
## 1                             2679
```

```r
# get the table of results for the first contrast
ann450ksub <- ann450k[match(rownames(M.norm),ann450k$Name),
c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=2, genelist=ann450ksub)
str(DMPs)
```

```
## 'data.frame':	338359 obs. of  28 variables:
##  $ chr                     : chr  "chr7" "chr1" "chr7" "chr11" ...
##  $ pos                     : int  155166290 246939249 155244149 493208 99736304 64711268 39666961 122535518 259396 44199597 ...
##  $ strand                  : chr  "+" "-" "-" "+" ...
##  $ Name                    : chr  "cg12820681" "cg06826710" "cg20987938" "cg13015350" ...
##  $ Probe_rs                : chr  NA "rs6426334" NA NA ...
##  $ Probe_maf               : num  NA 0.468 NA NA NA ...
##  $ CpG_rs                  : chr  NA NA NA NA ...
##  $ CpG_maf                 : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ SBE_rs                  : chr  NA NA NA NA ...
##  $ SBE_maf                 : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ Islands_Name            : chr  "chr7:155164557-155167854" "" "chr7:155244046-155244357" "chr11:494853-495156" ...
##  $ Relation_to_Island      : chr  "Island" "OpenSea" "Island" "N_Shore" ...
##  $ UCSC_RefGene_Name       : chr  "" "" "" "" ...
##  $ UCSC_RefGene_Accession  : chr  "" "" "" "" ...
##  $ UCSC_RefGene_Group      : chr  "" "" "" "" ...
##  $ Phantom                 : chr  "" "" "" "" ...
##  $ DMR                     : chr  "DMR" "" "" "" ...
##  $ Enhancer                : chr  "TRUE" "" "" "" ...
##  $ HMM_Island              : chr  "" "1:245005810-245005900" "7:154934147-154937118" "11:482969-483211" ...
##  $ Regulatory_Feature_Name : chr  "" "" "" "11:493116-493264" ...
##  $ Regulatory_Feature_Group: chr  "" "" "" "Unclassified_Cell_type_specific" ...
##  $ DHS                     : chr  "" "" "" "" ...
##  $ logFC                   : num  -0.568 0.756 -0.583 0.622 0.565 ...
##  $ AveExpr                 : num  -0.768 1.093 -0.816 0.976 0.949 ...
##  $ t                       : num  -10.6 10.4 -10.3 10.2 10.1 ...
##  $ P.Value                 : num  3.16e-22 1.77e-21 2.64e-21 7.39e-21 1.06e-20 ...
##  $ adj.P.Val               : num  1.07e-16 2.98e-16 2.98e-16 6.25e-16 7.16e-16 ...
##  $ B                       : num  39.7 38 37.6 36.6 36.3 ...
```

```r
# Saving the data.frame file as an Excel
write.table(DMPs, file="batch-DMPs.csv", sep=",", row.names=FALSE)

# As a quick check, plot the top 4 most significantly differentially methylated CpGs, Tissue
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
plotCpg(B.norm, cpg=cpg, pheno=meta_order_by_braakdf$broad_regions, type = "categorical", measure = "beta")
})
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

```
## $cg12820681
## NULL
## 
## $cg06826710
## NULL
## 
## $cg20987938
## NULL
## 
## $cg13015350
## NULL
```

```r
# As a quick check, plot the top 4 most significantly differentially methylated CpGs, Braak
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
plotCpg(B.norm, cpg=cpg, pheno=meta_order_by_braakdf$braak.stage, type = "continuous", measure = "beta")
})
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png)

```
## $cg12820681
## NULL
## 
## $cg06826710
## NULL
## 
## $cg20987938
## NULL
## 
## $cg13015350
## NULL
```

Differential methylation analysis of regions

Use dmrcate function to combine individual CpG statistics to identify differentially methylated regions. DMRs$results contains all of the regions found, with genomic annotations and p-values


```r
#DMRcate Contrast Copy
myannotation <- cpg.annotate(object = M.norm, datatype = "array", what = "M", arraytype = "450K", analysis.type="differential", design=design, coef=2) #113034 individually significant probes
```

```
## Your contrast returned 113154 individually significant probes. We recommend the default setting of pcutoff in dmrcate().
```

```r
str(myannotation)
```

```
## List of 6
##  $ ID    : Factor w/ 338359 levels "cg00000029","cg00000108",..: 180947 301360 201815 15178 158126 291516 2684 82058 283599 230470 ...
##  $ stat  : num [1:338359] 0.238 2.761 -2.719 -0.479 -0.217 ...
##  $ CHR   : Factor w/ 22 levels "chr1","chr10",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ pos   : int [1:338359] 15865 534242 710097 714177 758829 763119 790667 805102 805554 812539 ...
##  $ betafc: num [1:338359] 0.0011 0.013002 -0.014483 -0.000256 -0.001984 ...
##  $ indfdr: num [1:338359] 0.8968 0.0213 0.0237 0.7785 0.9065 ...
##  - attr(*, "row.names")= int [1:338359] 306465 97481 99227 274688 309116 209191 109491 69333 119238 29190 ...
##  - attr(*, "class")= chr "annot"
```

```r
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, pcutoff = 0.05)
```

```
## Fitting chr1...
```

```
## Fitting chr10...
```

```
## Fitting chr11...
```

```
## Fitting chr12...
```

```
## Fitting chr13...
```

```
## Fitting chr14...
```

```
## Fitting chr15...
```

```
## Fitting chr16...
```

```
## Fitting chr17...
```

```
## Fitting chr18...
```

```
## Fitting chr19...
```

```
## Fitting chr2...
```

```
## Fitting chr20...
```

```
## Fitting chr21...
```

```
## Fitting chr22...
```

```
## Fitting chr3...
```

```
## Fitting chr4...
```

```
## Fitting chr5...
```

```
## Fitting chr6...
```

```
## Fitting chr7...
```

```
## Fitting chr8...
```

```
## Fitting chr9...
```

```
## Demarcating regions...
```

```
## Done!
```

```r
#head(dmrcoutput$results)
#str(dmrcoutput$results)

# convert the regions to annotated genomic ranges
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
```

### Individual CpG Table Underlying DMRs


```r
load("Priest.RData")
### Individual CpG Table Underlying DMRs
fdat <- fData(Priest)
fdat$CHR <- as.character(fdat$CHR)

# Pulling out CpGs based on positions from ranges object.
sig.cpgs <- lapply(1:length(results.ranges), function(x){
  print(x)
  coords <- names(ranges(results.ranges))[x]
  chr <- sub(":.*", "", coords) 
  chr <- sub("chr", "", chr)
  bookends <- sub(".*:", "", coords)
  startcpg <- as.integer(sub("-.*", "", bookends))
  stopcpg <- as.integer(sub(".*-", "", bookends))
  cpgs <- rownames(fdat)[fdat$CHR %in% chr & fdat$MAPINFO >= 
                             startcpg & fdat$MAPINFO <= stopcpg]
})
sig.cpgs.list <- unlist(sig.cpgs)
```


```r
#Create a continuous variable from the categorical variable "brain region" using the "Cortex" value from the design matrix of 1 for cortex, 0 for cerebellum
colnames(design) <- c("Intercept", "Cortex", "Braak.stage", "Male", "Age")
design <- as.data.frame(design)

#Calculate delta betas using a linear model for brain region
delbeta<-sapply(1:nrow(B.norm), function(x) {
  z<-lm(unlist(B.norm[x,]) ~meta_order_by_braakdf$broad_regions+meta_order_by_braakdf$braak.stage+meta_order_by_braakdf$Sex+meta_order_by_braakdf$age.brain)
  intercept=z$coefficients[1]
  slope=z$coefficients[2]
  as.numeric((slope*max(design$Cortex, na.rm=T)+intercept)-intercept)
})
save(delbeta, file="GSE43414_batch_cor_delbeta.RData")
#load("GSE43414_batch_cor_delbeta.RData")
```


```r
#Making a dataframe based on CpGs from ranges object and attaching map info and DMRcate statistics. 
delbeta.ann <- cbind(rownames(B.norm),delbeta)
colnames(delbeta.ann) <- c("probe", "DB")
delbeta.ann <- as.data.frame(delbeta.ann)
delbeta.ann$DB <- as.numeric(as.character(delbeta.ann$DB))
delbeta.ann <- delbeta.ann[delbeta.ann$probe %in% sig.cpgs.list,]

dmrinput <- dmrcoutput$input
sub_dmrinput <- dmrinput[dmrinput$ID %in% sig.cpgs.list,]
sub_fdat <- fdat[fdat$ILMNID %in% sig.cpgs.list,]
head(sub_fdat)
```

```
##            MAPINF0-1 MAPINFO+1 Probe_start Probe_end SNPCpG n_SNPCpG
## cg00000029  53468111  53468113    53468112  53468162              NA
## cg00000108  37459205  37459207    37459206  37459256              NA
## cg00000165  91194673  91194675    91194624  91194674              NA
## cg00000292  28890099  28890101    28890100  28890150              NA
## cg00000321  41167801  41167803    41167752  41167802              NA
## cg00000363 230560792 230560794   230560793 230560843              NA
##              SNPprobe n_SNPprobe HIL_CpG_class
## cg00000029                    NA            HC
## cg00000108  rs9857774          1            LC
## cg00000165 rs76771611          1       ICshore
## cg00000292 rs62037371          1            IC
## cg00000321                    NA       ICshore
## cg00000363                    NA       ICshore
##                           HIL_CpG_Island_Name n_bp_repetitive AlleleA_Hits
## cg00000029  chr16_HCshore:53467967-53469412;.               0            1
## cg00000108                                .;.              48            1
## cg00000165   .;chr1_ICshore:91194238-91195462               0            1
## cg00000292       .;chr16_IC:28889913-28890373               0            1
## cg00000321   .;chr8_ICshore:41165660-41169277               0            1
## cg00000363 .;chr1_ICshore:230560253-230563060               0            1
##            AlleleB_Hits XY_Hits Autosomal_Hits Closest_TSS Closest_TSS_1
## cg00000029            0   XY_NO           A_NO    53468350      53468351
## cg00000108            0   XY_NO           A_NO    37458757      37458758
## cg00000165            0   XY_NO           A_NO    91182793      91182794
## cg00000292            0   XY_NO           A_NO    28889808      28889809
## cg00000321            0   XY_NO           A_NO    41166989      41166990
## cg00000363            0   XY_NO           A_NO   230513390     230513391
##            Distance_closest_TSS Closest_TSS_gene_name
## cg00000029                 -238                  RBL2
## cg00000108                  449               C3orf35
## cg00000165               -11880                BARHL2
## cg00000292                  292                ATP2A1
## cg00000321                 -812                 SFRP1
## cg00000363               -47402                 PGBD5
##            Closest_TSS_Transcript     ILMNID       NAME ADDRESSA_ID
## cg00000029              NM_005611 cg00000029 cg00000029    14782418
## cg00000108              CCDS46792 cg00000108 cg00000108    12709357
## cg00000165              NM_020063 cg00000165 cg00000165    12637463
## cg00000292              NM_004320 cg00000292 cg00000292    43764508
## cg00000321              NM_003012 cg00000321 cg00000321    62789509
## cg00000363              NM_024554 cg00000363 cg00000363    16661505
##                                              ALLELEA_PROBESEQ ADDRESSB_ID
## cg00000029 AACTATACTAACRAAAAAATATCCAAAAAACACTAACRTATAAAAATTTC          NA
## cg00000108 ATACAATAAAACAAACCTAAAATAATCCTAACTCCRCTATCATCCTAACC          NA
## cg00000165 CAAAATCTATTAATACAATAACTTTTAATAAAACAACTAAAACACACATC          NA
## cg00000292 AAAACATTAATTACCAACCRCTCTTCCAAAAAACACTTACCATTAAAACC          NA
## cg00000321 ATAAATACCCAATAAACCTAACTAAACTCCCTAAAAAACRAAACRAAAAC          NA
## cg00000363 RTCTTAACTTAACTTAATTTTCTCCTTAATCTAAAAAACTTTCCCTATCC          NA
##            ALLELEB_PROBESEQ INFINIUM_DESIGN_TYPE NEXT_BASE COLOR_CHANNEL
## cg00000029                                    II                        
## cg00000108                                    II                        
## cg00000165                                    II                        
## cg00000292                                    II                        
## cg00000321                                    II                        
## cg00000363                                    II                        
##                                                                                                                        FORWARD_SEQUENCE
## cg00000029 TTTTTTAGATAAGGATATCCAGGCGATGAGGAAGTTTTACTTCTGGGAACAGCCTGGATA[CG]AAACCTTCACACGTCAGTGTCTTTTGGACATTTTCTCGTCAGTACAGCCCTGTTGAATGT
## cg00000108 TCCATTTTGAAGGAAAAAAATGAAGGCTCTGAAAGTGTAAATCGCTTACTGAAGGGCACA[CG]GCCAGGATGACAGCGGAGCCAGGATCACCCCAGGTCTGTCTCATTGCATATGTCATGGCT
## cg00000165 CTAAGTGCAGTCAGGATCTGTTAGTACAGTGGCTTTTGATGGAACAGCTGAGGCACACAT[CG]CCCGTGGCATGGACTCCGGGGCCGAACGCTCACGACCAAGACTTTTGCCCTTTTGAAATG
## cg00000292 TGGGGTGAGTGAGACCACGGGCCTCACCCCGGACCAAGTTAAGCGGAATCTGGAGAAATA[CG]GCCTCAATGGTAAGTGTCCCTTGGAAGAGCGGCTGGTAATTAATGCCCTCCTGCACCCCC
## cg00000321 GAGGTCTGCTTGTAAATACCCAGTGGGCCTGGCTGGGCTCCCTGGAAGGCGAGGCGAAGG[CG]CAGTTGGAGCTGTTTGCTGTGAGCAGCACCTCTCCAGGTGGGGCCGCCCATGGTGGCCCT
## cg00000363 CTGCCCAATCGGTCCCTTCCTTCACTCCTCCCCATTCTTAACAAGAGATCTGGAATGGCG[CG]GACAGGGAAAGTTTCTCAGATTAAGGAGAAAACTAAGCCAAGTCAAGACGCCGCGGGTGG
##            GENOME_BUILD CHR   MAPINFO
## cg00000029           37  16  53468112
## cg00000108           37   3  37459206
## cg00000165           37   1  91194674
## cg00000292           37  16  28890100
## cg00000321           37   8  41167802
## cg00000363           37   1 230560793
##                                                     SOURCESEQ
## cg00000029 GCTGTACTGACGAGAAAATGTCCAAAAGACACTGACGTGTGAAGGTTTCG
## cg00000108 CGGCCAGGATGACAGCGGAGCCAGGATCACCCCAGGTCTGTCTCATTGCA
## cg00000165 AGGATCTGTTAGTACAGTGGCTTTTGATGGAACAGCTGAGGCACACATCG
## cg00000292 CGGCCTCAATGGTAAGTGTCCCTTGGAAGAGCGGCTGGTAATTAATGCCC
## cg00000321 CGCCTTCGCCTCGCCTTCCAGGGAGCCCAGCCAGGCCCACTGGGTATTTA
## cg00000363 TCTTGACTTGGCTTAGTTTTCTCCTTAATCTGAGAAACTTTCCCTGTCCG
##            CHROMOSOME_36 COORDINATE_36 STRAND PROBE_SNPS PROBE_SNPS_10
## cg00000029            16      52025613      F                         
## cg00000108             3      37434210      F  rs9857774              
## cg00000165             1      90967262      R                         
## cg00000292            16      28797601      F rs62037371              
## cg00000321             8      41286959      R                         
## cg00000363             1     228627416      F                         
##            RANDOM_LOCI METHYL27_LOCI UCSC_REFGENE_NAME
## cg00000029          NA            NA              RBL2
## cg00000108          NA            NA   C3orf35;C3orf35
## cg00000165          NA            NA                  
## cg00000292          NA          TRUE     ATP2A1;ATP2A1
## cg00000321          NA            NA             SFRP1
## cg00000363          NA            NA                  
##            UCSC_REFGENE_ACCESSION UCSC_REFGENE_GROUP
## cg00000029              NM_005611            TSS1500
## cg00000108    NM_178339;NM_178342         Body;3'UTR
## cg00000165                                          
## cg00000292    NM_004320;NM_173201    1stExon;1stExon
## cg00000321              NM_003012            TSS1500
## cg00000363                                          
##               UCSC_CPG_ISLANDS_NAME RELATION_TO_UCSC_CPG_ISLAND PHANTOM
## cg00000029  chr16:53468284-53469209                     N_Shore        
## cg00000108                                                             
## cg00000165   chr1:91190489-91192804                     S_Shore        
## cg00000292  chr16:28890954-28891868                     N_Shore        
## cg00000321   chr8:41165852-41167140                     S_Shore        
## cg00000363 chr1:230561103-230562702                     N_Shore        
##             DMR ENHANCER            HMM_ISLAND REGULATORY_FEATURE_NAME
## cg00000029            NA                          16:53467838-53469685
## cg00000108            NA                                              
## cg00000165 CDMR     TRUE   1:90967262-90967361                        
## cg00000292            NA                                              
## cg00000321            NA                                              
## cg00000363            NA 1:228627033-228629325                        
##            REGULATORY_FEATURE_GROUP  DHS Index   TargetID ProbeID_A
## cg00000029      Promoter_Associated TRUE     1 cg00000029  14782418
## cg00000108                            NA     2 cg00000108  12709357
## cg00000165                            NA     4 cg00000165  12637463
## cg00000292                            NA     7 cg00000292  43764508
## cg00000321                            NA     8 cg00000321  62789509
## cg00000363                            NA     9 cg00000363  16661505
##            ProbeID_B
## cg00000029  14782418
## cg00000108  12709357
## cg00000165  12637463
## cg00000292  43764508
## cg00000321  62789509
## cg00000363  16661505
```

```r
# Sanity check.
sub_fdat <- sub_fdat[order(sub_fdat$ILMNID),]
dim(sub_fdat)
```

```
## [1] 185494     57
```

```r
sub_dmrinput <- sub_dmrinput[order(sub_dmrinput$ID),]
dim(sub_dmrinput)
```

```
## [1] 158394     10
```

```r
#The dmr input list has fewer rows than sub_fdat, likely due to the probes containing NAs that were removed during filtering.

#Subset fdat probe list
sub_fdat <- sub_fdat[sub_fdat$ILMNID %in% sub_dmrinput$ID,]
identical(nrow(sub_fdat), nrow(sub_dmrinput)) #TRUE
```

```
## [1] TRUE
```

```r
dmrcate.in <- as.vector(sub_dmrinput[["ID"]])
fdat.in <- as.vector(sub_fdat[["ILMNID"]])
fdat.in <- fdat.in[fdat.in %in% dmrcate.in]
identical(dmrcate.in, fdat.in) 
```

```
## [1] FALSE
```

```r
sig.cpgs.fc <- data.frame(cpg=sub_dmrinput$ID, 
                          FC = sub_dmrinput$betafc,
                          FDR = sub_dmrinput$fdr,
                          Gene = sub_fdat$UCSC_REFGENE_NAME,
                          pvalue = sub_dmrinput$raw,
                          chr = sub_fdat$CHR,
                          pos = sub_fdat$MAPINFO,
                          DB = delbeta.ann$DB)
sig.cpgs.fc$chr <- gsub("^", "chr", sig.cpgs.fc$chr)
```


```r
# Assigning the correct DMR to each CpG (long running code)
sig.cpg.dmr <- lapply(1:nrow(sig.cpgs.fc), function(x){
  print(x)
  a<- which(sub(":.*", "", names(ranges(results.ranges))) == sig.cpgs.fc$chr[x] & 
              start(ranges(results.ranges)) <= sig.cpgs.fc$pos[x] & 
              end(ranges(results.ranges)) >= sig.cpgs.fc$pos[x])
  names(ranges(results.ranges))[a]
})
sig.cpgs.fc$dmr <- unlist(sig.cpg.dmr)
ELses.sig.cpgs.fc <- sig.cpgs.fc
write.table(ELses.sig.cpgs.fc, "significant.DMPs_batch.txt",sep="\t")
```

### Table Containing Annotation for DMRs
We will now produce a table containing information for the DMRs.

```r
# This code adds the mean, max and min delta beta values for each DMR.
head(ELses.sig.cpgs.fc)
```

```
##          cpg           FC          FDR                    Gene
## 1 cg00000029 -0.025958458 8.416575e-07                    RBL2
## 2 cg00000108  0.064001824 5.485042e-12         C3orf35;C3orf35
## 3 cg00000165 -0.036568643 7.416848e-16                        
## 4 cg00000292  0.009479659 4.486280e-03           ATP2A1;ATP2A1
## 5 cg00000363  0.011882490 1.365394e-02                        
## 6 cg00000622 -0.000845493 3.720137e-02 NIPA2;NIPA2;NIPA2;NIPA2
##         pvalue   chr       pos           DB                      dmr
## 1 1.826797e-07 chr16  53468112 -0.025958458  chr16:53467612-53469344
## 2 5.757560e-13  chr3  37459206  0.064001824   chr3:37458845-37459206
## 3 5.106496e-17  chr1  91194674 -0.036568643   chr1:91194674-91196488
## 4 1.943214e-03 chr16  28890100  0.009479659  chr16:28887830-28890100
## 5 6.658272e-03  chr1 230560793  0.011882490 chr1:230560318-230560793
## 6 2.034757e-02 chr15  23034447 -0.000845493  chr15:23034186-23034682
```

```r
mean_max_DB <- data.frame(DB=ELses.sig.cpgs.fc$DB, 
                          dmr = ELses.sig.cpgs.fc$dmr)
test <- aggregate(.~dmr, data=mean_max_DB, mean)
colnames(test) <- c("dmr", "mean_DB")
test.max <- aggregate(.~dmr, data=mean_max_DB, max)
colnames(test.max) <- c("dmr", "max_DB")
test.min <- aggregate(.~dmr, data=mean_max_DB, min)
colnames(test.min) <- c("dmr", "min_DB")
DB_info <- merge(test, test.max, by = "dmr")
DB_info <- merge(DB_info, test.min, by = "dmr")
# Making a second table with DMR info.
results <- dmrcoutput$results
#getting the mean for the DMR - you could do other manipulations here as well 
dmr.mean.fc <- lapply(1:length(unique(sig.cpgs.fc$dmr)), function(x) {
  mean(sig.cpgs.fc$FC[which(sig.cpgs.fc$dmr == unique(sig.cpgs.fc$dmr)[x])])
})
dmr.fc <- data.frame(dmr = unique(sig.cpgs.fc$dmr),
                     meanbetafc = unlist(dmr.mean.fc))
dmr.fc <- dmr.fc[order(match(dmr.fc$dmr, results$coord)),]
identical(as.character(dmr.fc$dmr), results$coord)
```

```
## [1] TRUE
```

```r
dmr.fc$dmr <- as.character(dmr.fc$dmr)
results <- merge(results, dmr.fc, by.x = "coord", by.y = "dmr") 
results$chr <- sub(":.*", "", results$coord)
bookends <- sub(".*:", "", results$coord)
results$start <- as.integer(sub("-.*", "", bookends))
results$end <- as.integer(sub(".*-", "", bookends))
results$width <- (results$end - results$start) + 1
results <- merge(results, DB_info, by.x = "coord",by.y = "dmr")
ELses_results <- results  
write.table(ELses_results, "significant.DMRs_batch.txt",sep="\t")

par(mfrow=c(1,2))
hist(ELses_results$max_DB, main="Max DB")
hist(ELses_results$mean_DB, main="Mean DB")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

We can see that the distribution of delta betas is centered around 0; this is more prominent for the mean DB score than the max DB score. Setting a DB cutoff of 0.05 will capture those that appear to have a real change between cortex and cerebellum.

# Look at top hits

```r
#Filter by mean DB >= 0.05
hits <- ELses_results %>% filter(mean_DB >= 0.05)
nrow(hits)
```

```
## [1] 1082
```

```r
#Arrange by mean DB
hits <- hits %>% arrange(desc(mean_DB))
hits$coord[1:2]
```

```
## [1] "chr11:674623-674687" "chr7:891044-891145"
```
Using the UCSC Genome Browser (hg19), the top two hits are located in DEAF1 and in an intergenic region upstream of MGMT.

#Volcano plot

```r
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# create the design for the model
design <- model.matrix(~meta_order_by_braakdf$broad_regions+meta_order_by_braakdf$braak.stage+meta_order_by_braakdf$Sex+meta_order_by_braakdf$age.brain, row.names=T)
fit <- lmFit(M.norm, design)
fit2 <- eBayes(fit)
summary(decideTests(fit2))
```

```
##    (Intercept) meta_order_by_braakdf$broad_regionscortex
## -1      150888                                     53038
## 0        27765                                    225205
## 1       159706                                     60116
##    meta_order_by_braakdf$braak.stage meta_order_by_braakdf$SexMALE
## -1                               504                          4437
## 0                             337330                        290621
## 1                                525                         43301
##    meta_order_by_braakdf$age.brain
## -1                             777
## 0                           334903
## 1                             2679
```

```r
ann450ksub <- ann450k[match(rownames(M.norm),ann450k$Name),
c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=2, genelist=ann450ksub)
str(DMPs)
```

```
## 'data.frame':	338359 obs. of  28 variables:
##  $ chr                     : chr  "chr7" "chr1" "chr7" "chr11" ...
##  $ pos                     : int  155166290 246939249 155244149 493208 99736304 64711268 39666961 122535518 259396 44199597 ...
##  $ strand                  : chr  "+" "-" "-" "+" ...
##  $ Name                    : chr  "cg12820681" "cg06826710" "cg20987938" "cg13015350" ...
##  $ Probe_rs                : chr  NA "rs6426334" NA NA ...
##  $ Probe_maf               : num  NA 0.468 NA NA NA ...
##  $ CpG_rs                  : chr  NA NA NA NA ...
##  $ CpG_maf                 : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ SBE_rs                  : chr  NA NA NA NA ...
##  $ SBE_maf                 : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ Islands_Name            : chr  "chr7:155164557-155167854" "" "chr7:155244046-155244357" "chr11:494853-495156" ...
##  $ Relation_to_Island      : chr  "Island" "OpenSea" "Island" "N_Shore" ...
##  $ UCSC_RefGene_Name       : chr  "" "" "" "" ...
##  $ UCSC_RefGene_Accession  : chr  "" "" "" "" ...
##  $ UCSC_RefGene_Group      : chr  "" "" "" "" ...
##  $ Phantom                 : chr  "" "" "" "" ...
##  $ DMR                     : chr  "DMR" "" "" "" ...
##  $ Enhancer                : chr  "TRUE" "" "" "" ...
##  $ HMM_Island              : chr  "" "1:245005810-245005900" "7:154934147-154937118" "11:482969-483211" ...
##  $ Regulatory_Feature_Name : chr  "" "" "" "11:493116-493264" ...
##  $ Regulatory_Feature_Group: chr  "" "" "" "Unclassified_Cell_type_specific" ...
##  $ DHS                     : chr  "" "" "" "" ...
##  $ logFC                   : num  -0.568 0.756 -0.583 0.622 0.565 ...
##  $ AveExpr                 : num  -0.768 1.093 -0.816 0.976 0.949 ...
##  $ t                       : num  -10.6 10.4 -10.3 10.2 10.1 ...
##  $ P.Value                 : num  3.16e-22 1.77e-21 2.64e-21 7.39e-21 1.06e-20 ...
##  $ adj.P.Val               : num  1.07e-16 2.98e-16 2.98e-16 6.25e-16 7.16e-16 ...
##  $ B                       : num  39.7 38 37.6 36.6 36.3 ...
```

```r
#Add DBs to topTable result
delbeta.ann2 <- cbind(rownames(B.norm),delbeta)
colnames(delbeta.ann2) <- c("probe", "DB")
delbeta.ann2 <- as.data.frame(delbeta.ann2)
delbeta.ann2$DB <- as.numeric(as.character(delbeta.ann2$DB))
DMPs <- DMPs %>% arrange(Name)
DMPs$DB <- delbeta.ann2$DB
gene_list <- DMPs

require(ggplot2)
##Highlight genes that have a DB >= 0.05 and a FDR <= 0.05
gene_list$threshold = as.factor(abs(gene_list$DB) >= 0.05 & gene_list$adj.P.Val <= 0.05)

##Construct the volcano plot
ggplot(data=gene_list, aes(x=DB, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  labs(legend.position = "none") +
  xlab("Delta Beta") + ylab("-log10 Adj P-Value")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)



