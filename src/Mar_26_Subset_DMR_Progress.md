Subset DMR
================

``` r
#Memory Cleaning - Necessary when data is large
#rm(GSE43414_cell_cor)
#rm(transpose_GSE43414_cell_cor)
#rm(GSE43414_cell_cor_sorted_by_brain_regions)
#rm(meta)
#rm(brain_regions)
#memory.limit()
#memory.size()
#gc()
```

``` r
source("http://bioconductor.org/biocLite.R")
```

    ## Bioconductor version 3.4 (BiocInstaller 1.24.0), ?biocLite for help

``` r
biocLite("wateRmelon")
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.2 (2016-10-31).

    ## Installing package(s) 'wateRmelon'

    ## package 'wateRmelon' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\GiLL\AppData\Local\Temp\RtmpwhjzGA\downloaded_packages

    ## installation path not writeable, unable to update packages: lattice

``` r
#Downloading 450k annotation data
source("https://bioconductor.org/biocLite.R")
```

    ## Bioconductor version 3.4 (BiocInstaller 1.24.0), ?biocLite for help

``` r
biocLite("IlluminaHumanMethylation450k.db")
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.2 (2016-10-31).

    ## Installing package(s) 'IlluminaHumanMethylation450k.db'

    ## installing the source package 'IlluminaHumanMethylation450k.db'

    ## Warning: running command '"C:/PROGRA~1/R/R-33~1.2/bin/x64/
    ## R" CMD INSTALL -l "C:\Users\GiLL\Documents\R\win-library\3.3" C:
    ## \Users\GiLL\AppData\Local\Temp\RtmpwhjzGA/downloaded_packages/
    ## IlluminaHumanMethylation450k.db_2.0.9.tar.gz' had status 1

    ## Warning in install.packages(pkgs = doing, lib = lib, ...): installation of
    ## package 'IlluminaHumanMethylation450k.db' had non-zero exit status

    ## installation path not writeable, unable to update packages: lattice

``` r
source("https://bioconductor.org/biocLite.R")
```

    ## Bioconductor version 3.4 (BiocInstaller 1.24.0), ?biocLite for help

``` r
biocLite("Sushi")
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.2 (2016-10-31).

    ## Installing package(s) 'Sushi'

    ## package 'Sushi' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\GiLL\AppData\Local\Temp\RtmpwhjzGA\downloaded_packages

    ## installation path not writeable, unable to update packages: lattice

``` r
library(GEOquery)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colnames,
    ##     do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, lengths, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff,
    ##     sort, table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
library(wateRmelon)
```

    ## Loading required package: limma

    ## Warning: package 'limma' was built under R version 3.3.3

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    ## Loading required package: matrixStats

    ## Warning: package 'matrixStats' was built under R version 3.3.3

    ## matrixStats v0.51.0 (2016-10-08) successfully loaded. See ?matrixStats for help.

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: methylumi

    ## Loading required package: scales

    ## Loading required package: reshape2

    ## Loading required package: ggplot2

    ## Loading required package: FDb.InfiniumMethylation.hg19

    ## Loading required package: GenomicFeatures

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 3.3.3

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     colMeans, colSums, expand.grid, rowMeans, rowSums

    ## Loading required package: IRanges

    ## Warning: package 'IRanges' was built under R version 3.3.3

    ## Loading required package: GenomeInfoDb

    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 3.3.3

    ## Loading required package: AnnotationDbi

    ## Loading required package: TxDb.Hsapiens.UCSC.hg19.knownGene

    ## Loading required package: org.Hs.eg.db

    ## 

    ## Loading required package: minfi

    ## Loading required package: SummarizedExperiment

    ## 
    ## Attaching package: 'SummarizedExperiment'

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     rowRanges

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## Warning: package 'XVector' was built under R version 3.3.3

    ## Loading required package: bumphunter

    ## Loading required package: foreach

    ## Warning: package 'foreach' was built under R version 3.3.3

    ## Loading required package: iterators

    ## Loading required package: locfit

    ## locfit 1.5-9.1    2013-03-22

    ## Loading required package: lumi

    ## 
    ## Attaching package: 'lumi'

    ## The following objects are masked from 'package:methylumi':
    ## 
    ##     estimateM, getHistory

    ## Loading required package: ROC

    ## Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19

    ## Loading required package: illuminaio

``` r
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
```

    ## 

``` r
library(matrixStats)
library(minfiData)
library(Gviz)
```

    ## Loading required package: grid

``` r
library(DMRcate)
```

    ## Warning: package 'DMRcate' was built under R version 3.3.3

    ## Loading required package: DSS

    ## Loading required package: bsseq

    ## 
    ## Attaching package: 'bsseq'

    ## The following object is masked from 'package:minfi':
    ## 
    ##     getMeth

    ## Loading required package: splines

    ## Loading required package: DMRcatedata

``` r
library(stringr)
library(tidyverse)
```

    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Warning: package 'readr' was built under R version 3.3.3

    ## Conflicts with tidy packages ----------------------------------------------

    ## accumulate(): purrr, foreach
    ## col_factor(): readr, scales
    ## collapse():   dplyr, Biostrings, IRanges
    ## combine():    dplyr, bsseq, lumi, methylumi, minfi, Biobase, BiocGenerics
    ## compact():    purrr, XVector
    ## count():      dplyr, matrixStats
    ## desc():       dplyr, IRanges
    ## discard():    purrr, scales
    ## expand():     tidyr, S4Vectors
    ## filter():     dplyr, stats
    ## first():      dplyr, S4Vectors
    ## lag():        dplyr, stats
    ## Position():   ggplot2, BiocGenerics, base
    ## reduce():     purrr, GenomicRanges, IRanges
    ## regroup():    dplyr, IRanges
    ## rename():     dplyr, S4Vectors
    ## select():     dplyr, AnnotationDbi
    ## simplify():   purrr, IRanges
    ## slice():      dplyr, XVector, IRanges
    ## when():       purrr, foreach

``` r
#Assigning annotation information for CpG probes on the 450k array
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
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

``` r
#Load data and meta files
load("C:/Users/GiLL/Desktop/Mar 22 Methyl/GSE43414_cell_cor.RData", verbose=TRUE) 
```

    ## Loading objects:
    ##   GSE43414_cell_cor

``` r
load("C:/Users/GiLL/Desktop/Mar 22 Methyl/Meta_batch_cor.RData", verbose=TRUE)
```

    ## Loading objects:
    ##   meta

``` r
## transpose data such that probe names are colnames, and rows are patient samples
transpose_GSE43414_cell_cor <- t(GSE43414_cell_cor)
## order metadata by brain region
meta_order_by_brain_regions <- meta %>% arrange(Tissue)
matches_GSE43414_cell_cor <- match(meta_order_by_brain_regions$barcode, rownames(transpose_GSE43414_cell_cor))
GSE43414_cell_cor_sorted_by_brain_regions <- t(transpose_GSE43414_cell_cor[matches_GSE43414_cell_cor,])

#Making tissue names syntactically valid for contrasts later (removes space, adds .)
make.names(meta_order_by_brain_regions$Tissue)
```

    ##   [1] "cerebellum"              "cerebellum"             
    ##   [3] "cerebellum"              "cerebellum"             
    ##   [5] "cerebellum"              "cerebellum"             
    ##   [7] "cerebellum"              "cerebellum"             
    ##   [9] "cerebellum"              "cerebellum"             
    ##  [11] "cerebellum"              "cerebellum"             
    ##  [13] "cerebellum"              "cerebellum"             
    ##  [15] "cerebellum"              "cerebellum"             
    ##  [17] "cerebellum"              "cerebellum"             
    ##  [19] "cerebellum"              "cerebellum"             
    ##  [21] "cerebellum"              "cerebellum"             
    ##  [23] "cerebellum"              "cerebellum"             
    ##  [25] "cerebellum"              "cerebellum"             
    ##  [27] "cerebellum"              "cerebellum"             
    ##  [29] "cerebellum"              "cerebellum"             
    ##  [31] "cerebellum"              "cerebellum"             
    ##  [33] "cerebellum"              "cerebellum"             
    ##  [35] "cerebellum"              "cerebellum"             
    ##  [37] "cerebellum"              "cerebellum"             
    ##  [39] "cerebellum"              "cerebellum"             
    ##  [41] "cerebellum"              "cerebellum"             
    ##  [43] "cerebellum"              "cerebellum"             
    ##  [45] "cerebellum"              "cerebellum"             
    ##  [47] "cerebellum"              "cerebellum"             
    ##  [49] "cerebellum"              "cerebellum"             
    ##  [51] "cerebellum"              "cerebellum"             
    ##  [53] "cerebellum"              "cerebellum"             
    ##  [55] "cerebellum"              "cerebellum"             
    ##  [57] "cerebellum"              "cerebellum"             
    ##  [59] "cerebellum"              "cerebellum"             
    ##  [61] "cerebellum"              "cerebellum"             
    ##  [63] "cerebellum"              "cerebellum"             
    ##  [65] "cerebellum"              "cerebellum"             
    ##  [67] "cerebellum"              "cerebellum"             
    ##  [69] "cerebellum"              "cerebellum"             
    ##  [71] "cerebellum"              "cerebellum"             
    ##  [73] "cerebellum"              "cerebellum"             
    ##  [75] "cerebellum"              "cerebellum"             
    ##  [77] "cerebellum"              "cerebellum"             
    ##  [79] "cerebellum"              "cerebellum"             
    ##  [81] "cerebellum"              "cerebellum"             
    ##  [83] "cerebellum"              "cerebellum"             
    ##  [85] "cerebellum"              "cerebellum"             
    ##  [87] "cerebellum"              "cerebellum"             
    ##  [89] "cerebellum"              "cerebellum"             
    ##  [91] "cerebellum"              "cerebellum"             
    ##  [93] "cerebellum"              "cerebellum"             
    ##  [95] "cerebellum"              "cerebellum"             
    ##  [97] "cerebellum"              "cerebellum"             
    ##  [99] "cerebellum"              "cerebellum"             
    ## [101] "cerebellum"              "cerebellum"             
    ## [103] "cerebellum"              "cerebellum"             
    ## [105] "cerebellum"              "cerebellum"             
    ## [107] "cerebellum"              "cerebellum"             
    ## [109] "cerebellum"              "cerebellum"             
    ## [111] "cerebellum"              "cerebellum"             
    ## [113] "cerebellum"              "cerebellum"             
    ## [115] "cerebellum"              "cerebellum"             
    ## [117] "cerebellum"              "cerebellum"             
    ## [119] "cerebellum"              "cerebellum"             
    ## [121] "cerebellum"              "cerebellum"             
    ## [123] "cerebellum"              "cerebellum"             
    ## [125] "cerebellum"              "cerebellum"             
    ## [127] "cerebellum"              "cerebellum"             
    ## [129] "cerebellum"              "cerebellum"             
    ## [131] "cerebellum"              "cerebellum"             
    ## [133] "cerebellum"              "cerebellum"             
    ## [135] "cerebellum"              "cerebellum"             
    ## [137] "cerebellum"              "cerebellum"             
    ## [139] "cerebellum"              "cerebellum"             
    ## [141] "cerebellum"              "cerebellum"             
    ## [143] "cerebellum"              "cerebellum"             
    ## [145] "cerebellum"              "cerebellum"             
    ## [147] "cerebellum"              "cerebellum"             
    ## [149] "cerebellum"              "cerebellum"             
    ## [151] "cerebellum"              "cerebellum"             
    ## [153] "cerebellum"              "cerebellum"             
    ## [155] "cerebellum"              "cerebellum"             
    ## [157] "cerebellum"              "cerebellum"             
    ## [159] "cerebellum"              "cerebellum"             
    ## [161] "cerebellum"              "cerebellum"             
    ## [163] "cerebellum"              "cerebellum"             
    ## [165] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [167] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [169] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [171] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [173] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [175] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [177] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [179] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [181] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [183] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [185] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [187] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [189] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [191] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [193] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [195] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [197] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [199] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [201] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [203] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [205] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [207] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [209] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [211] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [213] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [215] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [217] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [219] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [221] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [223] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [225] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [227] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [229] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [231] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [233] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [235] "entorhinal.cortex"       "entorhinal.cortex"      
    ## [237] "frontal.cortex"          "frontal.cortex"         
    ## [239] "frontal.cortex"          "frontal.cortex"         
    ## [241] "frontal.cortex"          "frontal.cortex"         
    ## [243] "frontal.cortex"          "frontal.cortex"         
    ## [245] "frontal.cortex"          "frontal.cortex"         
    ## [247] "frontal.cortex"          "frontal.cortex"         
    ## [249] "frontal.cortex"          "frontal.cortex"         
    ## [251] "frontal.cortex"          "frontal.cortex"         
    ## [253] "frontal.cortex"          "frontal.cortex"         
    ## [255] "frontal.cortex"          "frontal.cortex"         
    ## [257] "frontal.cortex"          "frontal.cortex"         
    ## [259] "frontal.cortex"          "frontal.cortex"         
    ## [261] "frontal.cortex"          "frontal.cortex"         
    ## [263] "frontal.cortex"          "frontal.cortex"         
    ## [265] "frontal.cortex"          "frontal.cortex"         
    ## [267] "frontal.cortex"          "frontal.cortex"         
    ## [269] "frontal.cortex"          "frontal.cortex"         
    ## [271] "frontal.cortex"          "frontal.cortex"         
    ## [273] "frontal.cortex"          "frontal.cortex"         
    ## [275] "frontal.cortex"          "frontal.cortex"         
    ## [277] "frontal.cortex"          "frontal.cortex"         
    ## [279] "frontal.cortex"          "frontal.cortex"         
    ## [281] "frontal.cortex"          "frontal.cortex"         
    ## [283] "frontal.cortex"          "frontal.cortex"         
    ## [285] "frontal.cortex"          "frontal.cortex"         
    ## [287] "frontal.cortex"          "frontal.cortex"         
    ## [289] "frontal.cortex"          "frontal.cortex"         
    ## [291] "frontal.cortex"          "frontal.cortex"         
    ## [293] "frontal.cortex"          "frontal.cortex"         
    ## [295] "frontal.cortex"          "frontal.cortex"         
    ## [297] "frontal.cortex"          "frontal.cortex"         
    ## [299] "frontal.cortex"          "frontal.cortex"         
    ## [301] "frontal.cortex"          "frontal.cortex"         
    ## [303] "frontal.cortex"          "frontal.cortex"         
    ## [305] "frontal.cortex"          "frontal.cortex"         
    ## [307] "frontal.cortex"          "frontal.cortex"         
    ## [309] "frontal.cortex"          "frontal.cortex"         
    ## [311] "frontal.cortex"          "frontal.cortex"         
    ## [313] "frontal.cortex"          "frontal.cortex"         
    ## [315] "frontal.cortex"          "frontal.cortex"         
    ## [317] "frontal.cortex"          "frontal.cortex"         
    ## [319] "frontal.cortex"          "frontal.cortex"         
    ## [321] "frontal.cortex"          "frontal.cortex"         
    ## [323] "frontal.cortex"          "frontal.cortex"         
    ## [325] "frontal.cortex"          "frontal.cortex"         
    ## [327] "frontal.cortex"          "frontal.cortex"         
    ## [329] "frontal.cortex"          "frontal.cortex"         
    ## [331] "frontal.cortex"          "frontal.cortex"         
    ## [333] "frontal.cortex"          "frontal.cortex"         
    ## [335] "frontal.cortex"          "frontal.cortex"         
    ## [337] "frontal.cortex"          "frontal.cortex"         
    ## [339] "frontal.cortex"          "frontal.cortex"         
    ## [341] "frontal.cortex"          "frontal.cortex"         
    ## [343] "frontal.cortex"          "frontal.cortex"         
    ## [345] "frontal.cortex"          "frontal.cortex"         
    ## [347] "frontal.cortex"          "frontal.cortex"         
    ## [349] "frontal.cortex"          "frontal.cortex"         
    ## [351] "frontal.cortex"          "frontal.cortex"         
    ## [353] "frontal.cortex"          "frontal.cortex"         
    ## [355] "frontal.cortex"          "frontal.cortex"         
    ## [357] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [359] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [361] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [363] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [365] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [367] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [369] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [371] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [373] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [375] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [377] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [379] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [381] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [383] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [385] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [387] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [389] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [391] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [393] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [395] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [397] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [399] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [401] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [403] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [405] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [407] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [409] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [411] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [413] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [415] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [417] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [419] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [421] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [423] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [425] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [427] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [429] "superior.temporal.gyrus" "superior.temporal.gyrus"
    ## [431] "superior.temporal.gyrus" "superior.temporal.gyrus"

``` r
#Creating beta and M value dataframes
B.norm <- GSE43414_cell_cor_sorted_by_brain_regions
B.normdf <- as.data.frame(B.norm)

M.norm <- beta2m(GSE43414_cell_cor_sorted_by_brain_regions)
M.normdf <- as.data.frame(M.norm)

#Subsetting columns so computer doesn't blow up
myvars <- c("6042316024_R01C01", "6042316031_R04C01", "6929718136_R04C02", "6969568087_R02C01", "6042316069_R04C02", "6055432060_R01C01", "6042316054_R03C01", "6042316103_R06C01")
subM.norm <- as.matrix(M.normdf[myvars])
View(subM.norm)

subB.norm <- as.matrix(B.normdf[myvars])
View(subB.norm)

submeta <- c("cerebellum", "cerebellum", "entorhinal.cortex", "entorhinal.cortex", "frontal.cortex", "frontal.cortex", "superior.temporal.gyrus", "superior.temporal.gyrus")

View(myvars)
myvarsdf <- as.data.frame(myvars)
View(myvarsdf)

submetadf <- as.data.frame(submeta)
View(submetadf)

submetabind <- cbind(myvarsdf,submetadf)
rownames(submetabind) <- submetabind$myvars

View(submetabind)
keep <- c("submeta")
submetabind <- submetabind[keep]
```

``` r
# Pre-DMRcate Setup
tissue <- submetabind$submeta
# create the design for the model
design <- model.matrix(~tissue)
View(design)

# fit the linear model
fit <- lmFit(subM.norm, design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(tissueentorhinal.cortex-tissuefrontal.cortex,
tissueentorhinal.cortex-tissuesuperior.temporal.gyrus,
levels=design)
```

    ## Warning in makeContrasts(tissueentorhinal.cortex - tissuefrontal.cortex, :
    ## Renaming (Intercept) to Intercept

``` r
contMatrix
```

    ##                                Contrasts
    ## Levels                          tissueentorhinal.cortex - tissuefrontal.cortex
    ##   Intercept                                                                  0
    ##   tissueentorhinal.cortex                                                    1
    ##   tissuefrontal.cortex                                                      -1
    ##   tissuesuperior.temporal.gyrus                                              0
    ##                                Contrasts
    ## Levels                          tissueentorhinal.cortex - tissuesuperior.temporal.gyrus
    ##   Intercept                                                                           0
    ##   tissueentorhinal.cortex                                                             1
    ##   tissuefrontal.cortex                                                                0
    ##   tissuesuperior.temporal.gyrus                                                      -1

``` r
# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
```

    ## Warning in contrasts.fit(fit, contMatrix): row names of contrasts don't
    ## match col names of coefficients

``` r
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))
```

    ##    tissueentorhinal.cortex - tissuefrontal.cortex
    ## -1                                              0
    ## 0                                          338359
    ## 1                                               0
    ##    tissueentorhinal.cortex - tissuesuperior.temporal.gyrus
    ## -1                                                       0
    ## 0                                                   338358
    ## 1                                                        1

``` r
# get the table of results for the first contrast
ann450kSub <- ann450k[match(rownames(subM.norm),ann450k$Name),
c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)
```

    ##              chr       pos strand       Name    Probe_rs Probe_maf
    ## cg18071071  chr4   1558451      + cg18071071        <NA>        NA
    ## cg12419862 chr22  24373484      + cg12419862        <NA>        NA
    ## cg03606954 chr11  43701349      - cg03606954        <NA>        NA
    ## cg04367395 chr17  27406907      - cg04367395        <NA>        NA
    ## cg14240646 chr10  27532245      - cg14240646        <NA>        NA
    ## cg12003941  chr6 168436019      - cg12003941 rs117759935  0.020399
    ##                 CpG_rs  CpG_maf      SBE_rs  SBE_maf
    ## cg18071071        <NA>       NA        <NA>       NA
    ## cg12419862        <NA>       NA        <NA>       NA
    ## cg03606954        <NA>       NA        <NA>       NA
    ## cg04367395 rs116625318 0.027663 rs116625318 0.027663
    ## cg14240646        <NA>       NA        <NA>       NA
    ## cg12003941        <NA>       NA        <NA>       NA
    ##                        Islands_Name Relation_to_Island UCSC_RefGene_Name
    ## cg18071071     chr4:1557561-1557772            S_Shore                  
    ## cg12419862  chr22:24372912-24373725             Island         LOC391322
    ## cg03606954  chr11:43702016-43702597            N_Shore          HSD17B12
    ## cg04367395  chr17:27406753-27408059             Island     MYO18A;MYO18A
    ## cg14240646  chr10:27530938-27531363            S_Shore       ACBD5;ACBD5
    ## cg12003941 chr6:168435835-168436086             Island       KIF25;KIF25
    ##            UCSC_RefGene_Accession UCSC_RefGene_Group Phantom DMR Enhancer
    ## cg18071071                                                               
    ## cg12419862           NM_001144931               Body                     
    ## cg03606954              NM_016142            TSS1500                     
    ## cg04367395    NM_078471;NM_203318          Body;Body                     
    ## cg14240646 NR_024150;NM_001042473    TSS1500;TSS1500                     
    ## cg12003941    NM_005355;NM_030615          Body;Body                 TRUE
    ##                       HMM_Island Regulatory_Feature_Name
    ## cg18071071     4:1528505-1528859       4:1558297-1558615
    ## cg12419862  22:22702803-22703712    22:24372239-24373584
    ## cg03606954                                              
    ## cg04367395  17:24430880-24432292                        
    ## cg14240646                          10:27530593-27532409
    ## cg12003941 6:168178688-168179417                        
    ##                   Regulatory_Feature_Group  DHS     logFC   AveExpr
    ## cg18071071 Unclassified_Cell_type_specific TRUE -2.060602 -1.607965
    ## cg12419862             Promoter_Associated      -1.898726 -3.844594
    ## cg03606954                                       1.786006 -2.858402
    ## cg04367395                                       2.967612  1.342130
    ## cg14240646             Promoter_Associated       2.619912 -1.512874
    ## cg12003941                                 TRUE  2.106232 -2.441058
    ##                    t      P.Value adj.P.Val         B
    ## cg18071071 -8.918197 2.568399e-06 0.8690409 -4.056561
    ## cg12419862 -7.845401 8.661157e-06 0.9999965 -4.077742
    ## cg03606954  7.150982 2.033107e-05 0.9999965 -4.095569
    ## cg04367395  7.088404 2.201743e-05 0.9999965 -4.097374
    ## cg14240646  6.829295 3.078092e-05 0.9999965 -4.105243
    ## cg12003941  6.631228 3.999314e-05 0.9999965 -4.111719

``` r
# Saving the data.frame file as an Excel
write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# As a quick check, plot the top 4 most significantly differentially methylated CpGs
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
plotCpg(subB.norm, cpg=cpg, pheno=submetabind$submeta, type = "categorical", measure = "beta")
})
```

![](Mar_26_Subset_DMR_files/figure-markdown_github/unnamed-chunk-3-1.png)

    ## $cg18071071
    ## NULL
    ## 
    ## $cg12419862
    ## NULL
    ## 
    ## $cg03606954
    ## NULL
    ## 
    ## $cg04367395
    ## NULL

Differential methylation analysis of regions

Use dmrcate function to combine individual CpG statistics to identify differentially methylated regions. DMRs$results contains all of the regions found, with genomic annotations and p-values

``` r
#DMRcate
design <- model.matrix(~submetabind$submeta)
myannotation <- cpg.annotate(datatype = "array", subM.norm, what = "M", arraytype = "450K", analysis.type="differential", design=design, contrasts = FALSE, cont.matrix = NULL, fdr = 0.05, coef=1)
```

    ## Your contrast returned 317295 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

``` r
dmrcoutput <- dmrcate(myannotation, lambda=1000)
```

    ## Fitting chr1...

    ## Fitting chr10...

    ## Fitting chr11...

    ## Fitting chr12...

    ## Fitting chr13...

    ## Fitting chr14...

    ## Fitting chr15...

    ## Fitting chr16...

    ## Fitting chr17...

    ## Fitting chr18...

    ## Fitting chr19...

    ## Fitting chr2...

    ## Fitting chr20...

    ## Fitting chr21...

    ## Fitting chr22...

    ## Fitting chr3...

    ## Fitting chr4...

    ## Fitting chr5...

    ## Fitting chr6...

    ## Fitting chr7...

    ## Fitting chr8...

    ## Fitting chr9...

    ## Demarcating regions...

    ## Done!

``` r
str(myannotation)
```

    ## List of 7
    ##  $ ID    : Factor w/ 338359 levels "cg00000029","cg00000108",..: 180947 301360 201815 15178 158126 291516 2684 82058 283599 230470 ...
    ##  $ stat  : num [1:338359] 8.45 11.01 8.59 -27.56 5.63 ...
    ##  $ CHR   : Factor w/ 22 levels "chr1","chr10",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ pos   : int [1:338359] 15865 534242 710097 714177 758829 763119 790667 805102 805554 812539 ...
    ##  $ betafc: num [1:338359] 0.8599 0.8247 0.7583 0.0306 0.6999 ...
    ##  $ indfdr: num [1:338359] 5.85e-06 5.09e-07 5.01e-06 3.95e-10 1.94e-04 ...
    ##  $ is.sig: logi [1:338359] TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  - attr(*, "row.names")= int [1:338359] 1 2 3 4 5 6 7 8 9 10 ...
    ##  - attr(*, "class")= chr "annot"

``` r
DMRs <- dmrcate(myannotation, lambda=1000, C=2)
```

    ## Fitting chr1...

    ## Fitting chr10...

    ## Fitting chr11...

    ## Fitting chr12...

    ## Fitting chr13...

    ## Fitting chr14...

    ## Fitting chr15...

    ## Fitting chr16...

    ## Fitting chr17...

    ## Fitting chr18...

    ## Fitting chr19...

    ## Fitting chr2...

    ## Fitting chr20...

    ## Fitting chr21...

    ## Fitting chr22...

    ## Fitting chr3...

    ## Fitting chr4...

    ## Fitting chr5...

    ## Fitting chr6...

    ## Fitting chr7...

    ## Fitting chr8...

    ## Fitting chr9...

    ## Demarcating regions...

    ## Done!

``` r
head(DMRs$results)
```

    ##                        coord no.cpgs minfdr Stouffer maxbetafc meanbetafc
    ## 90716 chr6:33156164-33181870     193      0        0 0.9657504  0.4229455
    ## 90732 chr6:33279563-33291947     154      0        0 0.9483754  0.3782921
    ## 90592 chr6:32036532-32059605     135      0        0 0.9625611  0.7573567
    ## 90724 chr6:33230625-33249087     126      0        0 0.9523675  0.5217330
    ## 90603 chr6:32114490-32123701     105      0        0 0.9266316  0.2943591
    ## 90670 chr6:32804628-32815831     101      0        0 0.9692180  0.5004982

``` r
# convert the regions to annotated genomic ranges
data(dmrcatedata)
results.ranges <- extractRanges(DMRs, genome = "hg19")

# Plotting the DMR
groups <- c(cerebellum="magenta", superior.temporal.gyrus="forestgreen", entorhinal.cortex="blue", frontal.cortex="yellow")
cols <- groups[as.character(submetabind$submeta)]
samps <- c(1:8)
DMR.plot(ranges=results.ranges, dmr=1, CpGs=subM.norm, what="Beta", arraytype = "450K",
phen.col=cols, genome="hg19", samps=samps, toscale=FALSE, plotmedians = TRUE)
```

![](Mar_26_Subset_DMR_files/figure-markdown_github/unnamed-chunk-4-1.png)
