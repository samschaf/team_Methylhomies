Data Normalization: BMIQ
========================================================
## Original author: Sumaiya Islam
## Edited by: Samantha Schaffner
## Date updated: March 31, 2017

Many inter- and intra-sample normalization methods are available for 450k methylation data; however, quantile normalization (QN) followed by beta-mixture quantile normalization (BMIQ) has been shown to be one of the most effective (Wang et al, 2015; PMC4623491). Quantile normalization is used to correct for broad differences in the distributions of methylation data, which can be due to factors such as tissue. BMIQ specifically corrects for differences in probe type distribution (further explained below), as the 450k array contains two types of probe chemistry. While each method can be used alone, we chose to increase our stringency and allow for between-sample analysis by applying both.

The beta values we downloaded from GEO (GSE43414) had already been quantile normalized, so we only applied BMIQ here.

### A. Set up working directory & packages

R version 3.2.3 (2015-12-10)

We will initially set our working directory and load our libraries.

```r
setwd("/home/sschaffner/team_Methylhomies")
unloadNamespace("minfi")
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```
## Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display
```

```
## Warning: 'rgl_init' failed, running with rgl.useNULL = TRUE
```

```r
unloadNamespace("mgcv")
unloadNamespace("lumi")
library(methylumi)
```

```
## Loading required package: Biobase
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist, unsplit
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: scales
```

```
## Loading required package: reshape2
```

```
## Loading required package: ggplot2
```

```
## Loading required package: matrixStats
```

```
## matrixStats v0.50.2 (2016-04-24) successfully loaded. See ?matrixStats for help.
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## Loading required package: FDb.InfiniumMethylation.hg19
```

```
## Loading required package: GenomicFeatures
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: GenomicRanges
```

```
## 
## Attaching package: 'GenomicRanges'
```

```
## The following object is masked from 'package:matrixStats':
## 
##     rowRanges
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: TxDb.Hsapiens.UCSC.hg19.knownGene
```

```
## Loading required package: org.Hs.eg.db
```

```
## Loading required package: DBI
```

```
## 
```

```
## Loading required package: minfi
```

```
## Loading required package: lattice
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## Loading required package: bumphunter
```

```
## Loading required package: foreach
```

```
## Loading required package: iterators
```

```
## Loading required package: locfit
```

```
## locfit 1.5-9.1 	 2013-03-22
```

```r
library(gplots)
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:IRanges':
## 
##     space
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
library(marray)
```

```
## Loading required package: limma
```

```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
library(lumi)
```

```
## 
## Attaching package: 'lumi'
```

```
## The following objects are masked from 'package:methylumi':
## 
##     estimateM, getHistory
```

```r
library(lattice)
library("RColorBrewer")
library(knitr)
library(xtable)
library(wateRmelon)
```

```
## Loading required package: ROC
```

```
## Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```r
library(limma)
library(sva)
```

```
## Loading required package: mgcv
```

```
## Loading required package: nlme
```

```
## 
## Attaching package: 'nlme'
```

```
## The following object is masked from 'package:Biostrings':
## 
##     collapse
```

```
## The following object is masked from 'package:IRanges':
## 
##     collapse
```

```
## This is mgcv 1.8-17. For overview type 'help("mgcv-package")'.
```

```
## Loading required package: genefilter
```

```
## 
## Attaching package: 'genefilter'
```

```
## The following objects are masked from 'package:ROC':
## 
##     AUC, pAUC
```

```
## The following objects are masked from 'package:methylumi':
## 
##     featureFilter, varFilter
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     rowSds, rowVars
```

```
## The following object is masked from 'package:base':
## 
##     anyNA
```

```r
library(RPMM)
```

```
## Loading required package: cluster
```

```
## 
## Attaching package: 'RPMM'
```

```
## The following object is masked from 'package:limma':
## 
##     ebayes
```

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:nlme':
## 
##     collapse
```

```
## The following object is masked from 'package:lumi':
## 
##     combine
```

```
## The following object is masked from 'package:methylumi':
## 
##     combine
```

```
## The following objects are masked from 'package:Biostrings':
## 
##     collapse, intersect, setdiff, setequal, union
```

```
## The following object is masked from 'package:XVector':
## 
##     slice
```

```
## The following object is masked from 'package:AnnotationDbi':
## 
##     select
```

```
## The following objects are masked from 'package:GenomicRanges':
## 
##     intersect, setdiff, union
```

```
## The following object is masked from 'package:GenomeInfoDb':
## 
##     intersect
```

```
## The following objects are masked from 'package:IRanges':
## 
##     collapse, desc, intersect, setdiff, slice, union
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     intersect, rename, setdiff, union
```

```
## The following object is masked from 'package:matrixStats':
## 
##     count
```

```
## The following object is masked from 'package:Biobase':
## 
##     combine
```

```
## The following objects are masked from 'package:BiocGenerics':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

### B. Load files



```r
load(file = "GSE43414_filtered.RData") # filtered data
load(file="Brain_meta_matched_GSE43414.RData")
load("Priest.RData") # Priest methylumi object for fData
dim(GSE43414_filtered)
```

```
## [1] 414904    432
```

```r
dim(Priest)
```

```
## Features  Samples 
##   485577        8
```

### C. Probe-type Normalization: BMIQ

The 450K Illumina Infinium Array has inherent variation associated with its methodologies which must be accounted for in our analyses. Much of this variation is attributed to the use of two types of probes used in the array, Type I and Type II.

Type I probes contain two bead types corresponding to an unmethylated (U) or methylated (M) status. Type I probes obtain methylation status as a result of fluoresence expressed after a single base pair extension occurs just after the target basepair, resulting in only one color channel being utilized (red). Type I probes also assume that any CpG sites underlying the probe are of the same status as the query site (methylated or unmethylated). The beta values for Type I probes are then determined by this formula b= M/(U + M). Importantly, Type I probes are enriched in regions of high CpG density (carry 3 or more CpG sites underlying the probe body), particularly those associated with promoters). Type II probes tend to occur in lower CpG density regions (carry 3 or less CpG sites underlying the probe body). Type II probes do not assume the methylation status of underlying CpG sites within the probe and so consist of a combination of degenerate probes with varying combinations of up to three underlying CpG sites. Type II probes also detect methylation status with a single base pair extension except that the site being extended is the CpG site of detection and so require two fluorescent colors green for methylated (M) and red for unmethylated (U) sites. Type II probe beta values are calculated using this formula b = Green (M)/(Red (U) + Green (M)). In terms of considerations for normalization, Type I probes have a much higher dynamic range than Type II probes, which may introduce an enrichment bias towards Type I probes when ranking probes in supervised analyses (Teschendorff et al 2013, Bioinformatics). 

Due to these inherent differences between Type I and Type II probes used in the Illumina Infinium 450K array several groups in the field have deveolped various normalization analyses to correct for the differences between these probes. We will be using an intra-sample normalization method that corrects for probe-type differences called BMIQ (Beta MIxture Quantile dilation) (Teschendorff et al 2013, Bioinformatics).

#### BMIQ (Beta MIxture Quantile dilation)

BMIQ is an intra-sample normalisation procedure, correcting the bias of type-2 probe values. BMIQ uses a 3-step procedure: (i) fitting of a 3-state beta mixture model, (ii) transformation of state-membership probabilities of type2 probes into quantiles of the type1 distribution, and (iii) a conformal transformation for the hemi-methylated probes. Exact details can be found in the reference (Teschendorff et al 2013, Bioinformatics).



```r
dim(Priest<-Priest[featureNames(Priest)%in%rownames(GSE43414_filtered),])
```

```
## Features  Samples 
##   414904        8
```

```r
head(probe_design<-as.character(fData(Priest)$INFINIUM_DESIGN_TYPE))
```

```
## [1] "II" "II" "II" "II" "II" "II"
```

```r
probe_design.v<- replace(probe_design, probe_design=="I", 1)
probe_design.cor<- replace(probe_design.v, probe_design.v=="II", 2)
probe_design.cor<-as.numeric(probe_design.cor)
identical(nrow(GSE43414_filtered), length(probe_design.cor))
```

```
## [1] TRUE
```


```r
# Run BMIQ across each sample
betas_normalized<-apply(GSE43414_filtered, 2, function(x) BMIQ(x,probe_design.cor)) # this code takes a long time---best to run it overnight esp if you have a lot of samples
```


```r
# extract normalized beta values and reshape
betas_normalized_betas<-lapply(1:length(betas_normalized), function(x) betas_normalized[[x]]$nbeta)
betas_normalized_betas<-do.call(rbind, betas_normalized_betas)
betas_norm.fin<-t(betas_normalized_betas)
colnames(betas_norm.fin)<-colnames(GSE43414_filtered)
head(betas_norm.fin)
```

```
##            6057825008_R02C02 6057825008_R03C01 6057825008_R04C01
## cg00000029         0.7155851         0.7473098         0.7150771
## cg00000108         0.4760463         0.5750774         0.5462654
## cg00000109         0.8311210         0.8553274         0.7011020
## cg00000165         0.4169461         0.2900396         0.3916608
## cg00000236         0.8090814         0.7763379         0.7936374
## cg00000289         0.4965480         0.4294358         0.4366022
##            6057825008_R04C02 6057825008_R05C01 6057825008_R05C02
## cg00000029         0.7230920         0.7412611         0.6856197
## cg00000108         0.5866269         0.5143037         0.4713554
## cg00000109         0.8590713         0.8218203         0.8010489
## cg00000165         0.2641833         0.4385457         0.4502902
## cg00000236         0.7268182         0.7838728         0.7124051
## cg00000289         0.3841312         0.3722781         0.4662913
##            6057825008_R06C01 6057825008_R06C02 6057825014_R01C02
## cg00000029         0.7335495         0.7423605         0.7065903
## cg00000108         0.4810017         0.5077457         0.4507651
## cg00000109         0.7192034         0.7345150         0.8393072
## cg00000165         0.4578980         0.2951594         0.3102756
## cg00000236         0.7264998         0.7903475         0.7582851
## cg00000289         0.5270033         0.5035221         0.4932998
##            6057825014_R02C02 6057825014_R03C01 6057825014_R03C02
## cg00000029         0.7398421         0.7325206         0.6650057
## cg00000108         0.6253629         0.4690692         0.3282851
## cg00000109         0.8105545         0.8432079         0.8155422
## cg00000165         0.5428855         0.5293544         0.4671881
## cg00000236         0.7627796         0.7998785         0.7990827
## cg00000289         0.4687571         0.4379764         0.5004206
##            6057825014_R04C01 6057825014_R04C02 6057825014_R06C01
## cg00000029         0.7503493         0.6692086         0.6281296
## cg00000108         0.5120762         0.4736127         0.4999381
## cg00000109         0.8388827         0.7698099         0.8278244
## cg00000165         0.3111092         0.5959031         0.4812945
## cg00000236         0.7122455         0.7335317         0.6380892
## cg00000289         0.4800653         0.4532797         0.5131484
##            6057825017_R01C01 6057825017_R01C02 6057825017_R02C02
## cg00000029         0.7547050         0.6736631         0.7283975
## cg00000108         0.3808441         0.5130942         0.5675379
## cg00000109         0.8118927         0.8546029         0.8509225
## cg00000165         0.2465659         0.3829607         0.3795256
## cg00000236         0.7849662         0.7672891         0.7864231
## cg00000289         0.5759853         0.4517695         0.5128083
##            6057825017_R03C01 6057825017_R03C02 6057825017_R04C01
## cg00000029         0.6980561         0.7013728         0.6924726
## cg00000108         0.5122192         0.3505210         0.5554337
## cg00000109         0.8394087         0.7966946         0.8265246
## cg00000165         0.3881700         0.6560528         0.4559561
## cg00000236         0.7554057         0.7566701         0.7471556
## cg00000289         0.5547272         0.4847385         0.4719426
##            6057825017_R04C02 6057825017_R05C01 6057825017_R05C02
## cg00000029         0.7003070         0.7313611         0.7223354
## cg00000108         0.6274366         0.5306652         0.5181015
## cg00000109         0.8585267         0.8240905         0.7548160
## cg00000165         0.2024706         0.3102421         0.2268066
## cg00000236         0.8002761         0.7401647         0.7648979
## cg00000289         0.4888861         0.4846905         0.5079131
##            6057825017_R06C01 6057825018_R01C01 6057825018_R02C02
## cg00000029         0.7633140         0.7151488         0.6582815
## cg00000108         0.6490924         0.3747713         0.4547296
## cg00000109         0.7891832         0.8294012         0.8127572
## cg00000165         0.3571132         0.3252707         0.3915431
## cg00000236         0.7653729         0.7708193         0.7846044
## cg00000289         0.4793284         0.4610800         0.4783734
##            6057825018_R03C02 6057825018_R04C01 6057825018_R04C02
## cg00000029         0.6213908         0.6501255         0.6527033
## cg00000108         0.4789395         0.4766254         0.4030058
## cg00000109         0.8525031         0.8244482         0.8632662
## cg00000165         0.4343551         0.2065450         0.3056248
## cg00000236         0.7646398         0.7311993         0.7903639
## cg00000289         0.4672966         0.4709675         0.5085050
##            6057825018_R05C01 6057825018_R05C02 6057825018_R06C01
## cg00000029         0.7527997         0.7243741         0.7152348
## cg00000108         0.5901884         0.5505611         0.6433144
## cg00000109         0.8062212         0.8271008         0.8015069
## cg00000165         0.2365263         0.1859365         0.1734938
## cg00000236         0.7495870         0.7718784         0.7475045
## cg00000289         0.4464124         0.5350191         0.5292506
##            6042316085_R01C01 6042316085_R01C02 6042316085_R02C02
## cg00000029         0.5863292         0.7546571         0.8570025
## cg00000108         0.3290140         0.5876477         0.6026739
## cg00000109         0.8072509         0.8099090         0.7256095
## cg00000165         0.4852637         0.2665095         0.2665629
## cg00000236         0.8073527         0.7786677         0.7076083
## cg00000289         0.5254314         0.5148180         0.2778758
##            6042316085_R03C01 6042316085_R03C02 6042316085_R05C01
## cg00000029         0.6703952         0.6796250         0.6333120
## cg00000108         0.5402998         0.5277018         0.5141434
## cg00000109         0.8550718         0.7921401         0.8237507
## cg00000165         0.2243528         0.2500278         0.2962147
## cg00000236         0.7307698         0.7474996         0.7095746
## cg00000289         0.4662759         0.4648340         0.4478230
##            6042316085_R05C02 6042316085_R06C01 6042316107_R01C02
## cg00000029         0.6682546         0.7318047         0.6793987
## cg00000108         0.4264600         0.6961005         0.4900278
## cg00000109         0.8049347         0.7459216         0.8205297
## cg00000165         0.2623667         0.4509597         0.2812345
## cg00000236         0.8074914         0.7258586         0.7758467
## cg00000289         0.4443292         0.4460273         0.4947861
##            6042316107_R03C01 6042316107_R04C02 6042316107_R05C01
## cg00000029         0.7123204         0.7011390         0.7099018
## cg00000108         0.5986989         0.5819851         0.6024593
## cg00000109         0.8267095         0.7974549         0.7033162
## cg00000165         0.5372126         0.3279275         0.3380993
## cg00000236         0.7748241         0.7878991         0.8439434
## cg00000289         0.5083868         0.4323291         0.4069517
##            6042316107_R05C02 6042316107_R06C01 6042316107_R06C02
## cg00000029         0.7039424         0.7521809         0.7489422
## cg00000108         0.5352970         0.3817610         0.4188992
## cg00000109         0.8456335         0.7275621         0.7635468
## cg00000165         0.2215484         0.2498215         0.2940681
## cg00000236         0.7381588         0.7340027         0.7297915
## cg00000289         0.4794181         0.4730627         0.4152364
##            6042316113_R01C01 6042316113_R01C02 6042316113_R02C01
## cg00000029         0.7135133         0.6766959         0.7389489
## cg00000108         0.6236640         0.4676754         0.7293088
## cg00000109         0.8060871         0.8159044         0.8295110
## cg00000165         0.4900234         0.3792105         0.3749096
## cg00000236         0.7754115         0.8034474         0.7815575
## cg00000289         0.5136255         0.4696242         0.5416748
##            6042316113_R02C02 6042316113_R03C02 6042316113_R04C02
## cg00000029         0.6514652         0.7395105         0.4410640
## cg00000108         0.6700423         0.5502467         0.9492810
## cg00000109         0.8535925         0.8243132         0.7472395
## cg00000165         0.6827235         0.5386003         0.1326731
## cg00000236         0.7866958         0.8012345         0.8526847
## cg00000289         0.3878583         0.4334214         0.3696050
##            6042316113_R05C01 6042316113_R05C02 6042316113_R06C01
## cg00000029         0.7705107         0.6714795         0.6716284
## cg00000108         0.5423803         0.5384886         0.5047308
## cg00000109         0.8575335         0.7706446         0.8292103
## cg00000165         0.2622692         0.5462160         0.3273333
## cg00000236         0.7981350         0.7515436         0.7510745
## cg00000289         0.4841544         0.5367255         0.5092499
##            6042316127_R01C01 6042316127_R01C02 6042316127_R02C02
## cg00000029         0.6807457         0.6997872         0.6726806
## cg00000108         0.4332197         0.3873681         0.4348633
## cg00000109         0.8006898         0.8534027         0.7839098
## cg00000165         0.2838326         0.7138943         0.4322258
## cg00000236         0.7310103         0.7676521         0.7975897
## cg00000289         0.4262618         0.5706739         0.4645224
##            6042316127_R03C01 6042316127_R03C02 6042316127_R04C01
## cg00000029         0.7689867         0.7237650         0.7464716
## cg00000108         0.5304383         0.4476618         0.5399321
## cg00000109         0.8571661         0.8535228         0.8491204
## cg00000165         0.5313372         0.3681693         0.3432446
## cg00000236         0.7906789         0.7647913         0.7595074
## cg00000289         0.4731158         0.5580238         0.4832075
##            6042316127_R04C02 6042316127_R05C02 6042316127_R06C01
## cg00000029         0.7103281         0.6540292         0.6775215
## cg00000108         0.6443829         0.3654771         0.5446684
## cg00000109         0.7970827         0.8538712         0.8644812
## cg00000165         0.5938369         0.1369145         0.3266228
## cg00000236         0.7728031         0.7448942         0.8047183
## cg00000289         0.4488339         0.5184027         0.5562286
##            6057825014_R01C02.1 6057825014_R02C02.1 6057825014_R03C01.1
## cg00000029           0.7491346           0.7862156           0.7912335
## cg00000108           0.4173299           0.6110330           0.4576655
## cg00000109           0.9375463           0.8853264           0.9395449
## cg00000165           0.2576824           0.5490669           0.5368872
## cg00000236           0.8029009           0.8436222           0.8613015
## cg00000289           0.5167327           0.5021914           0.5501859
##            6057825014_R03C02.1 6057825014_R04C01.1 6057825014_R04C02.1
## cg00000029           0.6947001           0.8042468           0.6533165
## cg00000108           0.2828165           0.4853211           0.4923968
## cg00000109           0.9405916           0.9448314           0.9104759
## cg00000165           0.4570909           0.2765419           0.6626910
## cg00000236           0.8507543           0.7561835           0.7614777
## cg00000289           0.4715845           0.4519967           0.5009036
##            6057825014_R06C01.1 6057825017_R01C01.1 6057825017_R01C02.1
## cg00000029           0.6555367           0.7881234           0.6759894
## cg00000108           0.4958942           0.3310383           0.5234797
## cg00000109           0.9323358           0.9200764           0.9372546
## cg00000165           0.4303525           0.2039439           0.3632490
## cg00000236           0.6835605           0.8135025           0.8352302
## cg00000289           0.4992227           0.5356037           0.4855475
##            6057825017_R02C01 6057825017_R02C02.1 6057825017_R03C01.1
## cg00000029         0.7588995           0.7515722           0.7368248
## cg00000108         0.4441467           0.5534015           0.5033450
## cg00000109         0.9723624           0.9730713           0.9412600
## cg00000165         0.5780590           0.3609646           0.3534531
## cg00000236         0.7925470           0.8116224           0.8200891
## cg00000289         0.4938019           0.5325529           0.5911767
##            6057825017_R03C02.1 6057825017_R04C01.1 6057825017_R04C02.1
## cg00000029           0.7016199           0.7108849           0.7361985
## cg00000108           0.3026039           0.5657915           0.6169407
## cg00000109           0.9582519           0.8937873           0.9346242
## cg00000165           0.6781932           0.4847591           0.1585663
## cg00000236           0.8307893           0.7939098           0.8389965
## cg00000289           0.4988779           0.5128353           0.4759829
##            6057825017_R05C01.1 6057825017_R05C02.1 6057825017_R06C01.1
## cg00000029           0.8100886           0.7707073           0.7743885
## cg00000108           0.5569919           0.5023789           0.6137444
## cg00000109           0.9245370           0.8734003           0.9073111
## cg00000165           0.2559633           0.2158548           0.3285839
## cg00000236           0.8287956           0.8208558           0.8292389
## cg00000289           0.5101939           0.4885447           0.4918013
##            6057825018_R01C01.1 6057825018_R02C02.1 6057825018_R03C01
## cg00000029           0.7353198           0.6539048        0.41023230
## cg00000108           0.3294345           0.4192477        0.99106536
## cg00000109           0.9170981           0.9158540        0.96203317
## cg00000165           0.2850746           0.4102094        0.09068355
## cg00000236           0.8392332           0.8612556        0.80160804
## cg00000289           0.4782483           0.4518788        0.48508924
##            6057825018_R03C02.1 6057825018_R04C01.1 6057825018_R04C02.1
## cg00000029           0.6483488           0.6797967           0.6721103
## cg00000108           0.4728200           0.4703839           0.3759843
## cg00000109           0.9383412           0.9405665           0.9406634
## cg00000165           0.4302484           0.1610697           0.2745874
## cg00000236           0.7931809           0.7483648           0.8576181
## cg00000289           0.5248759           0.4314067           0.5559808
##            6057825018_R05C01.1 6057825018_R05C02.1 6057825018_R06C01.1
## cg00000029           0.7880826           0.7549169           0.7431124
## cg00000108           0.5840685           0.5261700           0.6179170
## cg00000109           0.9175147           0.9506440           0.9189300
## cg00000165           0.1858399           0.1500332           0.1160653
## cg00000236           0.8097885           0.8198745           0.7943131
## cg00000289           0.4298458           0.5942159           0.5202997
##            6042316035_R01C01 6042316035_R02C02 6042316035_R03C01
## cg00000029         0.5648247         0.6442902         0.5090548
## cg00000108         0.9374204         0.9694230         0.9734507
## cg00000109         0.7850594         0.7880430         0.8158620
## cg00000165         0.2072560         0.2307040         0.2089230
## cg00000236         0.8777121         0.8534832         0.8700390
## cg00000289         0.4699768         0.5410782         0.4725598
##            6042316035_R03C02 6042316035_R04C01 6042316035_R05C01
## cg00000029         0.5439301         0.5709226         0.5664456
## cg00000108         0.9675922         0.9521103         0.9625666
## cg00000109         0.8271138         0.7801714         0.7553027
## cg00000165         0.2273942         0.1924357         0.2621559
## cg00000236         0.8615474         0.8475491         0.8557053
## cg00000289         0.5929521         0.4947272         0.4781464
##            6042316035_R05C02 6042316035_R06C02 6042316048_R01C01
## cg00000029         0.5837341         0.5611984         0.5590582
## cg00000108         0.9626112         0.9622780         0.9427216
## cg00000109         0.8127241         0.6908816         0.7609575
## cg00000165         0.1862612         0.2414585         0.2409154
## cg00000236         0.8198768         0.8741175         0.8163664
## cg00000289         0.5784699         0.4755858         0.4150420
##            6042316048_R01C02 6042316048_R02C01 6042316048_R03C01
## cg00000029         0.5412760         0.5641390         0.6184411
## cg00000108         0.9556588         0.9742960         0.9511156
## cg00000109         0.8768902         0.7426963         0.7684928
## cg00000165         0.2300810         0.2010412         0.1822481
## cg00000236         0.8809098         0.8290389         0.8908765
## cg00000289         0.4881680         0.5610942         0.4792138
##            6042316048_R03C02 6042316048_R04C01 6042316048_R04C02
## cg00000029         0.5371171         0.4573015         0.4940479
## cg00000108         0.9610825         0.9584332         0.9543361
## cg00000109         0.7685803         0.7595657         0.8353549
## cg00000165         0.1806652         0.1924786         0.1527446
## cg00000236         0.8821103         0.8481464         0.8595824
## cg00000289         0.5289176         0.5050688         0.4939040
##            6042316048_R05C01 6042316048_R05C02 6042316110_R01C01
## cg00000029         0.5235445         0.5690696         0.5718896
## cg00000108         0.9381575         0.9634064         0.9693222
## cg00000109         0.8194701         0.7840051         0.8221501
## cg00000165         0.2449164         0.2052229         0.1916265
## cg00000236         0.8511870         0.8572274         0.8665805
## cg00000289         0.5379502         0.5259747         0.5105639
##            6042316110_R01C02 6042316110_R03C02 6042316110_R04C01
## cg00000029         0.6775671         0.4963396         0.6195668
## cg00000108         0.9699721         0.9540233         0.9736557
## cg00000109         0.7491425         0.7939315         0.8063924
## cg00000165         0.2166596         0.2168971         0.2188742
## cg00000236         0.8981229         0.8796694         0.8439522
## cg00000289         0.5770005         0.5558024         0.5397282
##            6042316110_R04C02 6042316110_R05C01 6042316110_R06C01
## cg00000029         0.5306676         0.6518054         0.6585302
## cg00000108         0.9602342         0.9483648         0.9464948
## cg00000109         0.7659379         0.7829746         0.7124740
## cg00000165         0.1852846         0.2049528         0.1814644
## cg00000236         0.8577106         0.8420025         0.8818164
## cg00000289         0.4556471         0.5926411         0.5502970
##            6042316121_R02C02 6042316121_R03C01 6042316121_R03C02
## cg00000029         0.4654212         0.5630133         0.4623173
## cg00000108         0.9440368         0.9489927         0.9769402
## cg00000109         0.8177268         0.8101796         0.7950716
## cg00000165         0.2208312         0.2113132         0.1868420
## cg00000236         0.8409794         0.8976741         0.8445223
## cg00000289         0.5338225         0.5287386         0.4793433
##            6042316121_R04C02 6042316121_R05C01 6042316121_R05C02
## cg00000029         0.4904030         0.5304752         0.5738031
## cg00000108         0.9490174         0.9678450         0.9288039
## cg00000109         0.7259991         0.7865628         0.8048302
## cg00000165         0.1954731         0.1878726         0.2054240
## cg00000236         0.8630757         0.8401846         0.8734778
## cg00000289         0.4383439         0.4621512         0.5235788
##            6042316121_R06C01 6042316121_R06C02 6042316066_R01C01
## cg00000029         0.5373448         0.5027316         0.4753872
## cg00000108         0.9619446         0.9585639         0.9471931
## cg00000109         0.8458774         0.7173978         0.7628107
## cg00000165         0.2468511         0.2364343         0.1660496
## cg00000236         0.8896693         0.8665737         0.8546779
## cg00000289         0.5069966         0.5916035         0.4433447
##            6042316066_R02C01 6042316066_R02C02 6042316066_R03C01
## cg00000029         0.5051654         0.5805437         0.6092959
## cg00000108         0.9550010         0.9603668         0.9594025
## cg00000109         0.7443477         0.7664625         0.7813846
## cg00000165         0.1572296         0.1644106         0.2400123
## cg00000236         0.8788902         0.8596914         0.8194197
## cg00000289         0.3802116         0.5340450         0.4991248
##            6042316066_R04C01 6042316066_R04C02 6042316066_R05C01
## cg00000029         0.5165177         0.5167736         0.6151987
## cg00000108         0.9270888         0.9501668         0.9583326
## cg00000109         0.8502175         0.7879584         0.7527674
## cg00000165         0.2084616         0.1413274         0.2075420
## cg00000236         0.8650161         0.8589260         0.8593145
## cg00000289         0.6075691         0.5131851         0.4846585
##            6042316066_R06C01 6042316066_R06C02 6042316069_R01C01
## cg00000029         0.5439083         0.5099158         0.4216603
## cg00000108         0.9518331         0.9596469         0.9578153
## cg00000109         0.8414069         0.6914578         0.8349683
## cg00000165         0.1547783         0.2602677         0.1710548
## cg00000236         0.7866413         0.8512483         0.8902070
## cg00000289         0.4841570         0.5711923         0.5246502
##            6042316069_R01C02 6042316069_R02C01 6042316069_R03C01
## cg00000029         0.4536473         0.5513449         0.5545916
## cg00000108         0.9534153         0.9372124         0.9545802
## cg00000109         0.8103654         0.7600420         0.7828243
## cg00000165         0.1767534         0.1833675         0.1825023
## cg00000236         0.8672815         0.8227907         0.8464889
## cg00000289         0.4809465         0.4433966         0.4706374
##            6042316069_R03C02 6042316069_R04C02 6042316069_R05C01
## cg00000029         0.6124435         0.5335005         0.5516037
## cg00000108         0.9507139         0.9583435         0.9445112
## cg00000109         0.7734838         0.7984859         0.7194835
## cg00000165         0.2241170         0.1738689         0.2221011
## cg00000236         0.8584215         0.8818334         0.8504185
## cg00000289         0.5595562         0.5139040         0.5087503
##            6042316069_R06C02 6042316094_R01C02 6042316094_R02C01
## cg00000029         0.6060345         0.4640112         0.5599860
## cg00000108         0.9548677         0.9236432         0.9672835
## cg00000109         0.7566607         0.8459063         0.8146783
## cg00000165         0.2474629         0.1683365         0.1948722
## cg00000236         0.8913748         0.8370451         0.8251929
## cg00000289         0.5545140         0.5398459         0.4696633
##            6042316094_R03C02 6042316094_R04C01 6042316094_R04C02
## cg00000029         0.5499791         0.5292890         0.6694532
## cg00000108         0.9342578         0.9548269         0.9557720
## cg00000109         0.7886574         0.8209171         0.8084342
## cg00000165         0.1841903         0.2232893         0.1886650
## cg00000236         0.8541499         0.8755928         0.8713802
## cg00000289         0.6217677         0.5126326         0.5408360
##            6042316094_R05C01 6042316094_R05C02 6042316094_R06C01
## cg00000029         0.5196249         0.4852665         0.5885475
## cg00000108         0.9752753         0.9381301         0.9326208
## cg00000109         0.7228355         0.8282414         0.7864663
## cg00000165         0.2340116         0.1978492         0.1761000
## cg00000236         0.8690409         0.8726268         0.8555357
## cg00000289         0.5646950         0.5217890         0.5510647
##            6042316099_R01C01 6042316099_R01C02 6042316099_R02C01
## cg00000029         0.5122835         0.4805882         0.4906764
## cg00000108         0.9516297         0.9639121         0.9629007
## cg00000109                NA         0.7969902         0.7358754
## cg00000165         0.2062682         0.1739682         0.1810247
## cg00000236         0.8509288         0.8547022         0.8722068
## cg00000289         0.4405560         0.4906476         0.5123750
##            6042316099_R02C02 6042316099_R03C01 6042316099_R04C01
## cg00000029         0.5452895         0.6817959         0.5194849
## cg00000108         0.9672258         0.9655397         0.9416800
## cg00000109         0.8057926         0.7996384         0.7163574
## cg00000165         0.1638414         0.2459392         0.1940723
## cg00000236         0.8353157         0.8834347         0.8547508
## cg00000289         0.5124338         0.5740975         0.5159342
##            6042316099_R04C02 6042316099_R05C01 6969568082_R02C01
## cg00000029         0.4297615         0.4847706         0.5966640
## cg00000108         0.9392734         0.9506706         0.9719561
## cg00000109         0.8120894         0.7812931         0.8382244
## cg00000165         0.2000958         0.2011712         0.2201548
## cg00000236         0.8949947         0.8338669         0.8404914
## cg00000289         0.5234984         0.4912974         0.5056755
##            6969568082_R06C01 6969568082_R02C02 6969568082_R04C02
## cg00000029         0.4943082         0.6262704         0.5768440
## cg00000108         0.9637226         0.9604651         0.9501825
## cg00000109         0.6848135         0.7198529         0.7090254
## cg00000165         0.2282920         0.1556421         0.2314138
## cg00000236         0.8692345         0.8914621         0.8650333
## cg00000289         0.4178396         0.5251722         0.5161936
##            6969568082_R06C02 6969568084_R01C01 6969568084_R02C01
## cg00000029         0.6141492         0.5875177         0.5141288
## cg00000108         0.9579960         0.9451552         0.9539717
## cg00000109         0.7428723         0.7917714         0.6844736
## cg00000165         0.2258116         0.1785271         0.1808893
## cg00000236         0.8549799         0.8634997         0.8643876
## cg00000289         0.4226812         0.4682372         0.3985299
##            6969568084_R03C01 6969568084_R04C01 6969568084_R06C01
## cg00000029         0.5548962         0.5547428         0.5250814
## cg00000108         0.9621183         0.9437331         0.9649274
## cg00000109         0.6889236         0.7676988         0.7737852
## cg00000165         0.2364122         0.2001423         0.2294926
## cg00000236         0.8598610         0.8628635         0.8763880
## cg00000289         0.4491317         0.4998518         0.3819575
##            6969568084_R02C02 6969568084_R03C02 6969568084_R04C02
## cg00000029         0.5341048         0.5541846         0.6043138
## cg00000108         0.9567026         0.9513723         0.9636626
## cg00000109         0.6788877         0.7116606         0.7065202
## cg00000165         0.1841310         0.2100263         0.1862201
## cg00000236         0.8871925         0.8531816         0.8315817
## cg00000289         0.4902112         0.4470952         0.5062860
##            6969568084_R05C02 6969568087_R01C01 6969568087_R02C01
## cg00000029         0.5823729         0.5249667         0.6924179
## cg00000108         0.9736684         0.9510737         0.9652889
## cg00000109         0.7735198         0.6852531         0.7012814
## cg00000165         0.1714081         0.2553287         0.1793005
## cg00000236         0.8433863         0.8676442         0.8555818
## cg00000289         0.5272644         0.5353664         0.5438716
##            6969568087_R03C01 6969568087_R04C01 6969568087_R05C01
## cg00000029         0.5170160         0.5095897         0.4999419
## cg00000108         0.9538180         0.9419028         0.9415744
## cg00000109         0.7721032         0.7779493         0.7632831
## cg00000165         0.1773474         0.1669479         0.2111779
## cg00000236         0.8358303         0.9078334         0.8490148
## cg00000289         0.4336579         0.6272458         0.4239234
##            6969568087_R06C01 6969568087_R01C02 6969568087_R02C02
## cg00000029         0.5943415         0.4679759         0.5215820
## cg00000108         0.9703605         0.9643094         0.9477982
## cg00000109         0.7176204         0.7842813         0.8169286
## cg00000165         0.1853532         0.1596933         0.1717008
## cg00000236         0.9106414         0.8177979         0.8922681
## cg00000289         0.4718072         0.5228420         0.5031134
##            6969568087_R03C02 6969568087_R05C02 6969568087_R06C02
## cg00000029         0.6249249         0.5424356         0.4998601
## cg00000108         0.9644426         0.9439192         0.9292157
## cg00000109         0.7595445         0.7366149         0.7844510
## cg00000165         0.1816813         0.1902820         0.1872035
## cg00000236         0.8859874         0.8962023         0.8432548
## cg00000289         0.6167478         0.5764072         0.4332543
##            6969568118_R01C01 6969568118_R02C01 6969568118_R03C01
## cg00000029         0.5587847         0.5366289         0.5256195
## cg00000108         0.9590558         0.9672621         0.9298197
## cg00000109         0.6989730         0.7676475         0.8058072
## cg00000165         0.2363916         0.1705880         0.1938121
## cg00000236         0.8881869         0.8825145         0.8751198
## cg00000289         0.5200399         0.6284771         0.5472226
##            6969568118_R04C01 6969568118_R01C02 6969568118_R02C02
## cg00000029         0.4827513         0.5576476         0.5632094
## cg00000108         0.9524460         0.9581917         0.9599376
## cg00000109         0.7714995         0.7113394         0.7347177
## cg00000165         0.2231399         0.3643319         0.2332120
## cg00000236         0.8837956         0.8455748         0.8940032
## cg00000289         0.5058405         0.6668828         0.5502939
##            6969568118_R03C02 6969568118_R04C02 6969568118_R06C02
## cg00000029         0.6062242         0.6359101         0.6314785
## cg00000108         0.9755923         0.9507538         0.9609952
## cg00000109         0.7714065         0.7881279         0.7138720
## cg00000165         0.2347529         0.1639544         0.2205604
## cg00000236         0.8792152         0.8930244         0.8863521
## cg00000289         0.6016269         0.5600057         0.6144241
##            6929726046_R02C01 6929726046_R05C01 6929726046_R06C01
## cg00000029         0.5972168         0.6233957         0.6110711
## cg00000108         0.9588643         0.9594021         0.9622198
## cg00000109         0.8179758         0.7961111         0.7770863
## cg00000165         0.2054170         0.2254520         0.2392496
## cg00000236         0.8573842         0.8668127         0.8977855
## cg00000289         0.5603500         0.6418561         0.4746912
##            6929726046_R01C02 6929726046_R03C02 6929726046_R04C02
## cg00000029         0.6596143         0.6248318         0.4828038
## cg00000108         0.9575752         0.9679689         0.9459357
## cg00000109         0.7457175         0.7670225         0.7340656
## cg00000165         0.2593212         0.1872654         0.1532463
## cg00000236         0.8348067         0.8624587         0.8690241
## cg00000289         0.5019900         0.5918055         0.3927186
##            6929726046_R06C02 6929718123_R01C01 6929718123_R04C01
## cg00000029         0.6380119         0.5887385         0.6518362
## cg00000108         0.9570405         0.9557961         0.9598160
## cg00000109         0.7658559         0.7383195         0.8251310
## cg00000165         0.2459098         0.1848802         0.2198499
## cg00000236         0.8498176         0.8861981         0.8902118
## cg00000289         0.5126985         0.5433446         0.5496978
##            6929718123_R06C01 6929718123_R01C02 6929718123_R02C02
## cg00000029         0.6533115         0.6518642         0.6900728
## cg00000108         0.9475845         0.9671575         0.6022994
## cg00000109         0.7788388         0.7996852         0.8282015
## cg00000165         0.2258675         0.1684975         0.4686387
## cg00000236         0.8693348         0.9004155         0.7731752
## cg00000289         0.5333504         0.4975860         0.4582411
##            6929718123_R03C02 6929718123_R04C02 6929718123_R05C02
## cg00000029         0.6580946         0.5346194         0.5827275
## cg00000108         0.9449651         0.9603167         0.9649177
## cg00000109         0.7415941         0.7390978         0.7771280
## cg00000165         0.2005216         0.1654598         0.1400005
## cg00000236         0.8990295         0.8746821         0.8894425
## cg00000289         0.5411770         0.5342186         0.4454065
##            6929718123_R06C02 6929718136_R01C01 6929718136_R02C01
## cg00000029         0.6287792         0.6499048         0.5610762
## cg00000108         0.9609294         0.9526308         0.9588475
## cg00000109         0.8485557         0.8109103         0.7505746
## cg00000165         0.2642255         0.2137764         0.1317756
## cg00000236         0.8430929         0.9051429         0.8924989
## cg00000289         0.5285447         0.5997893         0.5984541
##            6929718136_R03C01 6929718136_R04C01 6929718136_R05C01
## cg00000029         0.6253596         0.7001711         0.5999848
## cg00000108         0.9562766         0.5989724         0.9666853
## cg00000109         0.8021144         0.7898402         0.7867059
## cg00000165         0.2652139         0.4329996         0.1573081
## cg00000236         0.8955675         0.7891602         0.8883816
## cg00000289         0.5733664         0.5291416         0.5651118
##            6929718136_R02C02 6929718136_R04C02 6929718136_R05C02
## cg00000029         0.6217419         0.6122813         0.6463810
## cg00000108         0.9573494         0.9296227         0.9703096
## cg00000109         0.7902620         0.7696537         0.7704401
## cg00000165         0.2184335         0.1981806         0.2064139
## cg00000236         0.8783428         0.8905697         0.8638230
## cg00000289         0.5219412         0.4933391         0.5533981
##            6929718136_R06C02 6929718138_R01C01 6929718138_R02C01
## cg00000029         0.6601443         0.5277092         0.5866813
## cg00000108         0.9805838         0.9552412         0.9386663
## cg00000109         0.7405599         0.7920923         0.7795591
## cg00000165         0.2440137         0.2367960         0.1990068
## cg00000236         0.8806126         0.8791655         0.8825385
## cg00000289         0.5665603         0.5293439         0.5424658
##            6929718138_R04C01 6929718138_R06C01 6929718138_R01C02
## cg00000029         0.6085816         0.6922041         0.5926276
## cg00000108         0.9482873         0.9579944         0.9251887
## cg00000109         0.7400660         0.7941928         0.7415514
## cg00000165         0.2550607         0.2442955         0.1993705
## cg00000236         0.9045073         0.8262563         0.8682785
## cg00000289         0.5809544         0.5831295         0.5913988
##            6929718138_R03C02 6929718138_R04C02 6929718138_R05C02
## cg00000029         0.6705911         0.6016980         0.6900630
## cg00000108         0.8324576         0.9568799         0.9559211
## cg00000109         0.8416751         0.7713785         0.7502481
## cg00000165         0.1859445         0.2143962         0.2071252
## cg00000236         0.8610728         0.8571997         0.8781304
## cg00000289         0.6249447         0.5873853         0.5180861
##            6929718138_R06C02 6042316054_R02C01 6042316054_R02C02
## cg00000029         0.5981074         0.4782172         0.6164615
## cg00000108         0.9504713         0.9487172         0.9681856
## cg00000109         0.7292595         0.7516800         0.7608622
## cg00000165         0.2463660         0.1911794         0.3016198
## cg00000236         0.8961420         0.8481909         0.8673007
## cg00000289         0.5650239         0.5101716         0.5281651
##            6042316054_R03C01 6042316054_R03C02 6042316054_R04C01
## cg00000029         0.5144383         0.6156955         0.5995615
## cg00000108         0.9645876         0.9718398         0.9478357
## cg00000109         0.7412285         0.7768048         0.7440045
## cg00000165         0.1968223         0.2136537         0.1847733
## cg00000236         0.8594819         0.8637425         0.8430882
## cg00000289         0.4316929         0.4711052         0.4460083
##            6042316054_R04C02 6042316054_R05C01 6042316054_R05C02
## cg00000029         0.4602862         0.5572565         0.4921731
## cg00000108         0.9569371         0.9574484         0.9656794
## cg00000109         0.7542408         0.7395082         0.7356372
## cg00000165         0.1587655         0.2005531         0.1646041
## cg00000236         0.8808484         0.8657415         0.8559763
## cg00000289         0.4341918         0.5362716         0.4836668
##            6042316054_R06C01 6042316063_R02C01 6042316063_R02C02
## cg00000029         0.5288727         0.5146489         0.5201955
## cg00000108         0.9429401         0.9504114         0.9634170
## cg00000109         0.7241723         0.7473021         0.7757578
## cg00000165         0.1689745         0.2096040         0.1710332
## cg00000236         0.8302525         0.8490026         0.8683969
## cg00000289         0.4177315         0.4682597         0.4952821
##            6042316063_R03C01 6042316063_R03C02 6042316063_R04C01
## cg00000029         0.5364068         0.5957896         0.5799662
## cg00000108         0.9603915         0.9654517         0.9625143
## cg00000109         0.7483685         0.7758874         0.7952112
## cg00000165         0.1996431         0.1522736         0.1716336
## cg00000236         0.8955455         0.8741475         0.8437540
## cg00000289         0.5589403         0.5676567         0.5835582
##            6042316063_R04C02 6042316063_R05C01 6042316063_R05C02
## cg00000029         0.5451724         0.5191032         0.5372289
## cg00000108         0.9671079         0.9495029         0.9577258
## cg00000109         0.7692675         0.6924889         0.7847343
## cg00000165         0.1642032         0.2450367         0.1751591
## cg00000236         0.8749714         0.8303970         0.8903241
## cg00000289         0.5265855         0.4408742         0.5466810
##            6042316063_R06C02 6042316065_R01C02 6042316065_R02C02
## cg00000029         0.5351374         0.5395369         0.4869192
## cg00000108         0.9642174         0.9634115         0.9618806
## cg00000109         0.7953915         0.8501019         0.7903546
## cg00000165         0.2314279         0.2583145         0.1791672
## cg00000236         0.8660145         0.8676971         0.8678378
## cg00000289         0.6182102         0.5883308         0.5962522
##            6042316065_R03C01 6042316065_R04C01 6042316065_R04C02
## cg00000029         0.5346656         0.4621745         0.5270289
## cg00000108         0.9546711         0.9608338         0.9747348
## cg00000109         0.8398410         0.7634881         0.7793421
## cg00000165         0.1685381         0.1661273         0.1799197
## cg00000236         0.8822457         0.8458234         0.7910773
## cg00000289         0.5234746         0.4766960         0.4409669
##            6042316065_R05C02 6042316065_R06C02 6042316103_R02C01
## cg00000029         0.5114927         0.5676406         0.5345802
## cg00000108         0.9661009         0.9312320         0.9662748
## cg00000109         0.7869129         0.7446451         0.8032611
## cg00000165         0.2318741         0.2212581         0.1626281
## cg00000236         0.8558516         0.8709915         0.8381887
## cg00000289         0.4487089         0.4532016         0.5019969
##            6042316103_R03C01 6042316103_R03C02 6042316103_R04C01
## cg00000029         0.4517364         0.4642570         0.5562916
## cg00000108         0.9494715         0.9484300         0.9642157
## cg00000109         0.8261364         0.7562284         0.7943939
## cg00000165         0.1408529         0.2029468         0.1993491
## cg00000236         0.8232630         0.8829268         0.8767461
## cg00000289         0.4123398         0.3777126         0.6171938
##            6042316103_R05C01 6042316103_R06C01 6042316103_R06C02
## cg00000029         0.5123803         0.4287862         0.5039373
## cg00000108         0.9580207         0.9490171         0.9667610
## cg00000109         0.8060757         0.7610099         0.7995670
## cg00000165         0.2156502         0.2191035         0.1239531
## cg00000236         0.8570358         0.8385848         0.8053751
## cg00000289         0.4886395         0.5028841         0.4762139
##            6042316036_R01C02 6042316036_R02C01 6042316036_R03C01
## cg00000029         0.5502360         0.4950063         0.5655792
## cg00000108         0.9416479         0.9662746         0.9721286
## cg00000109         0.7839472         0.7999386         0.8024621
## cg00000165         0.1549568         0.1789811         0.1822520
## cg00000236         0.8612282         0.8205982         0.8286890
## cg00000289         0.5316729         0.4180945         0.5558812
##            6042316036_R03C02 6042316036_R04C01 6042316036_R04C02
## cg00000029         0.5111778         0.4520282         0.5169098
## cg00000108         0.9653341         0.9567823         0.9593582
## cg00000109         0.8614218         0.7485247         0.8027678
## cg00000165         0.1853216         0.1719881         0.2167877
## cg00000236         0.8393892         0.8561293         0.8830092
## cg00000289         0.4810467         0.5422871         0.4923073
##            6042316036_R05C01 6042316036_R05C02 6042316036_R06C01
## cg00000029         0.5185402         0.5629312         0.5874069
## cg00000108         0.9750276         0.9529467         0.9633928
## cg00000109         0.7493501         0.7554330         0.7755352
## cg00000165         0.2474177         0.2073004         0.2218491
## cg00000236         0.8806050         0.8041446         0.8363458
## cg00000289         0.5817724         0.5281974         0.5076441
##            6042316036_R06C02 6042316050_R01C02 6042316050_R02C02
## cg00000029         0.5288324         0.5570775         0.6799940
## cg00000108         0.9809375         0.9653709         0.9492554
## cg00000109         0.7383075         0.7963882         0.7415052
## cg00000165         0.2496756         0.1884047         0.2210199
## cg00000236         0.8804731         0.8641942         0.8796539
## cg00000289         0.4280144         0.4700822         0.5319999
##            6042316050_R03C01 6042316050_R04C01 6042316050_R04C02
## cg00000029         0.5075827         0.5078638         0.4635262
## cg00000108         0.9433196         0.9561588         0.9530683
## cg00000109         0.7226100         0.7423379         0.8159256
## cg00000165         0.1800162         0.2071560         0.2113250
## cg00000236         0.8325195         0.8646320         0.8418305
## cg00000289         0.4948327         0.5530471         0.5460387
##            6042316050_R05C01 6042316050_R05C02 6042316050_R06C01
## cg00000029         0.5409550         0.4788625         0.5908822
## cg00000108         0.9540458         0.9757633         0.9574968
## cg00000109         0.7432095         0.7982669         0.7577749
## cg00000165         0.2465684         0.2343458         0.3013273
## cg00000236         0.8662031         0.7988073         0.8302720
## cg00000289         0.5213333         0.4415024         0.5246324
##            6042316050_R06C02 6042316053_R01C01 6042316053_R02C01
## cg00000029         0.5468009         0.5094040         0.5388631
## cg00000108         0.9472318         0.9419866         0.9597182
## cg00000109         0.7977494         0.8295560         0.7504758
## cg00000165         0.1766608         0.1977316         0.2281451
## cg00000236         0.8509583         0.8589864         0.8723527
## cg00000289         0.5384922         0.5345764         0.4574328
##            6042316053_R02C02 6042316053_R03C01 6042316053_R03C02
## cg00000029         0.5401414         0.5333585         0.5116654
## cg00000108         0.9704277         0.9684996         0.9492856
## cg00000109         0.8273585         0.7668999         0.7770455
## cg00000165         0.1900874         0.2175655         0.2112661
## cg00000236         0.8962432         0.8701645         0.8819420
## cg00000289         0.5573069         0.5000095         0.4440263
##            6042316053_R04C01 6042316053_R04C02 6042316053_R06C01
## cg00000029         0.5273325         0.4506555         0.5562778
## cg00000108         0.9555084         0.9648148         0.9793006
## cg00000109         0.8454360         0.7054634         0.7562836
## cg00000165         0.1984425         0.2222676         0.2382629
## cg00000236         0.8521733         0.8791608         0.8856776
## cg00000289         0.5195525         0.5161467         0.6034571
##            6042316053_R06C02 6042316061_R01C01 6042316061_R02C01
## cg00000029         0.5664858         0.5587900         0.5353188
## cg00000108         0.9305363         0.9522849         0.9559581
## cg00000109         0.8346305         0.7839407         0.7684399
## cg00000165         0.1954431         0.1900690         0.1963520
## cg00000236         0.8426782         0.8786027         0.9002134
## cg00000289         0.5343989         0.5083894         0.4877308
##            6042316061_R03C01 6042316061_R03C02 6042316061_R04C01
## cg00000029         0.4319818         0.5730776         0.5534999
## cg00000108         0.9524086         0.9372900         0.9535820
## cg00000109         0.7801626         0.8046207         0.8513558
## cg00000165         0.1718740         0.1693953         0.1909215
## cg00000236         0.8546690         0.8641682         0.8693821
## cg00000289         0.5192240         0.5581028         0.5981707
##            6042316061_R04C02 6042316061_R05C01 6042316061_R05C02
## cg00000029         0.5432400         0.6143826         0.5178285
## cg00000108         0.9620722         0.9568026         0.9534234
## cg00000109         0.7940719         0.8142096         0.7383457
## cg00000165         0.1568074         0.2107004         0.2006102
## cg00000236         0.8836336         0.8569796         0.8723986
## cg00000289         0.5650781         0.5593077         0.6097239
##            7796806022_R01C01 7796806022_R03C01 7796806022_R04C01
## cg00000029         0.6781872         0.4205138         0.4060192
## cg00000108         0.4817813         0.9583776         0.9532449
## cg00000109         0.8206792         0.8246873         0.7668464
## cg00000165         0.2792717         0.1776329         0.1429961
## cg00000236         0.6937727         0.8655539         0.8484270
## cg00000289         0.5211306         0.4258215         0.3893333
##            7796806022_R04C02 7796806022_R05C02 7796806022_R06C01
## cg00000029         0.6906002         0.3587955         0.6730440
## cg00000108         0.5725125         0.9532298         0.5012888
## cg00000109         0.8463561         0.7667181         0.7606820
## cg00000165         0.4739010         0.1614372         0.5356709
## cg00000236         0.7724523         0.8684135         0.7780402
## cg00000289         0.5251852         0.4439925         0.4793687
##            7786923046_R03C02 7786923046_R04C01 7786923046_R06C01
## cg00000029         0.4403790         0.4234597         0.4515836
## cg00000108         0.9501729         0.9609557         0.9735603
## cg00000109         0.7406954         0.7601529         0.7814367
## cg00000165         0.1286044         0.1279955         0.2057687
## cg00000236         0.9037331         0.7977105         0.8758453
## cg00000289         0.4655490         0.4306380         0.4194938
##            7796806038_R03C01 7796806038_R03C02 7796806038_R05C02
## cg00000029         0.4323942         0.4673821         0.6678296
## cg00000108         0.9598945         0.9633144         0.4103790
## cg00000109         0.7499998         0.7175238         0.8531827
## cg00000165         0.1288067         0.1342228         0.6289999
## cg00000236         0.8726983         0.8759394         0.7572251
## cg00000289         0.3597305         0.4319351         0.5169107
##            7786923063_R01C01 7786923063_R01C02 7786923063_R03C01
## cg00000029         0.7131747         0.4302477         0.7486647
## cg00000108         0.4643721         0.9617236         0.5920324
## cg00000109         0.8359494         0.7255942         0.8281899
## cg00000165         0.1309701         0.1302109         0.1666629
## cg00000236         0.7052862         0.8481570         0.7512689
## cg00000289         0.4495545         0.3955134         0.4620442
##            7786923063_R04C01 7786923063_R06C01 7786923107_R01C02
## cg00000029         0.7011278         0.6704281         0.5736766
## cg00000108         0.6190897         0.5243359         0.9571484
## cg00000109         0.8457495         0.8744288         0.7989618
## cg00000165         0.4392841         0.3983373         0.1897578
## cg00000236         0.7375130         0.7449580         0.8753346
## cg00000289         0.5339011         0.4333907         0.5992958
##            7786923107_R04C01 7786923107_R05C01 7786923107_R06C01
## cg00000029         0.5694866         0.4639869         0.5670369
## cg00000108         0.9480735         0.9450337         0.9575606
## cg00000109         0.7549413         0.8562149         0.7818841
## cg00000165         0.2454800         0.1594727         0.2151570
## cg00000236         0.8782595         0.8590555         0.8420168
## cg00000289         0.4939223         0.5168848         0.5674580
##            7796806016_R03C02 7796806016_R06C01 7796806002_R01C01
## cg00000029         0.5099681         0.5517690         0.5701424
## cg00000108         0.9619596         0.9536571         0.9512598
## cg00000109         0.7424323         0.7470648         0.7535935
## cg00000165         0.2155136         0.1952491         0.2707602
## cg00000236         0.8432236         0.8661246         0.9094962
## cg00000289         0.4744738         0.4816269         0.5327705
##            7796806002_R02C02 7796806002_R04C02 7796806002_R05C01
## cg00000029         0.5800413         0.5040908         0.5957194
## cg00000108         0.9508995         0.9533528         0.9653204
## cg00000109         0.7797671         0.7636010         0.7662014
## cg00000165         0.2307006         0.1951569         0.2356187
## cg00000236         0.8879652         0.8902108         0.8978864
## cg00000289         0.5414562         0.5509282         0.5780436
##            7796806002_R05C02 7796806002_R06C01 7796806002_R06C02
## cg00000029         0.4792577         0.4655978         0.6410201
## cg00000108         0.9566240         0.9426852         0.9617855
## cg00000109         0.7447697         0.7399154         0.7925691
## cg00000165         0.1849155         0.2404800         0.2521086
## cg00000236         0.8448357         0.8323868         0.8556126
## cg00000289         0.5122511         0.4970751         0.5793462
##            7796806029_R02C01 7796806029_R04C01 7796806029_R06C01
## cg00000029         0.5479653         0.5635585         0.5739276
## cg00000108         0.9532203         0.9409491         0.9591283
## cg00000109         0.7714734         0.8560301         0.7473326
## cg00000165         0.2112410         0.2024445         0.2385572
## cg00000236         0.8629999         0.8905994         0.8952651
## cg00000289         0.5610950         0.6201201         0.5146200
##            6042316024_R01C01 6042316024_R01C02 6042316024_R02C01
## cg00000029         0.8499977         0.6787593         0.7260670
## cg00000108         0.5743032         0.3883582         0.4875101
## cg00000109         0.8446022         0.8413743         0.8830643
## cg00000165         0.3555823         0.5457204         0.5620894
## cg00000236         0.4519425         0.7813738         0.7203926
## cg00000289         0.4736487         0.5501053         0.4740281
##            6042316024_R02C02 6042316024_R03C01 6042316024_R04C01
## cg00000029         0.7272392         0.5836817         0.6934428
## cg00000108         0.5210021         0.6013196         0.4900262
## cg00000109         0.8937810         0.8256995         0.8546968
## cg00000165         0.4039127         0.4699677         0.1105378
## cg00000236         0.7417594         0.7935990         0.7666092
## cg00000289         0.5500863         0.4143308         0.3921155
##            6042316024_R04C02 6042316024_R05C01 6042316024_R05C02
## cg00000029         0.6304856         0.5179943         0.7025902
## cg00000108         0.7405458         0.5098961         0.6666012
## cg00000109         0.8650148         0.8478986         0.8678955
## cg00000165         0.6228635         0.4550192         0.3822142
## cg00000236         0.8139429         0.6825037         0.7477922
## cg00000289         0.2612610         0.3556179         0.4007416
##            6042316024_R06C01 6042316024_R06C02 6042316031_R01C01
## cg00000029         0.6886309         0.5342506         0.7182587
## cg00000108         0.3389088         0.3916858         0.4188445
## cg00000109         0.7730849         0.7527215         0.8415704
## cg00000165         0.3803304         0.4277816         0.5931500
## cg00000236         0.8371392         0.7285412         0.7550779
## cg00000289         0.3655690         0.4317905         0.4223241
##            6042316031_R01C02 6042316031_R02C01 6042316031_R02C02
## cg00000029         0.7987426         0.6420701         0.8718785
## cg00000108         0.2313432         0.4662603         0.5825151
## cg00000109         0.8673390         0.8339884         0.7961466
## cg00000165         0.6455985         0.1709480         0.3468179
## cg00000236         0.7714645         0.7574329         0.6619781
## cg00000289         0.4068543         0.4346206         0.2467341
##            6042316031_R03C01 6042316031_R03C02 6042316031_R04C01
## cg00000029         0.5504552         0.7275835         0.6916727
## cg00000108         0.4741437         0.6224770         0.5512840
## cg00000109                NA         0.8352224         0.8900823
## cg00000165         0.6146838         0.4805424         0.3492074
## cg00000236         0.6692676         0.7368866         0.6707637
## cg00000289         0.3732900         0.5275620         0.4637646
##            6042316031_R04C02 6042316031_R05C01 6042316031_R05C02
## cg00000029         0.7298546         0.6237968         0.7684735
## cg00000108         0.5875540         0.6930464         0.4912222
## cg00000109         0.8444747         0.8380561         0.8493081
## cg00000165         0.5058626         0.3926057         0.4803445
## cg00000236         0.7635339         0.8333424         0.8065347
## cg00000289         0.5103992         0.4562873         0.4789690
##            6042316031_R06C01 6042316031_R06C02 6042316042_R01C01
## cg00000029         0.7428661         0.7009792         0.6935727
## cg00000108         0.3772110         0.5026222         0.5240654
## cg00000109         0.7898476         0.7210267         0.8693866
## cg00000165         0.3193151         0.6094232         0.6170143
## cg00000236         0.7138216         0.7326305         0.7036948
## cg00000289         0.4559711         0.3693611         0.4641241
##            6042316042_R01C02 6042316042_R02C01 6042316042_R02C02
## cg00000029         0.7029172         0.6953521         0.7670121
## cg00000108         0.6645697         0.5682274         0.4125199
## cg00000109         0.8411099         0.8753626         0.8610318
## cg00000165         0.2779230         0.2238846         0.5276414
## cg00000236         0.7858153         0.7094149         0.7636505
## cg00000289         0.5407614         0.5166911         0.5494303
##            6042316042_R03C01 6042316042_R03C02 6042316042_R04C01
## cg00000029         0.7424398         0.6966564         0.7566369
## cg00000108         0.4605814         0.5510655         0.6207285
## cg00000109         0.9005737         0.7545792         0.8446269
## cg00000165         0.4495984         0.4861589         0.4574978
## cg00000236         0.7737180         0.7208128         0.7911344
## cg00000289         0.3970302         0.4409761         0.5921708
##            6042316042_R04C02 6042316042_R05C01 6042316042_R05C02
## cg00000029         0.7276279         0.6841344         0.5751774
## cg00000108         0.4208097         0.5038505         0.3656196
## cg00000109         0.7978687         0.7364768         0.8475675
## cg00000165         0.3708217         0.5541429         0.3247886
## cg00000236         0.6677966         0.7947917         0.7118644
## cg00000289         0.4856095         0.4814133         0.4552446
##            6042316042_R06C01 6042316042_R06C02 6042316047_R01C01
## cg00000029         0.7518702         0.6956244         0.6472324
## cg00000108         0.4780018         0.5119136         0.3472212
## cg00000109         0.6835965         0.7942457         0.7608499
## cg00000165         0.3163628         0.1497471         0.3970281
## cg00000236         0.7576389         0.7564971         0.7261444
## cg00000289         0.4499905         0.4503959         0.5091997
##            6042316047_R02C01 6042316047_R02C02 6042316047_R04C01
## cg00000029         0.5579358         0.6299049         0.7152246
## cg00000108         0.5955491         0.3167344         0.5498319
## cg00000109         0.7918387         0.8005979         0.8147258
## cg00000165         0.4178293         0.4000462         0.2776512
## cg00000236         0.7669616         0.7622299         0.7864915
## cg00000289         0.4499647         0.4949483         0.4290076
##            6042316047_R04C02 6042316047_R05C01 6042316047_R05C02
## cg00000029         0.6014030         0.6239500         0.8290006
## cg00000108         0.4334625         0.6363119         0.6462174
## cg00000109         0.8474974         0.7937128         0.7759521
## cg00000165         0.3085229         0.5111976         0.4466503
## cg00000236         0.6381011         0.7477743         0.8061760
## cg00000289         0.4773341         0.4406010         0.5034268
##            6042316047_R06C01 6055432012_R01C01 6055432012_R01C02
## cg00000029         0.6841595         0.4516399         0.4201638
## cg00000108         0.5008996         0.9509667         0.9283984
## cg00000109         0.7999119         0.6662109         0.7277319
## cg00000165         0.3772891         0.2031367         0.1590228
## cg00000236         0.7601790         0.8170757         0.8863514
## cg00000289         0.4715811         0.3529905         0.3801086
##            6055432012_R02C02 6055432012_R03C01 6055432012_R03C02
## cg00000029         0.4755556         0.5911270         0.5531213
## cg00000108         0.9632966         0.9339734         0.9690483
## cg00000109         0.7741081         0.7047394         0.7620240
## cg00000165         0.2639605         0.1703867         0.1831686
## cg00000236         0.8444044         0.8676888         0.8360814
## cg00000289         0.5282824         0.4704994         0.5514046
##            6055432012_R04C01 6055432012_R04C02 6055432012_R05C01
## cg00000029         0.6290573         0.5617760         0.5911418
## cg00000108         0.9547374         0.9555764         0.9514664
## cg00000109         0.7989008         0.7511733         0.7384962
## cg00000165         0.2410444         0.1874488         0.2245689
## cg00000236         0.8211875         0.7735187         0.8242814
## cg00000289         0.5285452         0.4507246         0.5495218
##            6055432012_R05C02 6055432012_R06C01 6055432012_R06C02
## cg00000029         0.5053841         0.5006774         0.5211055
## cg00000108         0.9467769         0.9650805         0.9567450
## cg00000109         0.7489974         0.6416233         0.6409307
## cg00000165         0.1716574         0.2594423         0.2214115
## cg00000236         0.8320933         0.8352259         0.8326686
## cg00000289         0.4929175         0.3747036         0.4595919
##            6055432029_R01C01 6055432029_R01C02 6055432029_R02C01
## cg00000029         0.4799871         0.5373293         0.5305113
## cg00000108         0.9319006         0.9609018         0.9431126
## cg00000109         0.7809393         0.6855294         0.7515199
## cg00000165         0.2473582         0.2700898         0.1767576
## cg00000236         0.7899162         0.8745755         0.8789874
## cg00000289         0.4421692         0.3741054         0.4129008
##            6055432029_R02C02 6055432029_R03C01 6055432029_R03C02
## cg00000029         0.5425801         0.5179306         0.4821117
## cg00000108         0.9417624         0.9636772         0.9634733
## cg00000109         0.7499507         0.7611667         0.6982212
## cg00000165         0.2772991         0.1679264         0.1799652
## cg00000236         0.8642089         0.8372218         0.8648562
## cg00000289         0.4984650         0.3680239         0.4584586
##            6055432029_R04C01 6055432029_R04C02 6055432029_R05C01
## cg00000029         0.4936348         0.4566294         0.5163243
## cg00000108         0.9407324         0.9532810         0.9485083
## cg00000109         0.7676887         0.7421760         0.7575966
## cg00000165         0.2329490         0.2090710         0.1811596
## cg00000236         0.8319460         0.8448426         0.8691221
## cg00000289         0.3635469         0.4452745         0.3967234
##            6055432029_R05C02 6055432029_R06C01 6055432029_R06C02
## cg00000029         0.5033495         0.6024533         0.4628428
## cg00000108         0.9476676         0.9471408         0.9372163
## cg00000109         0.6308147         0.6511048         0.6940600
## cg00000165         0.1897467         0.2498455         0.2791103
## cg00000236         0.8681307         0.8169406         0.8380174
## cg00000289         0.4281337         0.4644772         0.4836611
##            6055432060_R01C01 6055432060_R01C02 6055432060_R02C01
## cg00000029         0.4829432         0.5413782         0.4801499
## cg00000108         0.9616275         0.9469582         0.9583114
## cg00000109         0.7531995         0.7171960         0.8030402
## cg00000165         0.1680447         0.2255815         0.1924589
## cg00000236         0.8345776         0.8669872         0.8851344
## cg00000289         0.4639929         0.3713283         0.4900167
##            6055432060_R02C02 6055432060_R03C01 6055432060_R03C02
## cg00000029         0.5292569         0.5599636         0.5734813
## cg00000108         0.9426623         0.9547719         0.9587247
## cg00000109         0.7353881         0.7403170         0.7256562
## cg00000165         0.2306229         0.1677258         0.1448243
## cg00000236         0.8384684         0.8408841         0.8672923
## cg00000289         0.3810561         0.4656779         0.4489619
##            6055432060_R04C01 6055432060_R04C02 6055432060_R05C01
## cg00000029         0.4913923         0.5675564         0.5522113
## cg00000108         0.9519172         0.9650075         0.9664229
## cg00000109         0.7932175         0.7790028         0.7353540
## cg00000165         0.2016985         0.1839972         0.2426929
## cg00000236         0.8476475         0.8676987         0.8115478
## cg00000289         0.4881656         0.4330823         0.4133958
##            6055432060_R05C02 6055432060_R06C01 6055432060_R06C02
## cg00000029         0.4407743         0.5686548         0.6285834
## cg00000108         0.9650732         0.9612120         0.9212871
## cg00000109         0.7436838         0.7328813         0.5602885
## cg00000165         0.1681255         0.2302866         0.2915162
## cg00000236         0.8437584         0.8232449         0.8001906
## cg00000289         0.4142157         0.4140550         0.3817463
##            6055432066_R01C01 6055432066_R01C02 6055432066_R02C01
## cg00000029         0.4233583         0.4980459         0.5482992
## cg00000108         0.9516946         0.9470463         0.9652811
## cg00000109         0.8447866         0.7646765         0.8025670
## cg00000165         0.1423678         0.1590431         0.1745046
## cg00000236         0.8444947         0.7966481         0.8389194
## cg00000289         0.4297797         0.3594095         0.5018162
##            6055432066_R02C02 6055432066_R03C02 6055432066_R04C01
## cg00000029         0.4485561         0.5095897         0.5459673
## cg00000108         0.9529337         0.9490563         0.9465242
## cg00000109         0.8104128         0.7220734         0.7356091
## cg00000165         0.2166607         0.1685572         0.1951692
## cg00000236         0.8721994         0.8491366         0.8594021
## cg00000289         0.4111646         0.3640349         0.4128150
##            6055432066_R04C02 6055432066_R05C01 6055432066_R05C02
## cg00000029         0.4936206         0.4964094         0.5259956
## cg00000108         0.9603824         0.9468254         0.9392589
## cg00000109         0.7354827         0.7123768         0.7729230
## cg00000165         0.2220140         0.1724654         0.1540467
## cg00000236         0.8428671         0.8329177         0.8290315
## cg00000289         0.4242383         0.3090899         0.3141723
##            6057825115_R01C01 6057825115_R03C01 6057825115_R05C01
## cg00000029         0.7189976         0.6897378         0.7709326
## cg00000108         0.7416019         0.5428493         0.5434197
## cg00000109         0.6459122         0.8064851         0.8597899
## cg00000165         0.2655265         0.3313507         0.3233863
## cg00000236         0.7805707         0.6889288         0.8391550
## cg00000289         0.3178251         0.4128081         0.6593576
##            6057825115_R01C02 6057825115_R03C02 6057825115_R05C02
## cg00000029         0.7076367         0.6399105         0.6230247
## cg00000108         0.4510971         0.5852294         0.5679963
## cg00000109         0.8227988         0.8492465         0.8256986
## cg00000165         0.2725287         0.2016239         0.4615946
## cg00000236         0.7498842         0.7322501         0.7660635
## cg00000289         0.5261637         0.5488899         0.4321527
##            6057825128_R01C01 6057825128_R03C01 6057825128_R05C01
## cg00000029         0.6717592         0.7522418         0.6347535
## cg00000108         0.4739439         0.5572308         0.2463536
## cg00000109         0.7746812         0.8341231         0.8013856
## cg00000165         0.4769035         0.2791574         0.2970243
## cg00000236         0.6723752         0.7882433         0.7398250
## cg00000289         0.4233790         0.5636732         0.4627684
##            6057825115_R02C01 6057825115_R04C01 6057825115_R06C01
## cg00000029         0.6501169         0.6825724         0.6243567
## cg00000108         0.5345400         0.5421511         0.5625968
## cg00000109         0.8547748         0.8693722         0.8103014
## cg00000165         0.6520707         0.2469365         0.4594300
## cg00000236         0.7317447         0.7450085         0.7317146
## cg00000289         0.4532717         0.5568660         0.5413773
##            6057825115_R02C02 6057825115_R04C02 6057825115_R06C02
## cg00000029         0.6396068         0.7773322         0.6525435
## cg00000108         0.4149925         0.3233349         0.4557251
## cg00000109         0.8611083         0.8151881         0.8655906
## cg00000165         0.5498851         0.5053686         0.2940515
## cg00000236         0.7517095         0.7503832         0.7570359
## cg00000289         0.4553110         0.5178411         0.5250577
##            6057825128_R01C02 6057825128_R03C02 6057825128_R05C02
## cg00000029         0.6787162         0.6785583         0.6508382
## cg00000108         0.5034493         0.4944187         0.4808232
## cg00000109         0.8463614         0.8409880         0.7539724
## cg00000165         0.1800716         0.2454531         0.6862180
## cg00000236         0.6689967         0.7462064         0.7260894
## cg00000289         0.3853612         0.5304513         0.3858073
```

```r
# save normalized data
GSE43414_BMIQ <- betas_norm.fin
save(GSE43414_BMIQ, file = "GSE43414_BMIQ.RData")
```

## Comparing raw and normalized datasets   
NOTE: The loop in this in this code to add probe type to betas.plot appears to run infinitely, or at least for days. After letting it run a long time then manually stopping it, it still created the correct annotations, allowing us to move forward with plotting. However it will not knit into an .md without getting stuck on this loop again, so the evaluation has been set to FALSE for the purpose of GitHub. The resultant plot can be found [here](https://github.com/STAT540-UBC/team_Methylhomies/blob/master/Images/raw_normalized_ggplot.png).

For all samples:

```r
## use ggplot2 to plot density plots of raw and normalized data
# extract beta matrix for raw and normalized datasets
betas.raw <- as.matrix(GSE43414_filtered)
betas.norm <- GSE43414_BMIQ

# randomly sample 10,000 probes from the datasets (same random probes in each dataset)
random.probes<-sample(1:nrow(betas.raw), 10000)
betas.raw.subset<-betas.raw[random.probes,]
betas.norm.subset<-betas.norm[random.probes,]

#Create data frame with probe type annotations
type <- as.character(fData(Priest)$INFINIUM_DESIGN_TYPE)
probe <- as.character(fData(Priest)$NAME)
probe.ann <- as.data.frame(cbind(probe, type))
probe.ann$type <- as.character(probe.ann$type)
probe.ann$probe <- as.character(probe.ann$probe)

# add probe type variable to raw betas
betas.raw.subset <- as.data.frame(betas.raw.subset)
for (i in 1:nrow(betas.raw.subset)){
  for (j in 1:nrow(probe.ann)){
    if (rownames(betas.raw.subset)[i] == #probe.ann$probe[j]){
      betas.raw.subset$type[i] <- probe.ann$type[j]
    }}}

# add probe type variable to normalized betas
betas.norm.subset <- as.data.frame(betas.norm.subset)
for (i in 1:nrow(betas.norm.subset)){
  for (j in 1:nrow(probe.ann)){
    if (rownames(betas.norm.subset)[i] == probe.ann$probe[j]){
      betas.norm.subset$type[i] <- probe.ann$type[j]
    }}}

# melt each dataset
head(betas.raw.melt<-melt(betas.raw.subset))
head(betas.norm.melt<-melt(betas.norm.subset))
# remove NAs
betas.raw.melt.clean<-betas.raw.melt[which(!(is.na(betas.raw.melt$value))),]
betas.norm.melt.clean<-betas.norm.melt[which(!(is.na(betas.norm.melt$value))),]
# add descriptor for each datatype and add meta data to each dataset
betas.raw.melt.clean$Data<-"Raw"
betas.norm.melt.clean$Data<-"Normalized"
head(betas.raw.melt.clean)

meta.raw <- Brain_matched
meta.norm <- Brain_matched
head(betas.raw.plot<-merge(betas.raw.melt.clean, meta.raw, by.x="variable", by.y="barcode"))
head(betas.norm.plot<-merge(betas.norm.melt.clean, meta.norm, by.x="variable", by.y="barcode"))

# combine both datasets
betas.plot<-rbind(betas.raw.plot, betas.norm.plot)
betas.plot$type<-factor(betas.plot$type, levels=c("I", "II"))
betas.plot$Data<-factor(betas.plot$Data, levels=c("Raw", "Normalized"))

# plot density plots
betas.plot$type <- as.character(betas.plot$type)
ggplot(betas.plot, aes(value, group=Subject, colour=type)) + geom_density() + theme_bw() + facet_wrap(~Data) + scale_color_manual(values = c("blue","red"))
```

We can see that there is improved peak-to-peak overlap of the Type 1 and Type 2 probes with the final normalized dataset as desired (ie Type 1 and 2 probes have more similar dynamic ranges for their beta values). 

We will continue on with our analysis using the BMIQ-normalized datasets.
