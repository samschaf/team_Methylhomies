---
output: 
  html_document: 
    keep_md: yes
---
Hannon et al. (2017) 450K Data Principal Component Analysis
========================================================
## Original author: Sumaiya Islam
## Updated by: Samantha Schaffner
## Date updated: March 13, 2017
  
### Script contents:
  - detection and correction for technical batch variation using PCA and ComBat, respectively, of post-mortem human brain samples analyzed by Illumina HM450K platform from Jonathan Mill's research group (PMC4844197). 
   
### A. Set up working directory & packages

R version 3.2.3 (2015-12-10)

We will initially set our working directory and load our libraries.



## Heat scree plot Function

```r
### Function of association meta variable with PC (ANOVA)
heat_scree_plot<-function(Loadings, Importance, Num, Order){
  adjust<-1-Importance[1]
  pca_adjusted<-Importance[2:length(Importance)]/adjust
  pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))
  
  scree<-ggplot(pca_df[which(pca_df$PC<Num),],aes(PC,adjusted_variance))+geom_bar(stat = "identity",color="black",fill="grey")+theme_bw()+
        theme(axis.text = element_text(size =12),
              axis.title = element_text(size =15),
              plot.margin=unit(c(1,1.5,0.2,2.25),"cm"))+ylab("Variance")+
    scale_x_continuous(breaks = seq(1,Num,1))
  
  #### Heat
  ## correlate meta with PCS
  ## Run anova of each PC on each meta data variable

  aov_PC_meta<-lapply(1:ncol(meta_categorical), function(covar) sapply(1:ncol(Loadings), function(PC) summary(aov(Loadings[,PC]~meta_categorical[,covar]))[[1]]$"Pr(>F)"[1]))
  cor_PC_meta<-lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Loadings), function(PC) (cor.test(Loadings[,PC],as.numeric(meta_continuous[,covar]),alternative = "two.sided", method="spearman", na.action=na.omit, exact=FALSE)$p.value)))
  names(aov_PC_meta)<-colnames(meta_categorical)
  names(cor_PC_meta)<-colnames(meta_continuous)
  aov_PC_meta<-do.call(rbind, aov_PC_meta)
  cor_PC_meta<-do.call(rbind, cor_PC_meta)
 aov_PC_meta<-rbind(aov_PC_meta, cor_PC_meta)
  aov_PC_meta<-as.data.frame(aov_PC_meta)
  #adjust
  aov_PC_meta_adjust<-aov_PC_meta[,2:ncol(aov_PC_meta)]
  
    
  #reshape
  avo<-aov_PC_meta_adjust[,1:(Num-1)]
  avo_heat_num<-apply(avo,2, as.numeric)
  avo_heat<-as.data.frame(avo_heat_num)
  colnames(avo_heat)<-sapply(1:(Num-1), function(x) paste("PC",x, sep=""))
  avo_heat$meta<-rownames(avo)
  avo_heat_melt<-melt(avo_heat, id=c("meta"))
  
  # cluster meta data
  ord <- Order
  meta_var_order<-unique(avo_heat_melt$meta)[rev(ord)]
  avo_heat_melt$meta <- factor(avo_heat_melt$meta, levels = meta_var_order)
  
  # color if sig
   avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]>=0.9){">=0.9"}else{
    if(avo_heat_melt$value[x]>=0.5){">=0.5"}else{
      if(avo_heat_melt$value[x]>=0.1){">=0.1"}else{"<0.1"}}})
  avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]<=0.001){"<=0.001"}else{
     if(avo_heat_melt$value[x]<=0.01){"<=0.01"}else{
       if(avo_heat_melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
  
  heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
  geom_tile(color = "black",size=0.5) +
  theme_gray(8)+scale_fill_manual(values=c("#084594","#4292c6","#9ecae1","#deebf7"))+
      theme(axis.text = element_text(size =10, color="black"),
            axis.text.x = element_text(),
          axis.title = element_text(size =15),
          legend.text = element_text(size =14),
          legend.title = element_text(size =12),
          legend.position = c(1, 0), legend.justification = c(1,0),
          plot.margin=unit(c(0,2.25,1,1),"cm"))+
    xlab("Principal Component")+ylab(NULL)
  
  grid.arrange(scree, heat, ncol=1)
}
```


### B. Load files

#### We will be analyzing the normalized and filtered Hannon et al. dataset
First, load required files and reshape meta data:

```r
load("GSE43414_BMIQ.RData") # normalized beta values
ncol(GSE43414_BMIQ) #432
```

```
## [1] 432
```

```r
load("Brain_meta_matched_GSE43414.RData") # associated meta data
nrow(Brain_matched) #441
```

```
## [1] 441
```

```r
#There are some individuals in the meta data not included in the BMIQ data. We need to re-match the betas and meta data.
for (j in 1:nrow(Brain_matched)){
  if (!(Brain_matched$barcode[j] %in% colnames(GSE43414_BMIQ))){
  Brain_matched <- Brain_matched[-j,]  
 }}
nrow(Brain_matched) #434
```

```
## [1] 434
```

```r
length(unique(Brain_matched$barcode)) #434
```

```
## [1] 434
```

```r
length(unique(colnames(GSE43414_BMIQ))) #432
```

```
## [1] 432
```

```r
#Two individuals are still in the meta data, not present in BMIQ data.
#Rearrange to look at barcodes of both datasets
Brain_matched <- Brain_matched %>% arrange(barcode)
GSE.t <- as.data.frame(t(GSE43414_BMIQ))
GSE.t$barcode <- rownames(GSE.t)
GSE.t <- GSE.t %>% arrange(barcode)
for (j in 1:nrow(Brain_matched)){
  if (!(Brain_matched$barcode[j] %in% GSE.t$barcode)){
  Brain_matched <- Brain_matched[-j,]  
  }}
nrow(Brain_matched) #432 - all good!
```

```
## [1] 432
```

```r
save(Brain_matched, file="Brain_matched.RData")

cell.proportions<-read.csv("cellprop_uncor.csv", header=TRUE, row.names=1) # predicted neuron and glial cell proportions based on CETS
head(cell.proportions)
```

```
##                      neuron      glia
## 6057825008_R02C02 0.4978760 0.5021240
## 6057825008_R03C01 0.4990904 0.5009096
## 6057825008_R04C01 0.4928983 0.5071017
## 6057825008_R04C02 0.4976021 0.5023979
## 6057825008_R05C01 0.4832074 0.5167926
## 6057825008_R05C02 0.5175410 0.4824590
```

```r
# check for NAs in data
ind<-is.row.na(GSE43414_BMIQ) # The function returns a vector of logical variables, one for each row of the matrix. The variable is TRUE if the row does not contain any missing values and FAlSE otherwise.
length(na.count<-which(ind=="FALSE")) # 76545 rows contain NAs
```

```
## [1] 76545
```

```r
GSE43414_BMIQ <- na.omit(GSE43414_BMIQ) #Remove probes with NAs
dim(GSE43414_BMIQ) #[1] 338359    432
```

```
## [1] 338359    432
```

```r
#338359 probes left
save(GSE43414_BMIQ, file="GSE43414_BMIQ_na.RData")

uncor.dat<-GSE43414_BMIQ
meta<-Brain_matched

#Restructure meta data and cell proportion data so sample order matches
meta<- meta %>% arrange(barcode)
cell.proportions$barcode <- rownames(cell.proportions)
cell.proportions<- cell.proportions %>% arrange(barcode)

#Add cell proportion information to meta data
identical(cell.proportions$barcode, meta$barcode) # TRUE
```

```
## [1] TRUE
```

```r
meta$Neuron<-as.numeric(cell.proportions$neuron)
meta$Glia<-as.numeric(cell.proportions$glia)

#Add row and chip information to metadata; this can be found in the sample barcodes
for (i in 1:nrow(meta)){
  meta$chip[i]<-paste(substr(meta$barcode[i], start=1, stop=10))
  meta$row[i]<-paste(substr(meta$barcode[i], start=13, stop=14))
}

meta$age.brain <- as.numeric(meta$age.brain)
```

```
## Warning: NAs introduced by coercion
```

```r
meta$braak.stage <- as.numeric(meta$braak.stage)
```

```
## Warning: NAs introduced by coercion
```

```r
str(meta)
```

```
## 'data.frame':	432 obs. of  16 variables:
##  $ series_id        : chr  "GSE43414" "GSE43414" "GSE43414" "GSE43414" ...
##  $ gsm              : chr  "GSM1069412" "GSM1069413" "GSM1069414" "GSM1069415" ...
##  $ Subject          : chr  "NA" "NA" "NA" "NA" ...
##  $ barcode          : chr  "6042316024_R01C01" "6042316024_R01C02" "6042316024_R02C01" "6042316024_R02C02" ...
##  $ lunnon.et.al     : chr  "FALSE" "FALSE" "FALSE" "FALSE" ...
##  $ tissue.code      : chr  "NA" "NA" "NA" "NA" ...
##  $ braak.stage      : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ Sex              : chr  "NA" "NA" "NA" "NA" ...
##  $ ad.disease.status: chr  "NA" "NA" "NA" "NA" ...
##  $ age.brain        : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ age.blood        : chr  "NA" "NA" "NA" "NA" ...
##  $ Tissue           : chr  "cerebellum" "cerebellum" "cerebellum" "cerebellum" ...
##  $ Neuron           : num  0.493 0.515 0.503 0.48 0.513 ...
##  $ Glia             : num  0.507 0.485 0.497 0.52 0.487 ...
##  $ chip             : chr  "6042316024" "6042316024" "6042316024" "6042316024" ...
##  $ row              : chr  "01" "01" "02" "02" ...
```

```r
save(meta, file="Meta_uncor.RData")
```


## PCA Scree Heatmap for uncorrected data


```r
## PCA
uncor.dat <- t(scale(t(as.matrix(uncor.dat))))
PCA_full<-princomp(uncor.dat[complete.cases(uncor.dat),])
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
adjust<-1-Importance[1]
pca_adjusted<-Importance[2:length(Importance)]/adjust
pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))

#Specify which covariates are categorical and/or categorical
colnames(meta)
```

```
##  [1] "series_id"         "gsm"               "Subject"          
##  [4] "barcode"           "lunnon.et.al"      "tissue.code"      
##  [7] "braak.stage"       "Sex"               "ad.disease.status"
## [10] "age.brain"         "age.blood"         "Tissue"           
## [13] "Neuron"            "Glia"              "chip"             
## [16] "row"
```

```r
meta_categorical<-meta[,c("ad.disease.status", "braak.stage", "Sex", "Tissue", "chip", "row")]  # input column numbers in meta that contain categorical variables
meta_continuous<-meta[,c("age.brain", "Neuron")] # input column numbers in meta that contain continuous variables
#meta_continuous<-data.frame(meta_continuous)

# Specify the number of PCs you want shown (usually # of samples in the dataset)
Num<-20

# Designate what order you want the variables to appear (continuous variables rbinded to categorical variables in function)
Order<-c(4,8,7,3,1,2,5,6)

#Apply function on PCA results, pulls in the meta data and beta values from above
heat_scree_plot(Loadings, Importance, Num, Order)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

The main contributors to variance in this data are chip, tissue, neuronal proportion, sex, and AD disease status. Brain region, AD status, and sex will be used as covariates in a linear regression model, while the large chip effect will be batch-corrected using ComBat.

## Batch correction using ComBat

ComBat is a function included in the SVA (surrogate variable analysis) package ((Johnson et al., 2007))[https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxj037]. It uses empirical Bayesian adjustment to correct for known sources of batch variation. Correction is usually performed first on either the variable which contributes more to overall variance or the variable with fewer batches -- row satisfies both of these requirements, so we will correct for row, followed by chip.


```r
#The following packages must be unloaded in order to use sva:
detach("package:wateRmelon", unload=TRUE)
detach("package:IlluminaHumanMethylation450kanno.ilmn12.hg19", unload=TRUE)
detach("package:lumi", unload=TRUE)
unloadNamespace("methylumi")
unloadNamespace("minfi")
unloadNamespace("mgcv")
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
## The following object is masked from 'package:dplyr':
## 
##     collapse
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
#Correction for row
row <- meta$row
modcombat <- model.matrix(~1,data=meta)
combat_edata <- ComBat(dat=GSE43414_BMIQ, batch=row, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
```

```
## Found 6 batches
## Adjusting for 0 covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r
#Correction for chip
chip <- meta$chip
GSE43414_batch_cor <- ComBat(dat=combat_edata, batch=chip, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
```

```
## Found 50 batches
## Adjusting for 0 covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r
save(GSE43414_batch_cor, file="GSE43414_batch_cor.RData")
```

## Predict cell proportions for batch-corrected data


```r
library(cets)
# load "brain dataset" from data file in cetsBrain
load("~/team_Methylhomies/cetsBrain.rda") # click on cetsBrain.rda file to place in workspace
dim(brain)
```

```
## [1] 10000   146
```

```r
brain[1:3, 1:4]
```

```
##            X1_5175.G.P1A1._7766130090_R01C01
## cg02689072                         0.9271406
## cg12093060                         0.8889666
## cg05940691                         0.9427095
##            X7_5175.N.P1A7._7766130090_R01C02
## cg02689072                         0.1190979
## cg12093060                         0.1356233
## cg05940691                         0.1757897
##            X2_813.N.P1A2._7766130090_R02C01
## cg02689072                        0.1258751
## cg12093060                        0.1268435
## cg05940691                        0.1476256
##            X8_1740.N.P1A8._7766130090_R02C02
## cg02689072                         0.1372099
## cg12093060                         0.1110450
## cg05940691                         0.1692726
```

```r
head(pdBrain)
```

```
##   ID1 celltype       diag    sex ethnicity age batch row array PMI
## 1   1        G Depression   Male Caucasian  47     1   A     1  22
## 2   2        N Depression   Male Caucasian  47     1   A     1  22
## 3   3        N    Control Female Caucasian  30     1   A     1  14
## 4   4        N    Control Female   African  13     1   A     1  17
## 5   5        N Depression Female   African  14     1   A     1  15
## 6   6        G Depression Female   African  14     1   A     1  15
```

Create the neuron and glia reference profiles:


```r
modelIdx <- list(neuron = pdBrain$celltype == "N", glia = pdBrain$celltype ==  "G")
 # getReference returns a 2-column matrix, representing reference profiles for the two cell types.
refProfile <- getReference(brain, modelIdx)
head(refProfile)
```

```
##               neuron       glia
## cg02689072 0.1256035 0.93511026
## cg12093060 0.1425625 0.88545915
## cg05940691 0.1739833 0.92097174
## cg05403655 0.1380149 0.91071764
## cg05699921 0.7125904 0.08176213
## cg00968638 0.8164890 0.07838638
```

#### For the brain datasets

Estimate the neuronal proportion:

The estProportion function returns an estimate of the percentage of cell type in the first column of its profile argument (neurons in this case). 

```r
prop <- estProportion(GSE43414_batch_cor, profile = refProfile)
prop<-as.data.frame(prop)
prop$glia<-apply(prop,1,function(x) 1-x)
colnames(prop)<- c("neuron", "glia")
head(prop)
```

```
##                      neuron      glia
## 6057825008_R02C02 0.4151500 0.5848500
## 6057825008_R03C01 0.4171250 0.5828750
## 6057825008_R04C01 0.3922126 0.6077874
## 6057825008_R04C02 0.3948113 0.6051887
## 6057825008_R05C01 0.3847930 0.6152070
## 6057825008_R05C02 0.4264143 0.5735857
```

```r
write.csv(prop, file = "cellprop_batch_cor.csv", row.names=T)
summary(prop)
```

```
##      neuron            glia       
##  Min.   :0.1223   Min.   :0.4437  
##  1st Qu.:0.3599   1st Qu.:0.5705  
##  Median :0.4014   Median :0.5986  
##  Mean   :0.3920   Mean   :0.6080  
##  3rd Qu.:0.4295   3rd Qu.:0.6401  
##  Max.   :0.5563   Max.   :0.8777
```

```r
plot(density(prop$neuron), main="Neuronal Proportion Density") 
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)


```r
#Restructure meta data and cell proportion data so sample order matches
cell.proportions <- prop
cell.proportions$barcode <- rownames(cell.proportions)
cell.proportions<- cell.proportions %>% arrange(barcode)

#Add cell proportion information to meta data
identical(cell.proportions$barcode, meta$barcode) # TRUE
```

```
## [1] TRUE
```

```r
meta$Neuron<-as.numeric(cell.proportions$neuron)
meta$Glia<-as.numeric(cell.proportions$glia)

#Save corrected meta data
save(meta, file="Meta_batch_cor.RData")
```

## PCA Scree Heatmap for batch-corrected data


```r
## PCA
load("Meta_batch_cor.RData")
load("GSE43414_batch_cor.RData")
batch.cor.dat <- t(scale(t(as.matrix(GSE43414_batch_cor))))
PCA_full<-princomp(batch.cor.dat[complete.cases(batch.cor.dat),])
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
adjust<-1-Importance[1]
pca_adjusted<-Importance[2:length(Importance)]/adjust
pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))

#Specify which covariates are categorical and/or categorical
colnames(meta)
```

```
##  [1] "series_id"         "gsm"               "Subject"          
##  [4] "barcode"           "lunnon.et.al"      "tissue.code"      
##  [7] "braak.stage"       "Sex"               "ad.disease.status"
## [10] "age.brain"         "age.blood"         "Tissue"           
## [13] "Neuron"            "Glia"              "chip"             
## [16] "row"
```

```r
meta_categorical<-meta[,c("ad.disease.status", "braak.stage", "Sex", "Tissue", "chip", "row")]  # input column numbers in meta that contain categorical variables
meta_continuous<-meta[,c("age.brain", "Neuron")] # input column numbers in meta that contain continuous variables
#meta_continuous<-data.frame(meta_continuous)

# Specify the number of PCs you want shown (usually # of samples in the dataset)
Num<-20

# Designate what order you want the variables to appear (continuous variables rbinded to categorical variables in function)
Order<-c(4,8,7,3,1,2,5,6)

#Apply function on PCA results, pulls in the meta data and beta values from above
heat_scree_plot(Loadings, Importance, Num, Order)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

Batch correction removed much of the initial variation, leaving only Neuron significantly associated. 

We will now perform cell-type correction based on the neuronal/glial proportions.


```r
all(rownames(prop)%in%colnames(GSE43414_batch_cor))
```

```
## [1] TRUE
```

```r
brain.cor.dat<- as.data.frame(GSE43414_batch_cor)

# fit methylation data for each probe in the dataset by the neuronal proportion
avebeta.lm<-apply(brain.cor.dat, 1, function(x){
  brain.sub<-prop[colnames(brain.cor.dat),]
  lm(x~neuron,data=brain.sub)
})

# obtain residuals for each probe across all samples (as a matrix)
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
head(residuals)
```

```
##            6057825008_R02C02 6057825008_R03C01 6057825008_R04C01
## cg00000029       0.009644362      6.904700e-02      -0.014971905
## cg00000108      -0.118312012      3.190453e-02       0.022485001
## cg00000165       0.042555314     -9.206013e-02      -0.008113457
## cg00000236       0.041324336      9.585216e-05       0.035291721
## cg00000289       0.046021847     -2.942238e-02      -0.011858657
## cg00000292      -0.027265001     -3.037839e-02       0.069229590
##            6057825008_R04C02 6057825008_R05C01 6057825008_R05C02
## cg00000029      -0.003094615        0.02417504       -0.06056986
## cg00000108       0.083495016       -0.02696131       -0.04599481
## cg00000165      -0.124498985        0.09227951        0.06749207
## cg00000236      -0.036138975        0.01935169       -0.06303527
## cg00000289      -0.061816231       -0.10519206        0.01980000
## cg00000292      -0.042933891        0.01030051       -0.06657292
##            6057825008_R06C01 6057825008_R06C02 6057825014_R01C02
## cg00000029        0.01936911       0.046413087       -0.02039818
## cg00000108       -0.03591174      -0.005362574       -0.08358714
## cg00000165        0.07375497      -0.098536744       -0.08323761
## cg00000236       -0.04819732       0.025184442       -0.01333449
## cg00000289        0.08501055       0.051056504        0.04014438
## cg00000292        0.01179781       0.101068180        0.01106595
##            6057825014_R02C02 6057825014_R03C01 6057825014_R03C02
## cg00000029      0.0057588048       0.006545043       -0.02959152
## cg00000108      0.1501256552      -0.050148103       -0.21363214
## cg00000165      0.1130681150       0.105881599        0.07603309
## cg00000236     -0.0136758379       0.036063378        0.03927211
## cg00000289      0.0000632328      -0.022120213        0.01059506
## cg00000292     -0.0175399881      -0.025224850        0.06201372
##            6057825014_R04C01 6057825014_R04C02 6057825014_R06C01
## cg00000029        0.10682675       -0.04161019      -0.099944249
## cg00000108       -0.01831239       -0.03169876      -0.002975556
## cg00000165       -0.05883330        0.14593913       0.059689078
## cg00000236       -0.05525200       -0.01306339      -0.096222774
## cg00000289       -0.02183546       -0.04565986       0.029602540
## cg00000292        0.02654973       -0.03994406       0.001104363
##            6057825017_R01C01 6057825017_R01C02 6057825017_R02C02
## cg00000029        0.08823746      -0.030182217       0.041259217
## cg00000108       -0.12867192       0.017279882       0.079964032
## cg00000165       -0.11640288       0.021967525      -0.009381011
## cg00000236        0.03307804       0.021514214       0.026074173
## cg00000289        0.10294602      -0.064520244       0.024166433
## cg00000292       -0.05171741       0.008951448       0.033122107
##            6057825017_R03C01 6057825017_R03C02 6057825017_R04C01
## cg00000029      -0.000800335      1.414351e-02       0.000959136
## cg00000108       0.028351324     -1.322579e-01       0.073050097
## cg00000165      -0.001144441      1.920127e-01       0.040261319
## cg00000236      -0.001028064     -5.279127e-05      -0.008949004
## cg00000289       0.085847368     -1.827926e-02      -0.034743735
## cg00000292      -0.027791082     -6.322804e-03      -0.091071211
##            6057825017_R04C02 6057825017_R05C01 6057825017_R05C02
## cg00000029      -0.007663867        0.04252071       0.073509606
## cg00000108       0.140879435        0.04994966      -0.016907338
## cg00000165      -0.156420506       -0.07548870      -0.088009116
## cg00000236       0.041570996       -0.01821487      -0.012976466
## cg00000289      -0.018261465       -0.01987637       0.066037834
## cg00000292       0.032744778        0.03113830       0.009316313
##            6057825017_R06C01 6057825018_R01C01 6057825018_R02C02
## cg00000029       0.083577944        0.04295547     -2.875974e-02
## cg00000108       0.146655917       -0.15172896     -5.956180e-02
## cg00000165       0.037679035        0.04323564      1.311768e-01
## cg00000236       0.005315443        0.01943276      3.705630e-02
## cg00000289       0.008963656       -0.03133623     -6.402540e-06
## cg00000292      -0.049964040        0.02712919      6.333376e-04
##            6057825018_R03C02 6057825018_R04C01 6057825018_R04C02
## cg00000029      -0.074592322       -0.04641627       -0.04441689
## cg00000108      -0.009740487       -0.02400726       -0.10789063
## cg00000165       0.142852974       -0.11743954       -0.01600545
## cg00000236      -0.001565922       -0.06412639        0.03425726
## cg00000289      -0.005288700       -0.01686963        0.05147007
## cg00000292      -0.013576908       -0.05582052       -0.03103015
##            6057825018_R05C01 6057825018_R05C02 6057825018_R06C01
## cg00000029        0.07394593       0.044094908       0.017186157
## cg00000108        0.10646689      -0.034680453       0.076602952
## cg00000165       -0.08642044      -0.083225155      -0.095732462
## cg00000236       -0.03039786       0.007698512      -0.002126411
## cg00000289       -0.05447005       0.058045948       0.049430821
## cg00000292        0.09382786      -0.016485152       0.020613985
##            6042316085_R01C01 6042316085_R01C02 6042316085_R02C02
## cg00000029       -0.09677134        0.05506534        0.15114168
## cg00000108       -0.22151114        0.02076482        0.05358784
## cg00000165        0.22253236       -0.01190461       -0.01341539
## cg00000236        0.05488996        0.02684334       -0.04945051
## cg00000289        0.04710567        0.03289837       -0.14382126
## cg00000292        0.01358481       -0.08189642        0.10691100
##            6042316085_R03C01 6042316085_R03C02 6042316085_R05C01
## cg00000029     -0.0257830177      -0.010792091      -0.048165834
## cg00000108     -0.0100634479      -0.019653020      -0.024135433
## cg00000165     -0.0521856243      -0.032787528       0.005388742
## cg00000236     -0.0273303906      -0.011544614      -0.049176341
## cg00000289     -0.0001425577      -0.002248039      -0.009564861
## cg00000292     -0.0084973151      -0.007218871      -0.033597518
##            6042316085_R05C02 6042316085_R06C01 6042316107_R01C02
## cg00000029      -0.020090081        0.03562160       -0.02351149
## cg00000108      -0.096403809        0.14765087       -0.09827285
## cg00000165      -0.022219003        0.12130073       -0.06363655
## cg00000236       0.053655333       -0.03202809       -0.01241425
## cg00000289      -0.009399374       -0.01103916        0.04011754
## cg00000292      -0.045160948        0.03551879        0.05036675
##            6042316107_R03C01 6042316107_R04C02 6042316107_R05C01
## cg00000029        0.02128861       -0.01122467      -0.003163051
## cg00000108        0.00387673        0.01363541       0.027543978
## cg00000165        0.13600029       -0.03927648      -0.034251420
## cg00000236       -0.01786251        0.01458619       0.074433497
## cg00000289        0.04840893       -0.02376664      -0.055641717
## cg00000292       -0.02898671       -0.03452141       0.024658061
##            6042316107_R05C02 6042316107_R06C01 6042316107_R06C02
## cg00000029      -0.005155409       0.075628692        0.06688841
## cg00000108      -0.030007982      -0.173990880       -0.11563363
## cg00000165      -0.115356058      -0.087403595       -0.05532614
## cg00000236      -0.035141032      -0.036596160       -0.05674156
## cg00000289       0.012691312       0.009169751       -0.04855525
## cg00000292      -0.024000769       0.029949721       -0.01619001
##            6042316113_R01C01 6042316113_R01C02 6042316113_R02C01
## cg00000029       0.008492966      -0.044422422        0.04927155
## cg00000108       0.068095230      -0.079734883        0.15560964
## cg00000165       0.090935086      -0.005390018       -0.01101530
## cg00000236      -0.005448572       0.028076279       -0.00117995
## cg00000289       0.063578412       0.002491951        0.08175850
## cg00000292      -0.032963036      -0.012396969       -0.07324259
##            6042316113_R02C02 6042316113_R03C02 6042316113_R04C02
## cg00000029      -0.088329496       0.043044377       -0.17769236
## cg00000108       0.117934201       0.005788785        0.20191673
## cg00000165       0.196121135       0.098061506       -0.14372126
## cg00000236       0.009778202       0.026071996        0.05874406
## cg00000289      -0.080405236      -0.032601742       -0.11882824
## cg00000292       0.112199291      -0.000819046        0.10191454
##            6042316113_R05C01 6042316113_R05C02 6042316113_R06C01
## cg00000029       0.079381880        0.01792407        0.02133714
## cg00000108      -0.037955176       -0.03560058       -0.03916770
## cg00000165      -0.057373438        0.10192038       -0.01267612
## cg00000236       0.017005060       -0.02438799       -0.03155621
## cg00000289       0.009223757        0.05467033        0.03690424
## cg00000292       0.010265346        0.03180345        0.03305237
##            6042316127_R01C01 6042316127_R01C02 6042316127_R02C02
## cg00000029        0.03189625        0.04294542       0.029858708
## cg00000108       -0.07916364       -0.11732431      -0.079247698
## cg00000165       -0.03671920        0.19735969       0.041588316
## cg00000236       -0.04914707       -0.01533972       0.017034895
## cg00000289       -0.03360330        0.08640729      -0.004151747
## cg00000292       -0.06480635       -0.05226464      -0.044768705
##            6042316127_R03C01 6042316127_R03C02 6042316127_R04C01
## cg00000029      0.0838819569       0.020465784        0.04642471
## cg00000108     -0.0224623985      -0.062682459        0.01977116
## cg00000165      0.0867482827       0.002725208       -0.01456110
## cg00000236      0.0084001310      -0.036581352       -0.04077012
## cg00000289     -0.0003779837       0.057316595       -0.02460193
## cg00000292     -0.0283767421       0.036172921       -0.02170709
##            6042316127_R04C02 6042316127_R05C02 6042316127_R06C01
## cg00000029      -0.016951086      -0.066539840      -0.042651361
## cg00000108       0.122371082      -0.104707496       0.044392190
## cg00000165       0.125680523      -0.171504901      -0.018121350
## cg00000236      -0.021150935      -0.036870311       0.006590522
## cg00000289      -0.059204631       0.008259659       0.049708158
## cg00000292       0.008676657      -0.029949917      -0.041136340
##            6057825014_R01C02.1 6057825014_R02C02.1 6057825014_R03C01.1
## cg00000029         0.032001527          0.05554062         0.079702396
## cg00000108        -0.042422817          0.09738808        -0.020170886
## cg00000165        -0.077940930          0.11268049         0.090782596
## cg00000236         0.002349947          0.03246560         0.051789357
## cg00000289         0.014535134         -0.01270106         0.040712163
## cg00000292         0.040175698         -0.04683235         0.007043674
##            6057825014_R03C02.1 6057825014_R04C01.1 6057825014_R04C02.1
## cg00000029         -0.02188395         0.116955306        -0.076784345
## cg00000108         -0.15514232        -0.023667772         0.021511036
## cg00000165          0.04352850        -0.091973428         0.150138964
## cg00000236          0.04336211        -0.042988563        -0.021600251
## cg00000289         -0.04149515        -0.077481122        -0.002275084
## cg00000292          0.03712512        -0.008500597        -0.012890504
##            6057825014_R06C01.1 6057825017_R01C01.1 6057825017_R01C02.1
## cg00000029         -0.07485323          0.06823091         -0.05109762
## cg00000108          0.03299839         -0.15475853          0.09188603
## cg00000165          0.04115494         -0.14776571         -0.03280392
## cg00000236         -0.08298449          0.01821688          0.04178411
## cg00000289         -0.01479691          0.03861281         -0.02201402
## cg00000292         -0.01308900         -0.02921541         -0.01050867
##            6057825017_R02C01 6057825017_R02C02.1 6057825017_R03C01.1
## cg00000029       0.043062844         0.040332763          0.01112531
## cg00000108      -0.024319421         0.120516092          0.07105877
## cg00000165       0.091410139        -0.044786953         -0.05409374
## cg00000236      -0.001594579         0.019592592          0.02873985
## cg00000289      -0.030529159         0.032904223          0.10555945
## cg00000292       0.035028801         0.006788925          0.01661545
##            6057825017_R03C02.1 6057825017_R04C01.1 6057825017_R04C02.1
## cg00000029        -0.025940424          0.03005037         0.031105299
## cg00000108        -0.169655783         -0.01627202         0.023583043
## cg00000165         0.146789320          0.14784504        -0.090936203
## cg00000236         0.040020962         -0.03675623         0.019849848
## cg00000289        -0.011859851          0.03259656        -0.006912628
## cg00000292         0.008130492         -0.05633568        -0.006356663
##            6057825017_R05C01.1 6057825017_R05C02.1 6057825017_R06C01.1
## cg00000029         0.071661533        0.0587050832         0.053457160
## cg00000108        -0.013543968       -0.0375840684         0.018893951
## cg00000165        -0.029800229       -0.0542438210         0.039272970
## cg00000236         0.008358139        0.0005958356         0.007193253
## cg00000289         0.028531207       -0.0022597906        -0.004401384
## cg00000292        -0.006802017       -0.0091894385        -0.047047797
##            6057825018_R01C01.1 6057825018_R02C02.1 6057825018_R03C01
## cg00000029         0.034948371        -0.009895674       -0.17943371
## cg00000108        -0.125108314        -0.072408506        0.21797938
## cg00000165        -0.004379867         0.086733633       -0.14192891
## cg00000236         0.019628481         0.045765207       -0.02789321
## cg00000289        -0.008184222        -0.037761372       -0.02344118
## cg00000292         0.003218503        -0.007993392        0.13911313
##            6057825018_R03C02.1 6057825018_R04C01.1 6057825018_R04C02.1
## cg00000029         -0.01589564          0.03754441          0.03134609
## cg00000108         -0.03992778         -0.10643467         -0.15170888
## cg00000165          0.08115366         -0.02318111          0.06289140
## cg00000236         -0.02970748         -0.04914377          0.02085616
## cg00000289          0.04430024         -0.03352143          0.04314305
## cg00000292         -0.01660227         -0.06716394         -0.03492891
##            6057825018_R05C01.1 6057825018_R05C02.1 6057825018_R06C01.1
## cg00000029          0.09937504         0.081390003          0.07289304
## cg00000108         -0.05736699        -0.080741210         -0.03419845
## cg00000165         -0.01319016        -0.047995678         -0.06241185
## cg00000236         -0.01138876        -0.004814477         -0.02297887
## cg00000289         -0.04888816         0.068479182          0.02160710
## cg00000292          0.02285773        -0.035014055         -0.02159311
##            6042316035_R01C01 6042316035_R02C02 6042316035_R03C01
## cg00000029      -0.048970456      -0.036373114      -0.086757920
## cg00000108       0.086375527       0.064852288       0.104133373
## cg00000165       0.007203867       0.006357423       0.006393753
## cg00000236       0.030166109       0.003971991       0.028539967
## cg00000289      -0.026610903      -0.011628566      -0.028058085
## cg00000292       0.001551847       0.054852972       0.046545952
##            6042316035_R03C02 6042316035_R04C01 6042316035_R05C01
## cg00000029       -0.05803014       0.019968850       -0.01505954
## cg00000108        0.10570267       0.013796492        0.02559819
## cg00000165        0.02452947      -0.038736399        0.06014500
## cg00000236        0.02274862       0.002821276        0.00220907
## cg00000289        0.05548152       0.013674307       -0.02835740
## cg00000292        0.03522538       0.011813358        0.01006673
##            6042316035_R05C02 6042316035_R06C02 6042316048_R01C01
## cg00000029      -0.006708644       -0.01629658      -0.005719858
## cg00000108      -0.003857607        0.02748118      -0.006868984
## cg00000165      -0.083765874        0.04095343       0.045866696
## cg00000236      -0.045049837        0.02406742      -0.037890218
## cg00000289       0.067509482       -0.02566225      -0.075890254
## cg00000292      -0.066636588        0.02322286       0.028695266
##            6042316048_R01C02 6042316048_R02C01 6042316048_R03C01
## cg00000029      -0.038307318       -0.02900953        0.05756239
## cg00000108       0.014424374        0.03770081       -0.02006521
## cg00000165       0.023442078       -0.03991822       -0.06827036
## cg00000236       0.043561868       -0.02734369        0.05191313
## cg00000289       0.005066534        0.06582683       -0.03797930
## cg00000292       0.008249420        0.05828668       -0.03355910
##            6042316048_R03C02 6042316048_R04C01 6042316048_R04C02
## cg00000029      -0.010304925       -0.06847880      -0.050578096
## cg00000108       0.035049544        0.04907161       0.012951364
## cg00000165      -0.041435091       -0.01130549      -0.111545998
## cg00000236       0.025641232       -0.01595664      -0.006998485
## cg00000289       0.003453812       -0.01242507      -0.047741522
## cg00000292       0.046571040        0.02635461       0.044716284
##            6042316048_R05C01 6042316048_R05C02 6042316110_R01C01
## cg00000029       -0.02222030       0.027209935       0.006248999
## cg00000108       -0.02398771       0.020748799       0.004882973
## cg00000165        0.08226064       0.006723145      -0.024733924
## cg00000236       -0.02075227      -0.005145125       0.002299315
## cg00000289        0.02085186       0.013330288      -0.031388125
## cg00000292        0.01432166      -0.094530166      -0.009250118
##            6042316110_R01C02 6042316110_R03C02 6042316110_R04C01
## cg00000029       0.054682611      -0.063068852       0.054170685
## cg00000108      -0.051668976       0.004359394       0.027315168
## cg00000165      -0.006826617       0.020328449       0.018159291
## cg00000236       0.034444756       0.035709200      -0.028438830
## cg00000289       0.013667862       0.049263513       0.004184905
## cg00000292       0.055934683      -0.070452425       0.003382441
##            6042316110_R04C02 6042316110_R05C01 6042316110_R06C01
## cg00000029      -0.007055549       0.025758172       0.079039989
## cg00000108       0.036583503      -0.039983301      -0.009434255
## cg00000165      -0.022340038      -0.023656034      -0.049859024
## cg00000236       0.005511324      -0.043540521       0.017505154
## cg00000289      -0.043408339       0.030059795       0.030966667
## cg00000292      -0.008179680      -0.005823909      -0.008148102
##            6042316121_R02C02 6042316121_R03C01 6042316121_R03C02
## cg00000029       -0.05088991      -0.014450200       -0.08214428
## cg00000108        0.02023539      -0.022113756        0.06179252
## cg00000165        0.04519375       0.006152607       -0.02285663
## cg00000236       -0.01669093       0.047619861       -0.01300075
## cg00000289        0.05326480       0.005852893       -0.02595497
## cg00000292        0.02875897      -0.019961883       -0.01905978
##            6042316121_R04C02 6042316121_R05C01 6042316121_R05C02
## cg00000029     -0.0317572615     -0.0009357763       0.080940001
## cg00000108      0.0392896848      0.0464490754       0.018815820
## cg00000165     -0.0006178874     -0.0203718322       0.005758404
## cg00000236      0.0222306585     -0.0257698072       0.031670416
## cg00000289     -0.0489449657     -0.0231911674       0.045460762
## cg00000292      0.0777549809      0.0230202411       0.081530826
##            6042316121_R06C01 6042316121_R06C02 6042316066_R01C01
## cg00000029      -0.009231385     -0.0761215395      -0.064913502
## cg00000108       0.035012864      0.0004704955       0.016029852
## cg00000165       0.061694539      0.0437912038      -0.052693683
## cg00000236       0.035791833     -0.0014415871       0.002985459
## cg00000289      -0.005468336      0.0571982952      -0.030770439
## cg00000292      -0.054346001     -0.0164205687      -0.004028796
##            6042316066_R02C01 6042316066_R02C02 6042316066_R03C01
## cg00000029       -0.04573521       0.058897114        0.04485883
## cg00000108        0.01152991       0.037044018       -0.02082535
## cg00000165       -0.08140122      -0.060995696        0.03979873
## cg00000236        0.03225710       0.009273412       -0.07043176
## cg00000289       -0.10408440       0.036107027       -0.04233301
## cg00000292        0.01038738      -0.010316536        0.03703241
##            6042316066_R04C01 6042316066_R04C02 6042316066_R05C01
## cg00000029      -0.057738211       0.022949170       0.062212790
## cg00000108      -0.052263612       0.039750123      -0.021052527
## cg00000165      -0.008350564      -0.064845164       0.010517116
## cg00000236       0.008360879       0.030079536       0.006480729
## cg00000289       0.066327414       0.055714486      -0.046596558
## cg00000292      -0.041113601       0.004458508      -0.025386211
##            6042316066_R06C01 6042316066_R06C02 6042316069_R01C01
## cg00000029      -0.001488203      -0.037017160       -0.09957635
## cg00000108       0.003925579       0.011489424        0.04169366
## cg00000165      -0.060835306       0.112281347       -0.04151961
## cg00000236      -0.076898855      -0.002631111        0.05100050
## cg00000289      -0.029132627       0.081300287        0.04909862
## cg00000292       0.076936400       0.026994503        0.01637416
##            6042316069_R01C02 6042316069_R02C01 6042316069_R03C01
## cg00000029      -0.078355310        0.01928522      0.0169584984
## cg00000108       0.007889516       -0.01231886      0.0090047830
## cg00000165      -0.020496933       -0.01149246     -0.0056358042
## cg00000236       0.028685311       -0.02656318      0.0001595185
## cg00000289      -0.013654003       -0.07096904     -0.0376055726
## cg00000292      -0.018828051       -0.06452832     -0.0005346722
##            6042316069_R03C02 6042316069_R04C02 6042316069_R05C01
## cg00000029       0.041930621       -0.02948333       0.003470199
## cg00000108      -0.032958562        0.01585219       0.022301566
## cg00000165       0.015107280       -0.05557115       0.021486230
## cg00000236      -0.009348200        0.04140572      -0.010273026
## cg00000289       0.002012162       -0.02092349      -0.022367655
## cg00000292       0.005081576        0.03247946       0.006541286
##            6042316069_R06C02 6042316094_R01C02 6042316094_R02C01
## cg00000029      3.605629e-02       -0.10115685        0.01755553
## cg00000108     -9.455634e-03       -0.02082222        0.06443412
## cg00000165      6.160533e-02       -0.08180593       -0.01771829
## cg00000236      2.547854e-02       -0.02666118       -0.03863971
## cg00000289     -5.114472e-05        0.02102734       -0.05880259
## cg00000292      7.908186e-02       -0.10608975       -0.02314262
##            6042316094_R03C02 6042316094_R04C01 6042316094_R04C02
## cg00000029      -0.002219126       -0.03383191       0.083990740
## cg00000108      -0.002873138        0.02408488      -0.053095266
## cg00000165      -0.036866112        0.02580358      -0.032326539
## cg00000236       0.002845042        0.02970981       0.008584349
## cg00000289       0.116512416       -0.02280947      -0.010749837
## cg00000292      -0.001487920        0.00929050      -0.037980943
##            6042316094_R05C01 6042316094_R05C02 6042316094_R06C01
## cg00000029      -0.038511289      -0.048095877        0.03157556
## cg00000108       0.048236892      -0.011232308       -0.05294410
## cg00000165       0.062692559       0.004813424       -0.04014223
## cg00000236       0.006302004       0.017870750       -0.01107892
## cg00000289       0.044210400       0.017632591        0.04156205
## cg00000292       0.014342408       0.034476619       -0.04315432
##            6042316099_R01C01 6042316099_R01C02 6042316099_R02C01
## cg00000029      -0.013183827      -0.063780164      -0.049529082
## cg00000108       0.016074480       0.040755435       0.036962360
## cg00000165       0.027959909      -0.034375083      -0.013054470
## cg00000236      -0.008419335       0.001489146       0.032716185
## cg00000289      -0.078218970      -0.024721657       0.008100698
## cg00000292       0.050346297      -0.030222018      -0.009803086
##            6042316099_R02C02 6042316099_R03C01 6042316099_R04C01
## cg00000029        0.01867524        0.05955723     -0.0021701537
## cg00000108        0.05622276       -0.02994785      0.0212107872
## cg00000165       -0.04016356        0.02857054     -0.0149821081
## cg00000236       -0.03101957        0.01205555      0.0004660098
## cg00000289        0.01181966        0.02951581      0.0388591999
## cg00000292        0.03338150        0.03148834      0.0604137254
##            6042316099_R04C02 6042316099_R05C01 6969568082_R02C01
## cg00000029      -0.070544718      -0.044788882      -0.004677607
## cg00000108       0.036327578       0.033781528       0.003258015
## cg00000165      -0.008756845      -0.012762685      -0.009632766
## cg00000236       0.045329548      -0.030961655      -0.042488496
## cg00000289       0.045650516      -0.002158985      -0.036296969
## cg00000292      -0.034406804       0.034547982       0.007355139
##            6969568082_R06C01 6969568082_R02C02 6969568082_R04C02
## cg00000029      -0.053327386       0.029950559      -0.017046193
## cg00000108       0.022803791      -0.013053329      -0.022232194
## cg00000165       0.021870714      -0.088889451       0.008735513
## cg00000236       0.004897495       0.026995113      -0.002580725
## cg00000289      -0.093823281       0.003090686      -0.015634705
## cg00000292      -0.022807293      -0.055528414       0.003929998
##            6969568082_R06C02 6969568084_R01C01 6969568084_R02C01
## cg00000029        0.06690351       0.014051919       -0.01573348
## cg00000108        0.01670513      -0.035051852        0.05077347
## cg00000165        0.02222211      -0.069232638       -0.03286932
## cg00000236       -0.02031604      -0.018627804        0.00425662
## cg00000289       -0.03179864      -0.006878855       -0.02941637
## cg00000292       -0.02820179      -0.079097067        0.03257560
##            6969568084_R03C01 6969568084_R04C01 6969568084_R06C01
## cg00000029       -0.03272153     -0.0191159913       -0.02538783
## cg00000108       -0.01829250     -0.0140895498        0.04973200
## cg00000165        0.03483726     -0.0211803191        0.03902259
## cg00000236       -0.01806050     -0.0005313369        0.02929129
## cg00000289       -0.02474035      0.0475475478       -0.05724456
## cg00000292        0.04059075     -0.0019386604       -0.01555619
##            6969568084_R02C02 6969568084_R03C02 6969568084_R04C02
## cg00000029      -0.045129557       0.016868899        0.04281913
## cg00000108       0.003994865       0.030919946        0.03190201
## cg00000165      -0.037484921       0.038220043       -0.02685203
## cg00000236       0.036629573       0.006970748       -0.02991920
## cg00000289       0.034942702      -0.035444805       -0.00288107
## cg00000292       0.034121746      -0.027954983       -0.02247140
##            6969568084_R05C02 6969568087_R01C01 6969568087_R02C01
## cg00000029      -0.013250679      -0.077946644      0.0986320158
## cg00000108       0.007726766      -0.017883630     -0.0071395684
## cg00000165      -0.050426272       0.078469675     -0.0417291962
## cg00000236      -0.026036591       0.011182239     -0.0076397970
## cg00000289      -0.004022972       0.002308149      0.0004805148
## cg00000292       0.012932202      -0.050620531      0.0531572165
##            6969568087_R03C01 6969568087_R04C01 6969568087_R05C01
## cg00000029       -0.07033743     -0.0566175835       -0.01394144
## cg00000108       -0.00348225      0.0008107322        0.02861689
## cg00000165       -0.02946870     -0.0379018057        0.04993027
## cg00000236       -0.02550410      0.0785697740       -0.00862133
## cg00000289       -0.08469936      0.1117014943       -0.04008420
## cg00000292        0.07262838     -0.0613951622        0.03286175
##            6969568087_R06C01 6969568087_R01C02 6969568087_R02C02
## cg00000029      0.0277609198      -0.096585645     -0.0223208540
## cg00000108      0.0006141955       0.023267204      0.0005760099
## cg00000165     -0.0167221287      -0.074969541     -0.0253989921
## cg00000236      0.0401562297      -0.058279341      0.0257302656
## cg00000289     -0.0609035890       0.001103639     -0.0005292531
## cg00000292     -0.0254577491      -0.020984464     -0.0623589081
##            6969568087_R03C02 6969568087_R05C02 6969568087_R06C02
## cg00000029       0.069628946      -0.033620322       -0.01216910
## cg00000108       0.001973082      -0.028021763        0.01953597
## cg00000165      -0.017800714      -0.008292071        0.01777206
## cg00000236       0.010239520       0.026693090       -0.01284403
## cg00000289       0.076942975       0.035076840       -0.02915578
## cg00000292       0.047967479       0.018000760       -0.03321252
##            6969568118_R01C01 6969568118_R02C01 6969568118_R03C01
## cg00000029      -0.037338334      -0.021526572       -0.03749075
## cg00000108      -0.026025177       0.041679037       -0.02282238
## cg00000165      -0.005052352      -0.076018516       -0.04493246
## cg00000236       0.013433279       0.018102876       -0.01097451
## cg00000289      -0.076304191       0.082225689       -0.01813444
## cg00000292       0.037190343      -0.005920771       -0.01608922
##            6969568118_R04C01 6969568118_R01C02 6969568118_R02C02
## cg00000029      -0.064269359       -0.04435358      -0.003722415
## cg00000108       0.053829357       -0.02252831       0.004816432
## cg00000165      -0.001225281        0.19437907       0.002568543
## cg00000236       0.011778616       -0.07096938       0.021954277
## cg00000289      -0.035465317        0.07776711      -0.014138917
## cg00000292      -0.071285316        0.07723081       0.043011237
##            6969568118_R03C02 6969568118_R04C02 6969568118_R06C02
## cg00000029       0.028787267        0.08417926       0.051556096
## cg00000108       0.021943740       -0.00963066      -0.013151690
## cg00000165      -0.007589860       -0.08668045      -0.025616240
## cg00000236      -0.002171499        0.02659214       0.008755928
## cg00000289       0.019347975       -0.01354645       0.015469226
## cg00000292       0.011847164       -0.01393619      -0.081640245
##            6929726046_R02C01 6929726046_R05C01 6929726046_R06C01
## cg00000029       -0.01644772       0.001640827     -4.926218e-03
## cg00000108        0.01038775      -0.012914749     -4.498309e-05
## cg00000165       -0.02235074       0.002216953      2.879400e-02
## cg00000236       -0.01628361      -0.009277111      3.122165e-02
## cg00000289        0.01582988       0.076060345     -5.439665e-02
## cg00000292       -0.01164375       0.012109070      4.544581e-02
##            6929726046_R01C02 6929726046_R03C02 6929726046_R04C02
## cg00000029        0.03558685       0.005942824       -0.10133106
## cg00000108       -0.03023766       0.004424879        0.03542042
## cg00000165        0.05302420      -0.052380443       -0.07579620
## cg00000236       -0.05059937      -0.002535402        0.02750112
## cg00000289       -0.03836783       0.040209018       -0.06872962
## cg00000292       -0.01729261      -0.081834274        0.03393768
##            6929726046_R06C02 6929718123_R01C01 6929718123_R04C01
## cg00000029       0.026141006        0.01729743       0.022304392
## cg00000108      -0.001523746        0.04425565       0.044186531
## cg00000165       0.030912478       -0.02998100      -0.024246356
## cg00000236      -0.019465211        0.04991722       0.023390596
## cg00000289      -0.019066375        0.04431602       0.026002911
## cg00000292       0.035076670       -0.01635131      -0.001925043
##            6929718123_R06C01 6929718123_R01C02 6929718123_R02C02
## cg00000029       0.004329684       0.013260086        0.09526864
## cg00000108       0.009269166       0.058319829       -0.21460463
## cg00000165      -0.024544970      -0.068344739        0.24054846
## cg00000236      -0.002904423       0.025472102       -0.06759236
## cg00000289      -0.007559549      -0.025803787       -0.03835756
## cg00000292      -0.034681951      -0.008278185       -0.01023279
##            6929718123_R03C02 6929718123_R04C02 6929718123_R05C02
## cg00000029       0.003633000       -0.08393219       -0.03047612
## cg00000108       0.013062767        0.07098336        0.07673064
## cg00000165      -0.051542699       -0.05876732       -0.07937070
## cg00000236       0.013824352        0.01437898        0.02995313
## cg00000289      -0.001109216        0.03745866       -0.05214824
## cg00000292      -0.004094120       -0.02093699        0.08182926
##            6929718123_R06C02 6929718136_R01C01 6929718136_R02C01
## cg00000029      -0.027055225       -0.01454914       -0.07327293
## cg00000108       0.030088243        0.01616009        0.02609431
## cg00000165       0.001677347       -0.04305192       -0.11904688
## cg00000236      -0.026314221        0.03218475        0.01825158
## cg00000289      -0.010721789        0.04840217        0.05337013
## cg00000292       0.039488902       -0.03080006        0.02877257
##            6929718136_R03C01 6929718136_R04C01 6929718136_R05C01
## cg00000029       0.003707032       0.141945052        0.01084711
## cg00000108       0.021595341      -0.243481690        0.07971104
## cg00000165       0.030636151       0.214694626       -0.08866109
## cg00000236       0.023061401      -0.073066169        0.01449140
## cg00000289       0.019950840       0.015864997        0.04284471
## cg00000292      -0.024537445      -0.003374063        0.03515945
##            6929718136_R02C02 6929718136_R04C02 6929718136_R05C02
## cg00000029       0.001743356      -0.005776088       0.022935519
## cg00000108       0.035266833       0.010142447       0.035618266
## cg00000165      -0.021412169      -0.039634291      -0.037411034
## cg00000236      -0.007188644       0.016690214      -0.015438295
## cg00000289      -0.043167095      -0.068233564       0.000662063
## cg00000292       0.006639525       0.048167890      -0.066090193
##            6929718136_R06C02 6929718138_R01C01 6929718138_R02C01
## cg00000029       0.032822927      -0.116145380       -0.04264571
## cg00000108       0.040939317       0.031874827        0.01602482
## cg00000165      -0.001523982      -0.006219180       -0.04029862
## cg00000236       0.003106411       0.005715255        0.01170597
## cg00000289       0.005903041      -0.033895144       -0.01889015
## cg00000292       0.058524196      -0.002019585       -0.06930508
##            6929718138_R04C01 6929718138_R06C01 6929718138_R01C02
## cg00000029       -0.01965923       0.071595482      -0.033859762
## cg00000108        0.02125885      -0.003748317      -0.068021887
## cg00000165        0.01528953       0.038785859      -0.043723690
## cg00000236        0.03507814      -0.065251321       0.004913748
## cg00000289        0.02645954       0.034455284       0.033328897
## cg00000292       -0.03934866      -0.011134764      -0.065586785
##            6929718138_R03C02 6929718138_R04C02 6929718138_R05C02
## cg00000029       0.089661215       -0.01804349      0.0617986432
## cg00000108      -0.137322229        0.01517420      0.0001309382
## cg00000165      -0.035403193       -0.01257305     -0.0342380996
## cg00000236       0.007913601       -0.01621569      0.0087242615
## cg00000289       0.103860548        0.03162578     -0.0265418229
## cg00000292       0.022697476       -0.03343145     -0.0066626634
##            6929718138_R06C02 6042316054_R02C01 6042316054_R02C02
## cg00000029       -0.02842293      -0.062976301      6.889188e-05
## cg00000108       -0.01597764       0.052999649      1.912159e-02
## cg00000165        0.04145893      -0.021760072      1.276173e-01
## cg00000236        0.02969823      -0.004521068      2.729477e-03
## cg00000289        0.01070902       0.030657659     -5.593795e-03
## cg00000292       -0.08239328       0.049257484      5.202113e-02
##            6042316054_R03C01 6042316054_R03C02 6042316054_R04C01
## cg00000029      -0.065841604       0.024734880        0.00628440
## cg00000108       0.056031409       0.059555339        0.02216390
## cg00000165      -0.027292932      -0.004727429       -0.04053245
## cg00000236       0.006207203       0.013780448       -0.02147427
## cg00000289      -0.068753396      -0.038918793       -0.05902741
## cg00000292      -0.007008479       0.005924709        0.02513701
##            6042316054_R04C02 6042316054_R05C01 6042316054_R05C02
## cg00000029       -0.09118299        0.06118746      -0.085249300
## cg00000108        0.06438194        0.02855368       0.005932598
## cg00000165       -0.06674148        0.03087454      -0.048195502
## cg00000236        0.04586734        0.02286315      -0.009068207
## cg00000289       -0.04602857        0.04509505      -0.038266721
## cg00000292        0.04487446        0.04614512      -0.037330566
##            6042316054_R06C01 6042316063_R02C01 6042316063_R02C02
## cg00000029        0.03175151      -0.028164993      -0.015394521
## cg00000108        0.03087354       0.004311434       0.035505315
## cg00000165       -0.01226921       0.037353898      -0.050789165
## cg00000236       -0.03383063      -0.020394693       0.008810408
## cg00000289       -0.05513587      -0.044668465      -0.010835464
## cg00000292        0.08955615       0.037300411      -0.019249265
##            6042316063_R03C01 6042316063_R03C02 6042316063_R04C01
## cg00000029     -0.0121808444       0.079905006        0.02024576
## cg00000108      0.0053798705       0.003827802       -0.04355415
## cg00000165     -0.0001022725      -0.069489161       -0.04543132
## cg00000236      0.0371158242       0.016256774       -0.03989263
## cg00000289      0.0328465794       0.045668429        0.03033313
## cg00000292      0.0038707393      -0.011662822        0.04354174
##            6042316063_R04C02 6042316063_R05C01 6042316063_R05C02
## cg00000029       0.007403209       -0.06747424      0.0064219196
## cg00000108       0.031769835       -0.03631444      0.0145688115
## cg00000165      -0.040306397        0.07838846     -0.0048089584
## cg00000236       0.027294264       -0.04795918      0.0526791898
## cg00000289       0.014438893       -0.09985097      0.0328427483
## cg00000292      -0.033818322       -0.07157251      0.0007256785
##            6042316063_R06C02 6042316065_R01C02 6042316065_R02C02
## cg00000029      -0.053337390       0.030886154       -0.02809092
## cg00000108      -0.028596618       0.014759244        0.02685445
## cg00000165       0.051557243       0.100403230       -0.01672746
## cg00000236      -0.001341383       0.028242329        0.03300033
## cg00000289       0.050724418       0.087212036        0.10845641
## cg00000292      -0.030222392       0.005683394       -0.04326878
##            6042316065_R03C01 6042316065_R04C01 6042316065_R04C02
## cg00000029        0.04377972      -0.078222070        0.01509872
## cg00000108        0.03217197       0.008765247        0.04220093
## cg00000165       -0.03499162      -0.040925367       -0.03816981
## cg00000236        0.04277546      -0.012844830       -0.06645660
## cg00000289        0.03664736      -0.027544030       -0.04943090
## cg00000292       -0.07755945       0.030331256        0.02943761
##            6042316065_R05C02 6042316065_R06C02 6042316103_R02C01
## cg00000029       -0.00222224        0.02978311       -0.01575728
## cg00000108        0.02066493       -0.08685874       -0.00964249
## cg00000165        0.05832922        0.02825846       -0.05995749
## cg00000236        0.01057512        0.01424094       -0.02117133
## cg00000289       -0.03428914       -0.06768299       -0.02103041
## cg00000292       -0.04638713        0.01475478        0.01431185
##            6042316103_R03C01 6042316103_R03C02 6042316103_R04C01
## cg00000029       -0.07213033      -0.043660279       0.034787158
## cg00000108        0.03248377       0.022563552       0.003726285
## cg00000165       -0.07194299       0.027982713      -0.005316469
## cg00000236       -0.01887901       0.049869616       0.028978239
## cg00000289       -0.04847780      -0.081979014       0.087534578
## cg00000292        0.10219542       0.005916879      -0.003655445
##            6042316103_R05C01 6042316103_R06C01 6042316103_R06C02
## cg00000029      -0.012638585      -0.105784541     -0.0041377865
## cg00000108      -0.008499276      -0.005022802      0.0526616930
## cg00000165       0.041860851       0.056628097     -0.1033157229
## cg00000236       0.015360202      -0.007297967     -0.0443600377
## cg00000289      -0.013090670       0.015427820     -0.0009227993
## cg00000292       0.015681991      -0.040175889      0.0423264707
##            6042316036_R01C02 6042316036_R02C01 6042316036_R03C01
## cg00000029        0.02121810       -0.01476756       0.035519362
## cg00000108       -0.03841845        0.04815954       0.006515987
## cg00000165       -0.06715237       -0.01141530      -0.015560394
## cg00000236        0.02270812       -0.02798921      -0.029511570
## cg00000289        0.01330893       -0.05349697       0.021082749
## cg00000292       -0.05898121        0.07358675      -0.024045097
##            6042316036_R03C02 6042316036_R04C01 6042316036_R04C02
## cg00000029      -0.020169658      -0.083400099      -0.033865300
## cg00000108       0.015569785       0.004324409       0.001288728
## cg00000165      -0.057357218      -0.075795397      -0.002968016
## cg00000236      -0.014187418       0.016889162       0.032851993
## cg00000289      -0.026192596       0.059133392      -0.034160767
## cg00000292       0.004438251      -0.048582687      -0.035370026
##            6042316036_R05C01 6042316036_R05C02 6042316036_R06C01
## cg00000029      -0.003059017       0.010803774       0.061840465
## cg00000108       0.067019626      -0.038265914       0.008681369
## cg00000165       0.067334149      -0.028574667       0.019808611
## cg00000236       0.040185559      -0.068297090      -0.021483548
## cg00000289       0.089692648       0.001043375      -0.004819482
## cg00000292      -0.077946905       0.113171840       0.005229846
##            6042316036_R06C02 6042316050_R01C02 6042316050_R02C02
## cg00000029       -0.02747549        0.03654060        0.16612500
## cg00000108        0.02929343        0.02552003       -0.03636615
## cg00000165        0.07411231       -0.06566733        0.01844118
## cg00000236        0.02463774        0.01512532        0.03862478
## cg00000289       -0.12158496       -0.04287669        0.03074710
## cg00000292       -0.04285982       -0.04465036        0.05365070
##            6042316050_R03C01 6042316050_R04C01 6042316050_R04C02
## cg00000029      -0.021306393      -0.019015995      -0.052624201
## cg00000108      -0.033184170      -0.003850467       0.011062054
## cg00000165      -0.068934807      -0.011733118       0.004249727
## cg00000236      -0.021645851       0.024551147      -0.002294912
## cg00000289      -0.005023793       0.069023533       0.075876834
## cg00000292       0.042623973      -0.039107386       0.027929830
##            6042316050_R05C01 6042316050_R05C02 6042316050_R06C01
## cg00000029       -0.02422391       -0.04929616       0.025616374
## cg00000108       -0.03986664        0.06146223      -0.006619355
## cg00000165        0.04871063        0.03859601       0.144316205
## cg00000236        0.01556205       -0.06560167      -0.037320484
## cg00000289       -0.01357754       -0.06486901      -0.030720055
## cg00000292        0.02135960        0.03474063       0.062938961
##            6042316050_R06C02 6042316053_R01C01 6042316053_R02C01
## cg00000029       0.004757818     -0.0450199057      -0.018813815
## cg00000108       0.017041487      0.0055172810       0.013950924
## cg00000165      -0.054384878     -0.0196291034       0.028875397
## cg00000236      -0.005915417     -0.0005276243       0.012278663
## cg00000289       0.012050355      0.0028729442      -0.090964643
## cg00000292      -0.008061355      0.0072924205      -0.007589214
##            6042316053_R02C02 6042316053_R03C01 6042316053_R03C02
## cg00000029        0.02210445       0.009729315      -0.016090746
## cg00000108        0.07373384       0.065612128       0.041686386
## cg00000165       -0.01816479       0.017620567       0.008388825
## cg00000236        0.03676097       0.013899334       0.023147760
## cg00000289        0.04593070      -0.015893855      -0.074446551
## cg00000292       -0.03972255      -0.009580992      -0.003038724
##            6042316053_R04C01 6042316053_R04C02 6042316053_R06C01
## cg00000029      6.819251e-03       -0.10812866      -0.016391577
## cg00000108      5.706497e-02        0.03636217       0.020237789
## cg00000165     -1.718582e-02        0.01601753       0.041849896
## cg00000236     -9.693715e-05        0.01294207       0.008164102
## cg00000289      9.106071e-03       -0.01975206       0.055242837
## cg00000292      2.776256e-02        0.08427268       0.065592869
##            6042316053_R06C02 6042316061_R01C01 6042316061_R02C01
## cg00000029      -0.002869923        0.02436390       0.002516371
## cg00000108      -0.027039262        0.02942921       0.035991508
## cg00000165      -0.044757296       -0.02721875      -0.012325216
## cg00000236      -0.023774086        0.02210199       0.042359031
## cg00000289      -0.020605387       -0.01204909      -0.031497707
## cg00000292      -0.003440287        0.01571647      -0.024614059
##            6042316061_R03C01 6042316061_R03C02 6042316061_R04C01
## cg00000029      -0.091758542       0.030738767      -0.021899694
## cg00000108       0.052486994       0.004411943      -0.004469086
## cg00000165      -0.044379894      -0.062946057      -0.033019630
## cg00000236       0.006289828       0.005424975       0.003318867
## cg00000289       0.020143192       0.035784943       0.054119449
## cg00000292      -0.072879315      -0.042403839      -0.084951771
##            6042316061_R04C02 6042316061_R05C01 6042316061_R05C02
## cg00000029      -0.005267181       0.056875260      -0.012557471
## cg00000108       0.035928202       0.009417450       0.050777539
## cg00000165      -0.077967883      -0.007466972      -0.002709187
## cg00000236       0.029490332      -0.005277500       0.025409193
## cg00000289       0.041962025       0.016565500       0.112016109
## cg00000292       0.014601011       0.057146411      -0.027244186
##            7796806022_R01C01 7796806022_R03C01 7796806022_R04C01
## cg00000029        0.18813456       -0.13387018       -0.08542779
## cg00000108       -0.38135616        0.03209525        0.08086202
## cg00000165        0.11249096       -0.03774922       -0.06534733
## cg00000236       -0.14346669        0.01072456        0.01841916
## cg00000289        0.02459381       -0.10728074       -0.06343301
## cg00000292       -0.03782085        0.03116156        0.06468448
##            7796806022_R04C02 7796806022_R05C02 7796806022_R06C01
## cg00000029        0.07720881      -0.096821870        0.06798454
## cg00000108       -0.10368844       0.086628374       -0.12252699
## cg00000165        0.13538551      -0.051369447        0.15144750
## cg00000236       -0.03101815       0.029927426       -0.01623913
## cg00000289        0.05926402      -0.007767647        0.02982745
## cg00000292       -0.08648967       0.013364678       -0.00499560
##            7786923046_R03C02 7786923046_R04C01 7786923046_R06C01
## cg00000029      -0.059579084       -0.10100534       -0.09029319
## cg00000108       0.085812363        0.05499180        0.05815749
## cg00000165      -0.081465928       -0.09370453       -0.04333321
## cg00000236       0.047493097       -0.01892043        0.02131681
## cg00000289       0.003018994       -0.05435531       -0.06834530
## cg00000292       0.040397106       -0.01028608       -0.04686433
##            7796806038_R03C01 7796806038_R03C02 7796806038_R05C02
## cg00000029       -0.07028833       -0.06011033       0.063965331
## cg00000108        0.07972570        0.07418552      -0.154073724
## cg00000165       -0.07263551       -0.07369056       0.218038211
## cg00000236        0.03162461        0.03195347      -0.033076627
## cg00000289       -0.08838488       -0.03409267       0.058823350
## cg00000292       -0.02417461        0.02298257      -0.009015551
##            7786923063_R01C01 7786923063_R01C02 7786923063_R03C01
## cg00000029       0.084120315       -0.07966619       0.107565693
## cg00000108      -0.132527424        0.08217896      -0.073618866
## cg00000165      -0.068157624       -0.07138852      -0.047655083
## cg00000236      -0.066878293        0.01821430      -0.038657236
## cg00000289      -0.002348422       -0.06473794       0.004273524
## cg00000292      -0.010316622        0.04358162      -0.068804910
##            7786923063_R04C01 7786923063_R06C01 7786923107_R01C02
## cg00000029        0.07975142       0.065235198       -0.02321370
## cg00000108       -0.06346170      -0.104621547        0.05607435
## cg00000165        0.09764119       0.076686649       -0.04702280
## cg00000236       -0.04706763      -0.042528528        0.03045289
## cg00000289        0.06546496      -0.019734406        0.08636131
## cg00000292       -0.07825407      -0.001615719        0.03867048
##            7786923107_R04C01 7786923107_R05C01 7786923107_R06C01
## cg00000029       -0.01550184     -0.0758852118       0.023190325
## cg00000108        0.06298056     -0.0022458174      -0.005766173
## cg00000165       -0.01514033     -0.0996448139      -0.007799828
## cg00000236        0.03454756      0.0007694774      -0.035524124
## cg00000289        0.01066947     -0.0052112529       0.044127226
## cg00000292        0.02381715     -0.0265738179       0.033480330
##            7796806016_R03C02 7796806016_R06C01 7796806002_R01C01
## cg00000029     -0.0117076383        0.04760094      -0.005657157
## cg00000108      0.0538313364        0.05270432      -0.042140634
## cg00000165     -0.0001215885       -0.02321470       0.085601875
## cg00000236     -0.0237272763        0.01187934       0.036496421
## cg00000289     -0.0554852032       -0.03283673      -0.033583579
## cg00000292     -0.0256462883        0.05908430       0.048082973
##            7796806002_R02C02 7796806002_R04C02 7796806002_R05C01
## cg00000029        0.06027932      -0.063366290       0.016721604
## cg00000108        0.01421109      -0.027569686      -0.026820021
## cg00000165        0.03319685      -0.051394728       0.008913765
## cg00000236        0.02837070       0.031059589       0.032584716
## cg00000289        0.02713931       0.013250506       0.031870837
## cg00000292        0.04924403      -0.009224635      -0.028560933
##            7796806002_R05C02 7796806002_R06C01 7796806002_R06C02
## cg00000029      -0.044320394      -0.076929697        0.05387500
## cg00000108       0.048438387       0.003275044       -0.03522353
## cg00000165      -0.046580412       0.037858194        0.02437887
## cg00000236      -0.012160706      -0.038022830       -0.02351239
## cg00000289       0.003231339      -0.032021573        0.01134227
## cg00000292       0.019046942      -0.039083623       -0.02730459
##            7796806029_R02C01 7796806029_R04C01 7796806029_R06C01
## cg00000029       0.011903133       -0.07215405       -0.06463801
## cg00000108       0.015235789        0.08683691        0.09526527
## cg00000165      -0.012766163       -0.09014602       -0.06715891
## cg00000236       0.005446803        0.02430563        0.02641312
## cg00000289       0.040148554        0.06224708       -0.01821791
## cg00000292      -0.058704905        0.02981101        0.02784194
##            6042316024_R01C01 6042316024_R01C02 6042316024_R02C01
## cg00000029        0.13160143      0.0309003908        0.05931944
## cg00000108       -0.05418732     -0.1378300117       -0.08682865
## cg00000165        0.02720036      0.1624529953        0.11771536
## cg00000236       -0.09140317      0.0008947312       -0.02249558
## cg00000289       -0.02663287      0.0364997804       -0.02375618
## cg00000292        0.06026970     -0.0563608721       -0.08780811
##            6042316024_R02C02 6042316024_R03C01 6042316024_R04C01
## cg00000029        0.05630710       -0.03307004       0.045953071
## cg00000108       -0.07556941        0.04214521      -0.041224227
## cg00000165        0.03091579        0.04203310      -0.195269868
## cg00000236       -0.01654299        0.01722181      -0.010365114
## cg00000289        0.03231427        0.02673120       0.003718176
## cg00000292       -0.02969992       -0.02721597      -0.050398658
##            6042316024_R04C02 6042316024_R05C01 6042316024_R05C02
## cg00000029       0.001134925      -0.083284867        0.04459715
## cg00000108       0.155770793      -0.002895653        0.10475599
## cg00000165       0.113019826       0.055131081       -0.02282359
## cg00000236       0.042350482      -0.054546509       -0.01889059
## cg00000289      -0.092821358      -0.033311623        0.01363040
## cg00000292       0.067111937       0.022895781        0.04824378
##            6042316024_R06C01 6042316024_R06C02 6042316031_R01C01
## cg00000029        0.03476491      -0.081917302       0.049862313
## cg00000108       -0.10438455      -0.076656682      -0.052429597
## cg00000165       -0.02357338      -0.005532757       0.081735406
## cg00000236        0.05601884      -0.037063729      -0.013678834
## cg00000289       -0.02028148       0.036296188       0.033690631
## cg00000292       -0.09059774       0.044646235      -0.001675932
##            6042316031_R01C02 6042316031_R02C01 6042316031_R02C02
## cg00000029        0.10613317      -0.026588608        0.13137355
## cg00000108       -0.17693377      -0.107156073        0.02829998
## cg00000165        0.11017344      -0.151969850       -0.03332942
## cg00000236        0.00109833       0.002985782       -0.04486064
## cg00000289        0.01489245      -0.004051470       -0.11878381
## cg00000292       -0.06248817       0.039582019       -0.11837306
##            6042316031_R03C01 6042316031_R03C02 6042316031_R04C01
## cg00000029       -0.09457093       0.028596707      -0.001650137
## cg00000108       -0.07483515       0.061760306       0.008202015
## cg00000165        0.14363027       0.105320943      -0.026303090
## cg00000236       -0.04058001       0.002694433      -0.051568889
## cg00000289       -0.03648718       0.059679964       0.014948006
## cg00000292        0.02678639      -0.028172520      -0.013549945
##            6042316031_R04C02 6042316031_R05C01 6042316031_R05C02
## cg00000029        0.02726591      -0.052111347        0.05959621
## cg00000108        0.04101154       0.130046805       -0.04815920
## cg00000165        0.08574028      -0.005669618        0.05131192
## cg00000236        0.01222532       0.059734677        0.04114773
## cg00000289        0.05063695       0.005421664        0.02059376
## cg00000292        0.01109067       0.011157734        0.02310104
##            6042316031_R06C01 6042316031_R06C02 6042316042_R01C01
## cg00000029        0.03905173       0.003761404       -0.01099087
## cg00000108       -0.12879800      -0.061924796       -0.04472029
## cg00000165       -0.05448958       0.144832094        0.14825684
## cg00000236       -0.02157438      -0.015563551       -0.04339151
## cg00000289        0.01296883      -0.078927869       -0.01076310
## cg00000292        0.03840321       0.050079615       -0.01132098
##            6042316042_R01C02 6042316042_R02C01 6042316042_R02C02
## cg00000029       -0.01639725       -0.02427202        0.06323928
## cg00000108        0.12608302        0.02833918       -0.10738712
## cg00000165       -0.09994209       -0.13085491        0.07291622
## cg00000236        0.04059687       -0.01343744        0.01647660
## cg00000289        0.04625213        0.02600557        0.05562698
## cg00000292       -0.01551951        0.01559204       -0.03696173
##            6042316042_R03C01 6042316042_R03C02 6042316042_R04C01
## cg00000029        0.03720825       -0.02091998       0.044028749
## cg00000108       -0.05503579        0.02197657       0.093907993
## cg00000165        0.02252452        0.02986686       0.002034427
## cg00000236        0.02630639       -0.02077790       0.039900155
## cg00000289       -0.05639211       -0.03406055       0.074806930
## cg00000292       -0.06999252        0.09849547      -0.026071462
##            6042316042_R04C02 6042316042_R05C01 6042316042_R05C02
## cg00000029       0.010939050       0.049058442      -0.056693697
## cg00000108      -0.092704577      -0.002244566      -0.111353855
## cg00000165      -0.047320786       0.176855798      -0.024269955
## cg00000236      -0.066097157       0.036926154      -0.041000484
## cg00000289       0.002923121       0.012593449      -0.006900669
## cg00000292      -0.013684580      -0.039110996      -0.030742850
##            6042316042_R06C01 6042316042_R06C02 6042316047_R01C01
## cg00000029       0.107907639       0.053039812       0.004223237
## cg00000108       0.006406907       0.055670274      -0.102056294
## cg00000165      -0.008285049      -0.157162403       0.027993417
## cg00000236       0.013236212      -0.001058689      -0.043318074
## cg00000289      -0.030012436      -0.021033399       0.049745855
## cg00000292       0.080821690       0.031784355       0.039859835
##            6042316047_R02C01 6042316047_R02C02 6042316047_R04C01
## cg00000029       -0.07848284      -0.014239080        0.03245215
## cg00000108        0.13925338      -0.117063707       -0.03146553
## cg00000165        0.03890555       0.025656519       -0.08323277
## cg00000236        0.01144864       0.006342663        0.02339155
## cg00000289       -0.02069707       0.033860074       -0.03458747
## cg00000292       -0.04837710      -0.027048478        0.06592116
##            6042316047_R04C02 6042316047_R05C01 6042316047_R05C02
## cg00000029       -0.04882195      -0.040543789        0.12211325
## cg00000108       -0.14148748       0.069881556        0.08187901
## cg00000165       -0.05320134       0.145331485        0.05618319
## cg00000236       -0.07570961      -0.009460791        0.03392479
## cg00000289        0.02695907      -0.029275857        0.04083855
## cg00000292       -0.02664182      -0.056525075        0.01151913
##            6042316047_R06C01 6055432012_R01C01 6055432012_R01C02
## cg00000029      0.0105043974       -0.05706777      -0.081807851
## cg00000108     -0.0465891624        0.03893135       0.010081923
## cg00000165      0.0006011391        0.01539467      -0.049496217
## cg00000236      0.0018219805       -0.01112808       0.066356197
## cg00000289      0.0101037084       -0.07847716      -0.052910055
## cg00000292     -0.0037612738       -0.03816983       0.004646135
##            6055432012_R02C02 6055432012_R03C01 6055432012_R03C02
## cg00000029      -0.074674379       0.067332871       0.014784883
## cg00000108       0.024433483      -0.003637894       0.041088201
## cg00000165       0.080288766      -0.059522376      -0.040606547
## cg00000236       0.007494356       0.039190771       0.003693136
## cg00000289       0.020179396       0.001497015       0.057282665
## cg00000292      -0.070652692       0.032568071       0.078483562
##            6055432012_R04C01 6055432012_R04C02 6055432012_R05C01
## cg00000029        0.01649184       0.042057685       0.044906172
## cg00000108       -0.06816652       0.030902790       0.006379901
## cg00000165        0.03327797      -0.014342591       0.032995495
## cg00000236       -0.03691982      -0.060473866      -0.008204133
## cg00000289       -0.01995989      -0.007091917       0.053027908
## cg00000292        0.03196323      -0.053647464       0.042291742
##            6055432012_R05C02 6055432012_R06C01 6055432012_R06C02
## cg00000029      -0.032805522       -0.01752553       0.020608457
## cg00000108       0.009386528        0.01945818       0.007117726
## cg00000165      -0.042942899        0.06201615       0.006711199
## cg00000236       0.005114657       -0.01916321      -0.020460544
## cg00000289       0.009806631       -0.05772295       0.044065739
## cg00000292       0.011437674        0.07419800       0.039002648
##            6055432029_R01C01 6055432029_R01C02 6055432029_R02C01
## cg00000029     -0.0649700920        0.05447062       0.050107315
## cg00000108     -0.0556355211        0.03919279       0.007696739
## cg00000165      0.0217803577        0.06535873      -0.075475073
## cg00000236     -0.0808419939        0.03431710       0.039236535
## cg00000289     -0.0005883725       -0.03812006      -0.001412742
## cg00000292      0.0366408100        0.01839698       0.030878630
##            6055432029_R02C02 6055432029_R03C01 6055432029_R03C02
## cg00000029        0.04719963        0.02183834       -0.03850940
## cg00000108       -0.01762682        0.02371325        0.04454174
## cg00000165        0.09772423       -0.07624939       -0.05790299
## cg00000236        0.01633175       -0.01090808        0.03428975
## cg00000289        0.07897550       -0.05492298        0.05088788
## cg00000292       -0.02089349       -0.05290145        0.05622128
##            6055432029_R04C01 6055432029_R04C02 6055432029_R05C01
## cg00000029      -0.026117764      -0.079696010        0.02546388
## cg00000108      -0.003369482       0.009788120        0.01702579
## cg00000165       0.015577644      -0.020448968       -0.04821398
## cg00000236      -0.012643014       0.002003399        0.03933155
## cg00000289      -0.065491154       0.026907127       -0.01088874
## cg00000292      -0.016611848      -0.065890712       -0.04762925
##            6055432029_R05C02 6055432029_R06C01 6055432029_R06C02
## cg00000029      -0.022250644       0.092836577       -0.08720197
## cg00000108       0.006885344      -0.009307008       -0.01645537
## cg00000165      -0.016827858       0.066510940        0.09102397
## cg00000236       0.033973061      -0.050703839       -0.01703021
## cg00000289      -0.010167752       0.017310744        0.03226558
## cg00000292       0.038722101       0.001061321        0.02348651
##            6055432060_R01C01 6055432060_R01C02 6055432060_R02C01
## cg00000029       -0.05263865       0.010121765       -0.07497713
## cg00000108        0.04093390      -0.001616758        0.01040243
## cg00000165       -0.05920586       0.025097092       -0.03715382
## cg00000236       -0.01778705       0.017225269        0.03591030
## cg00000289        0.02140306      -0.095041695        0.03458943
## cg00000292       -0.07436693      -0.043176160        0.01816218
##            6055432060_R02C02 6055432060_R03C01 6055432060_R03C02
## cg00000029      -0.008341070       0.051872378        0.07719475
## cg00000108      -0.021539670       0.024065897        0.05020271
## cg00000165       0.032093444      -0.052194372       -0.07013680
## cg00000236      -0.017870996      -0.004802947        0.03969545
## cg00000289      -0.079519213       0.036489063        0.02633210
## cg00000292       0.009448922       0.047695772        0.00553603
##            6055432060_R04C01 6055432060_R04C02 6055432060_R05C01
## cg00000029       -0.04037923        0.02924844       -0.01299644
## cg00000108        0.04815750        0.02830378        0.01688123
## cg00000165       -0.01219338       -0.04317743        0.03031256
## cg00000236        0.01510604        0.03845406       -0.03950375
## cg00000289        0.06740678       -0.02318877       -0.06359627
## cg00000292        0.03393979       -0.04448200        0.03068524
##            6055432060_R05C02 6055432060_R06C01 6055432060_R06C02
## cg00000029       -0.05780055       0.008648630       0.112959714
## cg00000108        0.04903105      -0.006478811       0.002019833
## cg00000165       -0.03865444       0.013405000       0.122231304
## cg00000236        0.03150675      -0.012384599      -0.014444333
## cg00000289        0.02186394      -0.021765345      -0.017942392
## cg00000292        0.04363831      -0.038695656       0.060628042
##            6055432066_R01C01 6055432066_R01C02 6055432066_R02C01
## cg00000029       -0.08197272       -0.01180591        0.03395896
## cg00000108        0.01312615        0.02562109        0.02895073
## cg00000165       -0.07400477       -0.04263296       -0.03445878
## cg00000236        0.02732357       -0.02713413       -0.02724923
## cg00000289        0.03156017       -0.04134360        0.06723065
## cg00000292       -0.02911341        0.02570434       -0.02952502
##            6055432066_R02C02 6055432066_R03C02 6055432066_R04C01
## cg00000029     -0.1005827792      -0.004897766       0.044388575
## cg00000108      0.0015713480       0.034774152       0.025461694
## cg00000165      0.0367625880      -0.036609823       0.007409343
## cg00000236      0.0262563800       0.004397718       0.017844174
## cg00000289     -0.0286168619      -0.062618899      -0.007760500
## cg00000292      0.0008279338       0.038128083       0.005849442
##            6055432066_R04C02 6055432066_R05C01 6055432066_R05C02
## cg00000029       -0.06996924       -0.05633214       -0.05141607
## cg00000108        0.10677810        0.11249514        0.09255569
## cg00000165       -0.01775536       -0.06249792       -0.08311094
## cg00000236        0.02500097        0.02244031        0.01604663
## cg00000289       -0.01313459       -0.05456316       -0.06860809
## cg00000292        0.02682502        0.06462833        0.03626712
##            6057825115_R01C01 6057825115_R03C01 6057825115_R05C01
## cg00000029       0.076089413        0.06210317        0.08963579
## cg00000108       0.002625882       -0.09512454       -0.10750391
## cg00000165       0.017665537        0.07323341        0.05697704
## cg00000236      -0.012359366       -0.07349746        0.02276646
## cg00000289      -0.062419910       -0.01186669        0.09104878
## cg00000292       0.017644968       -0.04327183       -0.04111961
##            6057825115_R01C02 6057825115_R03C02 6057825115_R05C02
## cg00000029        0.06496866       0.029213410      -0.004121194
## cg00000108       -0.14024724      -0.022019377      -0.019793064
## cg00000165        0.02191585      -0.086290639       0.118944154
## cg00000236       -0.03308526      -0.021246814       0.006532010
## cg00000289        0.04129460       0.065549924      -0.052523482
## cg00000292       -0.03288444       0.001458819      -0.010158508
##            6057825128_R01C01 6057825128_R03C01 6057825128_R05C01
## cg00000029       0.029363999        0.09458836      -0.030295568
## cg00000108      -0.058618219        0.01845351      -0.210129845
## cg00000165       0.082572034       -0.06804830      -0.050063427
## cg00000236      -0.073293825        0.04551846       0.001011198
## cg00000289      -0.062143826        0.04863256      -0.024841190
## cg00000292       0.004037172        0.08366896      -0.076226995
##            6057825115_R02C01 6057825115_R04C01 6057825115_R06C01
## cg00000029      -0.018093019       0.026027771      -0.042874777
## cg00000108       0.027622115       0.032471376       0.065702544
## cg00000165       0.188707569      -0.092212301       0.042099826
## cg00000236      -0.008417854       0.004470456      -0.006955403
## cg00000289      -0.038712473       0.056616600       0.045586314
## cg00000292      -0.007087620      -0.013230509      -0.006871899
##            6057825115_R02C02 6057825115_R04C02 6057825115_R06C02
## cg00000029      -0.031183979        0.11055844      -0.022352572
## cg00000108      -0.043180670       -0.12967317       0.060898660
## cg00000165       0.111230907        0.09249014      -0.108588300
## cg00000236       0.001738708       -0.01348703      -0.001402109
## cg00000289      -0.050626764        0.03223563       0.035948232
## cg00000292       0.005367511       -0.02664265       0.008032464
##            6057825128_R01C02 6057825128_R03C02 6057825128_R05C02
## cg00000029       0.027976500        0.02432124      -0.010126298
## cg00000108      -0.013023808       -0.03159508      -0.021826065
## cg00000165      -0.073492713       -0.03484038       0.180283104
## cg00000236      -0.047810724        0.02786902      -0.008821356
## cg00000289      -0.044334673        0.07721447      -0.042547102
## cg00000292       0.006020978       -0.00797084      -0.001349628
```

```r
colnames(residuals)<-colnames(brain.cor.dat)

# generate adjusted residuals by adding the mean beta of each probe to the residuals
adj.residuals<-residuals+matrix(apply(brain.cor.dat, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

r1<-as.data.frame(adj.residuals)
# check difference between corrected and uncorrected methylation data
all.equal(r1,brain.cor.dat)
```

```
##   [1] "Component \"6057825008_R02C02\": Mean relative difference: 0.00496866"    
##   [2] "Component \"6057825008_R03C01\": Mean relative difference: 0.005435411"   
##   [3] "Component \"6057825008_R04C01\": Mean relative difference: 4.633884e-05"  
##   [4] "Component \"6057825008_R04C02\": Mean relative difference: 0.0006000302"  
##   [5] "Component \"6057825008_R05C01\": Mean relative difference: 0.001557573"   
##   [6] "Component \"6057825008_R05C02\": Mean relative difference: 0.00732968"    
##   [7] "Component \"6057825008_R06C01\": Mean relative difference: 0.005541414"   
##   [8] "Component \"6057825008_R06C02\": Mean relative difference: 0.003643052"   
##   [9] "Component \"6057825014_R01C02\": Mean relative difference: 0.003580407"   
##  [10] "Component \"6057825014_R02C02\": Mean relative difference: 0.01135846"    
##  [11] "Component \"6057825014_R03C01\": Mean relative difference: 0.002925688"   
##  [12] "Component \"6057825014_R03C02\": Mean relative difference: 0.007995749"   
##  [13] "Component \"6057825014_R04C01\": Mean relative difference: 0.00274609"    
##  [14] "Component \"6057825014_R04C02\": Mean relative difference: 0.00209708"    
##  [15] "Component \"6057825014_R06C01\": Mean relative difference: 0.00296576"    
##  [16] "Component \"6057825017_R01C01\": Mean relative difference: 0.002302703"   
##  [17] "Component \"6057825017_R01C02\": Mean relative difference: 0.004800294"   
##  [18] "Component \"6057825017_R02C02\": Mean relative difference: 0.000226121"   
##  [19] "Component \"6057825017_R03C01\": Mean relative difference: 0.0002355218"  
##  [20] "Component \"6057825017_R03C02\": Mean relative difference: 0.003239404"   
##  [21] "Component \"6057825017_R04C01\": Mean relative difference: 0.003695475"   
##  [22] "Component \"6057825017_R04C02\": Mean relative difference: 0.005639822"   
##  [23] "Component \"6057825017_R05C01\": Mean relative difference: 0.003935768"   
##  [24] "Component \"6057825017_R05C02\": Mean relative difference: 0.009781422"   
##  [25] "Component \"6057825017_R06C01\": Mean relative difference: 0.007565584"   
##  [26] "Component \"6057825018_R01C01\": Mean relative difference: 0.00436801"    
##  [27] "Component \"6057825018_R02C02\": Mean relative difference: 0.006677067"   
##  [28] "Component \"6057825018_R03C02\": Mean relative difference: 0.005889199"   
##  [29] "Component \"6057825018_R04C01\": Mean relative difference: 0.0008703853"  
##  [30] "Component \"6057825018_R04C02\": Mean relative difference: 0.0008157072"  
##  [31] "Component \"6057825018_R05C01\": Mean relative difference: 0.005321763"   
##  [32] "Component \"6057825018_R05C02\": Mean relative difference: 0.003761469"   
##  [33] "Component \"6057825018_R06C01\": Mean relative difference: 0.0008528026"  
##  [34] "Component \"6042316085_R01C01\": Mean relative difference: 0.001362491"   
##  [35] "Component \"6042316085_R01C02\": Mean relative difference: 0.0009247082"  
##  [36] "Component \"6042316085_R02C02\": Mean relative difference: 0.001318078"   
##  [37] "Component \"6042316085_R03C01\": Mean relative difference: 0.003765813"   
##  [38] "Component \"6042316085_R03C02\": Mean relative difference: 4.259089e-05"  
##  [39] "Component \"6042316085_R05C01\": Mean relative difference: 0.007578581"   
##  [40] "Component \"6042316085_R05C02\": Mean relative difference: 0.005001492"   
##  [41] "Component \"6042316085_R06C01\": Mean relative difference: 0.002166955"   
##  [42] "Component \"6042316107_R01C02\": Mean relative difference: 0.007225208"   
##  [43] "Component \"6042316107_R03C01\": Mean relative difference: 0.001221392"   
##  [44] "Component \"6042316107_R04C02\": Mean relative difference: 0.001265529"   
##  [45] "Component \"6042316107_R05C01\": Mean relative difference: 0.004903604"   
##  [46] "Component \"6042316107_R05C02\": Mean relative difference: 0.002330012"   
##  [47] "Component \"6042316107_R06C01\": Mean relative difference: 0.004369221"   
##  [48] "Component \"6042316107_R06C02\": Mean relative difference: 0.00441306"    
##  [49] "Component \"6042316113_R01C01\": Mean relative difference: 0.002037917"   
##  [50] "Component \"6042316113_R01C02\": Mean relative difference: 0.002717013"   
##  [51] "Component \"6042316113_R02C01\": Mean relative difference: 0.005187567"   
##  [52] "Component \"6042316113_R02C02\": Mean relative difference: 0.001776807"   
##  [53] "Component \"6042316113_R03C02\": Mean relative difference: 0.0001076938"  
##  [54] "Component \"6042316113_R04C02\": Mean relative difference: 0.01827261"    
##  [55] "Component \"6042316113_R05C01\": Mean relative difference: 0.001021883"   
##  [56] "Component \"6042316113_R05C02\": Mean relative difference: 0.004076606"   
##  [57] "Component \"6042316113_R06C01\": Mean relative difference: 0.007470594"   
##  [58] "Component \"6042316127_R01C01\": Mean relative difference: 0.0102167"     
##  [59] "Component \"6042316127_R01C02\": Mean relative difference: 0.005968567"   
##  [60] "Component \"6042316127_R02C02\": Mean relative difference: 0.01228684"    
##  [61] "Component \"6042316127_R03C01\": Mean relative difference: 0.003051516"   
##  [62] "Component \"6042316127_R03C02\": Mean relative difference: 0.005390373"   
##  [63] "Component \"6042316127_R04C01\": Mean relative difference: 0.00522226"    
##  [64] "Component \"6042316127_R04C02\": Mean relative difference: 0.001390841"   
##  [65] "Component \"6042316127_R05C02\": Mean relative difference: 0.00424916"    
##  [66] "Component \"6042316127_R06C01\": Mean relative difference: 0.00193334"    
##  [67] "Component \"6057825014_R01C02.1\": Mean relative difference: 0.004843924" 
##  [68] "Component \"6057825014_R02C02.1\": Mean relative difference: 0.008445282" 
##  [69] "Component \"6057825014_R03C01.1\": Mean relative difference: 0.001511774" 
##  [70] "Component \"6057825014_R03C02.1\": Mean relative difference: 0.003289338" 
##  [71] "Component \"6057825014_R04C01.1\": Mean relative difference: 0.002898635" 
##  [72] "Component \"6057825014_R04C02.1\": Mean relative difference: 0.002645651" 
##  [73] "Component \"6057825014_R06C01.1\": Mean relative difference: 0.003700896" 
##  [74] "Component \"6057825017_R01C01.1\": Mean relative difference: 0.002867279" 
##  [75] "Component \"6057825017_R01C02.1\": Mean relative difference: 0.004731186" 
##  [76] "Component \"6057825017_R02C01\": Mean relative difference: 0.004461493"   
##  [77] "Component \"6057825017_R02C02.1\": Mean relative difference: 0.0004087513"
##  [78] "Component \"6057825017_R03C01.1\": Mean relative difference: 0.0006703869"
##  [79] "Component \"6057825017_R03C02.1\": Mean relative difference: 0.002687469" 
##  [80] "Component \"6057825017_R04C01.1\": Mean relative difference: 0.006535025" 
##  [81] "Component \"6057825017_R04C02.1\": Mean relative difference: 0.0008907846"
##  [82] "Component \"6057825017_R05C01.1\": Mean relative difference: 0.001306738" 
##  [83] "Component \"6057825017_R05C02.1\": Mean relative difference: 0.005775531" 
##  [84] "Component \"6057825017_R06C01.1\": Mean relative difference: 0.0007369412"
##  [85] "Component \"6057825018_R01C01.1\": Mean relative difference: 0.004626421" 
##  [86] "Component \"6057825018_R02C02.1\": Mean relative difference: 0.005822739" 
##  [87] "Component \"6057825018_R03C01\": Mean relative difference: 0.0160181"     
##  [88] "Component \"6057825018_R03C02.1\": Mean relative difference: 0.005088012" 
##  [89] "Component \"6057825018_R04C01.1\": Mean relative difference: 0.01014491"  
##  [90] "Component \"6057825018_R04C02.1\": Mean relative difference: 0.009084185" 
##  [91] "Component \"6057825018_R05C01.1\": Mean relative difference: 0.006965893" 
##  [92] "Component \"6057825018_R05C02.1\": Mean relative difference: 0.00913593"  
##  [93] "Component \"6057825018_R06C01.1\": Mean relative difference: 0.008635721" 
##  [94] "Component \"6042316035_R01C01\": Mean relative difference: 0.003025003"   
##  [95] "Component \"6042316035_R02C02\": Mean relative difference: 0.03019016"    
##  [96] "Component \"6042316035_R03C01\": Mean relative difference: 0.002640202"   
##  [97] "Component \"6042316035_R03C02\": Mean relative difference: 0.000383893"   
##  [98] "Component \"6042316035_R04C01\": Mean relative difference: 0.01626422"    
##  [99] "Component \"6042316035_R05C01\": Mean relative difference: 0.01438542"    
## [100] "Component \"6042316035_R05C02\": Mean relative difference: 0.01881193"    
## [101] "Component \"6042316035_R06C02\": Mean relative difference: 0.004159761"   
## [102] "Component \"6042316048_R01C01\": Mean relative difference: 0.01097784"    
## [103] "Component \"6042316048_R01C02\": Mean relative difference: 0.01096577"    
## [104] "Component \"6042316048_R02C01\": Mean relative difference: 0.00318383"    
## [105] "Component \"6042316048_R03C01\": Mean relative difference: 0.01466546"    
## [106] "Component \"6042316048_R03C02\": Mean relative difference: 0.002607726"   
## [107] "Component \"6042316048_R04C01\": Mean relative difference: 0.01744718"    
## [108] "Component \"6042316048_R04C02\": Mean relative difference: 0.0126941"     
## [109] "Component \"6042316048_R05C01\": Mean relative difference: 0.0107092"     
## [110] "Component \"6042316048_R05C02\": Mean relative difference: 0.0118656"     
## [111] "Component \"6042316110_R01C01\": Mean relative difference: 0.007723471"   
## [112] "Component \"6042316110_R01C02\": Mean relative difference: 0.05702117"    
## [113] "Component \"6042316110_R03C02\": Mean relative difference: 0.006406333"   
## [114] "Component \"6042316110_R04C01\": Mean relative difference: 0.006380365"   
## [115] "Component \"6042316110_R04C02\": Mean relative difference: 0.01298917"    
## [116] "Component \"6042316110_R05C01\": Mean relative difference: 0.04585777"    
## [117] "Component \"6042316110_R06C01\": Mean relative difference: 0.00837521"    
## [118] "Component \"6042316121_R02C02\": Mean relative difference: 0.02108877"    
## [119] "Component \"6042316121_R03C01\": Mean relative difference: 0.01563759"    
## [120] "Component \"6042316121_R03C02\": Mean relative difference: 0.007783756"   
## [121] "Component \"6042316121_R04C02\": Mean relative difference: 0.02096482"    
## [122] "Component \"6042316121_R05C01\": Mean relative difference: 0.01314197"    
## [123] "Component \"6042316121_R05C02\": Mean relative difference: 0.02587878"    
## [124] "Component \"6042316121_R06C01\": Mean relative difference: 0.009945733"   
## [125] "Component \"6042316121_R06C02\": Mean relative difference: 0.01538254"    
## [126] "Component \"6042316066_R01C01\": Mean relative difference: 0.01685687"    
## [127] "Component \"6042316066_R02C01\": Mean relative difference: 0.003834141"   
## [128] "Component \"6042316066_R02C02\": Mean relative difference: 0.01037131"    
## [129] "Component \"6042316066_R03C01\": Mean relative difference: 0.03529848"    
## [130] "Component \"6042316066_R04C01\": Mean relative difference: 0.01447792"    
## [131] "Component \"6042316066_R04C02\": Mean relative difference: 0.0317652"     
## [132] "Component \"6042316066_R05C01\": Mean relative difference: 0.02495469"    
## [133] "Component \"6042316066_R06C01\": Mean relative difference: 0.01841931"    
## [134] "Component \"6042316066_R06C02\": Mean relative difference: 0.009140217"   
## [135] "Component \"6042316069_R01C01\": Mean relative difference: 0.0159006"     
## [136] "Component \"6042316069_R01C02\": Mean relative difference: 0.002352151"   
## [137] "Component \"6042316069_R02C01\": Mean relative difference: 0.003051261"   
## [138] "Component \"6042316069_R03C01\": Mean relative difference: 0.006767432"   
## [139] "Component \"6042316069_R03C02\": Mean relative difference: 0.02841319"    
## [140] "Component \"6042316069_R04C02\": Mean relative difference: 0.00327766"    
## [141] "Component \"6042316069_R05C01\": Mean relative difference: 0.001069052"   
## [142] "Component \"6042316069_R06C02\": Mean relative difference: 0.02115219"    
## [143] "Component \"6042316094_R01C02\": Mean relative difference: 0.01314671"    
## [144] "Component \"6042316094_R02C01\": Mean relative difference: 0.01489143"    
## [145] "Component \"6042316094_R03C02\": Mean relative difference: 0.00961152"    
## [146] "Component \"6042316094_R04C01\": Mean relative difference: 4.90452e-05"   
## [147] "Component \"6042316094_R04C02\": Mean relative difference: 0.03977844"    
## [148] "Component \"6042316094_R05C01\": Mean relative difference: 0.01663161"    
## [149] "Component \"6042316094_R05C02\": Mean relative difference: 0.01036295"    
## [150] "Component \"6042316094_R06C01\": Mean relative difference: 0.009467891"   
## [151] "Component \"6042316099_R01C01\": Mean relative difference: 0.01643653"    
## [152] "Component \"6042316099_R01C02\": Mean relative difference: 0.01000982"    
## [153] "Component \"6042316099_R02C01\": Mean relative difference: 0.01078948"    
## [154] "Component \"6042316099_R02C02\": Mean relative difference: 0.0153154"     
## [155] "Component \"6042316099_R03C01\": Mean relative difference: 0.03791758"    
## [156] "Component \"6042316099_R04C01\": Mean relative difference: 0.02912613"    
## [157] "Component \"6042316099_R04C02\": Mean relative difference: 0.02415069"    
## [158] "Component \"6042316099_R05C01\": Mean relative difference: 0.009167749"   
## [159] "Component \"6969568082_R02C01\": Mean relative difference: 0.02513291"    
## [160] "Component \"6969568082_R06C01\": Mean relative difference: 0.004244719"   
## [161] "Component \"6969568082_R02C02\": Mean relative difference: 0.01925031"    
## [162] "Component \"6969568082_R04C02\": Mean relative difference: 0.01996969"    
## [163] "Component \"6969568082_R06C02\": Mean relative difference: 0.01162994"    
## [164] "Component \"6969568084_R01C01\": Mean relative difference: 0.01993514"    
## [165] "Component \"6969568084_R02C01\": Mean relative difference: 0.03387861"    
## [166] "Component \"6969568084_R03C01\": Mean relative difference: 0.02322136"    
## [167] "Component \"6969568084_R04C01\": Mean relative difference: 0.003819405"   
## [168] "Component \"6969568084_R06C01\": Mean relative difference: 0.01559838"    
## [169] "Component \"6969568084_R02C02\": Mean relative difference: 0.008827488"   
## [170] "Component \"6969568084_R03C02\": Mean relative difference: 0.02940078"    
## [171] "Component \"6969568084_R04C02\": Mean relative difference: 0.007539981"   
## [172] "Component \"6969568084_R05C02\": Mean relative difference: 0.01822607"    
## [173] "Component \"6969568087_R01C01\": Mean relative difference: 0.01264853"    
## [174] "Component \"6969568087_R02C01\": Mean relative difference: 0.02503985"    
## [175] "Component \"6969568087_R03C01\": Mean relative difference: 0.00924478"    
## [176] "Component \"6969568087_R04C01\": Mean relative difference: 0.008402063"   
## [177] "Component \"6969568087_R05C01\": Mean relative difference: 0.03081664"    
## [178] "Component \"6969568087_R06C01\": Mean relative difference: 0.02708594"    
## [179] "Component \"6969568087_R01C02\": Mean relative difference: 0.008607383"   
## [180] "Component \"6969568087_R02C02\": Mean relative difference: 0.0002210055"  
## [181] "Component \"6969568087_R03C02\": Mean relative difference: 0.01827404"    
## [182] "Component \"6969568087_R05C02\": Mean relative difference: 0.02154393"    
## [183] "Component \"6969568087_R06C02\": Mean relative difference: 0.02786363"    
## [184] "Component \"6969568118_R01C01\": Mean relative difference: 0.01907614"    
## [185] "Component \"6969568118_R02C01\": Mean relative difference: 0.01828736"    
## [186] "Component \"6969568118_R03C01\": Mean relative difference: 0.006366801"   
## [187] "Component \"6969568118_R04C01\": Mean relative difference: 0.0277587"     
## [188] "Component \"6969568118_R01C02\": Mean relative difference: 0.02139702"    
## [189] "Component \"6969568118_R02C02\": Mean relative difference: 0.003980756"   
## [190] "Component \"6969568118_R03C02\": Mean relative difference: 0.01287797"    
## [191] "Component \"6969568118_R04C02\": Mean relative difference: 0.001491991"   
## [192] "Component \"6969568118_R06C02\": Mean relative difference: 0.02228438"    
## [193] "Component \"6929726046_R02C01\": Mean relative difference: 0.01317862"    
## [194] "Component \"6929726046_R05C01\": Mean relative difference: 0.01848101"    
## [195] "Component \"6929726046_R06C01\": Mean relative difference: 0.01135326"    
## [196] "Component \"6929726046_R01C02\": Mean relative difference: 0.02310383"    
## [197] "Component \"6929726046_R03C02\": Mean relative difference: 0.01553694"    
## [198] "Component \"6929726046_R04C02\": Mean relative difference: 0.03531404"    
## [199] "Component \"6929726046_R06C02\": Mean relative difference: 0.009023896"   
## [200] "Component \"6929718123_R01C01\": Mean relative difference: 0.02675321"    
## [201] "Component \"6929718123_R04C01\": Mean relative difference: 0.004522195"   
## [202] "Component \"6929718123_R06C01\": Mean relative difference: 0.02319077"    
## [203] "Component \"6929718123_R01C02\": Mean relative difference: 0.009055981"   
## [204] "Component \"6929718123_R02C02\": Mean relative difference: 0.02493802"    
## [205] "Component \"6929718123_R03C02\": Mean relative difference: 0.02150643"    
## [206] "Component \"6929718123_R04C02\": Mean relative difference: 0.01509121"    
## [207] "Component \"6929718123_R05C02\": Mean relative difference: 0.01596926"    
## [208] "Component \"6929718123_R06C02\": Mean relative difference: 0.01879209"    
## [209] "Component \"6929718136_R01C01\": Mean relative difference: 0.0250551"     
## [210] "Component \"6929718136_R02C01\": Mean relative difference: 0.01160423"    
## [211] "Component \"6929718136_R03C01\": Mean relative difference: 0.01299015"    
## [212] "Component \"6929718136_R04C01\": Mean relative difference: 0.0321602"     
## [213] "Component \"6929718136_R05C01\": Mean relative difference: 0.02004094"    
## [214] "Component \"6929718136_R02C02\": Mean relative difference: 0.006430265"   
## [215] "Component \"6929718136_R04C02\": Mean relative difference: 0.004677076"   
## [216] "Component \"6929718136_R05C02\": Mean relative difference: 0.01199318"    
## [217] "Component \"6929718136_R06C02\": Mean relative difference: 0.017871"      
## [218] "Component \"6929718138_R01C01\": Mean relative difference: 0.008461659"   
## [219] "Component \"6929718138_R02C01\": Mean relative difference: 0.01075589"    
## [220] "Component \"6929718138_R04C01\": Mean relative difference: 0.01258192"    
## [221] "Component \"6929718138_R06C01\": Mean relative difference: 0.01856663"    
## [222] "Component \"6929718138_R01C02\": Mean relative difference: 0.02602545"    
## [223] "Component \"6929718138_R03C02\": Mean relative difference: 0.01620736"    
## [224] "Component \"6929718138_R04C02\": Mean relative difference: 0.02190776"    
## [225] "Component \"6929718138_R05C02\": Mean relative difference: 0.01872934"    
## [226] "Component \"6929718138_R06C02\": Mean relative difference: 0.02427671"    
## [227] "Component \"6042316054_R02C01\": Mean relative difference: 0.03351534"    
## [228] "Component \"6042316054_R02C02\": Mean relative difference: 0.01431873"    
## [229] "Component \"6042316054_R03C01\": Mean relative difference: 0.01267266"    
## [230] "Component \"6042316054_R03C02\": Mean relative difference: 0.006259984"   
## [231] "Component \"6042316054_R04C01\": Mean relative difference: 0.003426156"   
## [232] "Component \"6042316054_R04C02\": Mean relative difference: 0.02453277"    
## [233] "Component \"6042316054_R05C01\": Mean relative difference: 0.02579635"    
## [234] "Component \"6042316054_R05C02\": Mean relative difference: 0.001743104"   
## [235] "Component \"6042316054_R06C01\": Mean relative difference: 0.02774079"    
## [236] "Component \"6042316063_R02C01\": Mean relative difference: 0.003403901"   
## [237] "Component \"6042316063_R02C02\": Mean relative difference: 0.01639184"    
## [238] "Component \"6042316063_R03C01\": Mean relative difference: 0.002266325"   
## [239] "Component \"6042316063_R03C02\": Mean relative difference: 0.004204157"   
## [240] "Component \"6042316063_R04C01\": Mean relative difference: 0.03452113"    
## [241] "Component \"6042316063_R04C02\": Mean relative difference: 0.00926267"    
## [242] "Component \"6042316063_R05C01\": Mean relative difference: 0.01848149"    
## [243] "Component \"6042316063_R05C02\": Mean relative difference: 0.005896312"   
## [244] "Component \"6042316063_R06C02\": Mean relative difference: 0.03555742"    
## [245] "Component \"6042316065_R01C02\": Mean relative difference: 0.003940634"   
## [246] "Component \"6042316065_R02C02\": Mean relative difference: 0.01532043"    
## [247] "Component \"6042316065_R03C01\": Mean relative difference: 0.00916733"    
## [248] "Component \"6042316065_R04C01\": Mean relative difference: 0.01724593"    
## [249] "Component \"6042316065_R04C02\": Mean relative difference: 0.003717116"   
## [250] "Component \"6042316065_R05C02\": Mean relative difference: 0.005905245"   
## [251] "Component \"6042316065_R06C02\": Mean relative difference: 0.02955138"    
## [252] "Component \"6042316103_R02C01\": Mean relative difference: 0.02326612"    
## [253] "Component \"6042316103_R03C01\": Mean relative difference: 0.03004319"    
## [254] "Component \"6042316103_R03C02\": Mean relative difference: 0.009468147"   
## [255] "Component \"6042316103_R04C01\": Mean relative difference: 0.01248312"    
## [256] "Component \"6042316103_R05C01\": Mean relative difference: 0.006465546"   
## [257] "Component \"6042316103_R06C01\": Mean relative difference: 0.008472114"   
## [258] "Component \"6042316103_R06C02\": Mean relative difference: 0.02022558"    
## [259] "Component \"6042316036_R01C02\": Mean relative difference: 0.01312759"    
## [260] "Component \"6042316036_R02C01\": Mean relative difference: 0.01442225"    
## [261] "Component \"6042316036_R03C01\": Mean relative difference: 0.02424502"    
## [262] "Component \"6042316036_R03C02\": Mean relative difference: 0.008205954"   
## [263] "Component \"6042316036_R04C01\": Mean relative difference: 0.01417237"    
## [264] "Component \"6042316036_R04C02\": Mean relative difference: 0.0202197"     
## [265] "Component \"6042316036_R05C01\": Mean relative difference: 0.002319912"   
## [266] "Component \"6042316036_R05C02\": Mean relative difference: 0.01764863"    
## [267] "Component \"6042316036_R06C01\": Mean relative difference: 0.001607823"   
## [268] "Component \"6042316036_R06C02\": Mean relative difference: 0.01573343"    
## [269] "Component \"6042316050_R01C02\": Mean relative difference: 0.008432375"   
## [270] "Component \"6042316050_R02C02\": Mean relative difference: 0.003918071"   
## [271] "Component \"6042316050_R03C01\": Mean relative difference: 0.007804211"   
## [272] "Component \"6042316050_R04C01\": Mean relative difference: 0.00890045"    
## [273] "Component \"6042316050_R04C02\": Mean relative difference: 0.02345964"    
## [274] "Component \"6042316050_R05C01\": Mean relative difference: 0.0252386"     
## [275] "Component \"6042316050_R05C02\": Mean relative difference: 0.01389564"    
## [276] "Component \"6042316050_R06C01\": Mean relative difference: 0.0256979"     
## [277] "Component \"6042316050_R06C02\": Mean relative difference: 0.0002492784"  
## [278] "Component \"6042316053_R01C01\": Mean relative difference: 0.003858654"   
## [279] "Component \"6042316053_R02C01\": Mean relative difference: 0.01213877"    
## [280] "Component \"6042316053_R02C02\": Mean relative difference: 0.01231578"    
## [281] "Component \"6042316053_R03C01\": Mean relative difference: 0.0079572"     
## [282] "Component \"6042316053_R03C02\": Mean relative difference: 0.004950719"   
## [283] "Component \"6042316053_R04C01\": Mean relative difference: 0.01976087"    
## [284] "Component \"6042316053_R04C02\": Mean relative difference: 0.002886665"   
## [285] "Component \"6042316053_R06C01\": Mean relative difference: 0.02517602"    
## [286] "Component \"6042316053_R06C02\": Mean relative difference: 0.02427661"    
## [287] "Component \"6042316061_R01C01\": Mean relative difference: 0.004490717"   
## [288] "Component \"6042316061_R02C01\": Mean relative difference: 0.007566767"   
## [289] "Component \"6042316061_R03C01\": Mean relative difference: 0.02190479"    
## [290] "Component \"6042316061_R03C02\": Mean relative difference: 0.003363988"   
## [291] "Component \"6042316061_R04C01\": Mean relative difference: 0.02447568"    
## [292] "Component \"6042316061_R04C02\": Mean relative difference: 0.001489334"   
## [293] "Component \"6042316061_R05C01\": Mean relative difference: 0.016579"      
## [294] "Component \"6042316061_R05C02\": Mean relative difference: 0.01630221"    
## [295] "Component \"7796806022_R01C01\": Mean relative difference: 0.03125971"    
## [296] "Component \"7796806022_R03C01\": Mean relative difference: 0.002838725"   
## [297] "Component \"7796806022_R04C01\": Mean relative difference: 0.004310711"   
## [298] "Component \"7796806022_R04C02\": Mean relative difference: 0.01027752"    
## [299] "Component \"7796806022_R05C02\": Mean relative difference: 0.002717047"   
## [300] "Component \"7796806022_R06C01\": Mean relative difference: 0.01284044"    
## [301] "Component \"7786923046_R03C02\": Mean relative difference: 0.007326031"   
## [302] "Component \"7786923046_R04C01\": Mean relative difference: 0.01964353"    
## [303] "Component \"7786923046_R06C01\": Mean relative difference: 0.02157744"    
## [304] "Component \"7796806038_R03C01\": Mean relative difference: 0.0006789345"  
## [305] "Component \"7796806038_R03C02\": Mean relative difference: 0.00478102"    
## [306] "Component \"7796806038_R05C02\": Mean relative difference: 0.01247494"    
## [307] "Component \"7786923063_R01C01\": Mean relative difference: 0.01086864"    
## [308] "Component \"7786923063_R01C02\": Mean relative difference: 0.00038396"    
## [309] "Component \"7786923063_R03C01\": Mean relative difference: 0.01197688"    
## [310] "Component \"7786923063_R04C01\": Mean relative difference: 0.009441984"   
## [311] "Component \"7786923063_R06C01\": Mean relative difference: 0.01197567"    
## [312] "Component \"7786923107_R01C02\": Mean relative difference: 0.01824832"    
## [313] "Component \"7786923107_R04C01\": Mean relative difference: 0.01004376"    
## [314] "Component \"7786923107_R05C01\": Mean relative difference: 0.0157187"     
## [315] "Component \"7786923107_R06C01\": Mean relative difference: 0.003626232"   
## [316] "Component \"7796806016_R03C02\": Mean relative difference: 0.01022651"    
## [317] "Component \"7796806016_R06C01\": Mean relative difference: 0.02114362"    
## [318] "Component \"7796806002_R01C01\": Mean relative difference: 0.02818048"    
## [319] "Component \"7796806002_R02C02\": Mean relative difference: 0.01372534"    
## [320] "Component \"7796806002_R04C02\": Mean relative difference: 0.01610265"    
## [321] "Component \"7796806002_R05C01\": Mean relative difference: 0.03183213"    
## [322] "Component \"7796806002_R05C02\": Mean relative difference: 0.03116065"    
## [323] "Component \"7796806002_R06C01\": Mean relative difference: 0.01677425"    
## [324] "Component \"7796806002_R06C02\": Mean relative difference: 0.04224427"    
## [325] "Component \"7796806029_R02C01\": Mean relative difference: 0.006535722"   
## [326] "Component \"7796806029_R04C01\": Mean relative difference: 0.01299031"    
## [327] "Component \"7796806029_R06C01\": Mean relative difference: 0.01351426"    
## [328] "Component \"6042316024_R01C01\": Mean relative difference: 0.009453011"   
## [329] "Component \"6042316024_R01C02\": Mean relative difference: 0.01249914"    
## [330] "Component \"6042316024_R02C01\": Mean relative difference: 0.01015168"    
## [331] "Component \"6042316024_R02C02\": Mean relative difference: 0.006930546"   
## [332] "Component \"6042316024_R03C01\": Mean relative difference: 0.006195781"   
## [333] "Component \"6042316024_R04C01\": Mean relative difference: 0.00218195"    
## [334] "Component \"6042316024_R04C02\": Mean relative difference: 0.005207145"   
## [335] "Component \"6042316024_R05C01\": Mean relative difference: 0.004617321"   
## [336] "Component \"6042316024_R05C02\": Mean relative difference: 0.002342651"   
## [337] "Component \"6042316024_R06C01\": Mean relative difference: 0.002423493"   
## [338] "Component \"6042316024_R06C02\": Mean relative difference: 0.0002915535"  
## [339] "Component \"6042316031_R01C01\": Mean relative difference: 0.001993511"   
## [340] "Component \"6042316031_R01C02\": Mean relative difference: 0.003142209"   
## [341] "Component \"6042316031_R02C01\": Mean relative difference: 0.004723449"   
## [342] "Component \"6042316031_R02C02\": Mean relative difference: 0.005092183"   
## [343] "Component \"6042316031_R03C01\": Mean relative difference: 0.003185017"   
## [344] "Component \"6042316031_R03C02\": Mean relative difference: 0.00231878"    
## [345] "Component \"6042316031_R04C01\": Mean relative difference: 0.0008693159"  
## [346] "Component \"6042316031_R04C02\": Mean relative difference: 0.000359399"   
## [347] "Component \"6042316031_R05C01\": Mean relative difference: 0.0001541431"  
## [348] "Component \"6042316031_R05C02\": Mean relative difference: 0.001237776"   
## [349] "Component \"6042316031_R06C01\": Mean relative difference: 0.002955389"   
## [350] "Component \"6042316031_R06C02\": Mean relative difference: 0.007783328"   
## [351] "Component \"6042316042_R01C01\": Mean relative difference: 0.003730153"   
## [352] "Component \"6042316042_R01C02\": Mean relative difference: 0.000152295"   
## [353] "Component \"6042316042_R02C01\": Mean relative difference: 0.0002119113"  
## [354] "Component \"6042316042_R02C02\": Mean relative difference: 0.0009922626"  
## [355] "Component \"6042316042_R03C01\": Mean relative difference: 0.005320828"   
## [356] "Component \"6042316042_R03C02\": Mean relative difference: 5.033644e-05"  
## [357] "Component \"6042316042_R04C01\": Mean relative difference: 0.00299015"    
## [358] "Component \"6042316042_R04C02\": Mean relative difference: 0.0006143021"  
## [359] "Component \"6042316042_R05C01\": Mean relative difference: 0.002770974"   
## [360] "Component \"6042316042_R05C02\": Mean relative difference: 0.001934715"   
## [361] "Component \"6042316042_R06C01\": Mean relative difference: 0.002450904"   
## [362] "Component \"6042316042_R06C02\": Mean relative difference: 0.004098867"   
## [363] "Component \"6042316047_R01C01\": Mean relative difference: 0.001572654"   
## [364] "Component \"6042316047_R02C01\": Mean relative difference: 0.002677129"   
## [365] "Component \"6042316047_R02C02\": Mean relative difference: 0.0009943293"  
## [366] "Component \"6042316047_R04C01\": Mean relative difference: 0.0001393208"  
## [367] "Component \"6042316047_R04C02\": Mean relative difference: 0.007402547"   
## [368] "Component \"6042316047_R05C01\": Mean relative difference: 0.002563325"   
## [369] "Component \"6042316047_R05C02\": Mean relative difference: 0.000978776"   
## [370] "Component \"6042316047_R06C01\": Mean relative difference: 0.003609061"   
## [371] "Component \"6055432012_R01C01\": Mean relative difference: 0.02857318"    
## [372] "Component \"6055432012_R01C02\": Mean relative difference: 0.01536161"    
## [373] "Component \"6055432012_R02C02\": Mean relative difference: 0.0197553"     
## [374] "Component \"6055432012_R03C01\": Mean relative difference: 0.01311016"    
## [375] "Component \"6055432012_R03C02\": Mean relative difference: 0.0004541998"  
## [376] "Component \"6055432012_R04C01\": Mean relative difference: 0.05807818"    
## [377] "Component \"6055432012_R04C02\": Mean relative difference: 0.01385685"    
## [378] "Component \"6055432012_R05C01\": Mean relative difference: 0.004909099"   
## [379] "Component \"6055432012_R05C02\": Mean relative difference: 0.0005888904"  
## [380] "Component \"6055432012_R06C01\": Mean relative difference: 0.0003715455"  
## [381] "Component \"6055432012_R06C02\": Mean relative difference: 0.001980122"   
## [382] "Component \"6055432029_R01C01\": Mean relative difference: 0.03261063"    
## [383] "Component \"6055432029_R01C02\": Mean relative difference: 0.001072346"   
## [384] "Component \"6055432029_R02C01\": Mean relative difference: 0.01477134"    
## [385] "Component \"6055432029_R02C02\": Mean relative difference: 0.002791673"   
## [386] "Component \"6055432029_R03C01\": Mean relative difference: 0.005164535"   
## [387] "Component \"6055432029_R03C02\": Mean relative difference: 0.0141598"     
## [388] "Component \"6055432029_R04C01\": Mean relative difference: 0.007454842"   
## [389] "Component \"6055432029_R04C02\": Mean relative difference: 0.0006432537"  
## [390] "Component \"6055432029_R05C01\": Mean relative difference: 0.01067081"    
## [391] "Component \"6055432029_R05C02\": Mean relative difference: 0.01089964"    
## [392] "Component \"6055432029_R06C01\": Mean relative difference: 0.000142306"   
## [393] "Component \"6055432029_R06C02\": Mean relative difference: 0.01504659"    
## [394] "Component \"6055432060_R01C01\": Mean relative difference: 0.004724758"   
## [395] "Component \"6055432060_R01C02\": Mean relative difference: 0.001465363"   
## [396] "Component \"6055432060_R02C01\": Mean relative difference: 0.007281945"   
## [397] "Component \"6055432060_R02C02\": Mean relative difference: 0.005147039"   
## [398] "Component \"6055432060_R03C01\": Mean relative difference: 0.011944"      
## [399] "Component \"6055432060_R03C02\": Mean relative difference: 0.0199792"     
## [400] "Component \"6055432060_R04C01\": Mean relative difference: 0.02353726"    
## [401] "Component \"6055432060_R04C02\": Mean relative difference: 0.002521748"   
## [402] "Component \"6055432060_R05C01\": Mean relative difference: 0.01848032"    
## [403] "Component \"6055432060_R05C02\": Mean relative difference: 0.01206478"    
## [404] "Component \"6055432060_R06C01\": Mean relative difference: 0.02464545"    
## [405] "Component \"6055432060_R06C02\": Mean relative difference: 0.01806241"    
## [406] "Component \"6055432066_R01C01\": Mean relative difference: 0.0009121482"  
## [407] "Component \"6055432066_R01C02\": Mean relative difference: 0.005167764"   
## [408] "Component \"6055432066_R02C01\": Mean relative difference: 0.004141654"   
## [409] "Component \"6055432066_R02C02\": Mean relative difference: 0.008431437"   
## [410] "Component \"6055432066_R03C02\": Mean relative difference: 0.01526689"    
## [411] "Component \"6055432066_R04C01\": Mean relative difference: 0.008475197"   
## [412] "Component \"6055432066_R04C02\": Mean relative difference: 0.001524354"   
## [413] "Component \"6055432066_R05C01\": Mean relative difference: 0.004307322"   
## [414] "Component \"6055432066_R05C02\": Mean relative difference: 0.001269072"   
## [415] "Component \"6057825115_R01C01\": Mean relative difference: 0.005915681"   
## [416] "Component \"6057825115_R03C01\": Mean relative difference: 0.01032372"    
## [417] "Component \"6057825115_R05C01\": Mean relative difference: 0.003156536"   
## [418] "Component \"6057825115_R01C02\": Mean relative difference: 0.006668261"   
## [419] "Component \"6057825115_R03C02\": Mean relative difference: 0.006865827"   
## [420] "Component \"6057825115_R05C02\": Mean relative difference: 0.0005121951"  
## [421] "Component \"6057825128_R01C01\": Mean relative difference: 0.009882778"   
## [422] "Component \"6057825128_R03C01\": Mean relative difference: 0.009766438"   
## [423] "Component \"6057825128_R05C01\": Mean relative difference: 0.004559136"   
## [424] "Component \"6057825115_R02C01\": Mean relative difference: 0.0006340263"  
## [425] "Component \"6057825115_R04C01\": Mean relative difference: 0.001216475"   
## [426] "Component \"6057825115_R06C01\": Mean relative difference: 0.004890784"   
## [427] "Component \"6057825115_R02C02\": Mean relative difference: 0.004753345"   
## [428] "Component \"6057825115_R04C02\": Mean relative difference: 0.0006350344"  
## [429] "Component \"6057825115_R06C02\": Mean relative difference: 0.00410318"    
## [430] "Component \"6057825128_R01C02\": Mean relative difference: 0.004186846"   
## [431] "Component \"6057825128_R03C02\": Mean relative difference: 0.0007479567"  
## [432] "Component \"6057825128_R05C02\": Mean relative difference: 0.003928646"
```


To make sure we do not induce any NAs into the dataset when we convert the beta values back M-values (by log2 transformation), we need to ensure we do not have any corrected beta values that are greater or equal to zero or any beta values that are greater than 1.


```r
library(lumi)
adj.residuals[adj.residuals<=0]<-0.001 # convert any values that are less than or equal to zero to 0.001
adj.residuals[adj.residuals>1]<-0.999 # convert any values that are greater than 1 to 0.999
adj.M.values<-beta2m(adj.residuals)
any(is.na(adj.M.values)) # should be FALSE indicating there are no NAs
```

```
## [1] FALSE
```

Save corrected dataset:

```r
GSE43414_cell_cor<-adj.residuals
save(GSE43414_cell_cor, file="GSE43414_cell_cor.RData")
```

## PCA Scree Heatmap for cell-corrected data


```r
## PCA
load("GSE43414_cell_cor.RData")
cell.cor.dat <- t(scale(t(as.matrix(GSE43414_cell_cor))))
PCA_full<-princomp(cell.cor.dat[complete.cases(cell.cor.dat),])
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
adjust<-1-Importance[1]
pca_adjusted<-Importance[2:length(Importance)]/adjust
pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))

#Specify which covariates are categorical and/or categorical
colnames(meta)
```

```
##  [1] "series_id"         "gsm"               "Subject"          
##  [4] "barcode"           "lunnon.et.al"      "tissue.code"      
##  [7] "braak.stage"       "Sex"               "ad.disease.status"
## [10] "age.brain"         "age.blood"         "Tissue"           
## [13] "Neuron"            "Glia"              "chip"             
## [16] "row"
```

```r
meta_categorical<-meta[,c("ad.disease.status", "braak.stage", "Sex", "Tissue", "chip", "row")]  # input column numbers in meta that contain categorical variables
meta_continuous<-meta[,c("age.brain", "Neuron")] # input column numbers in meta that contain continuous variables
#meta_continuous<-data.frame(meta_continuous)

# Specify the number of PCs you want shown (usually # of samples in the dataset)
Num<-20

# Designate what order you want the variables to appear (continuous variables rbinded to categorical variables in function)
Order<-c(4,8,7,3,1,2,5,6)

#Apply function on PCA results, pulls in the meta data and beta values from above
heat_scree_plot(Loadings, Importance, Num, Order)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

Cell type correction reduced the effect of Neuron and brought back some effects of age, AD disease status, and braak stage, all of which we will correct for by using them as covariates in our linear model for differential methylation analysis.
