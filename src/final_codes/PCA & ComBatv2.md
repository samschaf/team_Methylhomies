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
```

```
## Error: package 'mgcv' is required by 'sva' so will not be detached
```

```r
library(sva)

#Correction for row
#row <- meta$row
modcombat <- model.matrix(~1,data=meta)
combat_edata <- ComBat(dat=GSE43414_BMIQ, batch=row, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
#Correction for chip
chip <- meta$chip
GSE43414_batch_cor <- ComBat(dat=combat_edata, batch=chip, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
```

```
## Found 50 batches
## Adjusting for 0 covariate(s) or covariate level(s)
```

```
## Error in ComBat(dat = combat_edata, batch = chip, mod = modcombat, par.prior = TRUE, : object 'combat_edata' not found
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
## 6057825008_R02C02 0.4147119 0.5852881
## 6057825008_R03C01 0.4171618 0.5828382
## 6057825008_R04C01 0.3917295 0.6082705
## 6057825008_R04C02 0.3951055 0.6048945
## 6057825008_R05C01 0.3844469 0.6155531
## 6057825008_R05C02 0.4269697 0.5730303
```

```r
write.csv(prop, file = "cellprop_batch_cor.csv", row.names=T)
summary(prop)
```

```
##      neuron            glia       
##  Min.   :0.1221   Min.   :0.4433  
##  1st Qu.:0.3594   1st Qu.:0.5706  
##  Median :0.4015   Median :0.5985  
##  Mean   :0.3920   Mean   :0.6080  
##  3rd Qu.:0.4294   3rd Qu.:0.6406  
##  Max.   :0.5567   Max.   :0.8779
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
## cg00000029       0.009525003       0.070592860       -0.02020091
## cg00000108      -0.117069798       0.033211730        0.01829410
## cg00000165       0.043429293      -0.090753919       -0.01035697
## cg00000236       0.041374332       0.001453763        0.03233354
## cg00000289       0.046886205      -0.028523841       -0.01437302
## cg00000292      -0.027014327      -0.029418274        0.06642361
##            6057825008_R04C02 6057825008_R05C01 6057825008_R05C02
## cg00000029      -0.003750301        0.02646673       -0.06008314
## cg00000108       0.082163474       -0.02339459       -0.04736922
## cg00000165      -0.125710780        0.09519386        0.06607353
## cg00000236      -0.036711703        0.02101008       -0.06365987
## cg00000289      -0.063298431       -0.10216312        0.01865110
## cg00000292      -0.044049634        0.01234043       -0.06749460
##            6057825008_R06C01 6057825008_R06C02 6057825014_R01C02
## cg00000029        0.02070551       0.045276210       -0.01950453
## cg00000108       -0.03476183      -0.004723149       -0.08412560
## cg00000165        0.07433291      -0.096821319       -0.08381397
## cg00000236       -0.04784718       0.024479148       -0.01295822
## cg00000289        0.08609512       0.051569793        0.03974945
## cg00000292        0.01291889       0.101145591        0.01136739
##            6057825014_R02C02 6057825014_R03C01 6057825014_R03C02
## cg00000029       0.008528736       0.004490618      -0.030596665
## cg00000108       0.151729631      -0.052246492      -0.214258908
## cg00000165       0.113409999       0.103990335       0.075096902
## cg00000236      -0.011680246       0.035230349       0.039064304
## cg00000289       0.001038367      -0.023755029       0.009434159
## cg00000292      -0.016428059      -0.026371647       0.061382498
##            6057825014_R04C01 6057825014_R04C02 6057825014_R06C01
## cg00000029        0.10539011       -0.04224952      -0.099975952
## cg00000108       -0.02014499       -0.03186148      -0.003000984
## cg00000165       -0.06114973        0.14566320       0.059481996
## cg00000236       -0.05634716       -0.01326680      -0.095904679
## cg00000289       -0.02435663       -0.04572934       0.029470718
## cg00000292        0.02430738       -0.04057204       0.001359941
##            6057825017_R01C01 6057825017_R01C02 6057825017_R02C02
## cg00000029        0.08804570      -0.031090323       0.040692226
## cg00000108       -0.12775602       0.017257773       0.080528130
## cg00000165       -0.11587112       0.022131097      -0.007135682
## cg00000236        0.03290744       0.020881644       0.025324012
## cg00000289        0.10290121      -0.064156472       0.026289021
## cg00000292       -0.05173063       0.008541188       0.034441848
##            6057825017_R03C01 6057825017_R03C02 6057825017_R04C01
## cg00000029     -0.0011821070      0.0148520477       0.002974855
## cg00000108      0.0282462288     -0.1322070423       0.073503851
## cg00000165     -0.0009758277      0.1920224731       0.040114299
## cg00000236     -0.0013513927      0.0003334903      -0.007640753
## cg00000289      0.0856283451     -0.0183886830      -0.034464116
## cg00000292     -0.0280134516     -0.0060917142      -0.090596510
##            6057825017_R04C02 6057825017_R05C01 6057825017_R05C02
## cg00000029      -0.006416495        0.04343863       0.071466975
## cg00000108       0.141095239        0.04995167      -0.017997835
## cg00000165      -0.156146594       -0.07557900      -0.088165059
## cg00000236       0.042466673       -0.01767145      -0.015586935
## cg00000289      -0.017240040       -0.01944490       0.064615498
## cg00000292       0.033921118        0.03184537       0.007823554
##            6057825017_R06C01 6057825018_R01C01 6057825018_R02C02
## cg00000029       0.085686419        0.04249866     -0.0288589312
## cg00000108       0.148088702       -0.15151608     -0.0599274790
## cg00000165       0.036376760        0.04396447      0.1304698704
## cg00000236       0.007228821        0.01845584      0.0370045029
## cg00000289       0.008153574       -0.03078267     -0.0005233982
## cg00000292      -0.050425370        0.02748958      0.0008736751
##            6057825018_R03C02 6057825018_R04C01 6057825018_R04C02
## cg00000029      -0.076570924       -0.04599775       -0.04157300
## cg00000108      -0.010148499       -0.02369679       -0.10708381
## cg00000165       0.143388604       -0.11767706       -0.01554385
## cg00000236      -0.004799001       -0.06239921        0.03858205
## cg00000289      -0.005547679       -0.01658293        0.05433695
## cg00000292      -0.014964945       -0.05567064       -0.02853222
##            6057825018_R05C01 6057825018_R05C02 6057825018_R06C01
## cg00000029        0.07304143       0.043160217        0.01703854
## cg00000108        0.10599255      -0.036241437        0.07644193
## cg00000165       -0.08634915      -0.085766504       -0.09542020
## cg00000236       -0.03154178       0.006836868       -0.00214790
## cg00000289       -0.05522182       0.056812723        0.04940359
## cg00000292        0.09398270      -0.017998829        0.02060347
##            6042316085_R01C01 6042316085_R01C02 6042316085_R02C02
## cg00000029       -0.09719197        0.05552285        0.15107259
## cg00000108       -0.22244582        0.02183267        0.05313474
## cg00000165        0.22119583       -0.01017307       -0.01415219
## cg00000236        0.05502685        0.02721486       -0.04982949
## cg00000289        0.04668555        0.03379190       -0.14463241
## cg00000292        0.01334482       -0.08146066        0.10716932
##            6042316085_R03C01 6042316085_R03C02 6042316085_R05C01
## cg00000029      -0.026631206      -0.011340090      -0.047415907
## cg00000108      -0.011168215      -0.020170853      -0.023236156
## cg00000165      -0.053778022      -0.033255135       0.006093443
## cg00000236      -0.028160989      -0.012160708      -0.048273005
## cg00000289      -0.001240425      -0.002853507      -0.008966356
## cg00000292      -0.009707560      -0.007880332      -0.032917376
##            6042316085_R05C02 6042316085_R06C01 6042316107_R01C02
## cg00000029      -0.019017608        0.03622185       -0.02225448
## cg00000108      -0.094166732        0.14812950       -0.09665097
## cg00000165      -0.019394564        0.12199836       -0.06170949
## cg00000236       0.054189873       -0.03126063       -0.01171074
## cg00000289      -0.007629746       -0.01017645        0.04228235
## cg00000292      -0.043622598        0.03627230        0.05155829
##            6042316107_R03C01 6042316107_R04C02 6042316107_R05C01
## cg00000029       0.022731110       -0.01052349      -0.003984181
## cg00000108       0.004247623        0.01421661       0.027103449
## cg00000165       0.136132162       -0.03936116      -0.035054320
## cg00000236      -0.016243650        0.01524351       0.073650870
## cg00000289       0.048008538       -0.02303742      -0.055974738
## cg00000292      -0.028735166       -0.03405659       0.024274448
##            6042316107_R05C02 6042316107_R06C01 6042316107_R06C02
## cg00000029      -0.004994045       0.073739097        0.06588930
## cg00000108      -0.029209466      -0.176092951       -0.11691818
## cg00000165      -0.114338451      -0.089675052       -0.05647480
## cg00000236      -0.034954865      -0.037851111       -0.05756172
## cg00000289       0.013857152       0.007359432       -0.04935289
## cg00000292      -0.023583906       0.028850040       -0.01709840
##            6042316113_R01C01 6042316113_R01C02 6042316113_R02C01
## cg00000029       0.007752629      -0.041975052       0.048453604
## cg00000108       0.067687028      -0.078000765       0.155545152
## cg00000165       0.090166443      -0.004031192      -0.010875829
## cg00000236      -0.005789062       0.029528444      -0.002316321
## cg00000289       0.061982840       0.004933197       0.081344763
## cg00000292      -0.034003346      -0.010732079      -0.073692043
##            6042316113_R02C02 6042316113_R03C02 6042316113_R04C02
## cg00000029      -0.090541269      0.0442373974       -0.17815450
## cg00000108       0.117299391      0.0059834651        0.20208603
## cg00000165       0.195798540      0.0983321850       -0.14405031
## cg00000236       0.008707109      0.0268728815        0.05772451
## cg00000289      -0.081895401     -0.0323873871       -0.11984625
## cg00000292       0.111691895     -0.0004275105        0.10063725
##            6042316113_R05C01 6042316113_R05C02 6042316113_R06C01
## cg00000029       0.079465849        0.01836898        0.02171104
## cg00000108      -0.037399212       -0.03464645       -0.03929138
## cg00000165      -0.056596094        0.10260072       -0.01317311
## cg00000236       0.018014504       -0.02441627       -0.03133315
## cg00000289       0.009909283        0.05552216        0.03652176
## cg00000292       0.010984339        0.03249368        0.03370825
##            6042316127_R01C01 6042316127_R01C02 6042316127_R02C02
## cg00000029        0.03274857        0.04283693        0.02957265
## cg00000108       -0.07870200       -0.11741419       -0.08008233
## cg00000165       -0.03616471        0.19691793        0.04083389
## cg00000236       -0.04870725       -0.01569955        0.01715022
## cg00000289       -0.03247032        0.08599833       -0.00509616
## cg00000292       -0.06428007       -0.05255017       -0.04558469
##            6042316127_R03C01 6042316127_R03C02 6042316127_R04C01
## cg00000029      0.0834048643        0.02065456        0.04453020
## cg00000108     -0.0226094498       -0.06258837        0.01758130
## cg00000165      0.0865758326        0.00283209       -0.01660666
## cg00000236      0.0083322795       -0.03629142       -0.04214904
## cg00000289     -0.0004878426        0.05748738       -0.02713189
## cg00000292     -0.0285695300        0.03648759       -0.02410976
##            6042316127_R04C02 6042316127_R05C02 6042316127_R06C01
## cg00000029      -0.015398413       -0.06681228      -0.043203320
## cg00000108       0.123316662       -0.10548952       0.044503149
## cg00000165       0.125812761       -0.17220985      -0.017047600
## cg00000236      -0.019588399       -0.03687227       0.005863133
## cg00000289      -0.058926982        0.00773777       0.049781665
## cg00000292       0.009867723       -0.03041628      -0.041224940
##            6057825014_R01C02.1 6057825014_R02C02.1 6057825014_R03C01.1
## cg00000029         0.031917761          0.05638226         0.080865080
## cg00000108        -0.041749655          0.09797152        -0.019518976
## cg00000165        -0.076684620          0.11318119         0.091010619
## cg00000236         0.002212422          0.03322435         0.052470035
## cg00000289         0.015415627         -0.01163340         0.041722861
## cg00000292         0.040590482         -0.04602570         0.008031219
##            6057825014_R03C02.1 6057825014_R04C01.1 6057825014_R04C02.1
## cg00000029         -0.02281221         0.117632350        -0.075908484
## cg00000108         -0.15540110        -0.023154622         0.021998323
## cg00000165          0.04304310        -0.091592966         0.150399250
## cg00000236          0.04230859        -0.042360012        -0.020906979
## cg00000289         -0.04186866        -0.077246825        -0.001903762
## cg00000292          0.03633316        -0.007218827        -0.011994785
##            6057825014_R06C01.1 6057825017_R01C01.1 6057825017_R01C02.1
## cg00000029         -0.07506793          0.06832150         -0.05149303
## cg00000108          0.03200831         -0.15479382          0.09158477
## cg00000165          0.04000881         -0.14688871         -0.03326245
## cg00000236         -0.08361028          0.01858319          0.04170527
## cg00000289         -0.01640581          0.03963063         -0.02317868
## cg00000292         -0.01476837         -0.02826531         -0.01211397
##            6057825017_R02C01 6057825017_R02C02.1 6057825017_R03C01.1
## cg00000029       0.043362693         0.039898293         0.009386645
## cg00000108      -0.022924534         0.120952478         0.069726937
## cg00000165       0.091824951        -0.044813829        -0.054675744
## cg00000236      -0.001492869         0.019228942         0.027348617
## cg00000289      -0.029278800         0.032703964         0.104526047
## cg00000292       0.036461802         0.006007878         0.014406056
##            6057825017_R03C02.1 6057825017_R04C01.1 6057825017_R04C02.1
## cg00000029        -0.025125048          0.02963278         0.030729296
## cg00000108        -0.169585599         -0.01674188         0.023115692
## cg00000165         0.147032627          0.14705201        -0.092260997
## cg00000236         0.040673512         -0.03787696         0.019374024
## cg00000289        -0.010756458          0.03145070        -0.008581487
## cg00000292         0.009855489         -0.05659835        -0.006909661
##            6057825017_R05C01.1 6057825017_R05C02.1 6057825017_R06C01.1
## cg00000029         0.071964359        0.0589237214         0.053040985
## cg00000108        -0.013062544       -0.0371919243         0.018603804
## cg00000165        -0.029031799       -0.0533190901         0.038485383
## cg00000236         0.008610090        0.0006341466         0.006363153
## cg00000289         0.029177879       -0.0015729210        -0.005699768
## cg00000292        -0.006385278       -0.0089354018        -0.047437859
##            6057825018_R01C01.1 6057825018_R02C02.1 6057825018_R03C01
## cg00000029         0.035146913        -0.010280036       -0.17936619
## cg00000108        -0.125522042        -0.072466137        0.21820914
## cg00000165        -0.004597825         0.086860927       -0.14174452
## cg00000236         0.020131141         0.044834190       -0.02606975
## cg00000289        -0.008343099        -0.037626024       -0.02234308
## cg00000292         0.003273169        -0.008203087        0.13913289
##            6057825018_R03C02.1 6057825018_R04C01.1 6057825018_R04C02.1
## cg00000029         -0.01521060          0.03707922          0.03158538
## cg00000108         -0.03904987         -0.10716449         -0.15162937
## cg00000165          0.08217242         -0.02466435          0.06347393
## cg00000236         -0.02890734         -0.04964353          0.02080019
## cg00000289          0.04613348         -0.03421331          0.04330305
## cg00000292         -0.01580886         -0.06809120         -0.03472898
##            6057825018_R05C01.1 6057825018_R05C02.1 6057825018_R06C01.1
## cg00000029          0.09900405         0.081827357          0.07329532
## cg00000108         -0.05748737        -0.080734887         -0.03371788
## cg00000165         -0.01254643        -0.048890378         -0.05919835
## cg00000236         -0.01187645        -0.004402863         -0.02291163
## cg00000289         -0.04869660         0.068604525          0.02292336
## cg00000292          0.02257524        -0.034820985         -0.02080653
##            6042316035_R01C01 6042316035_R02C02 6042316035_R03C01
## cg00000029      -0.049060682      -0.036440860      -0.087213760
## cg00000108       0.086490814       0.065209383       0.104080708
## cg00000165       0.006339355       0.006320546       0.005422119
## cg00000236       0.030391019       0.003846819       0.028632463
## cg00000289      -0.026943194      -0.011876653      -0.028578371
## cg00000292       0.001483787       0.054734040       0.046506698
##            6042316035_R03C02 6042316035_R04C01 6042316035_R05C01
## cg00000029       -0.05780442        0.02326641      -0.016514846
## cg00000108        0.10560050        0.01417517       0.024645606
## cg00000165        0.02444433       -0.03221666       0.057509739
## cg00000236        0.02318799        0.00323586       0.001942335
## cg00000289        0.05570432        0.01642082      -0.029376415
## cg00000292        0.03547259        0.01320370       0.009321514
##            6042316035_R05C02 6042316035_R06C02 6042316048_R01C01
## cg00000029      -0.006131913       -0.01561568      -0.005520239
## cg00000108      -0.004153000        0.02703421      -0.006234478
## cg00000165      -0.083264071        0.04240338       0.045165038
## cg00000236      -0.044733586        0.02419335      -0.037472195
## cg00000289       0.067579451       -0.02496540      -0.075663340
## cg00000292      -0.066072644        0.02344693       0.029328224
##            6042316048_R01C02 6042316048_R02C01 6042316048_R03C01
## cg00000029      -0.039011460       -0.02898324        0.05543540
## cg00000108       0.014251859        0.03850295       -0.01965111
## cg00000165       0.023360221       -0.03853800       -0.07377015
## cg00000236       0.042485744       -0.02807261        0.05230757
## cg00000289       0.004778309        0.06595653       -0.04020735
## cg00000292       0.007308333        0.05773239       -0.03458657
##            6042316048_R03C02 6042316048_R04C01 6042316048_R04C02
## cg00000029      -0.009872868      -0.068315753      -0.049785709
## cg00000108       0.034208437       0.048144427       0.013492525
## cg00000165      -0.040175409      -0.009288758      -0.109164338
## cg00000236       0.025449935      -0.017777147      -0.006556865
## cg00000289       0.004093321      -0.012243076      -0.046392201
## cg00000292       0.047032106       0.025451042       0.045298800
##            6042316048_R05C01 6042316048_R05C02 6042316110_R01C01
## cg00000029       -0.02271054       0.026638045      0.0046315551
## cg00000108       -0.02434734       0.021613326      0.0052940808
## cg00000165        0.08194382       0.005198342     -0.0286472939
## cg00000236       -0.02137920      -0.005519320      0.0007228938
## cg00000289        0.02008266       0.012418785     -0.0340581253
## cg00000292        0.01377061      -0.095070271     -0.0109272393
##            6042316110_R01C02 6042316110_R03C02 6042316110_R04C01
## cg00000029       0.054900957      -0.062339006       0.055303674
## cg00000108      -0.051143485       0.004524739       0.026604541
## cg00000165      -0.008450474       0.022417979       0.018766013
## cg00000236       0.035292910       0.036471946      -0.025887989
## cg00000289       0.013222321       0.050821073       0.005605406
## cg00000292       0.056118698      -0.069385652       0.004581783
##            6042316110_R04C02 6042316110_R05C01 6042316110_R06C01
## cg00000029       -0.01027883        0.02734340       0.078976209
## cg00000108        0.03395785       -0.03777229      -0.009109834
## cg00000165       -0.03076913       -0.01721035      -0.052670773
## cg00000236        0.00276445       -0.04125343       0.019430163
## cg00000289       -0.04680362        0.03196219       0.030843493
## cg00000292       -0.01100838       -0.00316421      -0.007188380
##            6042316121_R02C02 6042316121_R03C01 6042316121_R03C02
## cg00000029       -0.05120453      -0.015411142       -0.08141057
## cg00000108        0.01969382      -0.021577219        0.06211974
## cg00000165        0.04358720       0.004688276       -0.01991610
## cg00000236       -0.01671072       0.047577237       -0.01398279
## cg00000289        0.05263570       0.004755320       -0.02494410
## cg00000292        0.02910670      -0.021052680       -0.01933512
##            6042316121_R04C02 6042316121_R05C01 6042316121_R05C02
## cg00000029      -0.030332063      8.814893e-06       0.080502297
## cg00000108       0.039099859      4.689306e-02       0.018479040
## cg00000165       0.003556631     -1.937982e-02       0.005966361
## cg00000236       0.022187900     -2.427981e-02       0.030607682
## cg00000289      -0.046711232     -2.245195e-02       0.045391647
## cg00000292       0.078050060      2.417318e-02       0.081295331
##            6042316121_R06C01 6042316121_R06C02 6042316066_R01C01
## cg00000029      -0.007794764     -0.0758749470     -0.0683946235
## cg00000108       0.035573506      0.0001616167      0.0148435665
## cg00000165       0.063105089      0.0461163537     -0.0588040944
## cg00000236       0.037291477     -0.0034616597      0.0003335628
## cg00000289      -0.004512215      0.0568211658     -0.0332486203
## cg00000292      -0.053207751     -0.0177492230     -0.0067999439
##            6042316066_R02C01 6042316066_R02C02 6042316066_R03C01
## cg00000029       -0.04410021       0.059200999        0.04353232
## cg00000108        0.01280496       0.037307588       -0.02071126
## cg00000165       -0.07975278      -0.059859490        0.03715714
## cg00000236        0.03415538       0.009681529       -0.07039528
## cg00000289       -0.10301964       0.036569707       -0.04322132
## cg00000292        0.01190463      -0.009961152        0.03651826
##            6042316066_R04C01 6042316066_R04C02 6042316066_R05C01
## cg00000029      -0.056975792       0.022939295       0.063592867
## cg00000108      -0.052760430       0.040433871      -0.020446926
## cg00000165      -0.006762107      -0.065918480       0.017068729
## cg00000236       0.008305549       0.030809606       0.007348845
## cg00000289       0.066791377       0.055909350      -0.043235645
## cg00000292      -0.040931499       0.005083365      -0.023120329
##            6042316066_R06C01 6042316066_R06C02 6042316069_R01C01
## cg00000029     -0.0008697732      -0.037554124       -0.09922446
## cg00000108      0.0027382262       0.010707807        0.04112504
## cg00000165     -0.0574480602       0.112072260       -0.04421191
## cg00000236     -0.0765090492      -0.004351935        0.05177339
## cg00000289     -0.0278627908       0.080419410        0.04901919
## cg00000292      0.0773162777       0.025297960        0.01731084
##            6042316069_R01C02 6042316069_R02C01 6042316069_R03C01
## cg00000029      -0.079301616        0.01845019      0.0158010519
## cg00000108       0.007591255       -0.01262456      0.0104953027
## cg00000165      -0.022995319       -0.01107759     -0.0100117341
## cg00000236       0.027195684       -0.02770830      0.0017733727
## cg00000289      -0.015478116       -0.07155663     -0.0392783438
## cg00000292      -0.020619420       -0.06541270     -0.0002763331
##            6042316069_R03C02 6042316069_R04C02 6042316069_R05C01
## cg00000029       0.041479408       -0.03071540       0.005219646
## cg00000108      -0.034277675        0.01559627       0.021824287
## cg00000165       0.013501164       -0.05940247       0.025609491
## cg00000236      -0.009474948        0.04210027      -0.010210726
## cg00000289       0.001539494       -0.02197614      -0.020667219
## cg00000292       0.004611371        0.03237282       0.007638856
##            6042316069_R06C02 6042316094_R01C02 6042316094_R02C01
## cg00000029      0.0364078960       -0.10034367        0.01732613
## cg00000108     -0.0089452610       -0.02026066        0.06619795
## cg00000165      0.0634729623       -0.07928177       -0.01925757
## cg00000236      0.0249922801       -0.02695448       -0.03777191
## cg00000289      0.0003322975        0.02152556       -0.05920176
## cg00000292      0.0788700836       -0.10595913       -0.02280971
##            6042316094_R03C02 6042316094_R04C01 6042316094_R04C02
## cg00000029      -0.002811610      -0.034302477       0.083511958
## cg00000108      -0.003339336       0.023819745      -0.054654359
## cg00000165      -0.036686323       0.024009574      -0.027369375
## cg00000236       0.001662890       0.030280855       0.006077746
## cg00000289       0.116273543      -0.023196227      -0.009711256
## cg00000292      -0.002201441       0.009260314      -0.037642196
##            6042316094_R05C01 6042316094_R05C02 6042316094_R06C01
## cg00000029      -0.038626529      -0.049603375        0.03054594
## cg00000108       0.048217709      -0.010381090       -0.05335283
## cg00000165       0.060925219      -0.001458605       -0.04259399
## cg00000236       0.006718013       0.017200678       -0.01234153
## cg00000289       0.044408645       0.015242688        0.04069089
## cg00000292       0.014404850       0.032998799       -0.04388138
##            6042316099_R01C01 6042316099_R01C02 6042316099_R02C01
## cg00000029      -0.010332939      -0.064815331      -0.050362144
## cg00000108       0.016948450       0.040374604       0.037622835
## cg00000165       0.035943882      -0.037814725      -0.016144302
## cg00000236      -0.007837599       0.001716362       0.033412253
## cg00000289      -0.075339746      -0.026669607       0.006907667
## cg00000292       0.051376229      -0.030419242      -0.009953156
##            6042316099_R02C02 6042316099_R03C01 6042316099_R04C01
## cg00000029        0.02060274        0.06015589     -0.0021143736
## cg00000108        0.05587124       -0.02895343      0.0215609192
## cg00000165       -0.03588658        0.02792855     -0.0146016763
## cg00000236       -0.02891090        0.01411651      0.0004767018
## cg00000289        0.01408133        0.03048087      0.0391528743
## cg00000292        0.03446023        0.03268745      0.0603007041
##            6042316099_R04C02 6042316099_R05C01 6969568082_R02C01
## cg00000029      -0.071267937      -0.044964617      -0.003756796
## cg00000108       0.035061093       0.033776659       0.002763582
## cg00000165      -0.009887653      -0.013236280      -0.005303051
## cg00000236       0.043875176      -0.030983041      -0.042526827
## cg00000289       0.044750303      -0.002394082      -0.034778306
## cg00000292      -0.035378274       0.034532276       0.007959133
##            6969568082_R06C01 6969568082_R02C02 6969568082_R04C02
## cg00000029      -0.052449662      0.0280432012      -0.016906299
## cg00000108       0.023327013     -0.0132539616      -0.022028171
## cg00000165       0.023505949     -0.0950176513       0.010425237
## cg00000236       0.005817729      0.0266551291      -0.003336006
## cg00000289      -0.093155073      0.0002407292      -0.015128947
## cg00000292      -0.021662299     -0.0571219636       0.003836053
##            6969568082_R06C02 6969568084_R01C01 6969568084_R02C01
## cg00000029        0.06536794       0.014867075      -0.015385354
## cg00000108        0.01491589      -0.035225870       0.052310819
## cg00000165        0.02091153      -0.066411347      -0.033522297
## cg00000236       -0.02269027      -0.018889881       0.005225076
## cg00000289       -0.03289923      -0.006065534      -0.029160399
## cg00000292       -0.02984661      -0.078588444       0.032987494
##            6969568084_R03C01 6969568084_R04C01 6969568084_R06C01
## cg00000029       -0.03175264     -0.0186393886       -0.02557366
## cg00000108       -0.01804484     -0.0142591877        0.05032492
## cg00000165        0.03652135     -0.0198166357        0.03809651
## cg00000236       -0.01757326     -0.0007466135        0.02968469
## cg00000289       -0.02384056      0.0480044010       -0.05750709
## cg00000292        0.04092070     -0.0017707917       -0.01539359
##            6969568084_R02C02 6969568084_R03C02 6969568084_R04C02
## cg00000029       -0.04580745       0.017488235       0.043277312
## cg00000108        0.00355036       0.030037156       0.032521729
## cg00000165       -0.03993325       0.040669196      -0.027335392
## cg00000236        0.03776951       0.006519384      -0.029476818
## cg00000289        0.03408815      -0.034774486      -0.002714579
## cg00000292        0.03418772      -0.027550745      -0.021882609
##            6969568084_R05C02 6969568087_R01C01 6969568087_R02C01
## cg00000029      -0.012911043      -0.076352793      0.0976732203
## cg00000108       0.008398380      -0.016837447     -0.0075271454
## cg00000165      -0.052884080       0.081891659     -0.0445932933
## cg00000236      -0.023696326       0.012234783     -0.0085748179
## cg00000289      -0.003893991       0.004146153     -0.0008245228
## cg00000292       0.014238007      -0.049362878      0.0519326514
##            6969568087_R03C01 6969568087_R04C01 6969568087_R05C01
## cg00000029      -0.071206253     -0.0581830929      -0.014934962
## cg00000108      -0.003786984      0.0002280606       0.027278428
## cg00000165      -0.029928270     -0.0387045863       0.048005614
## cg00000236      -0.026730114      0.0775665728      -0.009611475
## cg00000289      -0.085722261      0.1110588333      -0.040586141
## cg00000292       0.072166383     -0.0631952747       0.032157232
##            6969568087_R06C01 6969568087_R01C02 6969568087_R02C02
## cg00000029       0.029070982     -0.0965677408     -2.151295e-02
## cg00000108       0.001067815      0.0235823904      3.967225e-04
## cg00000165      -0.011926182     -0.0748599469     -2.355422e-02
## cg00000236       0.040139659     -0.0590711660      2.579847e-02
## cg00000289      -0.059256434      0.0006260318      1.102858e-05
## cg00000292      -0.024225776     -0.0220102072     -6.225012e-02
##            6969568087_R03C02 6969568087_R05C02 6969568087_R06C02
## cg00000029        0.06913142      -0.033029347       -0.01351462
## cg00000108        0.00236962      -0.027591853        0.01966038
## cg00000165       -0.01955298      -0.005786944        0.01275705
## cg00000236        0.01092153       0.026427717       -0.01159112
## cg00000289        0.07646156       0.035479799       -0.03041406
## cg00000292        0.04825688       0.018567950       -0.03337892
##            6969568118_R01C01 6969568118_R02C01 6969568118_R03C01
## cg00000029     -3.626417e-02      -0.019291382       -0.03650825
## cg00000108     -2.726627e-02       0.041919622       -0.02106505
## cg00000165     -8.517114e-05      -0.072462588       -0.04358928
## cg00000236      9.231021e-03       0.018682463       -0.01021945
## cg00000289     -7.502991e-02       0.083885233       -0.01735682
## cg00000292      3.595174e-02      -0.004425351       -0.01503575
##            6969568118_R04C01 6969568118_R01C02 6969568118_R02C02
## cg00000029      -0.065625228       -0.04612299      -0.001704540
## cg00000108       0.054277629       -0.02302402       0.005835409
## cg00000165      -0.002916001        0.19292587       0.004311456
## cg00000236       0.012159579       -0.07167430       0.024584939
## cg00000289      -0.036154927        0.07657447      -0.012523736
## cg00000292      -0.071463632        0.07620954       0.045015985
##            6969568118_R03C02 6969568118_R04C02 6969568118_R06C02
## cg00000029       0.028210884        0.08134745       0.052220283
## cg00000108       0.021473194       -0.01042452      -0.013699993
## cg00000165      -0.011647032       -0.09020888      -0.025683790
## cg00000236       0.001106956        0.02326545       0.009295842
## cg00000289       0.018712376       -0.01626181       0.015624044
## cg00000292       0.012142197       -0.01681836      -0.081120267
##            6929726046_R02C01 6929726046_R05C01 6929726046_R06C01
## cg00000029       -0.01591791       0.002103127     -0.0055279128
## cg00000108        0.01004305      -0.012951007      0.0006939263
## cg00000165       -0.02163900       0.003972774      0.0261136526
## cg00000236       -0.01562664      -0.009076841      0.0318942686
## cg00000289        0.01619350       0.076466788     -0.0552103246
## cg00000292       -0.01096318       0.012376242      0.0452633377
##            6929726046_R01C02 6929726046_R03C02 6929726046_R04C02
## cg00000029        0.03516505       0.005037713       -0.10099852
## cg00000108       -0.03065375       0.004938669        0.03461101
## cg00000165        0.05103111      -0.052301000       -0.07399736
## cg00000236       -0.05028982      -0.004433956        0.02643122
## cg00000289       -0.03910164       0.039856294       -0.06814916
## cg00000292       -0.01763141      -0.082520022        0.03349092
##            6929726046_R06C02 6929718123_R01C01 6929718123_R04C01
## cg00000029      0.0273411123        0.01687089       0.021963524
## cg00000108     -0.0008893201        0.04402384       0.043894854
## cg00000165      0.0328506680       -0.03092592      -0.023547328
## cg00000236     -0.0181618510        0.04968344       0.022112465
## cg00000289     -0.0181351991        0.04414899       0.025999456
## cg00000292      0.0358201606       -0.01657393      -0.002487906
##            6929718123_R06C01 6929718123_R01C02 6929718123_R02C02
## cg00000029      0.0042568183       0.013656457        0.09505740
## cg00000108      0.0092669310       0.058455134       -0.21482880
## cg00000165     -0.0282756699      -0.067257273        0.24067795
## cg00000236     -0.0006359615       0.025255894       -0.06839388
## cg00000289     -0.0089900365      -0.025306447       -0.03882611
## cg00000292     -0.0346793240      -0.008113902       -0.01060699
##            6929718123_R03C02 6929718123_R04C02 6929718123_R05C02
## cg00000029       0.002632712       -0.08400590       -0.02954871
## cg00000108       0.012915441        0.07199200        0.07745217
## cg00000165      -0.051116176       -0.06019922       -0.07786390
## cg00000236       0.012732849        0.01607945        0.03023042
## cg00000289      -0.001791076        0.03777649       -0.05103605
## cg00000292      -0.005253108       -0.01999589        0.08234911
##            6929718123_R06C02 6929718136_R01C01 6929718136_R02C01
## cg00000029      -0.026874921       -0.01442799       -0.07294990
## cg00000108       0.030140793        0.01575447        0.02566432
## cg00000165       0.001600563       -0.04195798       -0.11873290
## cg00000236      -0.025849570        0.03148150        0.01824183
## cg00000289      -0.010751600        0.04891289        0.05347639
## cg00000292       0.039834577       -0.03062392        0.02903340
##            6929718136_R03C01 6929718136_R04C01 6929718136_R05C01
## cg00000029       0.005439189       0.141297935        0.01111295
## cg00000108       0.022096699      -0.243910307        0.08015618
## cg00000165       0.032267157       0.213659174       -0.08959802
## cg00000236       0.024651999      -0.073393656        0.01524148
## cg00000289       0.021648283       0.014860752        0.04297210
## cg00000292      -0.023250427      -0.003958495        0.03586106
##            6929718136_R02C02 6929718136_R04C02 6929718136_R05C02
## cg00000029       0.001027792      -0.007497916       0.024817753
## cg00000108       0.035506098       0.009306162       0.035675972
## cg00000165      -0.022615774      -0.041475259      -0.032453160
## cg00000236      -0.007489937       0.015139761      -0.015350423
## cg00000289      -0.043940148      -0.069967685       0.003492241
## cg00000292       0.006019994       0.046719699      -0.064713097
##            6929718136_R06C02 6929718138_R01C01 6929718138_R02C01
## cg00000029       0.033192978      -0.115378757       -0.04255914
## cg00000108       0.040820796       0.031993353        0.01574513
## cg00000165      -0.002459254      -0.004986770       -0.03846147
## cg00000236       0.004325150       0.006122217        0.01028543
## cg00000289       0.005832531      -0.032695282       -0.01846473
## cg00000292       0.059170675      -0.001499825       -0.06943810
##            6929718138_R04C01 6929718138_R06C01 6929718138_R01C02
## cg00000029       -0.02198186       0.071523301      -0.032656220
## cg00000108        0.02105635      -0.002482208      -0.066320129
## cg00000165        0.01224372       0.035426833      -0.041524312
## cg00000236        0.03430762      -0.064967517       0.006911663
## cg00000289        0.02380276       0.033731596       0.034653264
## cg00000292       -0.04113287      -0.011480593      -0.064632732
##            6929718138_R03C02 6929718138_R04C02 6929718138_R05C02
## cg00000029       0.089498607       -0.01780169      0.0627722262
## cg00000108      -0.138698226        0.01460660      0.0008109713
## cg00000165      -0.034957240       -0.01402490     -0.0340145041
## cg00000236       0.006355496       -0.01438964      0.0100731222
## cg00000289       0.103608694        0.03182887     -0.0261082072
## cg00000292       0.022239217       -0.03301156     -0.0060849846
##            6929718138_R06C02 6042316054_R02C01 6042316054_R02C02
## cg00000029      -0.029696028      -0.061705452     -0.0001869817
## cg00000108      -0.016339697       0.053734277      0.0186401175
## cg00000165       0.040929349      -0.017033996      0.1269579139
## cg00000236       0.027411691      -0.004797062      0.0034693734
## cg00000289       0.009596633       0.032275450     -0.0060829586
## cg00000292      -0.083626142       0.049646742      0.0522146212
##            6042316054_R03C01 6042316054_R03C02 6042316054_R04C01
## cg00000029      -0.065434931       0.023951053       0.005429512
## cg00000108       0.055649672       0.060415791       0.021129209
## cg00000165      -0.025486715      -0.006042480      -0.041790077
## cg00000236       0.006933794       0.012408974      -0.023289125
## cg00000289      -0.067960651      -0.039796423      -0.059832056
## cg00000292      -0.006406815       0.005273379       0.024414928
##            6042316054_R04C02 6042316054_R05C01 6042316054_R05C02
## cg00000029       -0.09141453        0.06066688      -0.086054049
## cg00000108        0.06376408        0.02927941       0.006148815
## cg00000165       -0.06681982        0.02957798      -0.049740360
## cg00000236        0.04653445        0.02257678      -0.008243795
## cg00000289       -0.04578799        0.04485957      -0.038887984
## cg00000292        0.04507350        0.04565290      -0.037174430
##            6042316054_R06C01 6042316063_R02C01 6042316063_R02C02
## cg00000029        0.03218854      -0.026733415      -0.014061848
## cg00000108        0.03164404       0.005888453       0.035340464
## cg00000165       -0.01479970       0.037873573      -0.047856213
## cg00000236       -0.03127086      -0.018150827       0.008906632
## cg00000289       -0.05601731      -0.044219755      -0.010264635
## cg00000292        0.09094012       0.038458968      -0.019017105
##            6042316063_R03C01 6042316063_R03C02 6042316063_R04C01
## cg00000029      -0.011703976       0.083252991        0.01482783
## cg00000108       0.003397503       0.002773884       -0.04350541
## cg00000165       0.003664608      -0.061962877       -0.05751139
## cg00000236       0.035421442       0.016812419       -0.04132439
## cg00000289       0.033157422       0.048086871        0.02701600
## cg00000292       0.002949498      -0.010607344        0.04146188
##            6042316063_R04C02 6042316063_R05C01 6042316063_R05C02
## cg00000029       0.008964868       -0.06580984      2.881161e-03
## cg00000108       0.032481569       -0.03660106      1.396992e-02
## cg00000165      -0.037454573        0.08490249     -1.246402e-02
## cg00000236       0.028135461       -0.05143377      5.359563e-02
## cg00000289       0.015536065       -0.09847670      3.094737e-02
## cg00000292      -0.033094410       -0.07144160     -3.196853e-07
##            6042316063_R06C02 6042316065_R01C02 6042316065_R02C02
## cg00000029      -0.053660251       0.029973920       -0.03035653
## cg00000108      -0.028119044       0.014905388        0.02673050
## cg00000165       0.051966650       0.098155454       -0.02042698
## cg00000236      -0.002328256       0.029205039        0.03266826
## cg00000289       0.051189583       0.087205911        0.10797053
## cg00000292      -0.030373680       0.005818307       -0.04391338
##            6042316065_R03C01 6042316065_R04C01 6042316065_R04C02
## cg00000029        0.04412606      -0.080346172        0.01342108
## cg00000108        0.03271747       0.008048537        0.04272314
## cg00000165       -0.03608750      -0.044160748       -0.04403907
## cg00000236        0.04420042      -0.013170533       -0.06717014
## cg00000289        0.03679828      -0.028923614       -0.05153409
## cg00000292       -0.07707353       0.029926124        0.02928421
##            6042316065_R05C02 6042316065_R06C02 6042316103_R02C01
## cg00000029      -0.001065635        0.03037769      -0.013475141
## cg00000108       0.020993911       -0.08853869      -0.009279845
## cg00000165       0.061023821        0.03368993      -0.055315920
## cg00000236       0.010487282        0.01141164      -0.020237191
## cg00000289      -0.033798345       -0.06701290      -0.019530557
## cg00000292      -0.045976376        0.01392841       0.015331839
##            6042316103_R03C01 6042316103_R03C02 6042316103_R04C01
## cg00000029       -0.07049675      -0.044401162       0.035914669
## cg00000108        0.03366428       0.023017967       0.004644858
## cg00000165       -0.06983302       0.026036529      -0.006741363
## cg00000236       -0.01740020       0.050295612       0.031102258
## cg00000289       -0.04767369      -0.082549853       0.088339674
## cg00000292        0.10300552       0.006568618      -0.001997972
##            6042316103_R05C01 6042316103_R06C01 6042316103_R06C02
## cg00000029      -0.012657028      -0.106441620     -0.0036168909
## cg00000108      -0.008146603      -0.005066321      0.0525197563
## cg00000165       0.041580015       0.056288355     -0.1004003924
## cg00000236       0.015994528      -0.007237126     -0.0459561650
## cg00000289      -0.013212392       0.015103861     -0.0002750954
## cg00000292       0.016178091      -0.040435090      0.0409785535
##            6042316036_R01C02 6042316036_R02C01 6042316036_R03C01
## cg00000029        0.02126428      -0.013245052       0.033876884
## cg00000108       -0.03949402       0.047636238       0.006266931
## cg00000165       -0.06537478      -0.008462758      -0.018742652
## cg00000236        0.02117234      -0.026980796      -0.030110216
## cg00000289        0.01343439      -0.052439339       0.019953454
## cg00000292       -0.05930910       0.073810556      -0.025158896
##            6042316036_R03C02 6042316036_R04C01 6042316036_R04C02
## cg00000029      -0.020907106      -0.084579466     -0.0346716754
## cg00000108       0.014643525       0.003761483     -0.0002073464
## cg00000165      -0.057185073      -0.076265868     -0.0051069691
## cg00000236      -0.015265501       0.016383883      0.0328586860
## cg00000289      -0.026832726       0.058926156     -0.0350564297
## cg00000292       0.003867463      -0.049061846     -0.0361321789
##            6042316036_R05C01 6042316036_R05C02 6042316036_R06C01
## cg00000029      -0.002485721       0.011684634       0.061304576
## cg00000108       0.066348105      -0.037259720       0.009468132
## cg00000165       0.067903419      -0.027008085       0.012222600
## cg00000236       0.040922971      -0.067750168      -0.020200215
## cg00000289       0.090683565       0.001836601      -0.006411179
## cg00000292      -0.077533496       0.114115212       0.005516430
##            6042316036_R06C02 6042316050_R01C02 6042316050_R02C02
## cg00000029       -0.02537055        0.03783734        0.16374437
## cg00000108        0.02797254        0.02747394       -0.03691701
## cg00000165        0.08276470       -0.06277006        0.00740945
## cg00000236        0.02394948        0.01653797        0.03695959
## cg00000289       -0.11923414       -0.04204276        0.02653764
## cg00000292       -0.04236289       -0.04366677        0.05120785
##            6042316050_R03C01 6042316050_R04C01 6042316050_R04C02
## cg00000029      -0.021378509      -0.017463451      -0.053544663
## cg00000108      -0.032120380      -0.002945187       0.010483758
## cg00000165      -0.068841975      -0.004804743       0.002061097
## cg00000236      -0.021349655       0.023279037      -0.001657772
## cg00000289      -0.005037868       0.070795491       0.075883694
## cg00000292       0.043026680      -0.038960488       0.028425544
##            6042316050_R05C01 6042316050_R05C02 6042316050_R06C01
## cg00000029       -0.02341107       -0.04916834       0.025108072
## cg00000108       -0.03890441        0.06058031      -0.006452124
## cg00000165        0.05029847        0.03962773       0.143665444
## cg00000236        0.01589811       -0.06616917      -0.037233191
## cg00000289       -0.01265088       -0.06485632      -0.031044573
## cg00000292        0.02143003        0.03455465       0.062639800
##            6042316050_R06C02 6042316053_R01C01 6042316053_R02C01
## cg00000029       0.004381966     -0.0442013313      -0.018033807
## cg00000108       0.016965989      0.0057743048       0.013164925
## cg00000165      -0.053631160     -0.0189638718       0.031690676
## cg00000236      -0.006833288      0.0003548099       0.011780609
## cg00000289       0.011784696      0.0035884921      -0.089940787
## cg00000292      -0.008834957      0.0082897312      -0.007364631
##            6042316053_R02C02 6042316053_R03C01 6042316053_R03C02
## cg00000029        0.02477670       0.009224736      -0.015057810
## cg00000108        0.07407507       0.065638445       0.041349334
## cg00000165       -0.01279638       0.016224954       0.011175614
## cg00000236        0.03763722       0.013779293       0.022651740
## cg00000289        0.04835605      -0.016458924      -0.073303556
## cg00000292       -0.03674999      -0.009947984      -0.002722211
##            6042316053_R04C01 6042316053_R04C02 6042316053_R06C01
## cg00000029      0.0081450461       -0.10812771      -0.015991923
## cg00000108      0.0566572406        0.03626709       0.020977529
## cg00000165     -0.0150681395        0.01556935       0.040119417
## cg00000236      0.0005674277        0.01194674       0.009606092
## cg00000289      0.0102389044       -0.02051692       0.055812989
## cg00000292      0.0288592745        0.08276688       0.066923716
##            6042316053_R06C02 6042316061_R01C01 6042316061_R02C01
## cg00000029      -0.002182810       0.026443207       0.001919673
## cg00000108      -0.027247740       0.029241302       0.036526760
## cg00000165      -0.038410911      -0.020602785      -0.015553771
## cg00000236      -0.024951964       0.021936429       0.042905461
## cg00000289      -0.019705908      -0.009599511      -0.032454289
## cg00000292      -0.004027709       0.016883055      -0.024656188
##            6042316061_R03C01 6042316061_R03C02 6042316061_R04C01
## cg00000029      -0.091972692       0.030945734      -0.023560082
## cg00000108       0.052439834       0.004542238      -0.004613775
## cg00000165      -0.045030428      -0.061675029      -0.036387606
## cg00000236       0.004937753       0.005306796       0.003171868
## cg00000289       0.019150147       0.036073264       0.052990022
## cg00000292      -0.074267156      -0.042238001      -0.086113111
##            6042316061_R04C02 6042316061_R05C01 6042316061_R05C02
## cg00000029      -0.007215331       0.056502872      -0.014509084
## cg00000108       0.036383453       0.010579149       0.050343296
## cg00000165      -0.084008717      -0.010016390      -0.006536696
## cg00000236       0.030346879      -0.004147245       0.024892604
## cg00000289       0.040474766       0.016285839       0.110716819
## cg00000292       0.014053171       0.057643620      -0.028630126
##            7796806022_R01C01 7796806022_R03C01 7796806022_R04C01
## cg00000029        0.18744771       -0.13614609       -0.08483534
## cg00000108       -0.38166418        0.03168692        0.08097925
## cg00000165        0.11263636       -0.04331662       -0.06444064
## cg00000236       -0.14332468        0.01101559        0.01904053
## cg00000289        0.02467328       -0.11000109       -0.06243557
## cg00000292       -0.03782123        0.03027121        0.06553806
##            7796806022_R04C02 7796806022_R05C02 7796806022_R06C01
## cg00000029        0.07745824      -0.097028135       0.067768094
## cg00000108       -0.10378706       0.086688480      -0.122690345
## cg00000165        0.13494711      -0.051787152       0.151072637
## cg00000236       -0.03057635       0.030481143      -0.016573462
## cg00000289        0.05933918      -0.007734859       0.029863629
## cg00000292       -0.08626546       0.013828718      -0.005375975
##            7786923046_R03C02 7786923046_R04C01 7786923046_R06C01
## cg00000029      -0.060143455       -0.10176546       -0.09058807
## cg00000108       0.085600483        0.05493026        0.05846638
## cg00000165      -0.081964408       -0.09505939       -0.04382538
## cg00000236       0.047087320       -0.01858320        0.02155304
## cg00000289       0.002139694       -0.05562929       -0.06901637
## cg00000292       0.039158647       -0.01048033       -0.04703010
##            7796806038_R03C01 7796806038_R03C02 7796806038_R05C02
## cg00000029       -0.07088822       -0.05941514       0.063566681
## cg00000108        0.07965419        0.07422844      -0.154996589
## cg00000165       -0.07335943       -0.07264574       0.217403520
## cg00000236        0.03135306        0.03246551      -0.033480293
## cg00000289       -0.08967890       -0.03286340       0.057863055
## cg00000292       -0.02503613        0.02392338      -0.009928595
##            7786923063_R01C01 7786923063_R01C02 7786923063_R03C01
## cg00000029       0.084279306       -0.07848811       0.107739176
## cg00000108      -0.132660994        0.08245633      -0.073735507
## cg00000165      -0.068599710       -0.06988094      -0.048222408
## cg00000236      -0.066949298        0.01859132      -0.038740654
## cg00000289      -0.002452441       -0.06269494       0.003887945
## cg00000292      -0.010147284        0.04461668      -0.068722671
##            7786923063_R04C01 7786923063_R06C01 7786923107_R01C02
## cg00000029        0.07926370      0.0657423087       -0.02357945
## cg00000108       -0.06370532     -0.1041888492        0.05553763
## cg00000165        0.09787381      0.0770659791       -0.04618282
## cg00000236       -0.04803286     -0.0421719014        0.02898696
## cg00000289        0.06593931     -0.0189475049        0.08621801
## cg00000292       -0.07805166     -0.0007801806        0.03733556
##            7786923107_R04C01 7786923107_R05C01 7786923107_R06C01
## cg00000029       -0.01536866     -0.0747450297       0.022861874
## cg00000108        0.06281888     -0.0022541661      -0.007413402
## cg00000165       -0.01450054     -0.0968441826      -0.006942708
## cg00000236        0.03434081      0.0007849476      -0.036948581
## cg00000289        0.01122122     -0.0036192708       0.043570376
## cg00000292        0.02374983     -0.0257661125       0.032549564
##            7796806016_R03C02 7796806016_R06C01 7796806002_R01C01
## cg00000029      -0.013299828        0.04765942      -0.004672541
## cg00000108       0.052989055        0.05171763      -0.040100028
## cg00000165      -0.004064604       -0.02266147       0.087604898
## cg00000236      -0.023926754        0.01096521       0.037578754
## cg00000289      -0.057234764       -0.03255953      -0.032714193
## cg00000292      -0.026745488        0.05877262       0.049079783
##            7796806002_R02C02 7796806002_R04C02 7796806002_R05C01
## cg00000029        0.05846698      -0.062534989       0.016182658
## cg00000108        0.01527527      -0.026855247      -0.027076343
## cg00000165        0.02622914      -0.048124486       0.006444171
## cg00000236        0.02819775       0.029794261       0.032718458
## cg00000289        0.02438061       0.013888863       0.030648451
## cg00000292        0.04833424      -0.009413213      -0.029286158
##            7796806002_R05C02 7796806002_R06C01 7796806002_R06C02
## cg00000029      -0.041845890      -0.079701404        0.05369435
## cg00000108       0.049554710       0.003021531       -0.03638587
## cg00000165      -0.041111950       0.034116503        0.02369035
## cg00000236      -0.011017481      -0.038769921       -0.02381843
## cg00000289       0.006572768      -0.034361526        0.01054564
## cg00000292       0.020744513      -0.040500970       -0.02773572
##            7796806029_R02C01 7796806029_R04C01 7796806029_R06C01
## cg00000029       0.014405617       -0.07104980       -0.06547997
## cg00000108       0.015335758        0.08693854        0.09538066
## cg00000165      -0.008890422       -0.08800809       -0.06826917
## cg00000236       0.007977211        0.02438712        0.02631607
## cg00000289       0.042987338        0.06343364       -0.01894643
## cg00000292      -0.056369859        0.03044089        0.02728267
##            6042316024_R01C01 6042316024_R01C02 6042316024_R02C01
## cg00000029        0.13174223       0.030880946        0.05973409
## cg00000108       -0.05392788      -0.138258645       -0.08667895
## cg00000165        0.02690851       0.162255717        0.11779361
## cg00000236       -0.09152106       0.001096319       -0.02230762
## cg00000289       -0.02640357       0.036194038       -0.02362847
## cg00000292        0.06049516      -0.056529083       -0.08759303
##            6042316024_R02C02 6042316024_R03C01 6042316024_R04C01
## cg00000029        0.05571753       -0.03223225       0.046106972
## cg00000108       -0.07584559        0.04305644      -0.041379037
## cg00000165        0.03022745        0.04310770      -0.195648197
## cg00000236       -0.01686378        0.01764229      -0.009773799
## cg00000289        0.03170064        0.02863046       0.003062593
## cg00000292       -0.03012374       -0.02658439      -0.050742922
##            6042316024_R04C02 6042316024_R05C01 6042316024_R05C02
## cg00000029      0.0008491659      -0.083521122        0.04398197
## cg00000108      0.1557431871      -0.002936145        0.10444680
## cg00000165      0.1125915464       0.054762054       -0.02333878
## cg00000236      0.0421729116      -0.054539766       -0.01943114
## cg00000289     -0.0929412939      -0.033693412        0.01263723
## cg00000292      0.0673631501       0.022897251        0.04811950
##            6042316024_R06C01 6042316024_R06C02 6042316031_R01C01
## cg00000029        0.03489566      -0.082134055       0.049799154
## cg00000108       -0.10354094      -0.076037734      -0.052550178
## cg00000165       -0.02268019      -0.005144283       0.081536027
## cg00000236        0.05589088      -0.037966831      -0.013586504
## cg00000289       -0.01904234       0.036738753       0.033426228
## cg00000292       -0.09040897       0.044413569      -0.001603192
##            6042316031_R01C02 6042316031_R02C01 6042316031_R02C02
## cg00000029       0.106547035      -0.027418918        0.13064585
## cg00000108      -0.178066386      -0.107855517        0.02850626
## cg00000165       0.109996805      -0.152215333       -0.03261629
## cg00000236       0.001603835       0.001986234       -0.04464024
## cg00000289       0.013953287      -0.004660693       -0.11814848
## cg00000292      -0.062754285       0.038853212       -0.11838602
##            6042316031_R03C01 6042316031_R03C02 6042316031_R04C01
## cg00000029       -0.09496067       0.029580111      -0.001522921
## cg00000108       -0.07552463       0.062614694       0.007773995
## cg00000165        0.14301151       0.105323048      -0.027164187
## cg00000236       -0.04110090       0.003498648      -0.051665129
## cg00000289       -0.03694021       0.060106818       0.014357925
## cg00000292        0.02609996      -0.027747655      -0.013781786
##            6042316031_R04C02 6042316031_R05C01 6042316031_R05C02
## cg00000029        0.02753569      -0.052034322        0.05993331
## cg00000108        0.04105850       0.129965985       -0.04774067
## cg00000165        0.08588239      -0.005640624        0.05158288
## cg00000236        0.01246445       0.059571964        0.04131724
## cg00000289        0.05075608       0.005325816        0.02088416
## cg00000292        0.01134005       0.011217845        0.02350899
##            6042316031_R06C01 6042316031_R06C02 6042316042_R01C01
## cg00000029        0.03927887       0.004841982       -0.01274311
## cg00000108       -0.12839488      -0.061122016       -0.04654393
## cg00000165       -0.05412422       0.145415596        0.14728701
## cg00000236       -0.02135674      -0.014840528       -0.04444729
## cg00000289        0.01331360      -0.078159870       -0.01212035
## cg00000292        0.03879734       0.050700025       -0.01257582
##            6042316042_R01C02 6042316042_R02C01 6042316042_R02C02
## cg00000029       -0.01778515       -0.02422059        0.06338158
## cg00000108        0.12501705        0.02867560       -0.10878414
## cg00000165       -0.09989340       -0.13127725        0.07230115
## cg00000236        0.03922316       -0.01347231        0.01649649
## cg00000289        0.04588707        0.02608303        0.05493241
## cg00000292       -0.01617864        0.01560616       -0.03758923
##            6042316042_R03C01 6042316042_R03C02 6042316042_R04C01
## cg00000029        0.03927605       -0.02102154       0.044857172
## cg00000108       -0.05479258        0.02147511       0.095166660
## cg00000165        0.02236321        0.02899115       0.002944548
## cg00000236        0.02763159       -0.02101441       0.040190591
## cg00000289       -0.05633123       -0.03485580       0.075774801
## cg00000292       -0.06957759        0.09836350      -0.024997334
##            6042316042_R04C02 6042316042_R05C01 6042316042_R05C02
## cg00000029       0.010070616       0.049033052      -0.056940775
## cg00000108      -0.090563865      -0.001807489      -0.111558599
## cg00000165      -0.046044887       0.178211384      -0.025412681
## cg00000236      -0.066033258       0.036748788      -0.041335935
## cg00000289       0.004155097       0.013230830      -0.007780935
## cg00000292      -0.013239367      -0.038793752      -0.031190071
##            6042316042_R06C01 6042316042_R06C02 6042316047_R01C01
## cg00000029       0.107922538      0.0536542979       0.004595428
## cg00000108       0.006239101      0.0575917843      -0.100994354
## cg00000165      -0.009414817     -0.1543405387       0.028213571
## cg00000236       0.013744025     -0.0007159008      -0.043095411
## cg00000289      -0.030798220     -0.0172056091       0.050391213
## cg00000292       0.080916866      0.0323787743       0.039882880
##            6042316047_R02C01 6042316047_R02C02 6042316047_R04C01
## cg00000029       -0.07830996      -0.015289871        0.03158604
## cg00000108        0.13844125      -0.119543256       -0.03224267
## cg00000165        0.03847529       0.023541052       -0.08344816
## cg00000236        0.01149913       0.005598539        0.02260636
## cg00000289       -0.02140977       0.031188566       -0.03495543
## cg00000292       -0.04816758      -0.027902013        0.06511913
##            6042316047_R04C02 6042316047_R05C01 6042316047_R05C02
## cg00000029       -0.04858612      -0.040338851        0.12307148
## cg00000108       -0.14139133       0.069490897        0.08309669
## cg00000165       -0.05337181       0.144818468        0.05724326
## cg00000236       -0.07554431      -0.009595584        0.03515006
## cg00000289        0.02661939      -0.029746908        0.04237534
## cg00000292       -0.02684201      -0.056894428        0.01305768
##            6042316047_R06C01 6055432012_R01C01 6055432012_R01C02
## cg00000029      0.0101262055      -0.055566536      -0.082003980
## cg00000108     -0.0469229043       0.040538509       0.008241057
## cg00000165      0.0004338743       0.017731467      -0.048407533
## cg00000236      0.0012212965      -0.008191588       0.065434810
## cg00000289      0.0097555779      -0.077386883      -0.052791967
## cg00000292     -0.0040051998      -0.035554608       0.004268112
##            6055432012_R02C02 6055432012_R03C01 6055432012_R03C02
## cg00000029      -0.076045214      0.0657775190       0.015904966
## cg00000108       0.025155484     -0.0020717468       0.040627191
## cg00000165       0.077817026     -0.0676854567      -0.035740755
## cg00000236       0.007337454      0.0392799059       0.003314289
## cg00000289       0.019512915     -0.0007527553       0.058381968
## cg00000292      -0.071506083      0.0319094908       0.078358811
##            6055432012_R04C01 6055432012_R04C02 6055432012_R05C01
## cg00000029        0.01686954       0.041820000       0.044922520
## cg00000108       -0.06732544       0.030190053       0.005014317
## cg00000165        0.03297988      -0.013567239       0.032871664
## cg00000236       -0.03672225      -0.061922108      -0.008523356
## cg00000289       -0.01996701      -0.007241382       0.053222520
## cg00000292        0.03198334      -0.054303070       0.042056721
##            6055432012_R05C02 6055432012_R06C01 6055432012_R06C02
## cg00000029      -0.032099534       -0.01578572       0.021109477
## cg00000108       0.009320197        0.02010022       0.006780074
## cg00000165      -0.040562795        0.06280885       0.006855183
## cg00000236       0.005009934       -0.01675915      -0.019738872
## cg00000289       0.010607533       -0.05606263       0.044084559
## cg00000292       0.011624588        0.07590858       0.039525940
##            6055432029_R01C01 6055432029_R01C02 6055432029_R02C01
## cg00000029      -0.066656294        0.05416508       0.049160301
## cg00000108      -0.058276323        0.04037802       0.008466849
## cg00000165       0.020188289        0.06478586      -0.076694466
## cg00000236      -0.082346911        0.03361143       0.038717233
## cg00000289      -0.002004906       -0.03776719      -0.002143195
## cg00000292       0.035165939        0.01793000       0.030209913
##            6055432029_R02C02 6055432029_R03C01 6055432029_R03C02
## cg00000029        0.04850399        0.02020419       -0.03796601
## cg00000108       -0.01702889        0.02406444        0.04403169
## cg00000165        0.09922660       -0.07906511       -0.05557837
## cg00000236        0.01668699       -0.01135535        0.03377880
## cg00000289        0.07913421       -0.05612406        0.05129789
## cg00000292       -0.02030854       -0.05349907        0.05573457
##            6055432029_R04C01 6055432029_R04C02 6055432029_R05C01
## cg00000029      -0.027299656      -0.076466626        0.02389190
## cg00000108      -0.003291552       0.010046185        0.01682947
## cg00000165       0.013346433      -0.015235784       -0.04991497
## cg00000236      -0.011601587       0.001799793        0.03840538
## cg00000289      -0.066216207       0.029131883       -0.01200215
## cg00000292      -0.016141834      -0.064776170       -0.04831045
##            6055432029_R05C02 6055432029_R06C01 6055432029_R06C02
## cg00000029      -0.023932615       0.092139797       -0.08774747
## cg00000108       0.007374485      -0.008777125       -0.01871640
## cg00000165      -0.019702542       0.062867813        0.09337595
## cg00000236       0.033838362      -0.049933924       -0.01941414
## cg00000289      -0.011137796       0.015396648        0.03154098
## cg00000292       0.038163146       0.001019793        0.02189382
##            6055432060_R01C01 6055432060_R01C02 6055432060_R02C01
## cg00000029       -0.05224327       0.007782106      -0.072899841
## cg00000108        0.04013449      -0.001826297       0.009821846
## cg00000165       -0.05694667       0.018628765      -0.032978107
## cg00000236       -0.01890670       0.017429599       0.036427264
## cg00000289        0.02187466      -0.097236640       0.036245320
## cg00000292       -0.07434584      -0.044361358       0.019920976
##            6055432060_R02C02 6055432060_R03C01 6055432060_R03C02
## cg00000029      -0.007427281       0.053570401       0.077528300
## cg00000108      -0.019085502       0.023493583       0.050786654
## cg00000165       0.032448400      -0.048002360      -0.070593794
## cg00000236      -0.015482981      -0.004998538       0.039857161
## cg00000289      -0.078194129       0.038288285       0.026522374
## cg00000292       0.011245118       0.048076199       0.005216072
##            6055432060_R04C01 6055432060_R04C02 6055432060_R05C01
## cg00000029       -0.04001888        0.02785771       -0.01196823
## cg00000108        0.04916380        0.02704269        0.01712823
## cg00000165       -0.01324103       -0.04289048        0.03148314
## cg00000236        0.01651502        0.03620444       -0.03933159
## cg00000289        0.06748040       -0.02413820       -0.06254339
## cg00000292        0.03457886       -0.04552373        0.03103822
##            6055432060_R05C02 6055432060_R06C01 6055432060_R06C02
## cg00000029       -0.05808184       0.009552082       0.112137105
## cg00000108        0.04824809      -0.005400281       0.001062355
## cg00000165       -0.03830549       0.014049179       0.121301409
## cg00000236        0.02994838      -0.010573605      -0.015797951
## cg00000289        0.02159809      -0.021020987      -0.018748522
## cg00000292        0.04339066      -0.038622901       0.060178476
##            6055432066_R01C01 6055432066_R01C02 6055432066_R02C01
## cg00000029       -0.08161919       -0.01183839        0.03379657
## cg00000108        0.01437174        0.02523331        0.02671996
## cg00000165       -0.07439629       -0.04217250       -0.03138775
## cg00000236        0.02883979       -0.02762452       -0.02927888
## cg00000289        0.03180259       -0.04119555        0.06789632
## cg00000292       -0.02875934        0.02524944       -0.03004168
##            6055432066_R02C02 6055432066_R03C02 6055432066_R04C01
## cg00000029     -0.1025742145      -0.003473520       0.044692550
## cg00000108      0.0024921841       0.035724852       0.026025098
## cg00000165      0.0317206539      -0.034497560       0.007379219
## cg00000236      0.0263725565       0.005325006       0.019032317
## cg00000289     -0.0309039199      -0.061469769      -0.007472100
## cg00000292     -0.0003636307       0.039213844       0.006562815
##            6055432066_R04C02 6055432066_R05C01 6055432066_R05C02
## cg00000029       -0.06977938       -0.05529461       -0.05190706
## cg00000108        0.10660640        0.11249734        0.09293135
## cg00000165       -0.01732375       -0.05949368       -0.08433079
## cg00000236        0.02566402        0.02275901        0.01582083
## cg00000289       -0.01279434       -0.05310926       -0.06900631
## cg00000292        0.02756338        0.06542107        0.03586143
##            6057825115_R01C01 6057825115_R03C01 6057825115_R05C01
## cg00000029        0.07649900        0.06194635        0.08980255
## cg00000108        0.00298607       -0.09518197       -0.10750285
## cg00000165        0.01593133        0.07262439        0.05686954
## cg00000236       -0.01196789       -0.07352060        0.02255461
## cg00000289       -0.06351472       -0.01242363        0.09093426
## cg00000292        0.01768188       -0.04347184       -0.04098556
##            6057825115_R01C02 6057825115_R03C02 6057825115_R05C02
## cg00000029        0.06433251       0.029932295      -0.004852332
## cg00000108       -0.14072408      -0.021183441      -0.020598517
## cg00000165        0.02083885      -0.086835204       0.118910428
## cg00000236       -0.03379950      -0.020805242       0.006315054
## cg00000289        0.04063348       0.065905165      -0.052967911
## cg00000292       -0.03347542       0.002124535      -0.010741469
##            6057825128_R01C01 6057825128_R03C01 6057825128_R05C01
## cg00000029       0.029861225        0.09498250      -0.029579462
## cg00000108      -0.057938486        0.01916358      -0.210114705
## cg00000165       0.082947297       -0.06737177      -0.049194524
## cg00000236      -0.072935922        0.04603360       0.001307227
## cg00000289      -0.061566472        0.04940591      -0.024346641
## cg00000292       0.004583374        0.08404911      -0.075541572
##            6057825115_R02C01 6057825115_R04C01 6057825115_R06C01
## cg00000029      -0.017578926       0.024390063      -0.043287578
## cg00000108       0.027813443       0.031290920       0.065432956
## cg00000165       0.188541633      -0.093163600       0.041552927
## cg00000236      -0.007996397       0.003239786      -0.007310528
## cg00000289      -0.039179159       0.055428550       0.045387341
## cg00000292      -0.006869605      -0.014813748      -0.007197507
##            6057825115_R02C02 6057825115_R04C02 6057825115_R06C02
## cg00000029      -0.030259298        0.11091701      -0.023363311
## cg00000108      -0.041914099       -0.13244187       0.060768727
## cg00000165       0.111906134        0.09127767      -0.108038614
## cg00000236       0.003835905       -0.01255700      -0.004495292
## cg00000289      -0.049677577        0.03048242       0.036905709
## cg00000292       0.007178196       -0.02787043       0.007433252
##            6057825128_R01C02 6057825128_R03C02 6057825128_R05C02
## cg00000029       0.027797866       0.024190104     -0.0095380391
## cg00000108      -0.013161955      -0.030108197     -0.0233524676
## cg00000165      -0.073982150      -0.033575623      0.1803399454
## cg00000236      -0.048365697       0.027395746     -0.0076681760
## cg00000289      -0.043858704       0.077474612     -0.0436154185
## cg00000292       0.005737034      -0.007986944     -0.0008568996
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
##   [1] "Component \"6057825008_R02C02\": Mean relative difference: 0.004872702"   
##   [2] "Component \"6057825008_R03C01\": Mean relative difference: 0.005437368"   
##   [3] "Component \"6057825008_R04C01\": Mean relative difference: 5.633047e-05"  
##   [4] "Component \"6057825008_R04C02\": Mean relative difference: 0.0006643696"  
##   [5] "Component \"6057825008_R05C01\": Mean relative difference: 0.001627737"   
##   [6] "Component \"6057825008_R05C02\": Mean relative difference: 0.007450999"   
##   [7] "Component \"6057825008_R06C01\": Mean relative difference: 0.005555386"   
##   [8] "Component \"6057825008_R06C02\": Mean relative difference: 0.003551623"   
##   [9] "Component \"6057825014_R01C02\": Mean relative difference: 0.003646135"   
##  [10] "Component \"6057825014_R02C02\": Mean relative difference: 0.01126144"    
##  [11] "Component \"6057825014_R03C01\": Mean relative difference: 0.002919782"   
##  [12] "Component \"6057825014_R03C02\": Mean relative difference: 0.008064902"   
##  [13] "Component \"6057825014_R04C01\": Mean relative difference: 0.002844363"   
##  [14] "Component \"6057825014_R04C02\": Mean relative difference: 0.002091153"   
##  [15] "Component \"6057825014_R06C01\": Mean relative difference: 0.003019723"   
##  [16] "Component \"6057825017_R01C01\": Mean relative difference: 0.00233843"    
##  [17] "Component \"6057825017_R01C02\": Mean relative difference: 0.004714361"   
##  [18] "Component \"6057825017_R02C02\": Mean relative difference: 0.0005545101"  
##  [19] "Component \"6057825017_R03C01\": Mean relative difference: 0.0002074783"  
##  [20] "Component \"6057825017_R03C02\": Mean relative difference: 0.003284349"   
##  [21] "Component \"6057825017_R04C01\": Mean relative difference: 0.003828654"   
##  [22] "Component \"6057825017_R04C02\": Mean relative difference: 0.005604286"   
##  [23] "Component \"6057825017_R05C01\": Mean relative difference: 0.003901175"   
##  [24] "Component \"6057825017_R05C02\": Mean relative difference: 0.009728503"   
##  [25] "Component \"6057825017_R06C01\": Mean relative difference: 0.00735733"    
##  [26] "Component \"6057825018_R01C01\": Mean relative difference: 0.004300583"   
##  [27] "Component \"6057825018_R02C02\": Mean relative difference: 0.006735824"   
##  [28] "Component \"6057825018_R03C02\": Mean relative difference: 0.005712804"   
##  [29] "Component \"6057825018_R04C01\": Mean relative difference: 0.0008979274"  
##  [30] "Component \"6057825018_R04C02\": Mean relative difference: 0.0006846581"  
##  [31] "Component \"6057825018_R05C01\": Mean relative difference: 0.005401661"   
##  [32] "Component \"6057825018_R05C02\": Mean relative difference: 0.003842435"   
##  [33] "Component \"6057825018_R06C01\": Mean relative difference: 0.0008588087"  
##  [34] "Component \"6042316085_R01C01\": Mean relative difference: 0.001455321"   
##  [35] "Component \"6042316085_R01C02\": Mean relative difference: 0.0009434173"  
##  [36] "Component \"6042316085_R02C02\": Mean relative difference: 0.001338882"   
##  [37] "Component \"6042316085_R03C01\": Mean relative difference: 0.003708021"   
##  [38] "Component \"6042316085_R03C02\": Mean relative difference: 4.521012e-05"  
##  [39] "Component \"6042316085_R05C01\": Mean relative difference: 0.007564031"   
##  [40] "Component \"6042316085_R05C02\": Mean relative difference: 0.004838072"   
##  [41] "Component \"6042316085_R06C01\": Mean relative difference: 0.002115199"   
##  [42] "Component \"6042316107_R01C02\": Mean relative difference: 0.007101057"   
##  [43] "Component \"6042316107_R03C01\": Mean relative difference: 0.001350662"   
##  [44] "Component \"6042316107_R04C02\": Mean relative difference: 0.001264206"   
##  [45] "Component \"6042316107_R05C01\": Mean relative difference: 0.004900669"   
##  [46] "Component \"6042316107_R05C02\": Mean relative difference: 0.002254598"   
##  [47] "Component \"6042316107_R06C01\": Mean relative difference: 0.004357517"   
##  [48] "Component \"6042316107_R06C02\": Mean relative difference: 0.004420433"   
##  [49] "Component \"6042316113_R01C01\": Mean relative difference: 0.002129347"   
##  [50] "Component \"6042316113_R01C02\": Mean relative difference: 0.002693894"   
##  [51] "Component \"6042316113_R02C01\": Mean relative difference: 0.005289272"   
##  [52] "Component \"6042316113_R02C02\": Mean relative difference: 0.00180388"    
##  [53] "Component \"6042316113_R03C02\": Mean relative difference: 0.0001535495"  
##  [54] "Component \"6042316113_R04C02\": Mean relative difference: 0.01827813"    
##  [55] "Component \"6042316113_R05C01\": Mean relative difference: 0.001035428"   
##  [56] "Component \"6042316113_R05C02\": Mean relative difference: 0.00403114"    
##  [57] "Component \"6042316113_R06C01\": Mean relative difference: 0.007511818"   
##  [58] "Component \"6042316127_R01C01\": Mean relative difference: 0.01022482"    
##  [59] "Component \"6042316127_R01C02\": Mean relative difference: 0.005937962"   
##  [60] "Component \"6042316127_R02C02\": Mean relative difference: 0.01230285"    
##  [61] "Component \"6042316127_R03C01\": Mean relative difference: 0.003051783"   
##  [62] "Component \"6042316127_R03C02\": Mean relative difference: 0.00540715"    
##  [63] "Component \"6042316127_R04C01\": Mean relative difference: 0.005236418"   
##  [64] "Component \"6042316127_R04C02\": Mean relative difference: 0.001225861"   
##  [65] "Component \"6042316127_R05C02\": Mean relative difference: 0.004274676"   
##  [66] "Component \"6042316127_R06C01\": Mean relative difference: 0.001821397"   
##  [67] "Component \"6057825014_R01C02.1\": Mean relative difference: 0.004811482" 
##  [68] "Component \"6057825014_R02C02.1\": Mean relative difference: 0.008434969" 
##  [69] "Component \"6057825014_R03C01.1\": Mean relative difference: 0.001517944" 
##  [70] "Component \"6057825014_R03C02.1\": Mean relative difference: 0.003234719" 
##  [71] "Component \"6057825014_R04C01.1\": Mean relative difference: 0.002911234" 
##  [72] "Component \"6057825014_R04C02.1\": Mean relative difference: 0.002700415" 
##  [73] "Component \"6057825014_R06C01.1\": Mean relative difference: 0.003810135" 
##  [74] "Component \"6057825017_R01C01.1\": Mean relative difference: 0.002889384" 
##  [75] "Component \"6057825017_R01C02.1\": Mean relative difference: 0.004760503" 
##  [76] "Component \"6057825017_R02C01\": Mean relative difference: 0.004550713"   
##  [77] "Component \"6057825017_R02C02.1\": Mean relative difference: 0.0004553485"
##  [78] "Component \"6057825017_R03C01.1\": Mean relative difference: 0.0007399756"
##  [79] "Component \"6057825017_R03C02.1\": Mean relative difference: 0.002674942" 
##  [80] "Component \"6057825017_R04C01.1\": Mean relative difference: 0.00652685"  
##  [81] "Component \"6057825017_R04C02.1\": Mean relative difference: 0.0009028381"
##  [82] "Component \"6057825017_R05C01.1\": Mean relative difference: 0.001304165" 
##  [83] "Component \"6057825017_R05C02.1\": Mean relative difference: 0.00577515"  
##  [84] "Component \"6057825017_R06C01.1\": Mean relative difference: 0.0007407852"
##  [85] "Component \"6057825018_R01C01.1\": Mean relative difference: 0.004654957" 
##  [86] "Component \"6057825018_R02C02.1\": Mean relative difference: 0.005804405" 
##  [87] "Component \"6057825018_R03C01\": Mean relative difference: 0.01599167"    
##  [88] "Component \"6057825018_R03C02.1\": Mean relative difference: 0.005075563" 
##  [89] "Component \"6057825018_R04C01.1\": Mean relative difference: 0.01015456"  
##  [90] "Component \"6057825018_R04C02.1\": Mean relative difference: 0.0090868"   
##  [91] "Component \"6057825018_R05C01.1\": Mean relative difference: 0.006965299" 
##  [92] "Component \"6057825018_R05C02.1\": Mean relative difference: 0.009145731" 
##  [93] "Component \"6057825018_R06C01.1\": Mean relative difference: 0.008648112" 
##  [94] "Component \"6042316035_R01C01\": Mean relative difference: 0.003027193"   
##  [95] "Component \"6042316035_R02C02\": Mean relative difference: 0.03023256"    
##  [96] "Component \"6042316035_R03C01\": Mean relative difference: 0.002645247"   
##  [97] "Component \"6042316035_R03C02\": Mean relative difference: 0.0003776781"  
##  [98] "Component \"6042316035_R04C01\": Mean relative difference: 0.01634874"    
##  [99] "Component \"6042316035_R05C01\": Mean relative difference: 0.01438354"    
## [100] "Component \"6042316035_R05C02\": Mean relative difference: 0.01879499"    
## [101] "Component \"6042316035_R06C02\": Mean relative difference: 0.004194196"   
## [102] "Component \"6042316048_R01C01\": Mean relative difference: 0.01093755"    
## [103] "Component \"6042316048_R01C02\": Mean relative difference: 0.01100715"    
## [104] "Component \"6042316048_R02C01\": Mean relative difference: 0.003109089"   
## [105] "Component \"6042316048_R03C01\": Mean relative difference: 0.0147153"     
## [106] "Component \"6042316048_R03C02\": Mean relative difference: 0.002647276"   
## [107] "Component \"6042316048_R04C01\": Mean relative difference: 0.01751674"    
## [108] "Component \"6042316048_R04C02\": Mean relative difference: 0.01271402"    
## [109] "Component \"6042316048_R05C01\": Mean relative difference: 0.01073759"    
## [110] "Component \"6042316048_R05C02\": Mean relative difference: 0.01188172"    
## [111] "Component \"6042316110_R01C01\": Mean relative difference: 0.0077139"     
## [112] "Component \"6042316110_R01C02\": Mean relative difference: 0.05707928"    
## [113] "Component \"6042316110_R03C02\": Mean relative difference: 0.00644626"    
## [114] "Component \"6042316110_R04C01\": Mean relative difference: 0.006377035"   
## [115] "Component \"6042316110_R04C02\": Mean relative difference: 0.01296536"    
## [116] "Component \"6042316110_R05C01\": Mean relative difference: 0.04567919"    
## [117] "Component \"6042316110_R06C01\": Mean relative difference: 0.008424955"   
## [118] "Component \"6042316121_R02C02\": Mean relative difference: 0.02111576"    
## [119] "Component \"6042316121_R03C01\": Mean relative difference: 0.01572267"    
## [120] "Component \"6042316121_R03C02\": Mean relative difference: 0.007822507"   
## [121] "Component \"6042316121_R04C02\": Mean relative difference: 0.02102033"    
## [122] "Component \"6042316121_R05C01\": Mean relative difference: 0.01312783"    
## [123] "Component \"6042316121_R05C02\": Mean relative difference: 0.02589649"    
## [124] "Component \"6042316121_R06C01\": Mean relative difference: 0.009924427"   
## [125] "Component \"6042316121_R06C02\": Mean relative difference: 0.01532472"    
## [126] "Component \"6042316066_R01C01\": Mean relative difference: 0.01686536"    
## [127] "Component \"6042316066_R02C01\": Mean relative difference: 0.003838417"   
## [128] "Component \"6042316066_R02C02\": Mean relative difference: 0.01037721"    
## [129] "Component \"6042316066_R03C01\": Mean relative difference: 0.0353738"     
## [130] "Component \"6042316066_R04C01\": Mean relative difference: 0.01443201"    
## [131] "Component \"6042316066_R04C02\": Mean relative difference: 0.03175767"    
## [132] "Component \"6042316066_R05C01\": Mean relative difference: 0.0248537"     
## [133] "Component \"6042316066_R06C01\": Mean relative difference: 0.0183975"     
## [134] "Component \"6042316066_R06C02\": Mean relative difference: 0.009160613"   
## [135] "Component \"6042316069_R01C01\": Mean relative difference: 0.01583073"    
## [136] "Component \"6042316069_R01C02\": Mean relative difference: 0.002316989"   
## [137] "Component \"6042316069_R02C01\": Mean relative difference: 0.003060593"   
## [138] "Component \"6042316069_R03C01\": Mean relative difference: 0.00690229"    
## [139] "Component \"6042316069_R03C02\": Mean relative difference: 0.02843519"    
## [140] "Component \"6042316069_R04C02\": Mean relative difference: 0.003362693"   
## [141] "Component \"6042316069_R05C01\": Mean relative difference: 0.000999397"   
## [142] "Component \"6042316069_R06C02\": Mean relative difference: 0.02108224"    
## [143] "Component \"6042316094_R01C02\": Mean relative difference: 0.01316861"    
## [144] "Component \"6042316094_R02C01\": Mean relative difference: 0.01480408"    
## [145] "Component \"6042316094_R03C02\": Mean relative difference: 0.009648902"   
## [146] "Component \"6042316094_R04C01\": Mean relative difference: 0.0001049265"  
## [147] "Component \"6042316094_R04C02\": Mean relative difference: 0.03967377"    
## [148] "Component \"6042316094_R05C01\": Mean relative difference: 0.01665309"    
## [149] "Component \"6042316094_R05C02\": Mean relative difference: 0.01027852"    
## [150] "Component \"6042316094_R06C01\": Mean relative difference: 0.009510692"   
## [151] "Component \"6042316099_R01C01\": Mean relative difference: 0.01646548"    
## [152] "Component \"6042316099_R01C02\": Mean relative difference: 0.00994956"    
## [153] "Component \"6042316099_R02C01\": Mean relative difference: 0.01072279"    
## [154] "Component \"6042316099_R02C02\": Mean relative difference: 0.01530076"    
## [155] "Component \"6042316099_R03C01\": Mean relative difference: 0.03798235"    
## [156] "Component \"6042316099_R04C01\": Mean relative difference: 0.02910005"    
## [157] "Component \"6042316099_R04C02\": Mean relative difference: 0.02417453"    
## [158] "Component \"6042316099_R05C01\": Mean relative difference: 0.009145505"   
## [159] "Component \"6969568082_R02C01\": Mean relative difference: 0.02497174"    
## [160] "Component \"6969568082_R06C01\": Mean relative difference: 0.004244367"   
## [161] "Component \"6969568082_R02C02\": Mean relative difference: 0.01949863"    
## [162] "Component \"6969568082_R04C02\": Mean relative difference: 0.01987194"    
## [163] "Component \"6969568082_R06C02\": Mean relative difference: 0.01164703"    
## [164] "Component \"6969568084_R01C01\": Mean relative difference: 0.01986453"    
## [165] "Component \"6969568084_R02C01\": Mean relative difference: 0.03384333"    
## [166] "Component \"6969568084_R03C01\": Mean relative difference: 0.02316603"    
## [167] "Component \"6969568084_R04C01\": Mean relative difference: 0.003767401"   
## [168] "Component \"6969568084_R06C01\": Mean relative difference: 0.01556818"    
## [169] "Component \"6969568084_R02C02\": Mean relative difference: 0.008951071"   
## [170] "Component \"6969568084_R03C02\": Mean relative difference: 0.02949538"    
## [171] "Component \"6969568084_R04C02\": Mean relative difference: 0.007451549"   
## [172] "Component \"6969568084_R05C02\": Mean relative difference: 0.01840231"    
## [173] "Component \"6969568087_R01C01\": Mean relative difference: 0.01247348"    
## [174] "Component \"6969568087_R02C01\": Mean relative difference: 0.02518891"    
## [175] "Component \"6969568087_R03C01\": Mean relative difference: 0.009249288"   
## [176] "Component \"6969568087_R04C01\": Mean relative difference: 0.008469972"   
## [177] "Component \"6969568087_R05C01\": Mean relative difference: 0.03088827"    
## [178] "Component \"6969568087_R06C01\": Mean relative difference: 0.0269358"     
## [179] "Component \"6969568087_R01C02\": Mean relative difference: 0.008633143"   
## [180] "Component \"6969568087_R02C02\": Mean relative difference: 0.0002420357"  
## [181] "Component \"6969568087_R03C02\": Mean relative difference: 0.01837536"    
## [182] "Component \"6969568087_R05C02\": Mean relative difference: 0.02148458"    
## [183] "Component \"6969568087_R06C02\": Mean relative difference: 0.02778382"    
## [184] "Component \"6969568118_R01C01\": Mean relative difference: 0.01855533"    
## [185] "Component \"6969568118_R02C01\": Mean relative difference: 0.01843133"    
## [186] "Component \"6969568118_R03C01\": Mean relative difference: 0.006398451"   
## [187] "Component \"6969568118_R04C01\": Mean relative difference: 0.02767989"    
## [188] "Component \"6969568118_R01C02\": Mean relative difference: 0.02158086"    
## [189] "Component \"6969568118_R02C02\": Mean relative difference: 0.003942197"   
## [190] "Component \"6969568118_R03C02\": Mean relative difference: 0.01327359"    
## [191] "Component \"6969568118_R04C02\": Mean relative difference: 0.001527506"   
## [192] "Component \"6969568118_R06C02\": Mean relative difference: 0.02226645"    
## [193] "Component \"6929726046_R02C01\": Mean relative difference: 0.01315447"    
## [194] "Component \"6929726046_R05C01\": Mean relative difference: 0.01841901"    
## [195] "Component \"6929726046_R06C01\": Mean relative difference: 0.01144575"    
## [196] "Component \"6929726046_R01C02\": Mean relative difference: 0.02320518"    
## [197] "Component \"6929726046_R03C02\": Mean relative difference: 0.01544129"    
## [198] "Component \"6929726046_R04C02\": Mean relative difference: 0.03540608"    
## [199] "Component \"6929726046_R06C02\": Mean relative difference: 0.008967212"   
## [200] "Component \"6929718123_R01C01\": Mean relative difference: 0.02675697"    
## [201] "Component \"6929718123_R04C01\": Mean relative difference: 0.004477023"   
## [202] "Component \"6929718123_R06C01\": Mean relative difference: 0.02349603"    
## [203] "Component \"6929718123_R01C02\": Mean relative difference: 0.009016129"   
## [204] "Component \"6929718123_R02C02\": Mean relative difference: 0.02496489"    
## [205] "Component \"6929718123_R03C02\": Mean relative difference: 0.0214809"     
## [206] "Component \"6929718123_R04C02\": Mean relative difference: 0.01504211"    
## [207] "Component \"6929718123_R05C02\": Mean relative difference: 0.0160026"     
## [208] "Component \"6929718123_R06C02\": Mean relative difference: 0.01886059"    
## [209] "Component \"6929718136_R01C01\": Mean relative difference: 0.02502001"    
## [210] "Component \"6929718136_R02C01\": Mean relative difference: 0.01157104"    
## [211] "Component \"6929718136_R03C01\": Mean relative difference: 0.01297259"    
## [212] "Component \"6929718136_R04C01\": Mean relative difference: 0.03218399"    
## [213] "Component \"6929718136_R05C01\": Mean relative difference: 0.01999106"    
## [214] "Component \"6929718136_R02C02\": Mean relative difference: 0.006470846"   
## [215] "Component \"6929718136_R04C02\": Mean relative difference: 0.00467562"    
## [216] "Component \"6929718136_R05C02\": Mean relative difference: 0.01177587"    
## [217] "Component \"6929718136_R06C02\": Mean relative difference: 0.01795778"    
## [218] "Component \"6929718138_R01C01\": Mean relative difference: 0.008411061"   
## [219] "Component \"6929718138_R02C01\": Mean relative difference: 0.0106125"     
## [220] "Component \"6929718138_R04C01\": Mean relative difference: 0.01269295"    
## [221] "Component \"6929718138_R06C01\": Mean relative difference: 0.0187008"     
## [222] "Component \"6929718138_R01C02\": Mean relative difference: 0.02594137"    
## [223] "Component \"6929718138_R03C02\": Mean relative difference: 0.01630087"    
## [224] "Component \"6929718138_R04C02\": Mean relative difference: 0.02199839"    
## [225] "Component \"6929718138_R05C02\": Mean relative difference: 0.01875333"    
## [226] "Component \"6929718138_R06C02\": Mean relative difference: 0.0242575"     
## [227] "Component \"6042316054_R02C01\": Mean relative difference: 0.03359797"    
## [228] "Component \"6042316054_R02C02\": Mean relative difference: 0.01439657"    
## [229] "Component \"6042316054_R03C01\": Mean relative difference: 0.01273794"    
## [230] "Component \"6042316054_R03C02\": Mean relative difference: 0.006297969"   
## [231] "Component \"6042316054_R04C01\": Mean relative difference: 0.003470436"   
## [232] "Component \"6042316054_R04C02\": Mean relative difference: 0.02455132"    
## [233] "Component \"6042316054_R05C01\": Mean relative difference: 0.02576663"    
## [234] "Component \"6042316054_R05C02\": Mean relative difference: 0.001818027"   
## [235] "Component \"6042316054_R06C01\": Mean relative difference: 0.02760376"    
## [236] "Component \"6042316063_R02C01\": Mean relative difference: 0.003437754"   
## [237] "Component \"6042316063_R02C02\": Mean relative difference: 0.01642774"    
## [238] "Component \"6042316063_R03C01\": Mean relative difference: 0.002151564"   
## [239] "Component \"6042316063_R03C02\": Mean relative difference: 0.004027741"   
## [240] "Component \"6042316063_R04C01\": Mean relative difference: 0.03504399"    
## [241] "Component \"6042316063_R04C02\": Mean relative difference: 0.009293138"   
## [242] "Component \"6042316063_R05C01\": Mean relative difference: 0.01821314"    
## [243] "Component \"6042316063_R05C02\": Mean relative difference: 0.005708954"   
## [244] "Component \"6042316063_R06C02\": Mean relative difference: 0.03552747"    
## [245] "Component \"6042316065_R01C02\": Mean relative difference: 0.003875907"   
## [246] "Component \"6042316065_R02C02\": Mean relative difference: 0.01532199"    
## [247] "Component \"6042316065_R03C01\": Mean relative difference: 0.00915521"    
## [248] "Component \"6042316065_R04C01\": Mean relative difference: 0.01737185"    
## [249] "Component \"6042316065_R04C02\": Mean relative difference: 0.00364056"    
## [250] "Component \"6042316065_R05C02\": Mean relative difference: 0.00592463"    
## [251] "Component \"6042316065_R06C02\": Mean relative difference: 0.02956687"    
## [252] "Component \"6042316103_R02C01\": Mean relative difference: 0.02327368"    
## [253] "Component \"6042316103_R03C01\": Mean relative difference: 0.03009274"    
## [254] "Component \"6042316103_R03C02\": Mean relative difference: 0.009467194"   
## [255] "Component \"6042316103_R04C01\": Mean relative difference: 0.01251461"    
## [256] "Component \"6042316103_R05C01\": Mean relative difference: 0.006469341"   
## [257] "Component \"6042316103_R06C01\": Mean relative difference: 0.008503333"   
## [258] "Component \"6042316103_R06C02\": Mean relative difference: 0.02037825"    
## [259] "Component \"6042316036_R01C02\": Mean relative difference: 0.01308609"    
## [260] "Component \"6042316036_R02C01\": Mean relative difference: 0.01447532"    
## [261] "Component \"6042316036_R03C01\": Mean relative difference: 0.02432675"    
## [262] "Component \"6042316036_R03C02\": Mean relative difference: 0.00818478"    
## [263] "Component \"6042316036_R04C01\": Mean relative difference: 0.01413332"    
## [264] "Component \"6042316036_R04C02\": Mean relative difference: 0.02022778"    
## [265] "Component \"6042316036_R05C01\": Mean relative difference: 0.002319252"   
## [266] "Component \"6042316036_R05C02\": Mean relative difference: 0.01761138"    
## [267] "Component \"6042316036_R06C01\": Mean relative difference: 0.001718769"   
## [268] "Component \"6042316036_R06C02\": Mean relative difference: 0.01558359"    
## [269] "Component \"6042316050_R01C02\": Mean relative difference: 0.008409702"   
## [270] "Component \"6042316050_R02C02\": Mean relative difference: 0.00404109"    
## [271] "Component \"6042316050_R03C01\": Mean relative difference: 0.007783286"   
## [272] "Component \"6042316050_R04C01\": Mean relative difference: 0.008987079"   
## [273] "Component \"6042316050_R04C02\": Mean relative difference: 0.02337768"    
## [274] "Component \"6042316050_R05C01\": Mean relative difference: 0.02513576"    
## [275] "Component \"6042316050_R05C02\": Mean relative difference: 0.01389189"    
## [276] "Component \"6042316050_R06C01\": Mean relative difference: 0.02570208"    
## [277] "Component \"6042316050_R06C02\": Mean relative difference: 0.0002740546"  
## [278] "Component \"6042316053_R01C01\": Mean relative difference: 0.003851543"   
## [279] "Component \"6042316053_R02C01\": Mean relative difference: 0.01206842"    
## [280] "Component \"6042316053_R02C02\": Mean relative difference: 0.01234632"    
## [281] "Component \"6042316053_R03C01\": Mean relative difference: 0.007941091"   
## [282] "Component \"6042316053_R03C02\": Mean relative difference: 0.004994915"   
## [283] "Component \"6042316053_R04C01\": Mean relative difference: 0.01976316"    
## [284] "Component \"6042316053_R04C02\": Mean relative difference: 0.002873379"   
## [285] "Component \"6042316053_R06C01\": Mean relative difference: 0.02522915"    
## [286] "Component \"6042316053_R06C02\": Mean relative difference: 0.02408628"    
## [287] "Component \"6042316061_R01C01\": Mean relative difference: 0.004598592"   
## [288] "Component \"6042316061_R02C01\": Mean relative difference: 0.007523021"   
## [289] "Component \"6042316061_R03C01\": Mean relative difference: 0.02193942"    
## [290] "Component \"6042316061_R03C02\": Mean relative difference: 0.003325828"   
## [291] "Component \"6042316061_R04C01\": Mean relative difference: 0.02457436"    
## [292] "Component \"6042316061_R04C02\": Mean relative difference: 0.001615841"   
## [293] "Component \"6042316061_R05C01\": Mean relative difference: 0.0166475"     
## [294] "Component \"6042316061_R05C02\": Mean relative difference: 0.01626925"    
## [295] "Component \"7796806022_R01C01\": Mean relative difference: 0.03126923"    
## [296] "Component \"7796806022_R03C01\": Mean relative difference: 0.002921855"   
## [297] "Component \"7796806022_R04C01\": Mean relative difference: 0.004296551"   
## [298] "Component \"7796806022_R04C02\": Mean relative difference: 0.01028703"    
## [299] "Component \"7796806022_R05C02\": Mean relative difference: 0.002706806"   
## [300] "Component \"7796806022_R06C01\": Mean relative difference: 0.01283297"    
## [301] "Component \"7786923046_R03C02\": Mean relative difference: 0.007313817"   
## [302] "Component \"7786923046_R04C01\": Mean relative difference: 0.01971045"    
## [303] "Component \"7786923046_R06C01\": Mean relative difference: 0.0216028"     
## [304] "Component \"7796806038_R03C01\": Mean relative difference: 0.0006774355"  
## [305] "Component \"7796806038_R03C02\": Mean relative difference: 0.004767744"   
## [306] "Component \"7796806038_R05C02\": Mean relative difference: 0.0124919"     
## [307] "Component \"7786923063_R01C01\": Mean relative difference: 0.01087149"    
## [308] "Component \"7786923063_R01C02\": Mean relative difference: 0.0004053854"  
## [309] "Component \"7786923063_R03C01\": Mean relative difference: 0.01197463"    
## [310] "Component \"7786923063_R04C01\": Mean relative difference: 0.009423631"   
## [311] "Component \"7786923063_R06C01\": Mean relative difference: 0.0119703"     
## [312] "Component \"7786923107_R01C02\": Mean relative difference: 0.01822742"    
## [313] "Component \"7786923107_R04C01\": Mean relative difference: 0.0100415"     
## [314] "Component \"7786923107_R05C01\": Mean relative difference: 0.01576623"    
## [315] "Component \"7786923107_R06C01\": Mean relative difference: 0.003573157"   
## [316] "Component \"7796806016_R03C02\": Mean relative difference: 0.01018704"    
## [317] "Component \"7796806016_R06C01\": Mean relative difference: 0.02116778"    
## [318] "Component \"7796806002_R01C01\": Mean relative difference: 0.02811815"    
## [319] "Component \"7796806002_R02C02\": Mean relative difference: 0.01364491"    
## [320] "Component \"7796806002_R04C02\": Mean relative difference: 0.0159756"     
## [321] "Component \"7796806002_R05C01\": Mean relative difference: 0.03190387"    
## [322] "Component \"7796806002_R05C02\": Mean relative difference: 0.03117594"    
## [323] "Component \"7796806002_R06C01\": Mean relative difference: 0.01673692"    
## [324] "Component \"7796806002_R06C02\": Mean relative difference: 0.04228396"    
## [325] "Component \"7796806029_R02C01\": Mean relative difference: 0.00656956"    
## [326] "Component \"7796806029_R04C01\": Mean relative difference: 0.01292663"    
## [327] "Component \"7796806029_R06C01\": Mean relative difference: 0.01353076"    
## [328] "Component \"6042316024_R01C01\": Mean relative difference: 0.009438669"   
## [329] "Component \"6042316024_R01C02\": Mean relative difference: 0.01249381"    
## [330] "Component \"6042316024_R02C01\": Mean relative difference: 0.01014034"    
## [331] "Component \"6042316024_R02C02\": Mean relative difference: 0.006919189"   
## [332] "Component \"6042316024_R03C01\": Mean relative difference: 0.006143275"   
## [333] "Component \"6042316024_R04C01\": Mean relative difference: 0.002234526"   
## [334] "Component \"6042316024_R04C02\": Mean relative difference: 0.005218839"   
## [335] "Component \"6042316024_R05C01\": Mean relative difference: 0.004627739"   
## [336] "Component \"6042316024_R05C02\": Mean relative difference: 0.002365027"   
## [337] "Component \"6042316024_R06C01\": Mean relative difference: 0.002380228"   
## [338] "Component \"6042316024_R06C02\": Mean relative difference: 0.0002030805"  
## [339] "Component \"6042316031_R01C01\": Mean relative difference: 0.001963721"   
## [340] "Component \"6042316031_R01C02\": Mean relative difference: 0.003089001"   
## [341] "Component \"6042316031_R02C01\": Mean relative difference: 0.004692797"   
## [342] "Component \"6042316031_R02C02\": Mean relative difference: 0.00503462"    
## [343] "Component \"6042316031_R03C01\": Mean relative difference: 0.003155243"   
## [344] "Component \"6042316031_R03C02\": Mean relative difference: 0.002389146"   
## [345] "Component \"6042316031_R04C01\": Mean relative difference: 0.0008178796"  
## [346] "Component \"6042316031_R04C02\": Mean relative difference: 0.0003864482"  
## [347] "Component \"6042316031_R05C01\": Mean relative difference: 0.0001678988"  
## [348] "Component \"6042316031_R05C02\": Mean relative difference: 0.001242318"   
## [349] "Component \"6042316031_R06C01\": Mean relative difference: 0.002957568"   
## [350] "Component \"6042316031_R06C02\": Mean relative difference: 0.007795855"   
## [351] "Component \"6042316042_R01C01\": Mean relative difference: 0.003740608"   
## [352] "Component \"6042316042_R01C02\": Mean relative difference: 0.0002513863"  
## [353] "Component \"6042316042_R02C01\": Mean relative difference: 0.000185125"   
## [354] "Component \"6042316042_R02C02\": Mean relative difference: 0.00108516"    
## [355] "Component \"6042316042_R03C01\": Mean relative difference: 0.005412411"   
## [356] "Component \"6042316042_R03C02\": Mean relative difference: 2.911581e-06"  
## [357] "Component \"6042316042_R04C01\": Mean relative difference: 0.003008291"   
## [358] "Component \"6042316042_R04C02\": Mean relative difference: 0.0004820865"  
## [359] "Component \"6042316042_R05C01\": Mean relative difference: 0.002736236"   
## [360] "Component \"6042316042_R05C02\": Mean relative difference: 0.001921139"   
## [361] "Component \"6042316042_R06C01\": Mean relative difference: 0.002488837"   
## [362] "Component \"6042316042_R06C02\": Mean relative difference: 0.003984859"   
## [363] "Component \"6042316047_R01C01\": Mean relative difference: 0.001540615"   
## [364] "Component \"6042316047_R02C01\": Mean relative difference: 0.002705884"   
## [365] "Component \"6042316047_R02C02\": Mean relative difference: 0.0008944688"  
## [366] "Component \"6042316047_R04C01\": Mean relative difference: 9.705169e-05"  
## [367] "Component \"6042316047_R04C02\": Mean relative difference: 0.007385043"   
## [368] "Component \"6042316047_R05C01\": Mean relative difference: 0.002602368"   
## [369] "Component \"6042316047_R05C02\": Mean relative difference: 0.0009315526"  
## [370] "Component \"6042316047_R06C01\": Mean relative difference: 0.003594905"   
## [371] "Component \"6055432012_R01C01\": Mean relative difference: 0.02849863"    
## [372] "Component \"6055432012_R01C02\": Mean relative difference: 0.01537963"    
## [373] "Component \"6055432012_R02C02\": Mean relative difference: 0.01978808"    
## [374] "Component \"6055432012_R03C01\": Mean relative difference: 0.01304361"    
## [375] "Component \"6055432012_R03C02\": Mean relative difference: 0.0005229822"  
## [376] "Component \"6055432012_R04C01\": Mean relative difference: 0.05810359"    
## [377] "Component \"6055432012_R04C02\": Mean relative difference: 0.01388945"    
## [378] "Component \"6055432012_R05C01\": Mean relative difference: 0.004895921"   
## [379] "Component \"6055432012_R05C02\": Mean relative difference: 0.000630913"   
## [380] "Component \"6055432012_R06C01\": Mean relative difference: 0.0003375166"  
## [381] "Component \"6055432012_R06C02\": Mean relative difference: 0.001982523"   
## [382] "Component \"6055432029_R01C01\": Mean relative difference: 0.03263762"    
## [383] "Component \"6055432029_R01C02\": Mean relative difference: 0.001096567"   
## [384] "Component \"6055432029_R02C01\": Mean relative difference: 0.01475497"    
## [385] "Component \"6055432029_R02C02\": Mean relative difference: 0.002788254"   
## [386] "Component \"6055432029_R03C01\": Mean relative difference: 0.005187445"   
## [387] "Component \"6055432029_R03C02\": Mean relative difference: 0.01411945"    
## [388] "Component \"6055432029_R04C01\": Mean relative difference: 0.007489459"   
## [389] "Component \"6055432029_R04C02\": Mean relative difference: 0.0007203887"  
## [390] "Component \"6055432029_R05C01\": Mean relative difference: 0.01067704"    
## [391] "Component \"6055432029_R05C02\": Mean relative difference: 0.01092016"    
## [392] "Component \"6055432029_R06C01\": Mean relative difference: 0.0001390944"  
## [393] "Component \"6055432029_R06C02\": Mean relative difference: 0.01510947"    
## [394] "Component \"6055432060_R01C01\": Mean relative difference: 0.004760981"   
## [395] "Component \"6055432060_R01C02\": Mean relative difference: 0.001484638"   
## [396] "Component \"6055432060_R02C01\": Mean relative difference: 0.007296348"   
## [397] "Component \"6055432060_R02C02\": Mean relative difference: 0.005127782"   
## [398] "Component \"6055432060_R03C01\": Mean relative difference: 0.01195787"    
## [399] "Component \"6055432060_R03C02\": Mean relative difference: 0.02002448"    
## [400] "Component \"6055432060_R04C01\": Mean relative difference: 0.02347892"    
## [401] "Component \"6055432060_R04C02\": Mean relative difference: 0.00254576"    
## [402] "Component \"6055432060_R05C01\": Mean relative difference: 0.01844558"    
## [403] "Component \"6055432060_R05C02\": Mean relative difference: 0.01204773"    
## [404] "Component \"6055432060_R06C01\": Mean relative difference: 0.02459126"    
## [405] "Component \"6055432060_R06C02\": Mean relative difference: 0.01808787"    
## [406] "Component \"6055432066_R01C01\": Mean relative difference: 0.0008818764"  
## [407] "Component \"6055432066_R01C02\": Mean relative difference: 0.005166243"   
## [408] "Component \"6055432066_R02C01\": Mean relative difference: 0.004224405"   
## [409] "Component \"6055432066_R02C02\": Mean relative difference: 0.008351799"   
## [410] "Component \"6055432066_R03C02\": Mean relative difference: 0.0152471"     
## [411] "Component \"6055432066_R04C01\": Mean relative difference: 0.008501154"   
## [412] "Component \"6055432066_R04C02\": Mean relative difference: 0.001558576"   
## [413] "Component \"6055432066_R05C01\": Mean relative difference: 0.004292709"   
## [414] "Component \"6055432066_R05C02\": Mean relative difference: 0.001228855"   
## [415] "Component \"6057825115_R01C01\": Mean relative difference: 0.005943066"   
## [416] "Component \"6057825115_R03C01\": Mean relative difference: 0.01032604"    
## [417] "Component \"6057825115_R05C01\": Mean relative difference: 0.003143275"   
## [418] "Component \"6057825115_R01C02\": Mean relative difference: 0.006664909"   
## [419] "Component \"6057825115_R03C02\": Mean relative difference: 0.006886484"   
## [420] "Component \"6057825115_R05C02\": Mean relative difference: 0.0004889234"  
## [421] "Component \"6057825128_R01C01\": Mean relative difference: 0.009831255"   
## [422] "Component \"6057825128_R03C01\": Mean relative difference: 0.009756067"   
## [423] "Component \"6057825128_R05C01\": Mean relative difference: 0.004545604"   
## [424] "Component \"6057825115_R02C01\": Mean relative difference: 0.0006895463"  
## [425] "Component \"6057825115_R04C01\": Mean relative difference: 0.001223777"   
## [426] "Component \"6057825115_R06C01\": Mean relative difference: 0.004916346"   
## [427] "Component \"6057825115_R02C02\": Mean relative difference: 0.00478114"    
## [428] "Component \"6057825115_R04C02\": Mean relative difference: 0.0004601983"  
## [429] "Component \"6057825115_R06C02\": Mean relative difference: 0.003932085"   
## [430] "Component \"6057825128_R01C02\": Mean relative difference: 0.00420022"    
## [431] "Component \"6057825128_R03C02\": Mean relative difference: 0.0006560607"  
## [432] "Component \"6057825128_R05C02\": Mean relative difference: 0.004039341"
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
