# correction Matters
Cassia Warren  
March 16, 2017  



**Load the data and check that the columns are in the same order**


```r
load("/Volumes/Lexar/New corrected/GSE43414_batch_cor.RData")

load("/Volumes/Lexar/New corrected/GSE43414_cell_cor.RData")

load("/Volumes/Lexar/New corrected/Meta_batch_cor.RData")
```

Remove Braak stage exludes and NAs

```r
meta2 <- na.omit(meta) #remove NA
meta <- meta2[!c(meta2$braak.stage=="Exclude"),] #remove exlucdes
```



Lisas code to rearrange

```r
## add another broad region column
meta$broad_regions <- ifelse(meta$Tissue == "cerebellum", "cerebellum","cortex")
meta$tissue_color <- lapply(meta$Tissue,function(x){
  if (x == "cerebellum") {y <- "red"}
  if (x == "frontal cortex") {y <- "blue"}
  if (x == "superior temporal gyrus") {y <- "orange"}
  if (x == "entorhinal cortex") {y <- "yellow"}
  y
  })

meta$broad_colors <- lapply(meta$broad_regions,function(x){
  if (x == "cerebellum") {y <- "red"}
  if (x == "cortex") {y <- "blue"}
  y
})

## transpose data such that probe names are colnames, and rows are patient samples
transpose_GSE43414_cell_cor <- t(GSE43414_cell_cor)
## order metadata by brain region
#data has barcode NOT gsm as column names
meta_order_by_brain_regions <- meta %>% arrange(Tissue)
matches_GSE43414_cell_cor <- match(meta_order_by_brain_regions$barcode, rownames(transpose_GSE43414_cell_cor))

GSE43414_cell_cor_sorted_by_brain_regions <- t(transpose_GSE43414_cell_cor[matches_GSE43414_cell_cor,])
#repeat for batch

transpose_GSE43414_batch_cor <- t(GSE43414_batch_cor)
## order metadata by brain region
GSE43414_batch_cor_sorted_by_brain_regions <- t(transpose_GSE43414_batch_cor[matches_GSE43414_cell_cor,])

#Check for identical 
identical(colnames(GSE43414_batch_cor_sorted_by_brain_regions),colnames(GSE43414_cell_cor_sorted_by_brain_regions)) # TRUE
```

```
## [1] TRUE
```

```r
identical(colnames(GSE43414_cell_cor_sorted_by_brain_regions),as.character(meta_order_by_brain_regions$barcode)) #TRUE
```

```
## [1] TRUE
```

**subset 100 Probes so my computer doesnt crash**

```r
brain_regions <- meta_order_by_brain_regions$broad_regions
tissue <- meta_order_by_brain_regions$Tissue
breaksList <- seq(0, 1, by = 0.2)
#sample of 100 probes
set.seed(1)
probes <- sample(rownames(GSE43414_cell_cor_sorted_by_brain_regions),100)
aheatmap(GSE43414_batch_cor_sorted_by_brain_regions[probes,], Colv = NA,annCol = list(Brain_Region = brain_regions, Tissue = tissue), main =  "Uncorrected heatmap of 100 probes")
```

![](Correction_matters.2_files/figure-html/sub-1.png)<!-- -->

```r
aheatmap(GSE43414_batch_cor_sorted_by_brain_regions[probes,], annCol = list(Brain_Region = brain_regions, Tissue = tissue), main =  "Uncorrected heatmap of 100 probes")
```

![](Correction_matters.2_files/figure-html/sub-2.png)<!-- -->

```r
aheatmap(GSE43414_cell_cor_sorted_by_brain_regions[probes,], color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(length(breaksList)), Colv = NA, annCol = list(Brain_Region = brain_regions, Tissue = tissue), main = "Corrected heatmap of 100 probes")
```

![](Correction_matters.2_files/figure-html/sub-3.png)<!-- -->

```r
aheatmap(GSE43414_cell_cor_sorted_by_brain_regions[probes,], annCol = list(Brain_Region = brain_regions, Tissue = tissue), main = "Corrected heatmap of 100 probes")
```

![](Correction_matters.2_files/figure-html/sub-4.png)<!-- -->

**Want to see what thier expression is and if it changes between corrected and non corercted data sets. **

```r
#compare expression of those 100 probes between data sets
# create a combined data set
batch_sub <- GSE43414_batch_cor_sorted_by_brain_regions[probes,]
t.batch_sub <- t(batch_sub)
colnames(t.batch_sub) <- gsub("cg", "batch", colnames(t.batch_sub))

cell_sub <- GSE43414_cell_cor_sorted_by_brain_regions[probes,]
t.cell_sub <- t(cell_sub)
colnames(t.cell_sub) <- gsub("cg", "cell", colnames(t.cell_sub))

identical(rownames(t.cell_sub),as.character(meta_order_by_brain_regions$barcode))
```

```
## [1] TRUE
```

```r
identical(rownames(t.batch_sub),as.character(meta_order_by_brain_regions$barcode))
```

```
## [1] TRUE
```

```r
Meta_probes <- data.frame(meta_order_by_brain_regions, t.cell_sub, t.batch_sub)

cell_met <- data.frame(Meta_probes[,c(2:4,12,17)], t.cell_sub)
# keeps "gms", "Subject", "barcode", "Tissue","broad_regions"
batch_met <- data.frame(Meta_probes[,c(2:4,12,17)], t.batch_sub)

cell_melt <- melt(cell_met, id.vars= c("gsm", "Tissue", "broad_regions", "Subject", "barcode"))
cell_melt$Data.set <- "cell_corrected"

batch_melt <- melt(batch_met, id.vars= c("gsm", "Tissue", "broad_regions", "Subject", "barcode"))
batch_melt$Data.set <- "not_corrected"

meta_probes_meled <- rbind(batch_melt,cell_melt)
meta_probes_meled$variable <- gsub("cell", "cg", meta_probes_meled$variable)
meta_probes_meled$variable <- gsub("batch", "cg", meta_probes_meled$variable)
```


```r
strip.plot <- function(t){
stripplot(t$variable ~ t$value | t$Data.set, t, groups = Tissue, layout = c(2, 1), auto.key = TRUE)
}

strip.plot(meta_probes_meled)
```

![](Correction_matters.2_files/figure-html/plots-1.png)<!-- -->

```r
#facet wrap based on brain region

strip.plot.facet <- function(t){
stripplot(t$variable ~ t$value | t$Data.set*t$Tissue , t, groups = Tissue, auto.key = TRUE)
}

strip.plot.facet(meta_probes_meled)
```

![](Correction_matters.2_files/figure-html/plots-2.png)<!-- -->

```r
#subset patients
subset_patient_6 <- meta_probes_meled[meta_probes_meled$Subject==6,]

strip.plot(subset_patient_6)
```

![](Correction_matters.2_files/figure-html/plots-3.png)<!-- -->

```r
strip.plot.facet(subset_patient_6)
```

![](Correction_matters.2_files/figure-html/plots-4.png)<!-- -->


```r
#sample to sample correlations based on brain region


aheatmap(cor(batch_sub),Colv = NA, Rowv = NA,  annCol = list(Brain_Region = brain_regions, Tissue = tissue), annRow = list(Brain_Region = brain_regions, Tissue = tissue),main = "uncorrected sample to sample correlation of 100 probes")
```

![](Correction_matters.2_files/figure-html/correlation-1.png)<!-- -->

```r
aheatmap(cor(cell_sub),Colv = NA, Rowv = NA,  annCol = list(Brain_Region = brain_regions, Tissue = tissue), annRow = list(Brain_Region = brain_regions, Tissue = tissue), main = "corrected sample to sample correlation of 100 probes")
```

![](Correction_matters.2_files/figure-html/correlation-2.png)<!-- -->
**Both correlated a lot with eachother**





```r
plotfunction2 <- function(n)
{
  ggplot(n, aes(factor(Tissue),value)) +
    geom_boxplot(aes(col = Tissue)) +
   facet_wrap(~Data.set) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tissue") +
  ylab("Methylation") +
  labs(title = paste0("Methylation values for probe ", n$variable))
  
}

plotfunction2(meta_probes_meled[meta_probes_meled$variable== "cg20935223",] )
```

![](Correction_matters.2_files/figure-html/boxes-1.png)<!-- -->

```r
plotfunction2(meta_probes_meled[meta_probes_meled$variable== "cg09766628",] )
```

![](Correction_matters.2_files/figure-html/boxes-2.png)<!-- -->

```r
plotfunction2(meta_probes_meled[meta_probes_meled$variable== "cg15821319",] )
```

![](Correction_matters.2_files/figure-html/boxes-3.png)<!-- -->

```r
plotfunction2(meta_probes_meled[meta_probes_meled$variable== "cg09654300",] )
```

![](Correction_matters.2_files/figure-html/boxes-4.png)<!-- -->

cg20745248
Box plots, combine cortexes but keep seperated


```r
plotfunction3<- function(n)
{
  ggplot(n, aes(factor(broad_regions),value)) +
    geom_boxplot(aes(col = Tissue)) +
   facet_wrap(~Data.set) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tissue") +
  ylab("Methylation") +
  labs(title = paste0("Methylation values for probe ", n$variable))
  
}



plotfunction3(meta_probes_meled[meta_probes_meled$variable== "cg20935223",] )
```

![](Correction_matters.2_files/figure-html/b2-1.png)<!-- -->

```r
plotfunction3(meta_probes_meled[meta_probes_meled$variable== "cg06528575",] )
```

![](Correction_matters.2_files/figure-html/b2-2.png)<!-- -->

combine 

```r
plotfunction4<- function(n)
{
  ggplot(n, aes(factor(broad_regions),value)) +
    geom_boxplot(aes(col =broad_regions)) +
   facet_wrap(~Data.set) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tissue") +
  ylab("Methylation") +
  labs(title = paste0("Methylation values for probe ", n$variable))
  
}

plotfunction4(meta_probes_meled[meta_probes_meled$variable== "cg20935223",] )
```

![](Correction_matters.2_files/figure-html/b3-1.png)<!-- -->

```r
plotfunction4(meta_probes_meled[meta_probes_meled$variable== "cg14813579",] )
```

![](Correction_matters.2_files/figure-html/b3-2.png)<!-- -->

one in mill paper cg02573091
sup cg01156747, cg08835221, cg06242242. cg20636526, cg12362118




```r
cell.probes <- read.csv("/Volumes/Lexar/lists for differental expression/cellprobes.csv")
colnames(cell.probes) <- "cell.probes"
  
batch.probes <- read.csv("/Volumes/Lexar/lists for differental expression/batchprobes.csv")
colnames(batch.probes) <- "batch.probes"
  
both.probes <- read.csv("/Volumes/Lexar/lists for differental expression/bothprobes.csv")
colnames(both.probes) <- "both.probes"

set.seed(8)
#cell.pr <- sample(t(cell.probes),1)
#batch.pr <- sample(t(batch.probes),1)
#both.pr <- sample(t(both.probes),1)
cell.pr <- "cg02707152" 
batch.pr <- "cg26156279"
both.pr <- "cg04427707" 
probes.2 <- c(cell.pr,batch.pr,both.pr)

#compare expression of those 100 probes between data sets
# create a combined data set
batch_sub <- GSE43414_batch_cor_sorted_by_brain_regions[probes.2,]
t.batch_sub <- t(batch_sub)
colnames(t.batch_sub) <- gsub("cg", "batch", colnames(t.batch_sub))

cell_sub <- GSE43414_cell_cor_sorted_by_brain_regions[probes.2,]
t.cell_sub <- t(cell_sub)
colnames(t.cell_sub) <- gsub("cg", "cell", colnames(t.cell_sub))

identical(rownames(t.cell_sub),as.character(meta_order_by_brain_regions$barcode))
```

```
## [1] TRUE
```

```r
identical(rownames(t.batch_sub),as.character(meta_order_by_brain_regions$barcode))
```

```
## [1] TRUE
```

```r
identical(as.character(cell_met$gsm),as.character(batch_met$gsm))
```

```
## [1] TRUE
```

```r
Meta_probes <- data.frame(meta_order_by_brain_regions, t.cell_sub, t.batch_sub)

cell_met <- data.frame(Meta_probes[,c(2:4,12,17)], t.cell_sub)
# keeps "gms", "Subject", "barcode", "Tissue","broad_regions"
batch_met <- data.frame(Meta_probes[,c(2:4,12,17)], t.batch_sub)

cell_melt <- melt(cell_met, id.vars= c("gsm", "Tissue", "broad_regions", "Subject", "barcode"))
cell_melt$Data.set <- "cell_corrected"

batch_melt <- melt(batch_met, id.vars= c("gsm", "Tissue", "broad_regions", "Subject", "barcode"))
batch_melt$Data.set <- "not_corrected"

meta_probes_meled <- rbind(batch_melt,cell_melt)
meta_probes_meled$variable <- gsub("cell", "cg", meta_probes_meled$variable)
meta_probes_meled$variable <- gsub("batch", "cg", meta_probes_meled$variable)
```

Top one only found in cell 

```r
plotfunction4(meta_probes_meled[meta_probes_meled$variable== cell.pr,] )
```

![](Correction_matters.2_files/figure-html/box2-1.png)<!-- -->

Top one only in batch

```r
plotfunction4(meta_probes_meled[meta_probes_meled$variable== batch.pr,] )
```

![](Correction_matters.2_files/figure-html/mo-1.png)<!-- -->

Top found in both

```r
plotfunction4(meta_probes_meled[meta_probes_meled$variable== both.pr,] )
```

![](Correction_matters.2_files/figure-html/mre-1.png)<!-- -->





