# Basic Heatmaps (dasen)




```r
load("~/Methylhomies/GSE59685_batch_cor.RData")

load("~/Methylhomies/GSE59685_cell_cor.RData")

load("~/Methylhomies/Meta_brain_cell_batch.RData")
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
transpose_GSE59685_cell_cor <- t(GSE59685_cell_cor)
## order metadata by brain region
#data has barcode NOT gsm as column names
meta_order_by_brain_regions <- meta %>% arrange(Tissue)
matches_GSE59685_cell_cor <- match(meta_order_by_brain_regions$gsm, rownames(transpose_GSE59685_cell_cor))

GSE59685_cell_cor_sorted_by_brain_regions <- t(transpose_GSE59685_cell_cor[matches_GSE59685_cell_cor,])
#repeat for batch

transpose_GSE59685_batch_cor <- t(GSE59685_batch_cor)
## order metadata by brain region
GSE59685_batch_cor_sorted_by_brain_regions <- t(transpose_GSE59685_batch_cor[matches_GSE59685_cell_cor,])

#Check for identical 
identical(colnames(GSE59685_batch_cor_sorted_by_brain_regions),colnames(GSE59685_cell_cor_sorted_by_brain_regions)) # TRUE
```

```
## [1] TRUE
```

```r
identical(colnames(GSE59685_cell_cor_sorted_by_brain_regions),as.character(meta_order_by_brain_regions$barcode)) #TRUE
```

```
## [1] FALSE
```



```r
brain_regions <- meta_order_by_brain_regions$broad_regions
tissue <- meta_order_by_brain_regions$Tissue

set.seed(1)
probes <- sample(rownames(GSE59685_cell_cor_sorted_by_brain_regions),1000)

aheatmap(GSE59685_batch_cor_sorted_by_brain_regions[probes,], Colv = TRUE,annCol = list(Brain_Region = brain_regions, Tissue = tissue), main =  "Uncorrected heatmap dasen", width = 2, height = 2)
```

![](Heatmaps__dasen__files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
aheatmap(GSE59685_cell_cor_sorted_by_brain_regions[probes,], Colv = TRUE, annCol = list(Brain_Region = brain_regions, Tissue = tissue), main = "Corrected heatmap dasen", width = 2, height = 2)
```

![](Heatmaps__dasen__files/figure-html/unnamed-chunk-1-2.png)<!-- -->


```r
# create a combined data set
batch_sub <- GSE59685_batch_cor_sorted_by_brain_regions

cell_sub <- GSE59685_cell_cor_sorted_by_brain_regions
```


```r
#sample to sample correlations based on brain region

aheatmap(cor(batch_sub),Colv = NA, Rowv = NA,  annCol = list(Brain_Region = brain_regions, Tissue = tissue), annRow = list(Brain_Region = brain_regions, Tissue = tissue),main = "Non-Cell Type Corrected dasen")
```

![](Heatmaps__dasen__files/figure-html/correlation-1.png)<!-- -->

```r
aheatmap(cor(cell_sub),Colv = NA, Rowv = NA,  annCol = list(Brain_Region = brain_regions, Tissue = tissue), annRow = list(Brain_Region = brain_regions, Tissue = tissue), main = "Cell Type Corrected dasen")
```

![](Heatmaps__dasen__files/figure-html/correlation-2.png)<!-- -->
