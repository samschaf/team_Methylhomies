library(lattice)
# source("https://bioconductor.org/biocLite.R")
# biocLite("ComplexHeatmap")

library(ggplot2)
library(plyr)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

load("/Volumes/Transcend/QN_BMIQ_cell_type_corrected_data/GSE43414_batch_cell_cor.RData", verbose=TRUE) ## GSE43414_cell_cor
load("/Volumes/Transcend/QN_BMIQ_cell_type_corrected_data/GSE43414_batch_cor.RData", verbose=TRUE) ## GSE43414_batch_cor
load("/Volumes/Transcend/QN_BMIQ_cell_type_corrected_data/Meta_batch_cor.RData", verbose=TRUE) ## meta

meta <- meta[!(meta$braak.stage=="NA" | meta$braak.stage=="Exclude"),]

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
meta_order_by_brain_regions <- meta %>% arrange(broad_regions) %>% arrange(Tissue)
matches_GSE43414_cell_cor <- match(meta_order_by_brain_regions$barcode, rownames(transpose_GSE43414_cell_cor))
GSE43414_cell_cor_sorted_by_brain_regions <- t(transpose_GSE43414_cell_cor[matches_GSE43414_cell_cor,])

#### ===========================================================================
## repeat for bach-corrected only data
## transpose data such that probe names are colnames, and rows are patient samples
transpose_GSE43414_batch_cor <- t(GSE43414_batch_cor)
## order metadata by brain region
# meta_order_by_brain_regions <- meta %>% arrange(Tissue)
matches_GSE43414_batch_cor <- match(meta_order_by_brain_regions$barcode, rownames(transpose_GSE43414_batch_cor))
GSE43414_batch_cor_sorted_by_brain_regions <- t(transpose_GSE43414_batch_cor[matches_GSE43414_batch_cor,])

### ============================================================================
## wilcoxon test
x <- meta[which(meta$broad_regions=="cerebellum"),]$barcode
y <- meta[which(meta$broad_regions=="cortex"),]$barcode
### list of DM probes: cg04197371, cg23112745, and cg06826710
probe.ids <- c("cg04197371", "cg23112745", "cg06826710")

### =======================================================================
## x=barcode names for cerebellum, y=barcode names for cortex, probe.id=list of DM probes, data.type=data of cell type corrected or non-cell type corrected
## outputs the table called probes which gives you a table of probe IDs, p.value, and adjusted p values for each probe
wilcox.test.probes <- function(x, y, probe.ids, data.type){
  probes <- data.frame(probe.ids)
  cerebellum <- data.type[as.character(probe.ids),x]
  cortex <- data.type[as.character(probe.ids),y] 
  for (i in probe.ids){ 
    a <- as.numeric(as.vector(cerebellum[as.character(i),]))
    b <- as.numeric(as.vector(cortex[as.character(i),]))
    probes[which(probes$probe.ids==as.character(i)),"p.value"] <- wilcox.test(a,b)$p.value
  }
  probes$adj.p <- p.adjust(probes$p.value,method="BH")
  probes <- data.frame(probes)
  probes
}
#### ====================================================================

wilcox.test.probes(x=x,y=y,probe.ids=probe.ids,data.type=GSE43414_cell_cor_sorted_by_brain_regions)  
wilcox.test.probes(x=x,y=y,probe.ids=probe.ids,data.type=GSE43414_batch_cor_sorted_by_brain_regions)  


