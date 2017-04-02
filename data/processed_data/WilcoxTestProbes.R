## x=barcode names for cerebellum, y=barcode names for cortex
## probe.id=list of DM probes, data.type=cell type corrected or non-cell type corrected data (the entire data frame)
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