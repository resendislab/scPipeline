


#### Clean memory
rm(list = ls())
library(stringr)

dr <- './Cytoscape/Clusters Subgroups/'

fls = dir(path=dr, pattern = ".[:.:]csv$")

grp <- "(BvsC_B)"
### Overexpressed genes using differential expression analysis
# genes.OE <- as.vector(t(
#   read.table(paste("./results/data/Intersected_genes/Genes_",
#                    grp,
#                    "_kmeans_Class_uMAP_6000_(logFC_1.0).txt",
#                    sep = ""),
#                    sep = "\t")))
genes <- read.table(paste("./results/data/Diff_genes/",
                          "Upregulated_kmeans_Class_uMAP_6000",
                          "(Grp_BvsGrp_C)_Grp_B_genes_(logFC_1.0).txt",
                sep = ""),
         sep = "\t", header = TRUE)
genes.OE <- rownames(genes)
genes.OE <- c("ANLN","ANP32E")
for (i in fls){
  tb= read.table(paste(dr,i,sep = ""),
                sep = ";", header = TRUE)
  #### Extract Genes from list
  a <- as.vector(tb$Genes)
  #### Erase spacial characters except ","
  #### ^ inside [] means negation
  a1 <- gsub("[^[:alnum:],-]","",a)
  #### Split string given "," separation
  b <- unlist(strsplit(a1,","))
  #### Genes that englobe the group
  Uniq.Gn <- unique(b)
  #### Intersection between differentially expresed genes and group genes
  Gs <- intersect(Uniq.Gn,genes.OE)
  print(c(i,":",Gs))
  #### save file
#  write.table(Gs,file = paste("./results/data/DEG_and_GSEA/",
#                              "Group_",grp,
#                              "_Genes_intersected_with_",
#                              i,sep = ""), 
#              row.names = FALSE, col.names = FALSE, sep = ";",
#              quote = FALSE)
#  #### intersected genes
}  
#  n <- length(b)
#  In.Gen <- 
#  Dta <- intersect(b[[1]],b[[2]])
#  j <- 2
#  while(j<n){
#    j <- j+1
#    Dta <- intersect(Dta,b[[j]])
#    print(j)
#  }
 

