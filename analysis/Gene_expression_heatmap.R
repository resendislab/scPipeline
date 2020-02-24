####### This script computes the heatmap of DGEA
#### Clean memory
rm(list = ls())

library(scde)

### Number of genes
ng <- c("6000")
logFC <-1.0
### Clustering algorithm
mth <- c("Kmeans",
         "ME")
### Dimensionality reduction space
spc <- c("uMAP","Var")

### Number of cluster for each method over the specific data space
N.Clust <- read.table(file = "./results/data/N_Clusters.csv",
                      sep = "\t", header = TRUE, row.names = 1)


  #### Number of genes id
  id.ng <- 1
  #### Clustering method id
  id.mth <- 1
  fl1 <- (i-1)%/%k.spc+1
  #### Dimensionality reduction space id
  id.spc <- 1

### number of clusters
n <- N.Clust[mth[id.mth],paste(spc[id.spc], ng[id.ng],sep = "_")]
### Clusters
comb <- paste(LETTERS[1:n],collapse = "vs")
###### Globally DEG for cluster A
markers.A <-read.table(paste("./results/data/Intersected_genes/Genes_",
                             "A_",tolower(mth[id.mth]),"_Class_",spc[id.spc],"_",
                             ng[id.ng],"_(logFC_",sprintf("%.1f",logFC)
                             ,")",".txt",sep = ""),
                         sep = "\t", header = FALSE)
###### Globally DEG for cluster B
markers.B <-read.table(paste("./results/data/Intersected_genes/Genes_",
                             "B_",tolower(mth[id.mth]),"_Class_",spc[id.spc],"_",
                             ng[id.ng],"_(logFC_",sprintf("%.1f",logFC)
                             ,")",".txt",sep = ""),
                       sep = "\t", header = FALSE)  
###### Globally DEG for cluster C
markers.C <-read.table(paste("./results/data/Intersected_genes/Genes_",
                                "C_",tolower(mth[id.mth]),"_Class_",spc[id.spc],"_",
                             ng[id.ng],"_(logFC_",sprintf("%.1f",logFC)
                             ,")",".txt",sep = ""),
                       sep = "\t", header = FALSE)
##### All markers
a <- c(as.vector(t(markers.A$V1)),
         as.vector(t(markers.C$V1)))
#a <- as.vector(t(markers.C$V1))
  write.csv2(t(a), file = "./Genetic_signatures/file.csv",quote = FALSE)
##### a few genes
#a <- c("MKI67", "TOP2A","FOXM1","ESR1","PGR","HER2")
######Count Matrix
  fl.data <- paste("./data/processed/Groups_",tolower(mth[id.mth]),"_Class_",
                   spc[id.spc],"/",tolower(mth[id.mth]),"_Class_",
                   spc[id.spc],"_",ng[id.ng],"_",comb,".csv",sep = "")
  MyData <- read.csv(file = fl.data,
                     header = TRUE, row.names = 1, sep = ";")
  ##### Clean count
  cd <- clean.counts(MyData, min.lib.size=100,
                     min.reads = 1, min.detected = 1)
##### Sum by sample  
sumcd <-colSums(as.matrix(cd))
##### Normalization
cd1 <- (cd*1e6)/sumcd
#### Extracting DGEA genes
mat.sub <- as.matrix(cd1[a,])
##### heatmap
heatmap(log(mat.sub+1), Rowv=NA, Colv = NA, scale = "none",
        ColSideColors=rainbow(n)[sg],
        labCol = c(rep("", length(a))),
        col<- colorRampPalette(c("red", "white", "blue"))(100) )
#### Other representation
library("pheatmap")
pheatmap(log(mat.sub+1), cutree_rows = n)
 

