##### This scripts filters differentially expressed genes (DEG) from every
##### cluster compared to rest of them, for example, given the Cluster A DEG
##### for AvsB and AvsC comparison, it gets the intersected genes from both
##### analyses, therefore, the global DEG are obtained.

### Clean memory
rm(list = ls())
library(stringr)

############ CUSTOM FUNCTIONS
Gns.Inter <- function (j,sp,logFC){
  ####### Gets all file names for Cluster j
  fls = dir(path="./results/data/Diff_genes/",
            pattern=paste(LETTERS[j],
                          "_genes_[:(:]logFC_",
                          toString(logFC),sep = ""))
  ###### trim the list of file names for sp analysis
  fls = fls[str_detect(fls, sp)]
  ##### Open all comparison files for Cluster j
  Dat <- lapply(fls, function(x) read.table(
    paste("./results/data/Diff_genes/",x,sep = ""),
    sep = "\t"))
 
  n = length(fls)
  ##### Lists intersections
  A <- intersect(rownames(Dat[[1]]),rownames(Dat[[2]]))
  while ((n-2)>0) {
    A <- intersect(rownames(Dat[[n]]),A)
    n = n-1
  }
  ##### saving file
  write.table(A,file = paste(path,"/Genes_",LETTERS[j],"_",
                             sp,"_(LogFC_",toString(logFC),").txt",sep = ""),
              sep = "\t", col.names = FALSE, quote = FALSE,
              row.names = FALSE)
}
#####################################

########### MAIN #######
### Number of considered genes
val=6000
### log Fold Change
logFC=1.0

fls.name <- c("kmeans_Class_uMAP",
              "kmeans_Class_Var",
              "ME_Class_uMAP",
              "ME_Class_Var")
fl<- 1
id.val <- 0
ng <- length(val)
path <- ("./results/data/Intersected_genes")
if (dir.exists(path)==FALSE){dir.create(path)}

#### For loop for every analysis
for (i in 1:(ng*4)){     
  id.fls <- (i-1)%/%ng+1
  ###  id.fls <- j
  if (id.fls==fl){
    id.val <-  id.val + 1
  }else{
    id.val <-  1
    fl <- fl + 1}
  ### name of root analysis
  sp = paste(fls.name[id.fls],"_",
             toString(val[id.val]),sep = "")
 
  print(paste("/Groups_",sp,".csv",sep ="" ))
  ###### Number of Clusters
  Clust <- read.csv(paste("./data/processed/",sp,".csv",sep ="" ),
                    sep = ",")
  N.clust = max(Clust[,2])
 
  ##### Gns.Inter computes the intersection from the DEG of every
  ##### comparison
  invisible(lapply(1:N.clust, function(x) Gns.Inter(x,sp,logFC)))
}
