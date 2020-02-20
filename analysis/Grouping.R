###### This Script uses the clustering information to re-arange the data into
###### subgroups. It generates two kind of files: first file is a .csv
###### containing count matrix of every Clusters, this file is used in the 
###### GSEA software. The second one is the files considering the multiple 
###### comparisson by pairs; it's used in the Differentially Gene Expression
###### Analysis.
rm(list=ls())
require(cluster)
library(gtools)


Grp.Assigment.ALL <- function (Clust,N.clust,MyData2,Gnames){
### Data for Cluster A
  SubData <- Ext.Data(Clust,1,MyData2)
  j <- 2
  while (j <= N.clust){
    ### Data for Cluster j
    SubData2 <- Ext.Data(Clust,j,MyData2)
  ### Merging Data of groups
    SubData <- merge(SubData,SubData2, by = "row.names", sort = FALSE)
    SubData <- SubData[,-1]
    rownames(SubData) <- Gnames
    j <- j+1
  }
  return(SubData)
} 
Grp.Assigment.Pairs <- function (path,N.clust, MyData2,Gnames,fls,ng){
  ## Possible combinations of two elements
  comb <- combinations(N.clust, 2, 1:N.clust,
                       set=TRUE, repeats.allowed=FALSE)
  n <- dim(comb)
  j <- 1
  while (j <= n[1]){
    ## Data Extraction for Clusters comb[j,]
    SubData1 <- Ext.Data(Clust,comb[j,1],MyData2)
    SubData2 <- Ext.Data(Clust,comb[j,2],MyData2)
    
    SubData <- merge(SubData1,SubData2, by = "row.names", sort = FALSE)
    SubData <- SubData[,-1]
    rownames(SubData) <- Gnames
    ##### Write comparation matrix by pairs
    write.csv2(SubData, file=paste(path,'/',fls,"_",
                                   toString(ng),"_",
                                   LETTERS[comb[j,1]],"vs",
                                   LETTERS[comb[j,2]],
                                   ".csv",sep ="" ),
               quote = FALSE, row.names = TRUE)
    j = j+1
  }
}

Ext.Data <- function(Clust,j,MyData2){
  ## Indexes where sample correpond to Cluster j
  b <- which(Clust[,1] %in% j)
  ## Count matrix for only Cluster j
  SubData1 <- MyData2[,b]
  n <- dim(SubData1)
  ## Name assigment to Group j
  nms<-unlist(lapply(1:n[2],FUN=function(x) paste("Grp",LETTERS[j]
                                                  ,x,sep="_")))
  colnames(SubData1) <- nms
  return(SubData1)
}

## Number of considered genes
ng <- 6000

fls <- c("kmeans_Class_uMAP",
         "kmeans_Class_Var",
         "ME_Class_uMAP",
         "ME_Class_Var")
#################### Load Data
#### Count Matrix
MyData <- read.csv(file = './data/processed/MCTS_b.csv',
                   header = TRUE, row.names = 1, sep = ";")
#### Overdispersion Matrix
Var <- read.csv("./data/processed/Var_b.csv",header = TRUE,
                row.names=1, sep = ",")

### Most Overdispersed genes using the variance, n=ng
rcmvar <- apply(Var, 1, var)
vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar),ng)]
Var.f <- Var[vi,]
idx <- match(rownames(Var.f),rownames(MyData))
### Count Matrix of the ng most dispersed genes
MyData2 <- MyData[idx,]
Gnames <- rownames(MyData2)
fl <- 1
id.val <- 0
####Load Data
for (i in 1:length(fls)){
## Saving path 
  path <- paste('./data/processed/Groups_',fls[i],sep="")
  print(path)
  if (dir.exists(path)==FALSE){dir.create(path)}
## Reading groups labels
  Clust <- read.csv(paste('./data/processed/',fls[i],'_',toString(ng),
                          ".csv",sep ="" ),sep = ",", row.names = 1)
## Getting the number of clusters
  N.clust = max(Clust[,1])
## Grp.Assigment function re-arange Count Matrix acording to the clusterins
## output
 SubData <- Grp.Assigment.ALL(Clust,N.clust,MyData2,Gnames)
#### File with All Groups 
  write.csv2(SubData, file=paste(path,'/',fls[i],'_',
                                 toString(ng),"_",
                                 paste(LETTERS[1:N.clust], collapse = "vs")
                                 ,".csv",sep ="" ),
             quote = FALSE, row.names = TRUE)
#### Files by pairs
  Grp.Assigment.Pairs(path,N.clust, MyData2,Gnames,fls[i],ng)
}
