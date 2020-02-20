###### This scripts computes Maximal Expectation (ME) cluster algorithm given 
###### two spaces as inputs:
######    1. Overdispersion multidimensional space
######    2. uMAP bidimensional space
###### For every space the optinum number of groups is computed.
###### Given the optinum number of groups, the asignation of every sample in
###### each group is done. Also, proportions are ploted im two bar plots.
rm(list = ls())
source("./analysis/Funn.R")
library(mclust)
library(umap)

#### Loading Overdispersion matrix File
Var <- read.csv("./data/processed/Var_b.csv",header = TRUE,
                row.names=1, sep = ",")
nm <- colnames(Var)

#### Number of considered genes
ng <- 6000

#### Variance sorted (Decreasing)
rcmvar <- apply(Var, 1, var)

n_uMAP <- vector (mode="integer", length = length(ng))
n_Var <- vector (mode="integer", length = 168)

### Only the top ng genes are considered
vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar),ng)]
Var.f <- Var[vi,]
### Correlation Matrix between cells
Var.Corr <-as.data.frame(cor(Var.f,method = "pearson"))
### Distance matrix
Dis.R <- 1-Var.Corr

##### uMAP computation
Data.umap = umap(Dis.R,random_state=1)

##### Evaluating ME ######
####### Using uMAP Space
  uMAP.BIC <- mclustBIC(Data.umap$layout,G=1:8)
  uMAPSummary <- summary(uMAP.BIC,data = Data.umap$layout)
  n_uMAP <- uMAPSummary[["G"]]
  G2.col <- (uMAPSummary$classification)
  names(G2.col) <- nm
  prop <- Proportion(G2.col,n_uMAP,nm)
   
  #### proportions
  MyplotBar(prop,'ME_uMAP',ng)
  
  write.csv(G2.col,file = paste("./data/processed/ME_Class_uMAP_",
                                toString(i),".csv",sep = ""),
            quote = FALSE, row.names = TRUE)

###### Using Overdipsertion matrix
  Var.BIC <- mclustBIC(t(Var.f),G=1:8)
  VarSummary <- summary(Var.BIC,data = t(Var.f))
  n_Var <- VarSummary[["G"]]
  G1.col <- (VarSummary$classification)
  nm <- names(G1.col)
  prop <- Proportion(G1.col,n_Var,nm)
  #### proportions
  MyplotBar(prop,'ME_Var',ng)
   
  write.csv(G1.col,file = paste("./data/processed/ME_Class_Var_",
                                 toString(i),".csv",sep = ""),
             quote = FALSE, row.names = TRUE)

#### Ploting clustering results into uMAP Space
MyplotScatter(ng,Data.umap$layout,G1.col,max(G1.col),
               G2.col,max(G2.col),"ME")
