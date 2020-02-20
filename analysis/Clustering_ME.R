###### This scripts computes Maximal Expectation (ME) cluster algorithm given 
###### two spaces as inputs:
######    1. Overdispersion multidimensional space
######    2. tSNE bidimensional space
###### For every space the optinum number of groups is computed.
###### Given the optinum number of groups, the asignation of every sample in
###### each group is done. Also, proportions are ploted im two bar plots.
setwd("~/Documents/MCF7 Single Cell/MCF7_SC_ddSeeker")
rm(list = ls())
source("./analysis/Funn.R")
library(mclust)
library(Rtsne)

#### Loading Overdispersion matrix File
Var <- read.csv("./data/processed/Var_b.csv",header = TRUE, sep = ",")
Gnames <- Var[,1]
rownames(Var) <- Gnames
Var$X <- NULL
nm <- colnames(Var)

#### Number of considered genes
ng <- 6000

#### Variance sorted (Decreasing)
rcmvar <- apply(Var, 1, var)

n_tsn <- vector (mode="integer", length = 168)
n_Var <- vector (mode="integer", length = 168)
for (i in c(12000)){

ng <- as.integer(i)
### Only the top ng genes are considered
vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar),ng)]
Var.f <- Var[vi,]
### Correlation Matrix between cells
Var.Corr <-as.data.frame(cor(Var.f,method = "pearson"))
### Distance matrix
Dis.R <- 1-Var.Corr

##### tSNE computation
set.seed(1)
#tSNE.pagoda <- Rtsne(Dis.R, perplexity=50,max_iter=3000, verbose = 0)
tSNE.pagoda <- Rtsne(Dis.R, is_distance=TRUE, perplexity=30,max_iter=3000,verbose=0)
#plot(tSNE.pagoda$Y,pch=c(21,24))

##### Evaluating ME ######
##### Using tSNE Space
tSNE.BIC <- mclustBIC(tSNE.pagoda$Y,G=1:8)
tSNESummary <- summary(tSNE.BIC,data = tSNE.pagoda$Y)
n_tSNE <- tSNESummary[["G"]]
G2.col <- (tSNESummary$classification)
names(G2.col) <- nm
#plot(tSNE.pagoda$Y,col="black",pch=c(21,24), 
#     bg=tSNESummary$classification, lwd=0.5,
#      xlab="tSNE[1]",
#      ylab="tSNE[2]")
# title("tSNE Space")
prop <- Proportion(G2.col,n_tSNE,nm)
 
#### proportions
MyplotBar(prop,'ME_tSNE',ng)

write.csv(G2.col,file = paste("./data/processed/ME_Class_tSNE_",
                              toString(i),".csv",sep = ""),
          quote = FALSE, row.names = TRUE)

##Using Overdipsertion matrix
tSNE.BIC <- mclustBIC(t(Var.f),G=1:8)
tSNESummary <- summary(tSNE.BIC,data = t(Var.f))
n_Var <- tSNESummary[["G"]]
G1.col <- (tSNESummary$classification)
nm <- names(G1.col)
#plot(tSNE.pagoda$Y,col="black",pch=c(21,24), 
#      bg=tSNESummary$classification, lwd=0.5,
#      xlab="tSNE[1]",
#      ylab="tSNE[2]")
#title("Overdispersion Space")
prop <- Proportion(G1.col,n_Var,nm)
#### proportions
MyplotBar(prop,'ME_Var',ng)

write.csv(G1.col,file = paste("./data/processed/ME_Class_Var_",
                              toString(i),".csv",sep = ""),
          quote = FALSE, row.names = TRUE)

MyplotScatter(ng,tSNE.pagoda$Y,G1.col,max(G1.col),
              G2.col,max(G2.col),"ME")

}
#write.csv(n_tsn,file="./n_tsn_b.csv",quote= FALSE, row.names=FALSE, col.names=FALSE)
#write.csv(n_Var,file="./n_Var_b.csv",quote= FALSE, row.names=FALSE, col.names=FALSE)
