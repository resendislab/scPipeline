###### This scripts computes kmeans cluster algorithm given two spaces as inputs:
######    1. Overdispersion multidimensional space
######    2. tSNE bidimensional space
###### For every space the optinum number of groups is computed by calculating
###### the SSE and taking the second derivative to get the most significative 
###### change in the SSE reduction
###### Given the optinum number of groups, the asignation of every sample in
###### each group is done. Also, proportions are ploted im two bar plots.
rm(list=ls())
#setwd("~/Documents/MCF7\ Single\ Cell/MCF7_SC_ddSeeker")
source("./analysis/Funn.R")
library(factoextra)
library(Rtsne)
library(NbClust)
library(umap)

#### Loading Overdispersion matrix File
Var <- read.csv("./data/processed/Var_b.csv",header = TRUE, sep = ",")
print("bbbbbbbbbbbbbb")
Gnames <- Var[,1]
rownames(Var) <- Gnames
Var$X <- NULL
nm <- colnames(Var)

#### Number of considered genes
ng <- seq(from = 1, to = 168, by = 1)
ng <- c(60,80,120)

#### Variance sorted (Decreasing)
rcmvar <- apply(Var, 1, var)

n <- length(ng)


k.var <- vector("integer",length = n)
k.tsne <- vector("integer",length = n)

label<- gsub("(D6|D19).*","\\1",colnames(Var))
label <- gsub("D6","1",label)
label <- gsub("D19","2",label)
table(label)

for (i in ng){

  sse_Var=matrix(data = 0,nrow = 17,ncol = 1)
  sse_tsne=matrix(data = 0,nrow = 17,ncol = 1)
  
### Only the top ng genes are considered
vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar),as.integer(i*100))]
Var.f <- Var[vi,]
### Correlation Matrix between cells
Var.Corr <-as.data.frame(cor(Var.f,method = "pearson"))
### Distance matrix
Dis.R <- 1-Var.Corr
##### t-SNE computation
set.seed(1)
tSNE.pagoda <- Rtsne(Dis.R,is_distance=TRUE, perplexity=50,max_iter=5000,verbose=0)
plot(tSNE.pagoda$Y, pch=19,col=as.integer(label))
title(paste("tSNE",toString(i*100), sep = " "))

Data.umap = umap(Dis.R,random_state=1)

#label<- grep("D6|D19",colnames(Var.f))
plot(Data.umap$layout, pch=19, col=as.integer(label))
title(paste("uMAP",toString(i*100), sep = " "))
# meth <- c( "ward.D","ward.D2", "single", "complete", "average", 
#            "mcquitty", "median", "centroid","kmeans")
# for (i in meth){
#   print(c(i,"***************************"))
# a <-NbClust(data = tSNE.pagoda$Y, diss = NULL, distance = "euclidean", min.nc = 2, 
#         max.nc = 15, method = i, index = "all", alphaBeale = 0.1)
# }
##### Evaluating Kmeans ######

# #### Finding the best number of groups based on SSE
# for (j in 2:18){
#   ##### Using multispace Overdispersion space 
#     set.seed(0)
#     k=kmeans(Dis.R, j, nstart = 25,iter.max = 100)
#     sse_Var[j-1]= sum(k$withinss)
#   ##### Using bidimensional tSNE space 
#     set.seed(0)
#     k=kmeans(tSNE.pagoda$Y, j, nstart = 25,iter.max = 100)
#     sse_tsne[j-1]=sum(k$withinss)
# }
# #### The function k.groups gets local maximun of the second derivative of SSE
# #### Also, it plots the SSE for different k values
# k.var[i] <- k.groups(sse_Var,'Var_Space',1)
# 
# k.tsne[i] <-k.groups(sse_tsne,'tSNE_Space',1)

}
#write.csv(k.tsne,file="./k_tsn.csv",quote= FALSE, row.names=FALSE, col.names=FALSE)
#write.csv(k.var,file="./k_Var.csv",quote= FALSE, row.names=FALSE, col.names=FALSE)

# 
# ####kmeans using the optinum value
# ## Overdispersion Space
# set.seed(0)
# k1=kmeans(Dis.R, k.var, nstart = 25,iter.max = 100)
# #### PCA space plot
# png(file=paste('./results/plots/Clusters_PCA_Kmeans_Var_',toString(ng),'.png',sep = ""))
#   fviz_cluster(k1, geom = "point",  data = Dis.R)+ggtitle("PCA_Kmeans_Var_8000")
#   dev.off()
# a <- k1$cluster
# G1.col <- a
# nm <- names(a)
# ## Saving classification
# write.csv(a, file = paste("./data/processed/kmeans_Class_Var_",toString(ng),
#                           ".csv",sep = ""),
#             quote = FALSE, row.names = TRUE)
# ###tSNE Space
# set.seed(0)
# k2=kmeans(tSNE.pagoda$Y,k.tsne, nstart = 25,iter.max = 100)
# #### PCA space plot
# png(file=paste('./results/plots/Clusters_PCA_Kmeans_tSNE_',toString(ng),'.png',sep = ""))
#   fviz_cluster(k2, geom = "point",  data = Dis.R)+ggtitle("PCA_Kmeans_tSNE_8000")
#   dev.off()
# a <- k2$cluster
# names(a) <- nm
# ## Saving classification
# write.csv(a, file = paste("./data/processed/kmeans_Class_tSNE_",toString(ng),
#                           ".csv",sep = ""),
#             quote = FALSE, row.names = TRUE)
# 
# G2.col <- a
# #### Ploting clustering into tSNE Space
# MyplotScatter(ng,tSNE.pagoda$Y,G2.col,max(G2.col),
#                 G1.col,max(G2.col))
# #### Computing proportions in each group given the sample day. 
# #### prop[[1]] <- it shows the proportion 
# prop.Var <- Proportion(k1$cluster, k.var, rownames(as.data.frame(k1$cluster)))
# prop.tSNE <- Proportion(k2$cluster, k.tsne, rownames(as.data.frame(k1$cluster)))
# ## Ploting Proportions of cells in each group. 
# MyplotBar(prop.Var,'Kmeans_Var',ng)
# MyplotBar(prop.tSNE,'Kmeans_tSNE',ng)
# ## Computing Dendograms
# # Dendograms must be taken carefully, because their computations requiere
# # to use hierachical clustering which is another clustering algorithm.
# Dendo.Km(k1$withinss,'Kmeans_Var',ng,k.var)
# Dendo.Km(k2$withinss,'Kmeans_tSNE',ng,k.tsne)
