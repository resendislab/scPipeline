###### This scripts computes kmeans cluster algorithm given two spaces as inputs:
######    1. Overdispersion multidimensional space
######    2. umap bidimensional space
###### For every space the optinum number of groups is computed by calculating
###### the SSE and taking the second derivative to get the most significative 
###### change in the SSE reduction
###### Given the optinum number of groups, the asignation of every sample in
###### each group is done. Also, proportions are ploted im two bar plots.
rm(list=ls())
source("./analysis/Funn.R")
library(factoextra)
library(NbClust)
library(umap)

#### Loading Overdispersion matrix File
Var <- read.csv("./data/processed/Var_b.csv",header = TRUE,
                row.names=1, sep = ",")
#### Number of considered genes
ng <- 6000

#### getting the labels based on time sample, 1 goes for D6 samples and
#### 2 for D19 sample
label<- gsub("(D6|D19).*","\\1",colnames(Var))
label <- gsub("D6","1",label)
label <- gsub("D19","2",label)
table(label)
#### vectors to save SSE given kmeans algorithm
sse_Var=matrix(data = 0,nrow = 17,ncol = 1)
sse_umap=matrix(data = 0,nrow = 17,ncol = 1)
  
#### Variance sorted (Decreasing)
rcmvar <- apply(Var, 1, var)
### Only the top ng genes are considered
vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar),ng)]
Var.f <- Var[vi,]
### Correlation Matrix between cells
Var.Corr <-as.data.frame(cor(Var.f,method = "pearson"))
### Distance matrix
Dis.R <- 1-Var.Corr
#  ##### t-SNE computation for comparations purposes
  #library(Rtsne)
  # set.seed(1)
  # tSNE <- Rtsne(Dis.R,is_distance=TRUE, perplexity=50,max_iter=5000,verbose=0)
  # plot(tSNE$Y, pch=19,col=as.integer(label))
  # title(paste("tSNE",toString(ng), sep = " "))
  
### uMAP computation
Data.umap = umap(Dis.R,random_state=1)
  
# uMAP space of data, colored by sample time
plot(Data.umap$layout, pch=19, col=as.integer(label))
title(paste("uMAP",toString(ng), sep = " "))
legend("bottomright",
    c('D6','D19'),
    pch = 21,
    pt.bg =c(1,2),
    bty = "n",
    cex = 1.5
)
##### Evaluating Kmeans ######
  #### Finding the best number of groups based on SSE
  for (j in 2:18){
    ##### Using multispace Overdispersion space
      set.seed(0)
      k=kmeans(Dis.R, j, nstart = 25,iter.max = 100)
      sse_Var[j-1]= sum(k$withinss)
    ##### Using bidimensional tSNE space
      set.seed(0)
      k=kmeans(Data.umap$layout, j, nstart = 25,iter.max = 100)
      sse_umap[j-1]=sum(k$withinss)
  }
#### The function k.groups gets local maximun of the second derivative of SSE
#### Also, it plots the SSE for different k values
k.var <- k.groups(sse_Var,'Var_Space',ng)
  
k.umap <-k.groups(sse_umap,'uMAP_Space',ng)

###uMAP Space with the optimal group number
  ### Overdispersion Space
      set.seed(0)
      k1=kmeans(Dis.R, k.var, nstart = 1,iter.max = 100)
      G1.col <- k1$cluster
      nm <- names(G1.col)
      ## Saving classification
      write.csv(G1.col, file = paste("./data/processed/kmeans_Class_Var_",toString(ng),
                                ".csv",sep = ""),
                quote = FALSE, row.names = TRUE)
  ### Using uMAP as input
      set.seed(0)
      k2=kmeans(Data.umap$layout,k.umap, nstart = 1,iter.max = 100)
      G2.col <- k2$cluster
      ### Saving classification
      write.csv(G2.col, file = paste("./data/processed/kmeans_Class_uMAP_",toString(ng[i]*100),
                               ".csv",sep = ""),
                 quote = FALSE, row.names = TRUE)
#### Ploting clustering results into uMAP Space
MyplotScatter(ng,Data.umap$layout,G1.col,max(G1.col),
              G2.col,max(G2.col),"kmeans")
      
### Computing the proportion of every group
prop.Var <- Proportion(k1$cluster, k.var, 
                       rownames(as.data.frame(k1$cluster)))
      
prop.uMAP <- Proportion(k2$cluster, k.umap, 
                        rownames(as.data.frame(k2$cluster)))


### Ploting Proportions of cells in each group.
### Sample percentage plot tells you the proportion of the total sample
### contained in every Cluster giveng the time
### Groups percentage plot tells you the composition of every cluster 
### giveng  the Cluster sample
MyplotBar(prop.Var,'Kmeans_Var',ng)
MyplotBar(prop.uMAP,'Kmeans_uMAP',ng)
