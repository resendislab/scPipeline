###setwd("~/Documents/MCF7 Single Cell")
library(gtools)
#library(VennDiagram)

rm(list = ls())

join <- function(dirGsea,fls,nm,tag){
  n <- length(fls)
  DaTa <- data.frame()
  while (n>0) {
    ### read files
    DaTa.A <- read.table(file = paste(dirGsea,fls[n],"/",
                                      nm,"_",tag[n],".xls",sep = ""),
                         sep = "\t", header = TRUE)
    ### Selec pathways with FDR<0.5
    DaTa.A <- DaTa.A[DaTa.A[,8]<=0.05,]
    ### Selec pathways with p<0.01
    DaTa.A <- DaTa.A[DaTa.A[,7]<=0.01,]
    DaTa <-rbind(DaTa,DaTa.A)
    n <- n-1
  }
  ### Valid if DaTa is not empty
  if (nrow(DaTa)!=0){
    ### Name re-structuration to be more readable
    DaTa$NAME <- gsub("_", " ", DaTa$NAME)
    ### Ordering Data
    DaTa <- with(DaTa,  DaTa[order(DaTa$NAME) , ])
    ### Erasing non informative columns
    DaTa[,c(2:3,10:14)] <-NULL
  }
  return(DaTa)
}


### Path to the folder that stores GSEA enriched results
dirGsea=("~/gsea_home/output/Project_SC/")

### Clustering medthod
mth=c("Kmeans","ME")
#mth="ME"
### Dimentional reduction method, Var goes for a non method used, 
### just the variance multispacial matrix
spc=c("tSNE","Var","uMAP")

### Number of used genes
ng=c("6000","8000","12000")
pathE <- "./results/data/Enrichment"
if (dir.exists(pathE)==FALSE){dir.create(pathE)}
####j <- space
####k <- number of genes 
####m <- 

### Number of cluster for each method over the specific data space
N.Clust <- read.table(file = "./results/data/N_Clusters.csv", 
                      sep = "\t", header = TRUE)
rownames(N.Clust) <- N.Clust[,1]
N.Clust <- N.Clust[-1]

### Number of total analysis
N.a <- length(mth)*length(spc)*length(ng)
fl<- 1
i.spc <- 0

for (i in 1:N.a){
  i.mth <- (i-1)%/%(length(spc)*length(ng))+1
  i.ng <- (i-(i.mth-1)*length(ng)*length(spc)-1)%/%(length(ng))+1
  ###  id.fls <- j
  if (i.spc == length(spc)){
    i.spc <-  0
  }
  i.spc <-  i.spc +1
#  print(c(mth[i.mth],spc[i.spc], ng[i.ng]))
  
  mylist <- list()
  n <- N.Clust[mth[i.mth],paste(spc[i.spc], ng[i.ng],sep = "_")]
  print(paste(mth[i.mth],spc[i.spc],ng[i.ng],n,"Groups",sep = "_"))

#    comb <- combinations(n, 2, LETTERS[1:n],
#                         set=TRUE, repeats.allowed=FALSE)
#    nc=choose(n,2)
  ### create results Folder for the enrichment
  pathT <- paste(pathE,"/",mth[i.mth],"_Class_",spc[i.spc],"_",ng[i.ng]
                 ,sep = "")

  if (dir.exists(pathT)==FALSE){dir.create(pathT)}
  
  for (p in 1:n){
      ### Folder name wheer data are stored
      fdr= paste(mth[i.mth],"_Class_",spc[i.spc],"_",ng[i.ng],"_",
                 LETTERS[p],sep = "")
      fls <- list.files(path=dirGsea,pattern =paste("^",fdr,
                                                    "vsREST",".",sep = "" ))

      ### Extract the number of analysis imputed by GSEA
      tag <- unlist(lapply(1:length(fls),FUN=function(x) 
        gsub(paste("^",fdr,".\\D*",sep = "" ),"",fls[x])))
      
      nm <- paste("gsea_report_for_Grp_",LETTERS[p],sep = "")
      
      ### Join all results form differents databases
      DaTa <- join(dirGsea,fls,nm,tag)

      write.table(DaTa,file = paste(pathT,"/",fdr,".csv",sep = ""), 
                  row.names = FALSE, col.names = TRUE, sep = ";",
                  quote = FALSE)

  }
}

