###### This Script generates the files to use GSEA software

library(scde)
library(tibble)

## Number of considered genes
ng=12000

fls <- c("kmeans_Class_tSNE",
         "kmeans_Class_Var",
         "ME_Class_tSNE",
         "ME_Class_Var")

fls <- c("kmeans_Class_uMAP",
         "ME_Class_uMAP")

## Number of cluster 
ngrp <- c(3,4,7,8)
###uMAP
ngrp <- c(3,7)
for (j in 1:2){
  
  grp <-ngrp[j]
#  base = paste("_Class_",spc,"_",ng[j],sep="")

  cmp<-paste(LETTERS[1:grp], collapse = "vs")
  ###base="Hallmarks" #### ("Hallmarks","KRH")
#for (i in 1:1){
  ## /Data_KRH_C1vsC2.csv
  a <- paste("./data/processed/Groups_",fls[j],'/',
             fls[j],'_',toString(ng),"_",cmp,".csv",sep ="" )
  MyData <- read.csv(file = a,
                   header = TRUE, row.names = 1, sep = ";")
  cd <- clean.counts(MyData, min.lib.size=1, min.reads = 1
                     , min.detected = 1)
  Gene.names <- rownames(cd)
  cd <- add_column(cd, Description = "na", .before = 1)
  cd <- add_column(cd, NAME = Gene.names, .before = 1)
  col <- colnames(cd)
  cd1 <- rbind(col,cd) ##Adding #NAME  #Description #Samples
  n <- dim(cd)
  nm<-unlist(lapply(3:n[2],FUN=function(x) ""))
  nwe <-c(n[1],n[2]-2,nm) ####Adding #genes #samples
  cd1 <- rbind(nwe,cd1)
  nwe <- c("#1.2","",nm) #### Adding "#1.2"
  cd1 <- rbind(nwe,cd1)
## Saving path for GSEA files
# .GCT file
  path <- paste('./data/processed/Gsea_files_',fls[j],sep="")
  if (dir.exists(path)==FALSE){dir.create(path)}
  a <- paste(path,"/",fls[j],'_',toString(ng),"_",cmp,".gct", sep = "")
  write.table(cd1, file = a, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
#}
# .cls file
  n <- dim(MyData)
  cls <- data.frame(matrix(NA,nrow = 1, ncol = n[2]))

  G1<-unlist(lapply(1:grp,FUN=function(x) paste("Grp_",LETTERS[x],sep="")))
  G2 <- paste(G1, collapse = "|")
  a <- paste("(",G2,").*", sep = "")
  
  ###a<-"(D6|D19).*"
  nm<-unlist(lapply(5:n[2],FUN=function(x) ""))
  sg <- gsub(a, "\\1", colnames(MyData))  ### adding Clases names of every sample
  cls <- rbind(sg,cls)
  nwe <- c("#",G1,nm) ##Clases names
  cls <- rbind(nwe,cls)
  sg <- factor(sg,levels = G1)
  sg1<-table(sg)
  nwe <- c(n[2],length(sg1),1,"",nm)####Adding #samples #Categories #1
  cls <- rbind(nwe,cls)
  cls <- cls[-c(4),]
  a <- paste(path,"/",fls[j],'_',toString(ng),"_",cmp,".cls", sep = "")
  write.table(cls, file = a, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}
