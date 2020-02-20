### This script computes the Venn diagram for the differentially 
### expresed genes

library(VennDiagram)
library(gtools)
rm(list=ls())
pth <- ("./results/data/Diff_genes")
#setwd("~/Documents/MCF7 Single Cell")

####
### Number of cluster for each method over the specific data space
N.Clust <- read.table(file = "./results/data/N_Clusters.csv", 
                      sep = "\t", header = TRUE)
rownames(N.Clust) <- N.Clust[,1]
N.Clust <- N.Clust[-1]

### Clustering medthod
mth=c("Kmeans")
### Dimentional reduction method 
spc=c("uMAP")
### number of genes
ng ="6000" 

### number of clusters
n <- N.Clust[mth,paste(spc, ng,sep = "_")]
### Combinations by pairs
comb <- combinations(n, 2, LETTERS[1:n],
                     set=TRUE, repeats.allowed=FALSE)
### number of comprations
n.comb <- choose(n,2)

leg <- unlist(lapply(1:n.comb, FUN=function(x) paste("Grp_",comb[x,1],"vs",
                                           "Grp_",comb[x,2],sep = "")))
cutOff = "2.0"
#for (j in 3:3){#length(ng)){
  
base = paste("Class_",spc,"_",ng,sep="")
  
### files corresponding to the 
fls <- list.files(path=dirGsea,pattern = paste("^",fdr,".",sep = "" ))
Data <- list()

for (i in 1:n.comb){
    
    # NvsR
    a <- paste(pth,"/regulated_",tolower(mth),"_",
               base,"(",leg[i],")",".txt", sep = "")
    deg <- read.table(file=a, header = TRUE, sep = "\t")
    Data[[leg[i]]] <- as.vector(rownames(deg))
}
  # NvsV
  #b <- paste("./",mth,"/Diff_Gen/Diff_Genes_SCDE_",base,"(",leg[2],")_CutOff_",cutOff,".txt", sep = "")
  #Data.b <- read.table(file=b,header = TRUE, col.names = TRUE)
  # RvsV
  #c <- paste("./",mth,"/Diff_Gen/Diff_Genes_SCDE_",base,"(",leg[3],")_CutOff_",cutOff,".txt", sep = "")
  #Data.c <- read.table(file=c,header = TRUE, col.names = TRUE)
D1 <- intersect(Data[[1]],Data[[2]])  #### A'
D2 <- intersect(Data[[1]],Data[[3]])  #### B'
D3 <- intersect(Data[[2]],Data[[3]])  #### C'
n5 <- length(intersect(D1,intersect(D2,D3)))  #### A(inter)B(inter)C
n1 <- length(D1)-n5 #### A
n2 <- length(D2)-n5 #### B
n3 <- length(D3)-n5 #### C
n4 <- (length(Data[[1]]) - length(intersect(Data[[1]],D1)) 
       - length(intersect(Data[[1]],D2)) - n5) #### A(inter) B
n6 <- (length(Data[[3]]) - length(intersect(Data[[3]],D1)) 
       - length(intersect(Data[[3]],D3)) - n5) #### A(inter) C
n7 <- (length(Data[[2]]) - length(intersect(Data[[2]],D1)) 
       - length(intersect(Data[[2]],D2)) - n5) #### B(inter) C

#############Venn Diagram
  
venn.plot <- draw.triple.venn(
    area1 = n1+n4+n5+n6,
    area2 = n2+n4+n5+n7,
    area3 = n3+n5+n6+n7,
    
    n12 = n4+n5,
    n23 = n7+n5,
    n13 = n6+n5,
    n123 = n5,
    
    category = c("Cluster A", "Cluster B", "Cluster C"),
    sep.dist = 0.1,
    fill = c("darkgrey", "tomato", "springgreen4"),
    lty = "blank",
    cex = 2,
    cat.cex = 2.,
    rotation.degree = 0,
    cat.pos = c(-5,8,178),
    set.dist = 0.5)
  fl=paste('./',mth,'/Diff_Gen/Venn_comp_',base,'.tiff',sep = "")
  tiff(filename = fl, compression = "lzw");
  grid.draw(venn.plot)
  dev.off()
  