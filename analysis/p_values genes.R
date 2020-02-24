#### Script that computes the p-value for the differentially expressed genes
#### genes analysis. Genes with p-value<0.05 and a fold change>2.0 are saved.
#### Also, gets overexpressed genes for every cluster given the actual
#### comparison.

#### Clean memory
rm(list = ls())
library(gtools)
############## Custom Functions
Gns.Fil <- function(comb,base,logFC){
 
  #### factor determining the clusters comparison
  G1 <- paste("Grp_",comb[1],sep = "")
  G2 <- paste("Grp_",comb[2],sep = "")
 
  #### Open saved space from SCDE analysis
  a <- paste("./results/data/Results_",base,
             "(",G1,"vs",G2,").RData",sep="")
  load(a)
 
  # 2-tailed p-value
  p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
  ##### Genes with a p-value < 0.05
  significant.genes <- which(p.values.adj<0.05)
  ord <- order(p.values.adj[significant.genes]) # order by p-value
  de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
  colnames(de) <- c("Lower bound","log2 fold change","Upper bound","p-value")
  ############ Upregulated Genes with a Log Fold Change of logFC
  UpG <- de[de$`log2 fold change`>=logFC,]
  a <- paste(path,"/Upregulated_",
             base,"(",G1,"vs",G2,")_",G1,"_genes_(logFC_",
             toString(logFC),").txt",sep="")
  write.table(UpG, file = a, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
  ############ Downregulated genes
  DwG <- de[de$`log2 fold change`<=-logFC,]
  a <- paste(path,"/Upregulated_",
             base,"(",G1,"vs",G2,")_",G2,"_genes_(logFC_",
             toString(logFC),").txt",sep="")
  write.table(DwG, file = a, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
  ########### Differentially expressed genes for both comparisons
  DG <- rbind(DwG,UpG)
  a <- paste(path,"/regulated_",
             base,"(",G1,"vs",G2,")_(logFC_",
             toString(logFC),").txt",sep="")
  write.table(DG, file = a, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
}
############################

#### Number of selected genes
val <- 6000
#### Log Fold Change CutOff
logFC <- 2
#### Space used to perform clustering
fls.name <- c("kmeans_Class_uMAP",
              "kmeans_Class_Var",
              "ME_Class_uMAP",
              "ME_Class_Var")
#### Saving results on the following folder
path <- ("./results/data/Diff_genes")
if (dir.exists(path)==FALSE){dir.create(path)}
fl<- 1
id.val <- 0
ng <- length(val)
## This loop perform computation for more than one gene number
for (j in 1:(ng*4)){   ###1:12  
  id.fls <- (j-1)%/%ng+1
  ###  id.fls <- j
  if (id.fls==fl){
    id.val <-  id.val + 1
  }else{
    id.val <-  1
    fl <- fl + 1}

##### Open classification
  Clust <- read.csv(paste("./data/processed/",fls.name[id.fls],
                          "_",toString(val[id.val]),".csv",sep ="" ),
                    sep = ",")
#### Cluster number   
  N.clust = max(Clust[,2])
### Possible combinations of two elements
  comb <- combinations(N.clust, 2, LETTERS[1:N.clust],
                       set=TRUE, repeats.allowed=FALSE)
  ### name of root analysis
  base <- paste(fls.name[id.fls],"_",
                toString(val[id.val]),sep = "")
#### Main structure to compute filtering
  invisible(lapply(1:choose(N.clust,2), function(x) Gns.Fil(comb[x,],
         base, logFC)))
}
