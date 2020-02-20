#### Script that computes the p-value for the differentially expresed genes
#### genes with p-value<0.05 and a fold change>2.0 are saved.

#### Clean memory
rm(list = ls())

library(gtools)
##################### Comparación por clustering

#### Number of selected genes
val <- 6000#c(6000,8000,12000)

#### Space used to perform clustering
# fls <- c("./kmeans",
#          "./kmeans",
#          "./ME",
#          "./ME")
fls <- c("./uMAP,
         ./uMAP")
# fls.name <- c("kmeans_Class_tSNE",
#               "kmeans_Class_Var",
#               "ME_Class_tSNE",
#               "ME_Class_Var")
fls.name <- c("kmeans_Class_uMAP",
              "ME_Class_uMAP")
#### Saving results on the following folder
path <- ("./results/data/Diff_genes_IPA")
if (dir.exists(path)==FALSE){dir.create(path)}
fl<- 1
id.val <- 0
ng <- length(val)
for (j in 1:(ng*2)){   ###1:12  
  id.fls <- (j-1)%/%ng+1
  ###  id.fls <- j
  if (id.fls==fl){
    id.val <-  id.val + 1
  }else{
    id.val <-  1
    fl <- fl + 1}

  print(paste("/Groups_",fls.name[id.fls],"_",
              toString(val[id.val]),".csv",sep ="" ))
  Clust <- read.csv(paste("./data/processed/",fls.name[id.fls],
                          "_",toString(val[id.val]),".csv",sep ="" ),
                    sep = ",")
  N.clust = max(Clust[,2])
  comb <- combinations(N.clust, 2, LETTERS[1:N.clust],
                       set=TRUE, repeats.allowed=FALSE)

  for(i in 1:choose(N.clust,2)){
#    f.name <- paste(fls[id.fls],"/Groups/",fls.name[id.fls],
#                    toString(val[id.val]),"_",
#                    paste(comb[i,], collapse = "vs"),
#                    ".csv",sep = "")
    ################scde
    # factor determining cell types
    G1 <- paste("Grp_",comb[i,1],sep = "")
    G2 <- paste("Grp_",comb[i,2],sep = "")
    base <- paste(fls.name[id.fls],"_",
                  toString(val[id.val]),sep = "")
    
    a <- paste("./results/data/Results_",base,
               "(",G1,"vs",G2,").RData",sep="")
    load(a)
    
    # 2-tailed p-value
    p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
    p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
#    significant.genes <- which(p.values.adj<0.05)
    ord <- order(p.values.adj) # order by p-value
    de <- cbind(ediff[2],p.values,p.values.adj)[ord,]
    colnames(de) <- c("log2 fold change","p-value","Adjusted p-value")
    ############ Upregulated Genes
#    UpG <- de[de$`log2 fold change`>=1.0,]
#    a <- paste(path,"/Upregulated_",
#               base,"(",G1,"vs",G2,")_",G1,"_genes_(logFC_1.0)",
#               ".txt",sep="")
#    write.table(UpG, file = a, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
    ############ Downregulated genes
#    DwG <- de[de$`log2 fold change`<=-1.0,]
#    a <- paste(path,"/Upregulated_",
#               base,"(",G1,"vs",G2,")_",G2,"_genes_(logFC_1.0)",
#               ".txt",sep="")
#    write.table(DwG, file = a, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
    ########### Differentially expresed gene for both comparisons
#    DG <- rbind(DwG,UpG)
    a <- paste(path,"/DE_",
               base,"(",G1,"vs",G2,")",
               ".txt",sep="")
    print(a)
    write.table(de, file = a, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
    if (i>=2){
      def <- merge(def,de, by="row.names")
      rownames(def) <- def[,1]
      def[,1] <- NULL
    }
    else{
      def <- de
      }
  }
  a <- paste(path,"/DE_",
             base,".txt",sep="")
  write.table(def, file = a, row.names = TRUE, col.names = TRUE, 
              sep = "\t", quote = FALSE)
}