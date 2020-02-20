#### This script find the differentially expresed through Karchenko et al. 
#### 2014 framework. Therefore, SCDE R package is used. Differential gene Expresion 
#### analysis is done by pairs considering all posible combinations. 
#### If there is three Cluster (A,B,C) for example,
#### three analises are performed and so on: AvsB, AvsC and BvsC.

### Clean memory
rm(list = ls())
library(scde)
library(gtools)
#####################

### Number of considered genes
val <- c(6000)

#### Evaluated Space
fls.name <- c("kmeans_Class_uMAP",
         "ME_Class_uMAP",
         "kmeans_Class_Var",
         "ME_Class_Var")

id.val <- 1
fl <- 0
ng <- length(val)
##### This loop is usefull when its considered more than one genes number 
for (j in 1:(ng*4)){  
  id.fls <- (j-1)%/%ng+1
###  id.fls <- j
   if (id.fls==fl){
     id.val <-  id.val + 1
  }else{
     id.val <-  1
     fl <- fl + 1}

  print(paste("/Groups_",fls.name[id.fls],"_",
        toString(val[id.val]),".csv",sep ="" ))
### Number of Cluster under that specific space and clustering method
    Clust <- read.csv(paste("./data/processed/",fls.name[id.fls],
                          "_",toString(val[id.val]),".csv",sep ="" ),
                    sep = ",", row.names = 1)

  N.clust = max(Clust[,1])
  ### Possible combinations of 2 elements given the cluster numbers 
  comb <- combinations(N.clust, 2, LETTERS[1:N.clust],
                       set=TRUE, repeats.allowed=FALSE)
  
 ### This loop specifically runs the DEA for every combitation 
  for(i in 1:choose(N.clust,2)){
    #### file name
    f.name <- paste("./data/processed/Groups_",fls.name[id.fls],
                    "/",fls.name[id.fls],"_",
                    toString(val[id.val]),"_",
                    paste(comb[i,], collapse = "vs"),
                    ".csv",sep = "")
    #### Load file 
    MyData <- read.csv(file = f.name,
                        header = TRUE, row.names = 1, sep = ";")
    
    ################scde
    ## factor determining cell types
    G1 <- paste("Grp_",comb[i,1],sep = "")
    G2 <- paste("Grp_",comb[i,2],sep = "")
    a <- paste("(",G1,"|",G2,").*", sep = "")
    ### Find the proportion of every compared cluster
    sg <- factor(gsub(a, "\\1", colnames(MyData)), levels = c(G1,G2))
    names(sg) <- colnames(MyData)  
    table(sg)
    ## clean up the dataset, minimun of reads per sample= 1, minimun reads
    ## per gene=1, mininum genes per sample = 1.
    cd <- clean.counts(MyData, min.lib.size=1, min.reads = 1, 
                       min.detected = 1)
    ###Fitting error model
    ## calculate models
    #### CAUTION #####
    #### The limiting factor isn't the numbers of cores, It is the RAM
    #### For this datasets, when used 2 cores the RAM demmanded is on
    #### average 28 GB.
    
    o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 2, 
                               threshold.segmentation = TRUE, 
                               save.crossfit.plots = FALSE, 
                               save.model.plots = FALSE, 
                               verbose = 1)
    head(o.ifm)
    
    ## filter out cells that don't show positive correlation with
    ## the expected expression magnitudes (very poor fits)
    valid.cells <- o.ifm$corr.a > 0
    table(valid.cells)
    o.ifm <- o.ifm[valid.cells, ]
    ## estimate gene expression prior
    o.prior <- scde.expression.prior(models = o.ifm, counts = cd, 
                                     length.out = 400, show.plot = FALSE)
    
    ## define two groups of cells
    groups <- factor(gsub(a, "\\1", rownames(o.ifm)), levels  =  c(G1,G2))
    names(groups) <- row.names(o.ifm)
    
    ## Run differential expression tests on all genes.
    ediff <- scde.expression.difference(o.ifm, cd, o.prior, 
                                        groups  =  groups, 
                                        n.randomizations  =  100, 
                                        n.cores  =  5, verbose  =  1)
    
    ## top upregulated genes (tail would show top downregulated ones)
#    head(ediff[order(ediff$Z, decreasing  =  TRUE), ])
    
#    tail(ediff[order(ediff$Z, decreasing  =  TRUE), ])
    ## table with all the results, showing most significantly different genes 
    ## (in both directions) on top
    base <- paste(fls.name[id.fls],"_",
                  toString(val[id.val]),sep = "")
    a <- paste("./results/data/Diff_Gen_Results_",base,"(",G1,"vs",G2,").txt",sep="")
    write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = a, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
    a <- paste("./results/data/Results_",base,"(",G1,"vs",G2,").RData",sep="")
    save.image(file = a)
    ## run the differential expression on a single gene, and visualize the results
    ## scde.test.gene.expression.difference("TPX2", models = o.ifm, counts = cd, prior = o.prior)
  }
}
