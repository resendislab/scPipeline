###### This script computes the gene Over-dispersion under the PAGODA
###### pipeline by Fan et al. 2016
###### NOTE: 
###### NOTE0: Till the date 6/01/2019 version of the flemix package must be 
###### <= 2.3-13 to be compatible with scde pakage.
###### NOTE1:Data des not need logarithmic transformation
rm(list=ls())

library(scde)

##### Loading Data
MyData <- read.csv(file = './data/processed/MCTS_b.csv',
                   header = TRUE, row.names = 1, sep = ";")

##### Cleaning data, minimun of genes per sample=1000, minimun of 
##### count per gene in all samples=1, Minimum number of cells a 
##### gene must be seen in=1

cd <- clean.counts(MyData, min.lib.size=1000, min.reads = 1, min.detected = 1)

##### 
##### It is recommended to use n-2 cores, n is the total number of cores 
##### in used computer, see NOTE0.
knn <- knn.error.models(cd, k = ncol(cd)/3, n.cores = 6, min.count.threshold = 2, 
                        min.nonfailed = 5, max.model.plots = 10)
##### Normalizing variance
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, 
                          n.cores = 6, plot = TRUE)

##### Controlling for sequencing depth
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

##### Over-dispersion matrix
Var <- as.data.frame(varinfo$mat)

##### Saving Over-dispersion matrix
write.csv(Var, file="./data/processed/Var_b.csv")

