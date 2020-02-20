rm(list=ls())

hall <- read.table("./h.all.v7.0.symbols.gmt", header = FALSE, sep = "\t",
                   fill = TRUE) 
C2 <- read.table("./c2.all.v7.0.symbols.gmt", header = FALSE, sep = "\t",
                 fill = TRUE)
genes <- c("ANLN","ARHGAP11A","ANP32E")
n.H <- dim(hall)
n.C <- dim(C2)
Genes.list <- list()
for (i in genes){
  datasets <- c()
#  i <- genes[1]
  paths.idx <- as.vector(unlist(lapply(1:n.H[1], function(x) gtlst(hall[x,],i,x))))
  paths <- gtpth(paths.idx, hall[,1:2])
  datasets <- c(datasets,paths)
  
  paths.idx <- as.vector(unlist(lapply(1:n.C[1], function(x) gtlst(C2[x,],i,x))))
  paths <- gtpth(paths.idx, C2[,1:2])
  datasets <- c(datasets,paths)
  
  Genes.list[[i]] <- datasets
}

n.obs <- sapply(Genes.list, length)
seq.max <- seq_len(max(n.obs))
mat <- sapply(Genes.list, "[", i = seq.max)

write.table(mat, file = "./results/data/DEG_and_GSEA/Genes&datasets.csv",
            sep = ";", quote = FALSE, col.names = TRUE, row.names = FALSE)

### get list
gtlst <- function(lst,i,x){
  idx <- which(lst==i)
  lg <- length(idx)
  if (lg >=1){return(x)}
}
### get paths
gtpth <- function(paths.idx, paths){
#  paths <- C2[,1:2]
  datasets <- c()
  for (i in paths.idx){
    i.temp <- i
    H0 <- 0
    while (H0==0) {
      lg <- grepl("^http.",as.vector(paths[i.temp,2]))
      if (lg==FALSE){
        i.temp = i.temp-1
      }else{
        H0 <- 1
      }
    }
    path.name <- as.vector(paths[i.temp,1])
#    print(c(i,path.name))
    datasets <- c(datasets,path.name)
  }
  return(datasets)
}
