#### Script to select and filter the  percentage 
#### of valid cells based on the count matrix

#rm(list = ls())
#pth <- "./Documents/ScPipe/6B_samples/Output"

Filter.Count.Matrix <- function(smp,b){
  #### smp <- sample name
  pth <- "./data/raw/ddSeeker_Output"
  library(splus2R)
  library(sjmisc)
  ### load raw count matrix
  GnsL1 <- read.table(gzfile(paste(pth,"/",smp,"_L1/",smp,
                                   "_L1.dge.txt.gz", sep = "")),
                      header = TRUE, quote ="", row.names=NULL) 
  GnsL1 <- Duplicated.Genes(GnsL1)
  
  GnsL2 <- read.table(gzfile(paste(pth,"/",smp,"_L2/",smp,
                                   "_L2.dge.txt.gz", sep = "")),
                      header = TRUE, quote ="", row.names=NULL) 
  GnsL2 <-Duplicated.Genes(GnsL2)
  
  GnsL3 <- read.table(gzfile(paste(pth,"/",smp,"_L3/",smp,
                                   "_L3.dge.txt.gz", sep = "")),
                      header = TRUE, quote ="", row.names=NULL) 
  GnsL3 <- Duplicated.Genes(GnsL3)
  
  GnsL4 <- read.table(gzfile(paste(pth,"/",smp,"_L4/",smp,
                                   "_L4.dge.txt.gz", sep = "")),
                      header = TRUE, quote ="", row.names=NULL) 
  GnsL4 <- Duplicated.Genes(GnsL4)
  
  DF<-list(GnsL1,GnsL2,GnsL3,GnsL4)
  COL<-unique(unlist(lapply(DF, colnames)))
  ROW<-unique(unlist(lapply(DF, rownames)))
  TOTAL<-matrix(data=0, nrow=length(ROW), ncol=length(COL), dimnames=list(ROW, COL))
  # Subsetting :
  for (df in DF) { 
    c <- matrix(data=0, nrow=length(ROW), ncol=length(COL), dimnames=list(ROW, COL))
    c[rownames(df), colnames(df)] <- as.matrix(df)
    #  print(c)
    TOTAL <- TOTAL+c
    #print( c)
  }
  Gns <- as.data.frame(TOTAL)
  write.table(Gns, file=paste("./data/raw/Count_Matrix_",
                              smp,".csv",sep = ""),
              quote= FALSE, row.names = TRUE, col.names = TRUE, sep = ",")
  rm(GnsL1,GnsL2,GnsL3,GnsL4)
  
  FCs <- inflec(Gns,smp,"1")
}

inflec <- function(Gns,smp,tit){
  Sm <- unlist(lapply(Gns,sum))
  Sm1 <- order(Sm,decreasing = TRUE)
  Sm <- Sm[Sm1]
  CSm <- cumsum(Sm)
  #####UMI count Plot
  png(file=paste('./results/plots/Sum_plot_',smp,'_',tit,'.png',sep = " "))
  mar.default <- c(5,4,4,2)
  par(cex.axis=1.8,mfrow=c(1, 1),oma=c(0, 0, 0, 0),
      mar = mar.default + c(0, 1, 0, 0))
  plot(Sm, log='xy', 
       xlab = "Bead barcode in descending order by genic UMI count", 
       ylab = "Genic UMI Count")
  dev.off()
  
  ### knee Plot
  library(inflection)
#  a <- check_curve(seq(1,length(CSm),1),CSm)
  b <- ede(seq(1,length(CSm),1),CSm, 1)
  b <- b[1]
  png(file=paste('./results/plots/Knee_plot_',smp,'_',tit,'.png',sep = " "))
  mar.default <- c(5,4,4,2)
  par(cex.axis=1.8,mfrow=c(1, 1),oma=c(0, 0, 0, 0),
      mar = mar.default + c(0, 1, 0, 0))
  plot(CSm, 
       ylab = "Cumulative sum of genic UMI Count", 
       xlab = "Bead barcode in descending order by genic UMI count")
  abline(v=b, col="red", lwd=3, lty=1)
  abline(h=CSm[b], col="red", lwd=3, lty=1)
  points(b,CSm[b], col="blue",bg='blue',pch=21)
  dev.off()
  # ###
  # 
  Gns <- Gns[,Sm1]
  # ## Cells passing knee filter
  FCs <- Gns[,1:b]
#  FCs <- Gns[,1:b[1,1]]
#  rm(Gns,TOTAL,DF,c,df)
  # 
  
  write.table(FCs, file=paste("./data/raw/Filtered_Count_Matrix_",
                              smp,".csv",sep = ""),
              quote= FALSE, row.names = TRUE, col.names = TRUE, sep = ",")
  return(FCs)
}

Duplicated.Genes <- function(Gns){
#  Gns <- MyData
  nms <- upperCase(Gns[,1])
  n <- dim(Gns)
  a <- nms[duplicated(nms)]
  for (i in a){
    idx <- which(nms == i)
    if (is_empty(idx)==FALSE){
      Gns[idx[1],2:n[2]] <- colSums(Gns[idx,2:n[2]])
      idx <- idx[-1]
      Gns <- Gns[-idx,]
      nms <- nms[-idx]
    }
  }
  rownames(Gns) <- upperCase(Gns[,1])
  Gns <- Gns[,-1]
  return(Gns)
}

