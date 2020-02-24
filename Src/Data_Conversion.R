###### This script takes the Data output files and convert them into a data
###### frame with genes as rows and samples as columns. It is used for
###### posterior analysis. The number 6 in labels marks for day sixth samples;
###### label 19 goes for day 19th samples.
rm(list=ls())

###### Loading raw Data
source("./Src/Filter_Count_Matrix.R")

Names_conv <- function(DaTa,lb){
  n <- dim(DaTa)
  names<-unlist(lapply(1:n[2],FUN=function(x) paste(lb,x,sep="_")))
  colnames(DaTa) <- names
  return(DaTa)
}


### load raw count matrix
Data.6B <- Filter.Count.Matrix("6B")

Data.19A <-Filter.Count.Matrix("19A")
Data.19B <-Filter.Count.Matrix("19B")

### Data Correction for 19A and 19B samples
Data.19A <- inflec(Data.19A,"19A","2")
Data.19B <- inflec(Data.19B,"19B","2")

library(scde)

###### Merging and transposing data
Data.6 <- as.data.frame(Data.6B)
rownames(Data.6) <- upperCase(rownames(Data.6))
n <- dim(Data.6)
names6<-unlist(lapply(1:n[2],FUN=function(x) paste("D6",x,sep="_")))
colnames(Data.6) <- names6

Gene.names <- c(rownames(Data.19A),rownames(Data.19B))
Data.19 <- merge(Data.19B,Data.19A, by = "row.names",
                 all = TRUE, sort = FALSE)
### Matriz count filled with Zeros for non matching cases
Data.19[is.na(Data.19)] <- 0

Gene.names <- Data.19$Row.names

Data.19 <- Data.19[,-1]
n <- dim(Data.19)
names19<-unlist(lapply(1:n[2],FUN=function(x) paste("D19",x,sep="_")))
colnames(Data.19) <- names19
rownames(Data.19) <- Gene.names

rm(Data.19B,Data.19A,Data.6B)

MyData <- merge(Data.6,Data.19, by = "row.names" , all = TRUE, sort = FALSE)
MyData[is.na(MyData)] <- 0
Gene.names <- MyData$Row.names
####MyData <- MyData[,-1]

##### Name conversion from hg19 to hg38
##### loading Names Function
source("./src/Names_function.R")

MyData[,1] <- Names.hg19_hg38(as.vector(Gene.names))

#### duplicated genes with equivalent names
MyData1 <- Duplicated.Genes(MyData)

##### Saving Data as csv File
##### MCTS goes for MultiCellular Tumor Spheroids
write.csv2(MyData1, file="./data/processed/MCTS_b.csv",
           quote = FALSE)
rm(Data.19,Data.6,names19,names6,n)


