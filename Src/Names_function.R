### This script gets the names and convert them from hg19 to hg38

Names.hg19_hg38 <- function(gns.names){
## gns.names <- as.vector(Gene.names)
  library(splus2R)
#  ROW <- c(upperCase(rownames(FCs)),upperCase(rownames(Ilu_6B)))
  #write.table(ROW, file="./Documents/MCF7 Single Cell/ddSeeker_pipeline_Count_Matrix/ROW_6B.csv",
  #            quote= FALSE, row.names = TRUE, col.names = TRUE, sep = ",")
  ROW <- upperCase(gns.names)
  ROW <- gsub("*[.]AS","-AS",ROW)
  ROW <- gsub("*[.]D","-D",ROW)
  ROW <- noquote(unique(ROW))
  # 
  ############# Genes WITh LOC name
  library(mygene)
  LOC.gns <-ROW[grep("^LOC", ROW)]
  LOC.ENTREZ <- gsub("LOC","",LOC.gns)
  LOC.Query <- queryMany(as.numeric(LOC.ENTREZ),scopes="entrezgene",
                         fields="entrezgene, symbol", species="human")
  LOC.GNS <- matrix(data = NA, nrow=length(LOC.ENTREZ), ncol=3)
  colnames(LOC.GNS) <- c("Symbol_Hg19","Entrez_ID","Symbol_Hg38")
  LOC.GNS[,1] <- LOC.gns
  LOC.GNS[,2] <- as.numeric(LOC.ENTREZ)
  LOC.GNS[,3] <- LOC.Query$symbol
  #### Not found LOC
  idx <- which(LOC.Query$notfound == TRUE)
  LOC.GNS[idx,3] <- LOC.gns[idx]
  
  #### Genes without LOC*named genes
  ROW <- ROW[! ROW %in% LOC.gns]
  
  GNS.query<-queryMany(ROW,scopes="symbol",
                       fields="entrezgene, symbol", species="human")
  aux <- GNS.query$notfound
  aux[is.na(GNS.query$notfound)] <- FALSE
  ## Not Found
  NF <- GNS.query$query[aux]
  ## Found genes
  aux <- which(GNS.query$entrezgene != "NA")
  FN <- matrix(data=NA, nrow = length(aux), ncol = 3)
  colnames(FN) <- c("Symbol_Hg19","Entrez_ID","Symbol_Hg38")
  FN[,1] <- GNS.query$query[aux]
  FN[,2] <- GNS.query$entrezgene[aux]
  FN[,3] <- GNS.query$query[aux]
  
  
  #### NOT Found and depreciated
  aux1 <- c(FN[,1],NF, ROW)
  NF1 <- aux1[ave(seq_along(aux1), aux1, FUN = length) == 1]
  NF <- c(NF,NF1)
  
  ##### Using the HGNC list, the not found genes were converted to hg38
  ##### https://www.genenames.org/download/custom/
  Nms.Con <- read.csv('~/Documents/MCF7 Single Cell/ddSeeker_pipeline_Count_Matrix/HGNC.csv',
                      header = TRUE ,quote= "", sep = ";")
  
  Genes.NF <- matrix(data=NA, nrow = length(NF), ncol = 3)
  colnames(Genes.NF) <- c("Symbol_Hg19","Entrez_ID","Symbol_Hg38")
  Genes.NF[,1] <- NF 
  
  Prv.names <- strsplit(as.vector(Nms.Con$Previous.symbols), ", ")
  Sin.names <- strsplit(as.vector(Nms.Con$Synonyms), ", ")
  
  for (i in 1:length(NF)){
    Prv <- sapply(Prv.names, FUN=function(x) NF[i] %in% x)
    idx <- which(Prv == TRUE)
    if (!isEmpty(idx==TRUE)){
      idx <- which(Prv==TRUE) 
      if (length(Nms.Con$Approved.symbol[idx]>1)){
        Genes.NF[i,3] <- as.character(Nms.Con$Approved.symbol[idx[1]])
      }else {
        Genes.NF[i,3] <- as.character(Nms.Con$Approved.symbol[idx])
      }
    } else{
      Sin <- sapply(Sin.names, FUN=function(x) NF[i] %in% x)
      idx <- which(Sin == TRUE)
      if (!isEmpty(idx==TRUE)){
        if (length(Nms.Con$Approved.symbol[idx]>1)){
          Genes.NF[i,3] <- as.character(Nms.Con$Approved.symbol[idx[1]])
        }else {
          Genes.NF[i,3] <- as.character(Nms.Con$Approved.symbol[idx])
        }
      }
    }
  }
  Genes.NF[is.na(Genes.NF[,3]),3] <- Genes.NF[is.na(Genes.NF[,3]),1] 
  
  ##### Merging LOC, not found and found genes. 
  Genes <- rbind(LOC.GNS, Genes.NF, FN)
  Genes <- Genes[!duplicated(Genes[,1]),]
  
  
  
  #### Save Names convertion
  
#  write.table(Genes,'./Documents/MCF7 Single Cell/ddSeeker_pipeline_Count_Matrix/hg19_hg38_genes.csv',
#              quote= FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  
  ### Converts names in the original list
  
#  Illu.hg38 <- Ilu_6B
  
  nms <- (gns.names)
  nms <- gsub("*[.]AS","-AS",nms)
  nms <- gsub("*[.]D","-D",nms)
  #rownames(Illu.hg38) <- nms
  idxs <- unlist(lapply(nms,FUN = function(x) which(x == Genes[,1])))
  NMS <- Genes[idxs,3]
#  idx <- NMS[duplicated(NMS)]
#  idxs <- unlist(lapply(idx,FUN = function(x) which(x == Genes[,3])))
#  Dnms <- Genes[idxs,1]
#  idxs <- unlist(lapply(Dnms,FUN = function(x) which(x == nms)))
#  NMS[idxs] <- Dnms
  return(NMS)
}