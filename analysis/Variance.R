###### This scripts computes the Following:
###### 1. By the Variance, the genes are sorted in order to take the most
######    overdispersed genes, in this study the number of chosen genes were
######    6000, we wiped out genes with variance less than
######    less than 16% of total range. Therefore, we only consider 53% of
######    the total variance.
rm(list=ls())

source("./analysis/Funn.R")
library(factoextra)
library(Rtsne)

#### Loading Overdispersion matrix File
Var <- read.csv("./data/processed/Var_b.csv",
                header = TRUE, sep = ",")
Gnames <- Var[,1]
rownames(Var) <- Gnames
Var$X <- NULL
nm <- colnames(Var)

#### Cut off Number of considered genes
ng <- 6000

#### Variance sorted (Decreasing)
rcmvar <- apply(Var, 1, var)
# Variance ordering
vi <- order(rcmvar, decreasing = TRUE)
Var.order <- rcmvar[vi]
# Law of total Variance
TotVar <- sum(Var.order)
PropVar <- lapply(1:length(Var.order),FUN= function(x) Var.order[x]/TotVar)
PropVar1 <- unlist(lapply(1:length(Var.order),
                   FUN = function(x) 100*(max(Var.order[1:x])-min(Var.order[1:x]))/max(Var.order)))
# Cumulative sum of variance as number of considered genes increased
CumVar <- cumsum(PropVar)*100
# Plotting Variance Cumulative Sum
png(file='./results/plots/Variance_Cumulative_Sum.png')
  par(cex.axis=1.8)
  plot(CumVar, xlab="Number of considered genes",
       ylab="Cumulative Variance",
       cex.lab=1.8,ylim=c(0,100), lty=1, lwd=1)
  abline(v=ng,lty=2,lwd=2 ,col="red")
  abline(h=CumVar[ng],lwd=2,lty=2,col="red")
  dev.off()
# Plotting Variance Range proportion given the genes number
png(file='./results/plots/Variance_Range_1.png')
  par(cex.axis=1.8)
  plot(100-PropVar1, xlab="Number of considered genes",
       ylab="Variance range",
       cex.lab=1.8,ylim=c(0,100), lty=1, lwd=1)
  abline(v=ng,lty=2,lwd=2 ,col="red")
  abline(h=100-PropVar1[ng],lwd=2,lty=2,col="red")
  dev.off()
#####