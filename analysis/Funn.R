
linePlot <- function(...){  ### x, y, tit, xlab, ylab, ng
  x <- list(...)  
  png(file=paste('./results/plots/',x[3],'.png',sep = ""))
  par(cex.axis=1.8)
  print(x[2])
  plot(x[1],x[2])
  plot(x[1],x[2], xlab=x[4],ylab=x[5],
       cex.lab=1.8,ylim=c(0,100), lty=1, lwd=1)
  if (length(x)==6){
    abline(v=x[6],lty=2,lwd=2 ,col="red") 
    abline(h=x[2][x[6]],lwd=2,lty=2,col="red") 
  }
  dev.off()
}

k.groups <- function(sse,tit,i){
 sse1 <- sse
 png(file=paste('./results/plots/SSE_',tit,'.png',sep = " "))
 mar.default <- c(5,4,4,2)
 par(cex.axis=1.8,mfrow=c(1, 1),oma=c(0, 0, 0, 0),
     mar = mar.default + c(0, 1, 0, 0))
 matplot(2:18, sse1, type="o",lty=1,pch = 19, lwd = 2.0,
         xlab = "Groups number (k)",
         ylab = "SSE",cex.lab=1.8)
 legend("topright",
        as.character(i),
        pch = 21,
        pt.bg = c(1:length(i)),
        bty = "n",
        pt.cex = 1.5,
        cex = 2
 )
 dev.off()
  ##### Slope of SSE
  m <-diff(sse)
  dm <- diff(m)
  maxi<-unlist(lapply(1:length(i),FUN=function(x) order(dm[,x],decreasing = TRUE)[1:3]))
  k <- c(3:18) 
  kmx <- k[maxi]
  kmx <- matrix(kmx, nrow=3,ncol=length(i))
  return(kmx[1,])
}

MyplotScatter <- function(i, Y, G1.col,G1.gr,G2.col,G2.gr,tit){
  png(file=paste('./results/plots/Clusters_',tit,
                 '_',toString(i),'.png',sep = "")
      ,width = 1200, height = 700, 
      units = "px", pointsize = 24)
  par(mfrow=c(1, 2), mar=c(4.1, 2.1, 1.5, 0.5),
      oma=c(0, 0, 0, 0))
  plot(Y,main=paste('Var_Space_',toString(i),sep = "")
       ,pch=c(21,24),bg = G1.col, col=1)
  legend("bottomright",
         c(LETTERS[1:G1.gr],'D6','D19'),
         pch = c(rep(22,G1.gr),19,17),
         pt.bg = c(1:G1.gr,1,1),
         bty = "n",
         pt.cex = 1,
         cex = 1
  )
  plot(Y,main=paste('uMAP_Space_',toString(i),sep = ""),pch=c(21,24),
       bg = G2.col, col=1)
  legend("bottomright",
         c(LETTERS[1:G2.gr],'D6','D19'),
         pch = c(rep(22,G2.gr),19,17),
         pt.bg = c(1:G2.gr,1,1),
         bty = "n",
         pt.cex = 1,
         cex = 1
  )
  dev.off()
}

Proportion <- function(clust,n,col){
  #  col <- colnames(col.cols)
  #  nm <- as.data.frame(clust)
  #  col <- rownames(nm)
  tot=c(134,230)
  nam=c("D6_","D19_")
  prop <- vector(length = n*2)
  prop1 <- vector(length = n*2)
  for(i in 0:(n*2-1)){
    id <- i%/%2+1
    id1 <- i%/%(2*id-1)+1
    b=which(clust %in% id)
    prop[i+1]<-length(grep(nam[id1],col[b],value=TRUE))/tot[id1]*100
    if (id1 == 2){
      prop1[i:(i+1)]=100*prop[i:(i+1)]/sum(prop[i:(i+1)])
    } 
  }
  prop <-array(prop,dim=c(2,n))
  prop1 <-array(prop1,dim=c(2,n))
  prop <- as.table(prop)
  prop1 <- as.table(prop1)
  rownames(prop) <- c('D6','D19')
  #  colnames(prop) <- c(toString(1:n))
  
  rownames(prop1) <- c('D6','D19')
  #  colnames(prop1) <- c(toString(1:n))
  return(list(prop,prop1))
}

MyplotBar <- function(prop,tit,i){
  png(file=paste('./results/plots/Proportions_',tit,'_',toString(i),'.png',
                 sep = " "),width = 1200, height = 700,
      units = "px", pointsize = 24)
  par(mfrow=c(1, 2), mar=c(4.1, 2.1, 1.5, 0.5),
      oma=c(0, 2, 0, 0))
  barplot(prop[[1]], ylim = c(0,100), col=c("darkblue","red"),
          legend = rownames(prop[[1]]),main = 'Sample percentage')
  barplot(prop[[2]], ylim = c(0,100),col=c("darkblue","red"),
          legend = rownames(prop[[2]]),main = "Groups percentage")
  dev.off()
}