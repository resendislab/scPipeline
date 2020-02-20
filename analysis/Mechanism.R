
rm(list = ls())

Reg_search <- function(a,Ups){ 
  df <- data.frame(Regulator=as.character(),
                   Tar.Num=as.integer())
  #for (i in Gns){
  #  a <- Gns.lst
  b <- lapply(Ups$Target.Molecules.in.Dataset, function(x)
    a[match(unlist(strsplit(as.vector(x),",")),a, nomatch = 0)])
  b.1 <- unlist(lapply(b, function(x) length(x)))
  
  #  b <- (lapply(Ups$Target.Molecules.in.Dataset, function(x)
  #    sum(match(unlist(strsplit(as.vector(x),",")),a, nomatch = 0)>0)))
  c <- which(b.1!=0)
  nms <- b[c]
  nms.1 <- unlist(lapply(nms, function(x) 
    paste(as.vector(x),collapse = ",")))
  df <- data.frame(Regulator=as.vector(Ups$Upstream.Regulator[c]),
                   Tar.Num=b.1[c],Regu.Genes=nms.1)
  df <- df[order(-df$Tar.Num),]
}

pth <- "./IPA/Upstream.txt"

Ups <- read.table(pth,header = TRUE, sep = "\t")
Ups <- Ups[,-9]

idx <- which(Ups$Molecule.Type == "chemical drug")
Ups <- Ups[-idx,]
idx <- which(Ups$Molecule.Type == "chemical reagent")
Ups <- Ups[-idx,]


Gns <- list("Gns_Met" = c("CXCL10", "S100A7", "S100A8", "S100A9"),
            "Gns_Prol" = c("MKI67","TOP2A","FOXM1"),
            "Gns_OXPHOS" = c("MT-CO1","MT-CO2","MT-CO3",
                             "MT-RNR1", "MT-ND4L","MT-ND5","MT-ND2",
                             "TMEM64"),
            "Gns_BvsC" = c("CASP4","H2BC21","RNF170"),
            "Gns_Erk" = c("ANLN","ANP32","ARHGAP11A"))
a <- unlist(Gns)
df <- Reg_search(a,Ups)
write.table(df,"./IPA/Regulators&targets.txt", quote = FALSE,
            sep = ";", col.names = TRUE, row.names = FALSE)

Regulators <- as.vector(df$Regulator)
df.Reg <- Reg_search(Regulators,Ups)
write.table(df.Reg,"./IPA/Regulators&Regulators.txt", quote = FALSE,
            sep = ";", col.names = TRUE, row.names = FALSE)


  
#}
