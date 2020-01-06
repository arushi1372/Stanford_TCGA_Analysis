library(dplyr)
library(data.table)

setwd("~/Desktop/Coursera_R/LUSC/CopyNumberGisticLUSC")
df <- fread("all_data_by_genes.txt", data.table=FALSE)
df2 <- subset(df, df$`Gene Symbol`=='CDKN2A'|df$`Gene Symbol`=='PIK3CA',)
df2 <- data.frame(t(df2))
setnames(df2,'X5264','PIK3CA')
setnames(df2,'X10928','CDKN2A')

df2 <- df2[4:504,]
df2$PIK3CA <- as.numeric(as.character(df2$PIK3CA))
df2$CDKN2A <- as.numeric(as.character(df2$CDKN2A))
df2$PI3K <- 0
df2$p16 <- 0

quantile(df2$PIK3CA,c(0.33,0.67))
df2$PI3K <- ifelse(df2$PIK3CA>1.459, 1,0)
df2$p16 <- ifelse(df2$CDKN2A< -0.372, 1,0)

setDT(df2, keep.rownames = TRUE)[]
setnames(df2, 'rn', 'Tumor_Sample_Barcode')

substring(df2$Tumor_Sample_Barcode, 1, 15) -> df2$Tumor_Sample_Barcode

df2$PIK3CA <- NULL
df2$CDKN2A <- NULL
