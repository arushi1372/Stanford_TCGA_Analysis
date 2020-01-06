library(data.table)
library(dplyr)
copydata <- fread("CCLE_copynumber_byGene_2013-12-03.txt", data.table=FALSE)

copydata <- data.frame(t(copydata))
setDT(copydata, keep.rownames = TRUE)[]

copydata$rn <- gsub( "_.*$", "", copydata$rn )

dfcell$V1 <- as.character(dfcell$V1)
dfcell[25, ] <- "SYMBOL"

cnvbreast <- merge(dfcell, copydata, by.x='V1', by.y='rn')

cnvbreast2 <- subset(copydata, ,(copydata[2,]='PIK3CA'))
                     
                     
                     , (copydata[2,]='CDKN2A'), rn))

