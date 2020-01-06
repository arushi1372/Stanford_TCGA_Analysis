library(data.table)
library(dplyr)

clin <- fread("LUSC.clin.merged.txt", data.table=FALSE)
clin <- clin[c(13, 17), ]
clin <- data.frame(t(clin))

rnaseq <- fread("LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", data.table=F)
rnaseq <- subset(rnaseq, grepl("GSK3", rnaseq$`Hybridization REF`), )
rnaseq <- data.frame(t(rnaseq))
setDT(rnaseq, keep.rownames=TRUE)[]

substring(rnaseq$rn, 6, 12) -> rnaseq$rn
substring(clin$X13, 6, 12) -> clin$X13

rnaseq <- rnaseq[!duplicated(rnaseq$rn),]
toupper(clin$X13) -> clin$X13

m1 <- merge(clin, rnaseq, by.x='X13', by.y='rn', all =TRUE)

plot(as.numeric(as.character(m1$X7485))  ~ as.numeric(as.character(m1$X17)))

setwd("~/Desktop/Coursera_R/LUSC/CopyNumberGisticLUSC")
copynum <- fread("all_data_by_genes.txt", data.table=F)

copynum <- subset(copynum, grepl("GSK3", copynum$`Gene Symbol`), )
copynum <- data.frame(t(copynum))
setDT(copynum, keep.rownames=TRUE)[]

substring(copynum$rn, 6, 12) -> copynum$rn
copynum <- copynum[!duplicated(copynum$rn),]
m2 <- merge(clin, copynum, by.x='X13', by.y='rn', all=TRUE)
plot(as.numeric(as.character(m2$X4871))  ~ as.numeric(as.character(m2$X17)), ylim = c(-2, 2))

cnv <- fread("all_thresholded.by_genes.txt", data.table=FALSE)
cnv <- subset(cnv, grepl("GSK3", cnv$`Gene Symbol`), )
cnv <- data.frame(t(cnv))
setDT(cnv, keep.rownames=TRUE)[]

substring(cnv$rn, 6, 12) -> cnv$rn
cnv <- cnv[!duplicated(cnv$rn),]
m3 <- merge(clin, cnv, by.x='X13', by.y='rn', all=TRUE)
plot(as.numeric(as.character(m3$X4871))  ~ as.numeric(as.character(m3$X17)), ylim = c(-2, 2))


