#set wd
setwd("~/Desktop/Coursera_R/LUSC/ClinicalRaw")
clin <- fread("LUSC.clin.merged.txt", data.table=FALSE)
clin <- clin[c(13, 17), ]
toupper(clin[1,]) -> clin[1,]
colnames(clin) <- clin[1,]
clin <- t(clin)
colnames(clin) <- clin[1,]
clin <- clin[-1,]
clin <- data.frame(clin)
clin$patient.days_to_death <- as.numeric(as.character(clin$patient.days_to_death))

rppa2 <- data.frame(t(rppa2))
substring(rppa2$Tumor_Sample_Barcode, 1, 12) -> rppa2$Tumor_Sample_Barcode
rppa6 = merge(rppa2, clin, by.x=colnames(rppa2)[1], by.y= colnames(clin)[1])
rppa6 = na.omit(rppa6)
rppa6[ ,2:225] <- as.data.frame(lapply(rppa6[,2:225],function(x) as.numeric(as.character(x))))

#correlate protein expression with survival in patient group
proteinsforclin <- subset(rppa6, grepl("1", rppa6$Group), )
plot(proteinsforclin$SQSTM1.P62.LCK.LIGAND, proteinsforclin$patient.days_to_death)

cor.test(proteinsforclin$AR.AR, proteinsforclin$patient.days_to_death)

#make a table of findings
library(Hmisc)
res <- rcorr(as.matrix(proteinsforclin[,-(1:2)]))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res2 <- flattenCorrMatrix(res$r, res$P)
res2 <- subset(res2, res2$column=="patient.days_to_death", )
res2$column <- NULL

#crossreference limma and clinical
setDT(out.table, keep.rownames = TRUE)[]
clinandlimma <- merge(res2, out.table, by.x='row', by.y='rn', all=TRUE)

crossreference <- subset(clinandlimma, clinandlimma$cor<0 & clinandlimma$logFC>0 & clinandlimma$P.Value<0.1, )

#correlate phosphoratios
patientsforphospho <- rppa[1,]
colnames(patientsforphospho) <- colnames(phosphoratios)
phosphoratios <- rbind(phosphoratios, patientsforphospho)
phosphoratios[28, ] <- colnames(phosphoratios)
colnames(phosphoratios) <- phosphoratios2[27, ]
substring(colnames(phosphoratios), 1, 12) -> colnames(phosphoratios)
phosphoratios <- data.frame(t(phosphoratios))
phosphoratios <- phosphoratios[-164, ]
rownames(phosphoratios) <- phosphoratios$Tumor_Sample_Barcode
substring(phosphoratios$Tumor_Sample_Barcode, 1, 12) -> phosphoratios$Tumor_Sample_Barcode
clinandphospho = merge(phosphoratios, clin, by.x=colnames(phosphoratios)[27], by.y= colnames(clin)[1])
clinandphospho <- clinandphospho[grepl("X",clinandphospho$X28),]

res <- rcorr(as.matrix(clinandphospho[,-c(1,28)]))
res3 <- flattenCorrMatrix(res$r, res$P)
res3 <- subset(res3, res3$column=="patient.days_to_death", )
res3$column <- NULL

#merge limma and cor for phosphoratios
phosphotable$rn <- gsub("\\ ", ".", phosphotable$rn)
mergephospho <- merge(res3, phosphotable, by.x=colnames(res3)[1], by.y=colnames(phosphotable)[1])
