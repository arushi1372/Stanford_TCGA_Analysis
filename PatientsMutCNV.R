cnvchange <- cnvdata

cnvchange$PIK3CA <- ifelse(cnvdata$PIK3CA!=0, 1, 0)
cnvchange$CDKN2A <- ifelse(cnvdata$CDKN2A!=0, 1, 0)
cnvchange$common <- NULL

cnvdata2 <- cnvchange

substring(cnvdata2$Tumor_Sample_Barcode, 1, 15) -> cnvdata2$Tumor_Sample_Barcode
substring(mutationgenes$Tumor_Sample_Barcode, 1, 15) -> mutationgenes$Tumor_Sample_Barcode
intersect(mutationgenes$Tumor_Sample_Barcode, cnvdata2$Tumor_Sample_Barcode)

mutationgenes$PIK3CAMut <- 0
mutationgenes$PIK3CAMut <- ifelse(mutationgenes$Hugo_Symbol=='PIK3CA'& mutationgenes$Variant_Classification!='Silent', 1, 0)
mutationgenes$CDKN2AMut <- 0

mutationgenes$CDKN2AMut <- ifelse(mutationgenes$Hugo_Symbol=='CDKN2A'& mutationgenes$Variant_Classification!='Silent', 1, 0)

mutationgenes$Hugo_Symbol <- NULL
mutationgenes$Variant_Classification <- NULL
mutationgenes$Variant_Type <- NULL
mutationgenes$Matched_Norm_Sample_Barcode<- NULL

m1 <- merge(mutationgenes, cnvdata2, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode", all=TRUE)
m1[is.na(m1)] <- 0

#setwd
load("methdata.RData")
setnames(patientmeth, "X11","CDKN2AMethyl")
setDT(patientmeth, keep.rownames = TRUE)[]
setnames(patientmeth, "rn", "Tumor_Sample_Barcode")
substring(patientmeth$Tumor_Sample_Barcode, 1, 15) -> patientmeth$Tumor_Sample_Barcode
gsub(".", "-", patientmeth$Tumor_Sample_Barcode, fixed = TRUE) -> patientmeth$Tumor_Sample_Barcode


patientmeth2 <- patientmeth
patientmeth2$column_3 <- 0
patientmeth2 <- patientmeth2 %>% mutate(column_3 = Tumor_Sample_Barcode) 
substring(patientmeth2$column_3, 14, 15) -> patientmeth2$column_3
patientmeth2 <-patientmeth2[!(patientmeth2$column_3=="11"),]
patientmeth2$column_3 <- NULL
m3 <- merge(m1, patientmeth2, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode", all=TRUE)
m3[is.na(m3)] <- 0