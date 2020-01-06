load("methdata.RData")
setnames(patientmeth, "X11","CDKN2AMethyl")
setDT(patientmeth, keep.rownames = TRUE)[]
setnames(patientmeth, "rn", "Tumor_Sample_Barcode")
substring(patientmeth$Tumor_Sample_Barcode, 1, 15) -> patientmeth$Tumor_Sample_Barcode
gsub(".", "-", patientmeth$Tumor_Sample_Barcode, fixed = TRUE) -> patientmeth$Tumor_Sample_Barcode
m2 <- merge(m1, patientmeth, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode", all=TRUE)
m2[is.na(m2)] <- 0