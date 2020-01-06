#which cell lines aren't present in this data
df2.sub = df2[grepl("_BREAST",df2$Tumor_Sample_Barcode),]

df2.sub$Tumor_Sample_Barcode <- gsub( "_.*$", "", df2.sub$Tumor_Sample_Barcode)
df2.sub2 = subset(df2.sub, , select = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification'))

df2.sub2 <- df2.sub2[order(df2.sub2$Tumor_Sample_Barcode),]

setDT(df2.sub2, keep.rownames = TRUE)[]
df2.sub2$rn <- NULL

df2.sub2$Tumor_Sample_Barcode <- as.character(df2.sub2$Tumor_Sample_Barcode)

intersect(dfcell$V1, df2.sub2$Tumor_Sample_Barcode)

dfcell$Line <- paste("NCI", dfcell$Line, sep="")
intersect(dfcell$Line, df2.sub2$Tumor_Sample_Barcode)
