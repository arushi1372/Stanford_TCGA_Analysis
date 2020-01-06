library(data.table)
library(dplyr)

mutdata <- fread("BRCA-TP.final_analysis_set.maf", data.table = F)
mutationgenes <- mutdata %>% filter(Hugo_Symbol=='CDKN2A'|Hugo_Symbol == 'PIK3CA')
mutationgenes <- subset(mutationgenes, Hugo_Symbol =='CDKN2A'|Hugo_Symbol== 'PIK3CA', select=c(Hugo_Symbol, Variant_Classification, Variant_Type, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode))
mutationgenes <- subset(mutationgenes, Variant_Classification!='Silent', select=Hugo_Symbol:Matched_Norm_Sample_Barcode)

mutationgenes <- subset(mutationgenes, mutationgenes$Tumor_Sample_Barcode == 'TCGA-LQ-A4E4-01A-11D-A25Q-09', )
mutationgenes <- subset(mutationgenes, , select=c(Tumor_Sample_Barcode))
mutationgenes <- data.frame(mutationgenes[1, ])
substring(mutationgenes, 1, 15) -> mutationgenes

#------------

cnvraw <- fread("all_data_by_genes.txt", data.table=F)
cnvraw <- subset(cnvraw, cnvraw$`Gene Symbol`=='CDKN2A'|cnvraw$`Gene Symbol`=='PIK3CA',)
cnvraw <- data.frame(t(cnvraw))
setnames(cnvraw,'X5264','PIK3CA')
setnames(cnvraw,'X10928','CDKN2A')

cnvraw <- cnvraw[4:504,]
cnvraw$PIK3CA <- as.numeric(as.character(cnvraw$PIK3CA))
cnvraw$CDKN2A <- as.numeric(as.character(cnvraw$CDKN2A))

quantile(cnvraw$PIK3CA,c(1-0.0867,0.0867))

cnvraw$PI3K <- 0
cnvraw$p16 <- 0
cnvraw$PI3K <- ifelse(cnvraw$PIK3CA>0.7342, 1,0)
quantile(cnvraw$CDKN2A,c(.75,0.25))
cnvraw$p16 <- ifelse(cnvraw$CDKN2A < (-0.315), 1,0)

setDT(cnvraw, keep.rownames = TRUE)[]
setnames(cnvraw, 'rn', 'Tumor_Sample_Barcode')

substring(cnvraw$Tumor_Sample_Barcode, 1, 15) -> cnvraw$Tumor_Sample_Barcode

cnvraw$PIK3CA <- NULL
cnvraw$CDKN2A <- NULL

#--------
intersect(mutationgenes, cnvraw$Tumor_Sample_Barcode)


