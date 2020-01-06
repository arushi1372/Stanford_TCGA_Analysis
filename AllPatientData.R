#adjustedcnv

#set wd ------------------------
setwd("~/Desktop/Coursera_R/LUSC/MutSig2CVLUSC")
#use Mutation Analyses data named MutSig2CV for data
mutationdata <- fread("LUSC-TP.final_analysis_set.maf", data.table =FALSE)

#filter out the data pertaining to CDKN2A and PIK3CA
mutationgenes <- mutationdata %>% filter(Hugo_Symbol=='CDKN2A'|Hugo_Symbol == 'PIK3CA')

#select only the relevant columns
mutationgenes <- subset(mutationgenes, Hugo_Symbol =='CDKN2A'|Hugo_Symbol== 'PIK3CA', select=c(Hugo_Symbol, Variant_Classification, Variant_Type, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode))

#delete any data that had silent mutations which would not have had any effect
mutationgenes <- subset(mutationgenes, Variant_Classification!='Silent', select=Hugo_Symbol:Matched_Norm_Sample_Barcode)

#check for duplicates within each dataset separately
duplCD <- mutationgenes %>% filter(Hugo_Symbol=='CDKN2A')
duplPI <- mutationgenes %>% filter(Hugo_Symbol == 'PIK3CA')
duplicated(duplCD$Tumor_Sample_Barcode)
duplicated(duplPI$Tumor_Sample_Barcode)

#delete the duplicate rows in the PIK3CA dataset
duplPI <- duplPI[!duplicated(duplPI$Tumor_Sample_Barcode),]

#combine the unique data from each set into one dataframe
mutgenes <- rbind(duplCD, duplPI)

#create a table that documents the frequency of a patient ID
#n_occur <- data.frame(table(mutgenes$Tumor_Sample_Barcode))

#the patient IDs in which a mutation in PIK3CA and a mutation in CDKN2A occurred
#doublemut <- subset(n_occur, Freq==2, select=Var1)

#-------------------

substring(df2$Tumor_Sample_Barcode, 1, 15) -> df2$Tumor_Sample_Barcode
substring(mutgenes$Tumor_Sample_Barcode, 1, 15) -> mutgenes$Tumor_Sample_Barcode
intersect(mutgenes$Tumor_Sample_Barcode, df2$Tumor_Sample_Barcode)

mutgenes$PIK3CAMut <- 0
mutgenes$PIK3CAMut <- ifelse(mutgenes$Hugo_Symbol=='PIK3CA'& mutgenes$Variant_Classification!='Silent', 1, 0)
mutgenes$CDKN2AMut <- 0
mutgenes$CDKN2AMut <- ifelse(mutgenes$Hugo_Symbol=='CDKN2A'& mutgenes$Variant_Classification!='Silent', 1, 0)

mutgenes$Hugo_Symbol <- NULL
mutgenes$Variant_Classification <- NULL
mutgenes$Variant_Type <- NULL
mutgenes$Matched_Norm_Sample_Barcode<- NULL

m1 <- merge(mutgenes, df2, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode", all=TRUE)
m1[is.na(m1)] <- 0

#################################stop here
#setwd
setwd("~/Desktop/Coursera_R/LUSC")
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

rm(mutationdata)
rm(duplCD, duplPI, temp, m1,df)
