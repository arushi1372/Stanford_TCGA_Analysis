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
n_occur <- data.frame(table(mutgenes$Tumor_Sample_Barcode))

#the patient IDs in which a mutation in PIK3CA and a mutation in CDKN2A occurred
doublemut <- subset(n_occur, Freq==2, select=Var1)