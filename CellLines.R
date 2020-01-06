library(data.table)

df <- fread("CCLE_cellline_mapping.txt", data.table=F)
df2 = fread("cell_line_mut_ann.maf", data.table=F)

df2.sub.genes = df2[grepl("PIK3CA|CDKN2A",df2$Hugo_Symbol),]
df2.sub.genes = df2.sub.genes[grepl("_BREAST",df2.sub.genes$Tumor_Sample_Barcode),]
df2.sub.genes = subset(df2.sub.genes, , select = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification'))

df = df[grepl("_LUNG",df$`CCLE name`),]
df = subset(df, , select = c('CCLE name', 'Hist Subtype1'))

dfboth <- merge(df, df2.sub.genes, by.x='CCLE name', by.y='Tumor_Sample_Barcode', all=T)

dfcell <- read.csv("LUSC Protein Analysis - CDK1 Cell Lines.csv", header = F)
dfcell <- dfcell[,-1]

dfboth$`CCLE name` <- gsub( "_.*$", "", dfboth$`CCLE name` )

dfboth[70:168, 1] <- substring((dfboth[70:168, 1]), 4)

dff <- merge(dfcell, dfboth, by.y='CCLE name', by.x='CCLE')

setdiff(dfcell$CCLE, dff$CCLE)
