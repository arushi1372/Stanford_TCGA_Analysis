#open CNV file for all data with a threshold of 0.1 used
#0 represents no loss of duplicate, + and - represent amp. and del.

cnvalldatathreshold <- fread("all_thresholded.by_genes.txt", data.table=FALSE)

#filter out genes
cnvdata <- cnvalldatathreshold %>% filter(`Gene Symbol`=='CDKN2A'|`Gene Symbol`== 'PIK3CA')

#transpose data, create new column, and set appropriate names
cnvdata <- data.frame(t(cnvdata))
cnvdata$X3 <- 0
setnames(cnvdata, 'X1', 'PIK3CA')
setnames(cnvdata, 'X2', 'CDKN2A')
setnames(cnvdata, 'X3', 'common')

#label in the third column, the patients who had a copy number
#change in both the genes
cnvdata$common <- ifelse(cnvdata$PIK3CA!=0 & cnvdata$CDKN2A!=0, 1, 0)

#separate out only wanted data
cnvdata <- cnvdata[4:504,]
