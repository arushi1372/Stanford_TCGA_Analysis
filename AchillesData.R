setwd("~/Desktop/Coursera_R/LUSC")
library(data.table)
library(dplyr)
achillesall <- fread("Achilles_All.csv", data.table=F)
achillessample <- fread("Achilles_v2.4_SampleInfo_small.txt", data.table=F)

achilleslung <- subset(achillessample, grepl("LUNG", achillessample$Name), )
achilleslung <- subset(achilleslung, ,select=c("Name","Type", "Subtype","PIK3CA","Mut_source"))

rownames(achillesall) <- achillesall[,1]
achillesall <- subset(achillesall, , grepl("LUNG", colnames(achillesall)))
achillesall$medians <- 0

achillesall <- data.matrix(achillesall, rownames.force=T)
achillesall <- data.frame(achillesall)

for(x in 1:nrow(achillesall)) {achillesall[x,22] <- median(as.numeric(achillesall[x,1:21]))}

#achillesall$shane_median = apply(t(achillesall[,-22]),2,median)

#remove small cell lung cancer
rownames(achilleslung) <- achilleslung[,1]
achilleslung$Name <- NULL
achilleslung <- as.data.frame(t(achilleslung))
achilleslung$NCIH524_LUNG <- NULL
achilleslung$NCIH82_LUNG <- NULL
achilleslung$medians = 'Lung NSCLC'
achillesall <- rbind(achillesall, achilleslung)
achillesnonsmall <- achillesall[,grepl("NSCLC",achillesall[5712,])]
achillesnonsmall <- achillesnonsmall[-(5712:5715),]
for(x in 1:nrow(achillesnonsmall)) {achillesnonsmall[x,20] <- median(as.numeric(achillesall[x,1:19]))}


achillesnonsmall$Gene = rownames(achillesnonsmall)
achillesnonsmall[,!grepl("Gene",colnames(achillesnonsmall))] = sapply(achillesnonsmall[,!grepl("Gene",colnames(achillesnonsmall))],as.numeric)

#do the same with squamous cell
achillessquamous <- achillesall[,grepl("Squamous",achillesall[5713,])]
achillessquamous <- achillessquamous[-(5712:5715),]
achillessquamous$medians <- 0
for(x in 1:nrow(achillessquamous)) {achillessquamous[x,3] <- median(as.numeric(achillessquamous[x,1:2]))}

#gene expression correlations
setwd("~/Desktop/Coursera_R/LUSC/rnaseq_genes_normalized")
rnadata <- fread("LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", data.table=FALSE)
substring(colnames(rnadata), 1, 16) -> colnames(rnadata)
rownames(rnadata) <- rnadata$`Hybridization RE`
rnadata <- rnadata[,grepl("01A",colnames(rnadata))]
substring(colnames(rnadata), 1, 12) -> colnames(rnadata)
rnadata <- data.frame(t(rnadata))
bothgroups <- bothgroups[-3, ]
bothgroups <- bothgroups[-261, ]
rownames(bothgroups) <- bothgroups$Tumor_Sample_Barcode
rnadata$gene_id <- NULL
setDT(rnadata, keep.rownames = TRUE)[]
rnadata$rn <- as.character(rnadata$rn)
bothgroups$Tumor_Sample_Barcode <- as.character(bothgroups$Tumor_Sample_Barcode)
substring(bothgroups$Tumor_Sample_Barcode, 1, 12) -> bothgroups$Tumor_Sample_Barcode
rnaanalysis = merge(rnadata, bothgroups, by.x='rn', by.y='Tumor_Sample_Barcode')

rownames(rnaanalysis) <- rnaanalysis$rn
rnaanalysis$rn <- NULL

rnaanalysis <- lapply(rnaanalysis[1:496,], function(x) as.numeric(as.character(x)))
rnaanalysis <- data.frame(rnaanalysis)
rnaanalysis <- rnaanalysis[complete.cases(rnaanalysis),]
rnaanalysis <- t(rnaanalysis)
rnaanalysis <- data.frame(rnaanalysis)
colnames(rnaanalysis) <- rnaanalysis[20532,]
rnaanalysis <- rnaanalysis[-20532, ]
design3 <- cbind(Grp1=1, Grp2vs1=as.numeric(grepl('1', colnames(rnaanalysis))))
fit2 <- lmFit(rnaanalysis,design3)
fit2 <- eBayes(fit2)
rnalimma = topTable(fit2,coef=2,number = nrow(rnaanalysis))

logdata <- log2(rnaanalysis + 1)
design3 <- cbind(Grp1=1, Grp2vs1=
                   as.numeric(colnames(logdata))
                 # round(runif(ncol(logdata)))
                 )
fit2 <- lmFit(logdata,design3)
fit2 <- eBayes(fit2)
loglimma = topTable(fit2,coef=2,number = nrow(logdata))

qqplot(loglimma$P.Value,seq(0,1,by = 1/nrow(loglimma)))

#aggregate for duplicate ones

achillessquamous$medians <- NULL
achillessquamous <- data.frame(t(achillessquamous))
achillessquamous <- lapply(achillessquamous[-3,], function(x) as.numeric(as.character(x)))
achillessquamous <- data.frame(achillessquamous)
achillessquamous[3, ] <- colnames(achillessquamous)
achillessquamous[3,] <- gsub("\\_.*","",achillessquamous[3,])
achillessquamous <- data.frame(t(achillessquamous))
namesofrows <- rownames(achillessquamous)
achillessquamous <- lapply(achillessquamous[,-3], function(x) as.numeric(as.character(x)))
achillessquamous <- data.frame(achillessquamous)
rownames(achillessquamous) <- namesofrows
achillessquamous$Gene <- rownames(achillessquamous)
achillessquamous$Gene <- gsub("\\_.*","",achillessquamous$Gene)
achillessquamous3 <- aggregate(achillessquamous$X1, list(Gene = achillessquamous$Gene),mean)
achillessquamous4 <- aggregate(achillessquamous$X2, list(Gene = achillessquamous$Gene),mean)
achillessquamous <- merge(achillessquamous3, achillessquamous4, by.x='Gene', by.y='Gene')
rm(achillessquamous3, achillessquamous4)
achillessquamous$median <- 0
for(x in 1:nrow(achillessquamous)) {achillessquamous[x,4] <- median(as.numeric(achillessquamous[x,2:3]))}

achillesnonsmall$medians <- NULL
namesofrows <- rownames(achillesnonsmall)
achillesnonsmall <- lapply(achillesnonsmall[,-20], function(x) as.numeric(as.character(x)))
achillesnonsmall <- data.frame(achillesnonsmall)
rownames(achillesnonsmall) <- namesofrows
achillesnonsmall$Gene <- rownames(achillesnonsmall)
achillesnonsmall$Gene <- gsub("\\_.*","",achillesnonsmall$Gene)
achillesnonsmall2 <- data.frame(achillesnonsmall$Gene)
achillesnonsmall2 <- achillesnonsmall2[!duplicated(achillesnonsmall2$achillesnonsmall.Gene),]
achillesnonsmall2 <- data.frame(achillesnonsmall2)
for (i in 1:19) {
  achillesnonsmall3 <- aggregate(as.numeric(as.character(achillesnonsmall[,i])), list(Gene = achillesnonsmall$Gene),mean)
  achillesnonsmall2 <- merge(achillesnonsmall2, achillesnonsmall3, by.x='achillesnonsmall2', by.y='Gene')
}
rm(achillesnonsmall3)
statsachilles <- achillesnonsmall2
achillesnonsmall2$median <- 0
for(x in 1:nrow(achillesnonsmall2)) {achillesnonsmall2[x,21] <- median(as.numeric(achillesnonsmall2[x,2:20]))}

#achilles and rna

setDT(loglimma, keep.rownames = TRUE)[]
loglimma$rn <- gsub("\\..*","",loglimma$rn)
achillesnonsmall2 <- subset(achillesnonsmall2, ,select=c('median', 'achillesnonsmall2'))
achillesandrna <- merge(achillesnonsmall2, loglimma, by.x=colnames(achillesnonsmall2)[2], by.y=colnames(loglimma)[1])

achillessquamous <- subset(achillessquamous, ,select=c('Gene', 'median'))
squamousandrna <- merge(achillessquamous, loglimma, by.x=colnames(achillessquamous)[1], by.y=colnames(loglimma)[1])

#write.csv(achillesandrna, file = "AchillesNonSmallComputed.csv", fileEncoding = "macroman")

achillesandrna$sum <- 0
for (i in 1:nrow(achillesandrna)) {
  achillesandrna[i,9] <- ((as.numeric(as.character(achillesandrna[i,3]))) - (as.numeric(as.character(achillesandrna[i,2]))))
}

achillesandrna2 <- subset(achillesandrna, achillesandrna$median < 0, )

#scale
achillesandrna2$median <- scale(achillesandrna2$median)
achillesandrna2$logFC <- scale(achillesandrna2$logFC)
achillesandrna2$sum <- achillesandrna2$logFC - achillesandrna2$median

#AAR
plot(achillesandrna2$median, achillesandrna2$logFC, main="logFC and Medians", 
     ylab="logFC ", xlab="median ", ylim = c(0,2), xlim=c(-6,-3))
install.packages("calibrate")
library("calibrate", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
textxy(achillesandrna2$median, achillesandrna2$logFC, achillesandrna2$achillesnonsmall2)

squamousandrna$sum <- 0
squamousandrna$sum <- squamousandrna$logFC - squamousandrna$median
plot(squamousandrna$median, squamousandrna$logFC, main="logFC and Medians: squamous cell lines", 
     ylab="logFC ", xlab="median ", ylim = c(0,3), xlim=c(-3,0),col=ifelse(squamousandrna$Gene=='PRDX1', "red", "black"),
     pch=ifelse(squamousandrna$Gene=='PRDX1', 19, 1), cex=ifelse(squamousandrna$Gene=='PRDX1', 2, 1))
textxy(squamousandrna$median, squamousandrna$logFC, squamousandrna$Gene)

dotchart(squamousandrna$median, main="Medians of Achilles Squamous Lines", xlab="median ", xlim=c(-3,0), col=ifelse(squamousandrna$Gene=='PRDX1', "red", "black"),pch=ifelse(squamousandrna$Gene=='PRDX1', 19, 1))
dotchart(squamousandrna$logFC, main="LogFC Values of RNA", xlab="logFC ", xlim=c(0,4), col=ifelse(squamousandrna$Gene=='PRDX1', "red", "black"),pch=ifelse(squamousandrna$Gene=='PRDX1', 19, 1))
dotchart(clinandlimma$cor, main="Correlation Values of Clinical Data", xlab="cor", col=ifelse(clinandlimma$row=='PRDX1.PRDX1', "red", "black"),pch=ifelse(clinandlimma$row=='PRDX1.PRDX1', 19, 1))
dotchart(clinandlimma$logFC, main="LogFC Values of Proteins", xlab="logFC", col=ifelse(clinandlimma$row=='PRDX1.PRDX1', "red", "black"),pch=ifelse(clinandlimma$row=='PRDX1.PRDX1', 19, 1))
dotchart(achillesandrna$median, main="Medians of Achilles Non-small Cell Lines", xlab="medians", col=ifelse(achillesandrna$achillesnonsmall2=='PRDX1', "red", "black"),pch=ifelse(achillesandrna$achillesnonsmall2=='PRDX1', 19, 1))
dotchart(tvalues$tvalues, main="p-values for T-tests on Achilles Non-small Data", xlab="p-values", xlim=c(0,0.05), col=ifelse(tvalues$V2=='PRDX1', "red", "black"),pch=ifelse(tvalues$V2=='PRDX1', 19, 1))

#statistical tests achilles
statsachilles <- t(statsachilles)
colnames(statsachilles) <- statsachilles[1,]
statsachilles <- data.frame(statsachilles)
statsachilles <- statsachilles[-1,]
statsachilles <- lapply(statsachilles, function(x) as.numeric(as.character(x)))
statsachilles <- data.frame(statsachilles)
tvalues <- 0
tvalues <- data_frame(tvalues)
names = list(colnames(statsachilles))
names <- unlist(names)
for(i in 1:ncol(statsachilles)) {
x = t.test(statsachilles[,i], mu=0, alternative = "less")$p.value
tvalues[i,1] <- x
tvalues[i,2] <- names[i]
}

mean(squamousandrna$median>=-0.79099367)
mean(achillesandrna2$median>=-5.1244666)
mean(achillesandrna2$logF<=1.303121)
mean(clinandlimma$logFC<=0.02465076)
         
#clinandlimma2 <- subset(clinandlimma, clinandlimma$cor <0, )
#rppa7 <- data.frame(t(rppa6))
#setDT(rppa7, keep.rownames = TRUE)[]
#rppa7 <- subset(rppa7, rppa7=rppa7[c(1,2,146,226), ], )
