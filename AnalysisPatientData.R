#analyzing the patient data LUSC

m3t <- m3
m3t$Group <- 0
m3t$Group <- ifelse(((m3t$PI3K==1|m3t$PIK3CAMut==1) & (m3t$p16==1|m3t$CDKN2AMethyl==1|m3t$CDKN2AMut==1)), 1, 0)
patientgroup <- subset(m3t, m3t$Group==1, )

#m3t$Group <- ifelse(((m3t$PIK3CA==1|m3t$PIK3CAMut==1) & (m3t$CDKN2A==1|m3t$CDKN2AMut==1)), 1, 0)
#patientgroup2 <- subset(m3t, m3t$Group==1, )

bothgroups <- subset(m3t, , select=c('Tumor_Sample_Barcode','Group'))

#protein analysis
setwd("~/Desktop/Coursera_R/LUSC/RPPAAnalyses")
rpparaw <- fread("LUSC-TP.normalized.gct", data.table=FALSE)
rpparaw <- data.frame(rpparaw[,-1], row.names=rpparaw[,1])
rpparaw <- data.frame(t(rpparaw))
setDT(rpparaw, keep.rownames = TRUE)[]
gsub(".", "-", rpparaw$rn, fixed = TRUE) -> rpparaw$rn
substring(rpparaw$rn, 1, 15) -> rpparaw$rn
rpparaw <- rpparaw[-1,]
setnames(rpparaw, 'rn', 'Tumor_Sample_Barcode')
rppa <- merge(bothgroups, rpparaw, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode", all=TRUE)
na.omit(rppa) -> rppa
rppa <- data.frame(t(rppa))

library(limma)

rppa2 <- rppa
colnames(rppa2) <- ifelse(rppa2[2,]=='1', 'X','Y')
design <- cbind(Grp1=1, Grp2vs1=as.numeric(grepl('X', colnames(rppa2))))
lapply(rppa2[3:225,], function(x) as.numeric(as.character(x))) -> rppa3
rppa3 <- data.frame(rppa3)
rownames(rppa3) <- rownames(rppa2[3:225,])
fit <- lmFit(rppa3,design)
fit <- eBayes(fit)
out.table = topTable(fit,coef=2,number = nrow(rppa3))


#                         logFC AveExpr     t  P.Value adj.P.Val     B
#PIK3CA..PI3K.P110.ALPHA  0.197   0.512  7.12 6.60e-12  1.47e-09 16.46
#ERRFI1.MIG.6            -0.141   0.510 -6.29 9.83e-10  1.10e-07 11.59
#PARK7.DJ.1              -0.162   0.481 -5.64 3.53e-08  2.62e-06  8.12
#SQSTM1.P62.LCK.LIGAND    0.206   0.503  4.77 2.70e-06  1.51e-04  3.95
#KEAP1.KEAP1              0.148   0.505  4.40 1.48e-05  6.59e-04  2.33
#TFRC.TFRC                0.198   0.504  4.32 2.03e-05  7.56e-04  2.03
#PDK1.PDK1               -0.108   0.506 -4.21 3.27e-05  1.01e-03  1.57
#EIF4E.EIF4E             -0.124   0.500 -4.19 3.61e-05  1.01e-03  1.48
#CDH3.P.CADHERIN         -0.138   0.505 -4.09 5.48e-05  1.32e-03  1.09
#BAX.BAX                 -0.134   0.513 -4.07 5.90e-05  1.32e-03  1.02

rppa3[(duplicated(substr(rownames(rppa3), start=1, stop=4)) | duplicated(substr(rownames(rppa3), start=1, stop=4), fromLast=TRUE)), ] -> phosporatio2
phosporatio2  <- phosporatio2[4:115,]
phosporatio2  <- phosporatio2[-(13:17), ]
setDT(phosporatio2, keep.rownames = TRUE)[]
phosporatio2  <- phosporatio2[-(15:18), ]
phosporatio2  <- phosporatio2[-(19:20), ]
phosporatio2  <- phosporatio2[-(90), ]
phosporatio2  <- phosporatio2[-(91), ]
phosporatio2  <- phosporatio2[-(84:87), ]
phosporatio2  <- phosporatio2[-(80:81), ]
phosporatio2  <- phosporatio2[-(80:81), ]
phosporatio2  <- phosporatio2[-(6:7), ]
phosporatio2  <- phosporatio2[-(22), ]
phosporatio2  <- phosporatio2[-(26), ]
phosporatio2  <- phosporatio2[-(30:31), ]
phosporatio2  <- phosporatio2[-(32:36), ]
phosporatio2  <- phosporatio2[-(5), ]
phosporatio2  <- phosporatio2[-(29:32), ]
phosporatio2  <- phosporatio2[-(33:34), ]
phosporatio2  <- phosporatio2[-(35:36), ]
phosporatio2  <- phosporatio2[-(35:37), ]
phosporatio2  <- phosporatio2[-(42:46), ]
phosporatio2  <- phosporatio2[-(47), ]
phosporatio2  <- phosporatio2[-(53:55), ]
phosporatio2  <- phosporatio2[-(55), ]

phosporatio2 <- data.frame(phosporatio2)
phosporatio2[59:84,] <- 0
phosporatio2[59,1] <- "EIF4EBP1.4E.BP1 Ratio"
phosporatio2[59, 2:330] <- ((phosporatio2[2, 2:330] + phosporatio2[3, 2:330] + phosporatio2[4, 2:330]) - phosporatio2[1, 2:330])
phosporatio2[60,1] <- "PRKAA1.AMPK Ratio"
phosporatio2[60, 2:330] <- (phosporatio2[6, 2:330] - phosporatio2[5, 2:330])
phosporatio2[61,1] <- "AKT1.AKT2.AKT3.AKT Ratio"
phosporatio2[61, 2:330] <- (phosporatio2[8, 2:330] - (phosporatio2[7, 2:330] + phosporatio2[9, 2:330]))
phosporatio2[62,1] <- "RAF1.C.RAF Ratio"
phosporatio2[62, 2:330] <- (phosporatio2[11, 2:330] - phosporatio2[10, 2:330])
phosporatio2[63,1] <- "CHEK1.CHK1 Ratio"
phosporatio2[63, 2:330] <- (phosporatio2[13, 2:330] - phosporatio2[12, 2:330])
phosporatio2[64,1] <- "CHEK2.CHK2 Ratio"
phosporatio2[64, 2:330] <- (phosporatio2[15, 2:330] - phosporatio2[14, 2:330])
phosporatio2[65,1] <- "EGFR.EGFR Ratio"
phosporatio2[65, 2:330] <- ((phosporatio2[18, 2:330] + phosporatio2[17, 2:330]) - phosporatio2[16, 2:330])
phosporatio2[66,1] <- "ESR1.ER.ALPHA Ratio"
phosporatio2[66, 2:330] <- (phosporatio2[20, 2:330] - phosporatio2[19, 2:330])
phosporatio2[67,1] <- "FOXO3.FOXO3A Ratio"
phosporatio2[67, 2:330] <- (phosporatio2[22, 2:330] - phosporatio2[21, 2:330])
phosporatio2[68,1] <- "GSK3A.GSK3B.GSK3.ALPHA.BETA Ratio"
phosporatio2[68, 2:330] <- (phosporatio2[24, 2:330] - phosporatio2[23, 2:330])
phosporatio2[69,1] <- "ERBB2.HER2 Ratio"
phosporatio2[69, 2:330] <- (phosporatio2[26, 2:330] - phosporatio2[25, 2:330])
phosporatio2[70,1] <- "ERBB3.HER3 Ratio"
phosporatio2[70, 2:330] <- (phosporatio2[28, 2:330] - phosporatio2[27, 2:330])
phosporatio2[71,1] <- "PDK1.PDK1 Ratio"
phosporatio2[71, 2:330] <- (phosporatio2[30, 2:330] - phosporatio2[29, 2:330])
phosporatio2[72,1] <- "PEA15.PEA15 Ratio"
phosporatio2[72, 2:330] <- (phosporatio2[32, 2:330] - phosporatio2[31, 2:330])
phosporatio2[73,1] <- "PRKCA..PKC.ALPHA Ratio"
phosporatio2[73, 2:330] <- (phosporatio2[34, 2:330] - phosporatio2[33, 2:330])
phosporatio2[74,1] <- "RB1.RB Ratio"
phosporatio2[74, 2:330] <- (phosporatio2[36, 2:330] - phosporatio2[35, 2:330])
phosporatio2[75,1] <- "RICTOR.RICTOR Ratio"
phosporatio2[75, 2:330] <- (phosporatio2[38, 2:330] - phosporatio2[37, 2:330])
phosporatio2[76,1] <- "RPS6.S6 Ratio"
phosporatio2[76, 2:330] <- ((phosporatio2[41, 2:330] + phosporatio2[40, 2:330]) - phosporatio2[39, 2:330])
phosporatio2[77,1] <- "SRC.SRC Ratio"
phosporatio2[77, 2:330] <- ((phosporatio2[44, 2:330] + phosporatio2[43, 2:330]) - phosporatio2[42, 2:330])
phosporatio2[78,1] <- "TSC2.TUBERIN Ratio"
phosporatio2[78, 2:330] <- (phosporatio2[46, 2:330] - phosporatio2[45, 2:330])
phosporatio2[79,1] <- "YAP1.YAP Ratio"
phosporatio2[79, 2:330] <- (phosporatio2[48, 2:330] - phosporatio2[47, 2:330])
phosporatio2[80,1] <- "YBX1.YB.1 Ratio"
phosporatio2[80, 2:330] <- (phosporatio2[50, 2:330] - phosporatio2[49, 2:330])
phosporatio2[81,1] <- "MTOR.MTOR Ratio"
phosporatio2[81, 2:330] <- (phosporatio2[52, 2:330] - phosporatio2[51, 2:330])
phosporatio2[82,1] <- "MAPK14.P38_MAPK Ratio"
phosporatio2[82, 2:330] <- (phosporatio2[54, 2:330] - phosporatio2[53, 2:330])
phosporatio2[83,1] <- "RPS6KB1.P70S6K Ratio"
phosporatio2[83, 2:330] <- (phosporatio2[56, 2:330] - phosporatio2[55, 2:330])
phosporatio2[84,1] <- "RPS6KA1.P90RSK Ratio"
phosporatio2[84, 2:330] <- (phosporatio2[58, 2:330] - phosporatio2[57, 2:330])


rownames(phosporatio2) <- phosporatio2$rn
phosporatio2$rn <- NULL

#limma of phosphoratios
design2 <- cbind(Grp1=1, Grp2vs1=as.numeric(grepl('X', colnames(phosporatio2))))
fit2 <- lmFit(phosporatio2, design2)
fit2 <- eBayes(fit2)
phosphotable = topTable(fit2,coef=2,number = nrow(phosporatio2))
setDT(phosphotable, keep.rownames = TRUE)[]
phosphotable <- subset(phosphotable, grepl("Ratio", phosphotable$rn), )

phosphoratios <- subset(phosporatio2, grepl("Ratio", rownames(phosporatio2), ))
phosphoratios[27, ] <- colnames(phosphoratios)

#rppa2[1,] <- as.data.frame(lapply(rppa2[1, ],function(x) as.character(x)))

#n <- rppa2["Tumor_Sample_Barcode", ]
#colnames(phosphoratios) <- n

#colnames(phosphoratios) <- phosphoratios[27, ]
#phosphoratios[27, ] <- phosphoratios[28, ]