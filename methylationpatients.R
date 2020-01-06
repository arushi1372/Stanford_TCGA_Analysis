
df <- fread("xx01",data.table=FALSE)
df <- data.frame(df[1:5,])
df$V1 <- NULL
ind <- seq(1, ncol(df), by=2)
df <- df[,ind ]
ind <- seq(1, ncol(df), by=2)
df <- df[,ind ]
         
patients <- fread("patients",data.table=FALSE)
patients$`Hybridization REF`<- NULL
ind <- seq(1, ncol(patients), by=2)
patients <- patients[,ind]
ind <- seq(1, ncol(patients), by=2)
patients <- patients[,ind]

temp <- df
colnames(temp) <- colnames(patients)

temp[6:10,]<- 0
temp[6, ] <- ifelse(temp[1,]>0.200000, 1, 0)
temp[7, ] <- ifelse(temp[2,]>0.200000, 1, 0)
temp[8, ] <- ifelse(temp[3,]>0.200000, 1, 0)
temp[9, ] <- ifelse(temp[4,]>0.200000, 1, 0)
temp[10, ] <- ifelse(temp[5,]>0.200000, 1, 0)

temp[11,] <- 0
temp[11,] <- ifelse((temp[6,]=='1'|temp[7,]=='1'|temp[8,]=='1'|temp[9,]=='1'|temp[10,]=='1'),1,0)
patientmeth <- data.frame(temp[11,])
patientmeth <- data.frame(t(patientmeth))