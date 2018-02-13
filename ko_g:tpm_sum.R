####import data####
ko_metadata <- read.csv("ko_metadata.csv")
head(ko_metadata)
####only want unique KO#s####
ko_num <- subset(ko_metadata, select = c(KO_number))
ko_num <- unique(ko_num)
ko_num$KO_number <- sub("^", "KO:", ko_num$KO_number)
head(ko_num)

####lets see if we can't take TPM from a timepoint and map it to the correct KO# ####
MT_01 <- read.csv("071401_MT_TPM_KO.csv")
head(MT_01)

MT_01$KO_number <- MT_01$KO_num
KO_TPM_sum <- aggregate(MT_01$TPM, by=list(KO_number=MT_01$KO_num), FUN=sum)
head(KO_TPM_sum)
KO_number <- merge(ko_num, KO_TPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[2] <- "MT_01"
head(KO_number)

####do it again do it again like 15 more times###
MT_02 <- read.csv("071402_MT_TPM_KO.csv")
MT_02$KO_number <- MT_02$KO_num
KO_TPM_sum2 <- aggregate(MT_02$TPM, by=list(KO_number=MT_02$KO_num), FUN=sum)
head(KO_TPM_sum2)
KO_number <- merge(KO_number, KO_TPM_sum2, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[3] <- "MT_02"
head(KO_number)

MT_03 <- read.csv("071403_MT_TPM_KO.csv")
MT_03$KO_number <- MT_03$KO_num
KO_TPM_sum3 <- aggregate(MT_03$TPM, by=list(KO_number=MT_03$KO_num), FUN=sum)
head(KO_TPM_sum3)
KO_number <- merge(KO_number, KO_TPM_sum3, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[4] <- "MT_03"
head(KO_number)

MT_05B <- read.csv("071405B_MT_TPM_KO.csv")
MT_05B$KO_number <- MT_05B$KO_num
KO_TPM_sum5B <- aggregate(MT_05B$TPM, by=list(KO_number=MT_05B$KO_num), FUN=sum)
head(KO_TPM_sum5B)
KO_number <- merge(KO_number, KO_TPM_sum5B, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[5] <- "MT_05B"
head(KO_number)

MT_05C <- read.csv("071405C_MT_TPM_KO.csv")
MT_05C$KO_number <- MT_05C$KO_num
KO_TPM_sum5C <- aggregate(MT_05C$TPM, by=list(KO_number=MT_05C$KO_num), FUN=sum)
head(KO_TPM_sum5C)
KO_number <- merge(KO_number, KO_TPM_sum5C, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[6] <- "MT_05C"
head(KO_number)

MT_06 <- read.csv("071406_MT_TPM_KO.csv")
MT_06$KO_number <- MT_06$KO_num
KO_TPM_sum6 <- aggregate(MT_06$TPM, by=list(KO_number=MT_06$KO_num), FUN=sum)
head(KO_TPM_sum6)
KO_number <- merge(KO_number, KO_TPM_sum6, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[7] <- "MT_06"
head(KO_number)

MT_07 <- read.csv("071407_MT_TPM_KO.csv")
MT_07$KO_number <- MT_07$KO_num
KO_TPM_sum7 <- aggregate(MT_07$TPM, by=list(KO_number=MT_07$KO_num), FUN=sum)
head(KO_TPM_sum7)
KO_number <- merge(KO_number, KO_TPM_sum7, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[8] <- "MT_07"
head(KO_number)

MT_08 <- read.csv("071408_MT_TPM_KO.csv")
MT_08$KO_number <- MT_08$KO_num
KO_TPM_sum8 <- aggregate(MT_08$TPM, by=list(KO_number=MT_08$KO_num), FUN=sum)
head(KO_TPM_sum8)
KO_number <- merge(KO_number, KO_TPM_sum8, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[9] <- "MT_08"
head(KO_number)

MT_09A <- read.csv("071409A_MT_TPM_KO.csv")
MT_09A$KO_number <- MT_09A$KO_num
KO_TPM_sum9A <- aggregate(MT_09A$TPM, by=list(KO_number=MT_09A$KO_num), FUN=sum)
head(KO_TPM_sum9A)
KO_number <- merge(KO_number, KO_TPM_sum9A, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[10] <- "MT_09A"
head(KO_number)

MT_09B <- read.csv("071409B_MT_TPM_KO.csv")
MT_09B$KO_number <- MT_09B$KO_num
KO_TPM_sum9B <- aggregate(MT_09B$TPM, by=list(KO_number=MT_09B$KO_num), FUN=sum)
head(KO_TPM_sum9B)
KO_number <- merge(KO_number, KO_TPM_sum9B, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[11] <- "MT_09B"
head(KO_number)

MT_10 <- read.csv("071410_MT_TPM_KO.csv")
MT_10$KO_number <- MT_10$KO_num
KO_TPM_sum10 <- aggregate(MT_10$TPM, by=list(KO_number=MT_10$KO_num), FUN=sum)
head(KO_TPM_sum10)
KO_number <- merge(KO_number, KO_TPM_sum10, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[12] <- "MT_10"
head(KO_number)

MT_11 <- read.csv("071411_MT_TPM_KO.csv")
MT_11$KO_number <- MT_11$KO_num
KO_TPM_sum11 <- aggregate(MT_11$TPM, by=list(KO_number=MT_11$KO_num), FUN=sum)
head(KO_TPM_sum11)
KO_number <- merge(KO_number, KO_TPM_sum11, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[13] <- "MT_11"
head(KO_number)

MT_12A <- read.csv("071412A_MT_TPM_KO.csv")
MT_12A$KO_number <- MT_12A$KO_num
KO_TPM_sum12A <- aggregate(MT_12A$TPM, by=list(KO_number=MT_12A$KO_num), FUN=sum)
head(KO_TPM_sum12A)
KO_number <- merge(KO_number, KO_TPM_sum12A, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[14] <- "MT_12A"
head(KO_number)

MT_12B <- read.csv("071412B_MT_TPM_KO.csv")
MT_12B$KO_number <- MT_12B$KO_num
KO_TPM_sum12B <- aggregate(MT_12B$TPM, by=list(KO_number=MT_12B$KO_num), FUN=sum)
head(KO_TPM_sum12B)
KO_number <- merge(KO_number, KO_TPM_sum12B, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[15] <- "MT_12B"
head(KO_number)

MT_13A <- read.csv("071413A_MT_TPM_KO.csv")
MT_13A$KO_number <- MT_13A$KO_num
KO_TPM_sum13A <- aggregate(MT_13A$TPM, by=list(KO_number=MT_13A$KO_num), FUN=sum)
head(KO_TPM_sum13A)
KO_number <- merge(KO_number, KO_TPM_sum13A, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[16] <- "MT_13A"
head(KO_number)

MT_13C <- read.csv("071413C_MT_TPM_KO.csv")
MT_13C$KO_number <- MT_13C$KO_num
KO_TPM_sum13C <- aggregate(MT_13C$TPM, by=list(KO_number=MT_13C$KO_num), FUN=sum)
head(KO_TPM_sum13C)
KO_number <- merge(KO_number, KO_TPM_sum13C, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[17] <- "MT_13C"
head(KO_number)

MG_01 <- read.csv("071401_MG_TPM_KO.csv")
head(MG_01)

MG_01$KO_number <- MG_01$KO_num
KO_GPM_sum <- aggregate(MG_01$GPM, by=list(KO_number=MG_01$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[18] <- "MG_01"
head(KO_number)

MG_09B <- read.csv("071409B_MG_GPM_KO.csv")
head(MG_09B)

MG_09B$KO_number <- MG_09B$KO_num
KO_GPM_sum <- aggregate(MG_09B$GPM, by=list(KO_number=MG_09B$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[19] <- "MG_09B"
head(KO_number)

MG_13B <- read.csv("071413B_MG_GPM_KO.csv")
head(MG_13B)

MG_13B$KO_number <- MG_13B$KO_num
KO_GPM_sum <- aggregate(MG_13B$GPM, by=list(KO_number=MG_13B$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[20] <- "MG_13B"
head(KO_number)

KO_number[is.na(KO_number)] <- 0

write.csv(KO_number, "MGMT_KOnum_GPMTPM.csv", row.names=FALSE)
