####import data####
ko_metadata <- read.csv("ko_metadata.csv")
head(ko_metadata)
####only want unique KO#s####
ko_metadata$KO_number <- sub("^", "KO:", ko_metadata$KO_number)
head(ko_metadata)

ko_num <- subset(ko_metadata, select = c(KO_number))
ko_num <- unique(ko_num)
head(ko_num)
####lets see if we can't take TPM from a timepoint and map it to the correct KO# ####
MT_01 <- read.csv("metatranscriptome_tpm/01_MT_TPM.csv")
head(MT_01)

MT_01$KO_number <- MT_01$KO_num
KO_TPM_sum <- aggregate(MT_01$rg, by=list(KO_number=MT_01$KO_num), FUN=sum)
head(KO_TPM_sum)
KO_number <- merge(ko_num, KO_TPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[2] <- "MT_01"
head(KO_number)

####do it again do it again like 15 more times###
MT_02 <- read.csv("metatranscriptome_tpm/02_MT_TPM.csv")
MT_02$KO_number <- MT_02$KO_num
KO_TPM_sum2 <- aggregate(MT_02$rg, by=list(KO_number=MT_02$KO_num), FUN=sum)
head(KO_TPM_sum2)
KO_number <- merge(KO_number, KO_TPM_sum2, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[3] <- "MT_02"
head(KO_number)

MT_03 <- read.csv("metatranscriptome_tpm/03_MT_TPM.csv")
MT_03$KO_number <- MT_03$KO_num
KO_TPM_sum3 <- aggregate(MT_03$rg, by=list(KO_number=MT_03$KO_num), FUN=sum)
head(KO_TPM_sum3)
KO_number <- merge(KO_number, KO_TPM_sum3, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[4] <- "MT_03"
head(KO_number)

MT_05B <- read.csv("metatranscriptome_tpm/05B_MT_TPM.csv")
MT_05B$KO_number <- MT_05B$KO_num
KO_TPM_sum5B <- aggregate(MT_05B$rg, by=list(KO_number=MT_05B$KO_num), FUN=sum)
head(KO_TPM_sum5B)
KO_number <- merge(KO_number, KO_TPM_sum5B, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[5] <- "MT_05B"
head(KO_number)

MT_05C <- read.csv("metatranscriptome_tpm/05C_MT_TPM.csv")
MT_05C$KO_number <- MT_05C$KO_num
KO_TPM_sum5C <- aggregate(MT_05C$rg, by=list(KO_number=MT_05C$KO_num), FUN=sum)
head(KO_TPM_sum5C)
KO_number <- merge(KO_number, KO_TPM_sum5C, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[6] <- "MT_05C"
head(KO_number)

MT_06 <- read.csv("metatranscriptome_tpm/06_MT_TPM.csv")
MT_06$KO_number <- MT_06$KO_num
KO_TPM_sum6 <- aggregate(MT_06$rg, by=list(KO_number=MT_06$KO_num), FUN=sum)
head(KO_TPM_sum6)
KO_number <- merge(KO_number, KO_TPM_sum6, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[7] <- "MT_06"
head(KO_number)

MT_07 <- read.csv("metatranscriptome_tpm/07_MT_TPM.csv")
MT_07$KO_number <- MT_07$KO_num
KO_TPM_sum7 <- aggregate(MT_07$rg, by=list(KO_number=MT_07$KO_num), FUN=sum)
head(KO_TPM_sum7)
KO_number <- merge(KO_number, KO_TPM_sum7, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[8] <- "MT_07"
head(KO_number)

MT_08 <- read.csv("metatranscriptome_tpm/08_MT_TPM.csv")
MT_08$KO_number <- MT_08$KO_num
KO_TPM_sum8 <- aggregate(MT_08$rg, by=list(KO_number=MT_08$KO_num), FUN=sum)
head(KO_TPM_sum8)
KO_number <- merge(KO_number, KO_TPM_sum8, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[9] <- "MT_08"
head(KO_number)

MT_09A <- read.csv("metatranscriptome_tpm/09A_MT_TPM.csv")
MT_09A$KO_number <- MT_09A$KO_num
KO_TPM_sum9A <- aggregate(MT_09A$rg, by=list(KO_number=MT_09A$KO_num), FUN=sum)
head(KO_TPM_sum9A)
KO_number <- merge(KO_number, KO_TPM_sum9A, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[10] <- "MT_09A"
head(KO_number)

MT_09B <- read.csv("metatranscriptome_tpm/09B_MT_TPM.csv")
MT_09B$KO_number <- MT_09B$KO_num
KO_TPM_sum9B <- aggregate(MT_09B$rg, by=list(KO_number=MT_09B$KO_num), FUN=sum)
head(KO_TPM_sum9B)
KO_number <- merge(KO_number, KO_TPM_sum9B, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[11] <- "MT_09B"
head(KO_number)

MT_10 <- read.csv("metatranscriptome_tpm/10_MT_TPM.csv")
MT_10$KO_number <- MT_10$KO_num
KO_TPM_sum10 <- aggregate(MT_10$rg, by=list(KO_number=MT_10$KO_num), FUN=sum)
head(KO_TPM_sum10)
KO_number <- merge(KO_number, KO_TPM_sum10, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[12] <- "MT_10"
head(KO_number)

MT_11 <- read.csv("metatranscriptome_tpm/11_MT_TPM.csv")
MT_11$KO_number <- MT_11$KO_num
KO_TPM_sum11 <- aggregate(MT_11$rg, by=list(KO_number=MT_11$KO_num), FUN=sum)
head(KO_TPM_sum11)
KO_number <- merge(KO_number, KO_TPM_sum11, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[13] <- "MT_11"
head(KO_number)

MT_12A <- read.csv("metatranscriptome_tpm/12A_MT_TPM.csv")
MT_12A$KO_number <- MT_12A$KO_num
KO_TPM_sum12A <- aggregate(MT_12A$rg, by=list(KO_number=MT_12A$KO_num), FUN=sum)
head(KO_TPM_sum12A)
KO_number <- merge(KO_number, KO_TPM_sum12A, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[14] <- "MT_12A"
head(KO_number)

MT_12B <- read.csv("metatranscriptome_tpm/12B_MT_TPM.csv")
MT_12B$KO_number <- MT_12B$KO_num
KO_TPM_sum12B <- aggregate(MT_12B$rg, by=list(KO_number=MT_12B$KO_num), FUN=sum)
head(KO_TPM_sum12B)
KO_number <- merge(KO_number, KO_TPM_sum12B, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[15] <- "MT_12B"
head(KO_number)

MT_13A <- read.csv("metatranscriptome_tpm/13A_MT_TPM.csv")
MT_13A$KO_number <- MT_13A$KO_num
KO_TPM_sum13A <- aggregate(MT_13A$rg, by=list(KO_number=MT_13A$KO_num), FUN=sum)
head(KO_TPM_sum13A)
KO_number <- merge(KO_number, KO_TPM_sum13A, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[16] <- "MT_13A"
head(KO_number)

MT_13C <- read.csv("metatranscriptome_tpm/13C_MT_TPM.csv")
MT_13C$KO_number <- MT_13C$KO_num
KO_TPM_sum13C <- aggregate(MT_13C$rg, by=list(KO_number=MT_13C$KO_num), FUN=sum)
head(KO_TPM_sum13C)
KO_number <- merge(KO_number, KO_TPM_sum13C, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[17] <- "MT_13C"
head(KO_number)

MG_01 <- read.csv("metagenome_gpm/01_MG_GPM.csv")
head(MG_01)

MG_01$KO_number <- MG_01$KO_num
KO_GPM_sum <- aggregate(MG_01$rg, by=list(KO_number=MG_01$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[18] <- "MG_01"
head(KO_number)

MG_09B <- read.csv("metagenome_gpm/09B_MG_GPM.csv")
head(MG_09B)

MG_09B$KO_number <- MG_09B$KO_num
KO_GPM_sum <- aggregate(MG_09B$rg, by=list(KO_number=MG_09B$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[19] <- "MG_09B"
head(KO_number)

MG_13B <- read.csv("metagenome_gpm/13B_MG_GPM.csv")
head(MG_13B)

MG_13B$KO_number <- MG_13B$KO_num
KO_GPM_sum <- aggregate(MG_13B$rg, by=list(KO_number=MG_13B$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[20] <- "MG_13B"
head(KO_number)

MG_02 <- read.csv("metagenome_gpm/02_MG_GPM.csv")
head(MG_02)

MG_02$KO_number <- MG_02$KO_num
KO_GPM_sum <- aggregate(MG_02$rg, by=list(KO_number=MG_02$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[21] <- "MG_02"
head(KO_number)

MG_04 <- read.csv("metagenome_gpm/04_MG_GPM.csv")
head(MG_04)

MG_04$KO_number <- MG_04$KO_num
KO_GPM_sum <- aggregate(MG_04$rg, by=list(KO_number=MG_04$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[22] <- "MG_04"
head(KO_number)

MG_05A <- read.csv("metagenome_gpm/05A_MG_GPM.csv")
head(MG_05A)

MG_05A$KO_number <- MG_05A$KO_num
KO_GPM_sum <- aggregate(MG_05A$rg, by=list(KO_number=MG_05A$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[23] <- "MG_05A"
head(KO_number)

KO_number[is.na(KO_number)] <- 0

MG_03 <- read.csv("metagenome_gpm/03_MG_GPM.csv")
head(MG_03)

MG_03$KO_number <- MG_03$KO_num
KO_GPM_sum <- aggregate(MG_03$rg, by=list(KO_number=MG_03$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[24] <- "MG_03"
head(KO_number)

MG_05C <- read.csv("metagenome_gpm/05C_MG_GPM.csv")
head(MG_05C)

MG_05C$KO_number <- MG_05C$KO_num
KO_GPM_sum <- aggregate(MG_05C$rg, by=list(KO_number=MG_05C$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[25] <- "MG_05C"
head(KO_number)

MG_06 <- read.csv("metagenome_gpm/06_MG_GPM.csv")
head(MG_06)

MG_06$KO_number <- MG_06$KO_num
KO_GPM_sum <- aggregate(MG_06$rg, by=list(KO_number=MG_06$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[26] <- "MG_06"
head(KO_number)

MG_07 <- read.csv("metagenome_gpm/07_MG_GPM.csv")
head(MG_07)

MG_07$KO_number <- MG_07$KO_num
KO_GPM_sum <- aggregate(MG_07$rg, by=list(KO_number=MG_07$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[27] <- "MG_07"
head(KO_number)

MG_09A <- read.csv("metagenome_gpm/09A_MG_GPM.csv")
head(MG_09A)

MG_09A$KO_number <- MG_09A$KO_num
KO_GPM_sum <- aggregate(MG_09A$rg, by=list(KO_number=MG_09A$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[28] <- "MG_09A"
head(KO_number)

MG_10 <- read.csv("metagenome_gpm/10_MG_GPM.csv")
head(MG_10)

MG_10$KO_number <- MG_10$KO_num
KO_GPM_sum <- aggregate(MG_10$rg, by=list(KO_number=MG_10$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[29] <- "MG_10"
head(KO_number)

MG_11 <- read.csv("metagenome_gpm/11_MG_GPM.csv")
head(MG_11)

MG_11$KO_number <- MG_11$KO_num
KO_GPM_sum <- aggregate(MG_11$rg, by=list(KO_number=MG_11$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[30] <- "MG_11"
head(KO_number)

MG_12A <- read.csv("metagenome_gpm/12A_MG_GPM.csv")
head(MG_12A)

MG_12A$KO_number <- MG_12A$KO_num
KO_GPM_sum <- aggregate(MG_12A$rg, by=list(KO_number=MG_12A$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[31] <- "MG_12A"
head(KO_number)

MG_12B <- read.csv("metagenome_gpm/12B_MG_GPM.csv")
head(MG_12B)

MG_12B$KO_number <- MG_12B$KO_num
KO_GPM_sum <- aggregate(MG_12B$rg, by=list(KO_number=MG_12B$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[32] <- "MG_12B"
head(KO_number)

MG_13A <- read.csv("metagenome_gpm/13A_MG_GPM.csv")
head(MG_13A)

MG_13A$KO_number <- MG_13A$KO_num
KO_GPM_sum <- aggregate(MG_13A$rg, by=list(KO_number=MG_13A$KO_num), FUN=sum)
head(KO_GPM_sum)
KO_number <- merge(KO_number, KO_GPM_sum, by="KO_number", all.x=TRUE)
head(KO_number)
tail(KO_number)
colnames(KO_number)[33] <- "MG_13A"
head(KO_number)

KO_number[is.na(KO_number)] <- 0

write.csv(KO_number, "MGMT_KOnum_reads.csv", row.names=FALSE)
