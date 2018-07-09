require(dplyr)
####import data####
ec <- read.table("enzymelist.txt", sep="\t", quote="")
ec <- rename(ec, EC=V1, name=V2)

####lets see if we can't take TPM from a timepoint and map it to the correct KO# ####
MT_01 <- read.csv("metatranscriptome_tpm/01_MT_TPM.csv")
head(MT_01)

ec_TPM_sum1 <- aggregate(MT_01$rg, by=list(EC=MT_01$EC), FUN=sum)
head(ec_TPM_sum1)
ec_TPM <- merge(ec, ec_TPM_sum1, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[3] <- "MT_01"
head(ec_TPM)
####do it again do it again like 15 more times###
MT_02 <- read.csv("metatranscriptome_tpm/02_MT_TPM.csv")
head(MT_02)

ec_TPM_sum2 <- aggregate(MT_02$rg, by=list(EC=MT_02$EC), FUN=sum)
head(ec_TPM_sum2)
ec_TPM <- merge(ec_TPM, ec_TPM_sum2, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[4] <- "MT_02"
head(ec_TPM)

MT_03 <- read.csv("metatranscriptome_tpm/03_MT_TPM.csv")
head(MT_03)
ec_TPM_sum3 <- aggregate(MT_03$rg, by=list(EC=MT_03$EC), FUN=sum)
head(ec_TPM_sum3)
ec_TPM <- merge(ec_TPM, ec_TPM_sum3, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[5] <- "MT_03"
head(ec_TPM)

MT_05B <- read.csv("metatranscriptome_tpm/05B_MT_TPM.csv")
head(MT_05B)
ec_TPM_sum5B <- aggregate(MT_05B$rg, by=list(EC=MT_05B$EC), FUN=sum)
head(ec_TPM_sum5B)
ec_TPM <- merge(ec_TPM, ec_TPM_sum5B, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[6] <- "MT_05B"
head(ec_TPM, n=50)

MT_05C <- read.csv("metatranscriptome_tpm/05C_MT_TPM.csv")
head(MT_05C)
ec_TPM_sum5C <- aggregate(MT_05C$rg, by=list(EC=MT_05C$EC), FUN=sum)
head(ec_TPM_sum5C)
ec_TPM <- merge(ec_TPM, ec_TPM_sum5C, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[7] <- "MT_05C"
head(ec_TPM)

MT_06 <- read.csv("metatranscriptome_tpm/06_MT_TPM.csv")
head(MT_06)
ec_TPM_sum6 <- aggregate(MT_06$rg, by=list(EC=MT_06$EC), FUN=sum)
head(ec_TPM_sum6)
ec_TPM <- merge(ec_TPM, ec_TPM_sum6, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[8] <- "MT_06"
head(ec_TPM)

MT_07 <- read.csv("metatranscriptome_tpm/07_MT_TPM.csv")
head(MT_07)
ec_TPM_sum7 <- aggregate(MT_07$rg, by=list(EC=MT_07$EC), FUN=sum)
head(ec_TPM_sum7)
ec_TPM <- merge(ec_TPM, ec_TPM_sum7, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[9] <- "MT_07"
head(ec_TPM, n=20)

MT_08 <- read.csv("metatranscriptome_tpm/08_MT_TPM.csv")
head(MT_08)
ec_TPM_sum8 <- aggregate(MT_08$rg, by=list(EC=MT_08$EC), FUN=sum)
head(ec_TPM_sum8)
ec_TPM <- merge(ec_TPM, ec_TPM_sum8, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[10] <- "MT_08"
head(ec_TPM)

MT_09A <- read.csv("metatranscriptome_tpm/09A_MT_TPM.csv")
head(MT_09A)
ec_TPM_sum9A <- aggregate(MT_09A$rg, by=list(EC=MT_09A$EC), FUN=sum)
head(ec_TPM_sum9A)
ec_TPM <- merge(ec_TPM, ec_TPM_sum9A, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[11] <- "MT_09A"
head(ec_TPM)

MT_09B <- read.csv("metatranscriptome_tpm/09B_MT_TPM.csv")
head(MT_09B)
ec_TPM_sum9B <- aggregate(MT_09B$rg, by=list(EC=MT_09B$EC), FUN=sum)
head(ec_TPM_sum9B)
ec_TPM <- merge(ec_TPM, ec_TPM_sum9B, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[12] <- "MT_09B"
head(ec_TPM)

MT_10 <- read.csv("metatranscriptome_tpm/10_MT_TPM.csv")
head(MT_10)
ec_TPM_sum10 <- aggregate(MT_10$rg, by=list(EC=MT_10$EC), FUN=sum)
head(ec_TPM_sum10)
ec_TPM <- merge(ec_TPM, ec_TPM_sum10, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[13] <- "MT_10"
head(ec_TPM)

MT_11 <- read.csv("metatranscriptome_tpm/11_MT_TPM.csv")
head(MT_11)
ec_TPM_sum11 <- aggregate(MT_11$rg, by=list(EC=MT_11$EC), FUN=sum)
head(ec_TPM_sum11)
ec_TPM <- merge(ec_TPM, ec_TPM_sum11, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[14] <- "MT_11"
head(ec_TPM)

MT_12A <- read.csv("metatranscriptome_tpm/12A_MT_TPM.csv")
head(MT_12A)
ec_TPM_sum12A <- aggregate(MT_12A$rg, by=list(EC=MT_12A$EC), FUN=sum)
head(ec_TPM_sum12A)
ec_TPM <- merge(ec_TPM, ec_TPM_sum12A, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[15] <- "MT_12A"
head(ec_TPM)

MT_12B <- read.csv("metatranscriptome_tpm/12B_MT_TPM.csv")
head(MT_12B)
ec_TPM_sum12B <- aggregate(MT_12B$rg, by=list(EC=MT_12B$EC), FUN=sum)
head(ec_TPM_sum12B)
ec_TPM <- merge(ec_TPM, ec_TPM_sum12B, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[16] <- "MT_12B"
head(ec_TPM)

MT_13A <- read.csv("metatranscriptome_tpm/13A_MT_TPM.csv")
head(MT_13A)
ec_TPM_sum13A <- aggregate(MT_13A$rg, by=list(EC=MT_13A$EC), FUN=sum)
head(ec_TPM_sum13A)
ec_TPM <- merge(ec_TPM, ec_TPM_sum13A, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[17] <- "MT_13A"
head(ec_TPM)

MT_13C <- read.csv("metatranscriptome_tpm/13C_MT_TPM.csv")
head(MT_13C)
ec_TPM_sum13C <- aggregate(MT_13C$rg, by=list(EC=MT_13C$EC), FUN=sum)
head(ec_TPM_sum13C)
ec_TPM <- merge(ec_TPM, ec_TPM_sum13C, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[18] <- "MT_13C"
head(ec_TPM)

MG_01 <- read.csv("metagenome_gpm/01_MG_GPM.csv")
head(MG_01)
ec_GPM_sum01 <- aggregate(MG_01$rg, by=list(EC=MG_01$EC), FUN=sum)
head(ec_GPM_sum01)
ec_TPM <- merge(ec_TPM, ec_GPM_sum01, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[19] <- "MG_01"
head(ec_TPM)

MG_02 <- read.csv("metagenome_gpm/02_MG_GPM.csv")
head(MG_02)
ec_GPM_sum02 <- aggregate(MG_02$rg, by=list(EC=MG_02$EC), FUN=sum)
head(ec_GPM_sum02)
ec_TPM <- merge(ec_TPM, ec_GPM_sum02, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[20] <- "MG_02"
head(ec_TPM)

MG_03 <- read.csv("metagenome_gpm/03_MG_GPM.csv")
head(MG_03)
ec_GPM_sum03 <- aggregate(MG_03$rg, by=list(EC=MG_03$EC), FUN=sum)
head(ec_GPM_sum03)
ec_TPM <- merge(ec_TPM, ec_GPM_sum03, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[21] <- "MG_03"
head(ec_TPM)

MG_04 <- read.csv("metagenome_gpm/04_MG_GPM.csv")
head(MG_04)
ec_GPM_sum04 <- aggregate(MG_04$rg, by=list(EC=MG_04$EC), FUN=sum)
head(ec_GPM_sum04)
ec_TPM <- merge(ec_TPM, ec_GPM_sum04, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[22] <- "MG_04"
head(ec_TPM)

MG_05A <- read.csv("metagenome_gpm/05A_MG_GPM.csv")
head(MG_05A)
ec_GPM_sum05A <- aggregate(MG_05A$rg, by=list(EC=MG_05A$EC), FUN=sum)
head(ec_GPM_sum05A)
ec_TPM <- merge(ec_TPM, ec_GPM_sum05A, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[23] <- "MG_05A"
head(ec_TPM)

MG_05C <- read.csv("metagenome_gpm/05C_MG_GPM.csv")
head(MG_05C)
ec_GPM_sum05C <- aggregate(MG_05C$rg, by=list(EC=MG_05C$EC), FUN=sum)
head(ec_GPM_sum05C)
ec_TPM <- merge(ec_TPM, ec_GPM_sum05C, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[24] <- "MG_05C"
head(ec_TPM)

MG_06 <- read.csv("metagenome_gpm/06_MG_GPM.csv")
head(MG_06)
ec_GPM_sum06 <- aggregate(MG_06$rg, by=list(EC=MG_06$EC), FUN=sum)
head(ec_GPM_sum06)
ec_TPM <- merge(ec_TPM, ec_GPM_sum06, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[25] <- "MG_06"
head(ec_TPM)

MG_07 <- read.csv("metagenome_gpm/07_MG_GPM.csv")
head(MG_07)
ec_GPM_sum07 <- aggregate(MG_07$rg, by=list(EC=MG_07$EC), FUN=sum)
head(ec_GPM_sum07)
ec_TPM <- merge(ec_TPM, ec_GPM_sum07, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[26] <- "MG_07"
head(ec_TPM)

MG_09A <- read.csv("metagenome_gpm/09A_MG_GPM.csv")
head(MG_09A)
ec_GPM_sum09A <- aggregate(MG_09A$rg, by=list(EC=MG_09A$EC), FUN=sum)
head(ec_GPM_sum09A)
ec_TPM <- merge(ec_TPM, ec_GPM_sum09A, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[27] <- "MG_09A"
head(ec_TPM)

MG_09B <- read.csv("metagenome_gpm/09B_MG_GPM.csv")
head(MG_09B)
ec_GPM_sum09B <- aggregate(MG_09B$rg, by=list(EC=MG_09B$EC), FUN=sum)
head(ec_GPM_sum09B)
ec_TPM <- merge(ec_TPM, ec_GPM_sum09B, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[28] <- "MG_09B"
head(ec_TPM)

MG_10 <- read.csv("metagenome_gpm/10_MG_GPM.csv")
head(MG_10)
ec_GPM_sum10 <- aggregate(MG_10$rg, by=list(EC=MG_10$EC), FUN=sum)
head(ec_GPM_sum10)
ec_TPM <- merge(ec_TPM, ec_GPM_sum10, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[29] <- "MG_10"
head(ec_TPM)

MG_11 <- read.csv("metagenome_gpm/11_MG_GPM.csv")
head(MG_11)
ec_GPM_sum11 <- aggregate(MG_11$rg, by=list(EC=MG_11$EC), FUN=sum)
head(ec_GPM_sum11)
ec_TPM <- merge(ec_TPM, ec_GPM_sum11, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[30] <- "MG_11"
head(ec_TPM)

MG_12A <- read.csv("metagenome_gpm/12A_MG_GPM.csv")
head(MG_12A)
ec_GPM_sum12A <- aggregate(MG_12A$rg, by=list(EC=MG_12A$EC), FUN=sum)
head(ec_GPM_sum12A)
ec_TPM <- merge(ec_TPM, ec_GPM_sum12A, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[31] <- "MG_12A"
head(ec_TPM)

MG_12B <- read.csv("metagenome_gpm/12B_MG_GPM.csv")
head(MG_12B)
ec_GPM_sum12B <- aggregate(MG_12B$rg, by=list(EC=MG_12B$EC), FUN=sum)
head(ec_GPM_sum12B)
ec_TPM <- merge(ec_TPM, ec_GPM_sum12B, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[32] <- "MG_12B"
head(ec_TPM)

MG_13A <- read.csv("metagenome_gpm/13A_MG_GPM.csv")
head(MG_13A)
ec_GPM_sum13A <- aggregate(MG_13A$rg, by=list(EC=MG_13A$EC), FUN=sum)
head(ec_GPM_sum13A)
ec_TPM <- merge(ec_TPM, ec_GPM_sum13A, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[33] <- "MG_13A"
head(ec_TPM)

MG_13B <- read.csv("metagenome_gpm/13B_MG_GPM.csv")
head(MG_13B)
ec_GPM_sum13B <- aggregate(MG_13B$rg, by=list(EC=MG_13B$EC), FUN=sum)
head(ec_GPM_sum13B)
ec_TPM <- merge(ec_TPM, ec_GPM_sum13B, by="EC", all.x=TRUE)
head(ec_TPM)
colnames(ec_TPM)[34] <- "MG_13B"
head(ec_TPM)

ec_TPM[is.na(ec_TPM)] <- 0
write.csv(ec_TPM, "ec_reads.csv", row.names=FALSE)

