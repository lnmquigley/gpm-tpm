require(dplyr)
require(splitstackshape)
MG_01 <- read.table("CAZy_results/01_MG_CAZy.tab", sep="\t")
head(MG_01)
MG_01 <- cSplit(MG_01, "V1", ":")
MG_01 <- rename(MG_01, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_01 <- subset(MG_01, select=-c(V6, V7, V8, V9, V10, V12))
MG_01 <- subset(MG_01, align_length >=30)
MG_01$type <- "metagenome"
MG_01$time <- 1
MG_01_GPM <- read.csv("metagenome_gpm/01_MG_GPM.csv")
head(MG_01_GPM)
MG_01 <- merge(MG_01, MG_01_GPM, all=FALSE)

MG_02 <- read.table("CAZy_results/02_MG_CAZy.tab", sep="\t")
head(MG_02)
MG_02 <- cSplit(MG_02, "V1", ":")
MG_02 <- rename(MG_02, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_02 <- subset(MG_02, select=-c(V6, V7, V8, V9, V10, V12))
MG_02 <- subset(MG_02, align_length >=30)
MG_02$type <- "metagenome"
MG_02$time <- 2
MG_02_GPM <- read.csv("metagenome_gpm/02_MG_GPM.csv")
head(MG_02_GPM)
MG_02 <- merge(MG_02, MG_02_GPM, all=FALSE)

MG_03 <- read.table("CAZy_results/03_MG_CAZy.tab", sep="\t")
head(MG_03)
MG_03 <- cSplit(MG_03, "V1", ":")
MG_03 <- rename(MG_03, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_03 <- subset(MG_03, select=-c(V6, V7, V8, V9, V10, V12))
MG_03 <- subset(MG_03, align_length >=30)
MG_03$type <- "metagenome"
MG_03$time <- 3
MG_03_GPM <- read.csv("metagenome_gpm/03_MG_GPM.csv")
head(MG_03_GPM)
MG_03 <- merge(MG_03, MG_03_GPM, all=FALSE)

MG_04 <- read.table("CAZy_results/04_MG_CAZy.tab", sep="\t")
head(MG_04)
MG_04 <- cSplit(MG_04, "V1", ":")
MG_04 <- rename(MG_04, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_04 <- subset(MG_04, select=-c(V6, V7, V8, V9, V10, V12))
MG_04 <- subset(MG_04, align_length >=30)
MG_04$type <- "metagenome"
MG_04$time <- 4
MG_04_GPM <- read.csv("metagenome_gpm/04_MG_GPM.csv")
head(MG_04_GPM)
MG_04 <- merge(MG_04, MG_04_GPM, all=FALSE)

MG_05A <- read.table("CAZy_results/05A_MG_CAZy.tab", sep="\t")
head(MG_05A)
MG_05A <- cSplit(MG_05A, "V1", ":")
MG_05A <- rename(MG_05A, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_05A <- subset(MG_05A, select=-c(V6, V7, V8, V9, V10, V12))
MG_05A <- subset(MG_05A, align_length >=30)
MG_05A$type <- "metagenome"
MG_05A$time <- 5
MG_05A_GPM <- read.csv("metagenome_gpm/05A_MG_GPM.csv")
head(MG_05A_GPM)
MG_05A <- merge(MG_05A, MG_05A_GPM, all=FALSE)

MG_05C <- read.table("CAZy_results/05C_MG_CAZy.tab", sep="\t")
head(MG_05C)
MG_05C <- cSplit(MG_05C, "V1", ":")
MG_05C <- rename(MG_05C, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_05C <- subset(MG_05C, select=-c(V6, V7, V8, V9, V10, V12))
MG_05C <- subset(MG_05C, align_length >=30)
MG_05C$type <- "metagenome"
MG_05C$time <- 5.5
MG_05C_GPM <- read.csv("metagenome_gpm/05C_MG_GPM.csv")
head(MG_05C_GPM)
MG_05C <- merge(MG_05C, MG_05C_GPM, all=FALSE)

MG_06 <- read.table("CAZy_results/06_MG_CAZy.tab", sep="\t")
head(MG_06)
MG_06 <- cSplit(MG_06, "V1", ":")
MG_06 <- rename(MG_06, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_06 <- subset(MG_06, select=-c(V6, V7, V8, V9, V10, V12))
MG_06 <- subset(MG_06, align_length >=30)
MG_06$type <- "metagenome"
MG_06$time <- 6
MG_06_GPM <- read.csv("metagenome_gpm/06_MG_GPM.csv")
head(MG_06_GPM)
MG_06 <- merge(MG_06, MG_06_GPM, all=FALSE)

MG_07 <- read.table("CAZy_results/07_MG_CAZy.tab", sep="\t")
head(MG_07)
MG_07 <- cSplit(MG_07, "V1", ":")
MG_07 <- rename(MG_07, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_07 <- subset(MG_07, select=-c(V6, V7, V8, V9, V10, V12))
MG_07 <- subset(MG_07, align_length >=30)
MG_07$type <- "metagenome"
MG_07$time <- 7
MG_07_GPM <- read.csv("metagenome_gpm/07_MG_GPM.csv")
head(MG_07_GPM)
MG_07 <- merge(MG_07, MG_07_GPM, all=FALSE)

MG_09A <- read.table("CAZy_results/09A_MG_CAZy.tab", sep="\t")
head(MG_09A)
MG_09A <- cSplit(MG_09A, "V1", ":")
MG_09A <- rename(MG_09A, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_09A <- subset(MG_09A, select=-c(V6, V7, V8, V9, V10, V12))
MG_09A <- subset(MG_09A, align_length >=30)
MG_09A$type <- "metagenome"
MG_09A$time <- 9
MG_09A_GPM <- read.csv("metagenome_gpm/09A_MG_GPM.csv")
head(MG_09A_GPM)
MG_09A <- merge(MG_09A, MG_09A_GPM, all=FALSE)

MG_09B <- read.table("CAZy_results/09B_MG_CAZy.tab", sep="\t")
head(MG_09B)
MG_09B <- cSplit(MG_09B, "V1", ":")
MG_09B <- rename(MG_09B, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_09B <- subset(MG_09B, select=-c(V6, V7, V8, V9, V10, V12))
MG_09B <- subset(MG_09B, align_length >=30)
MG_09B$type <- "metagenome"
MG_09B$time <- 9.5
MG_09B_GPM <- read.csv("metagenome_gpm/09B_MG_GPM.csv")
head(MG_09B_GPM)
MG_09B <- merge(MG_09B, MG_09B_GPM, all=FALSE)

MG_10 <- read.table("CAZy_results/10_MG_CAZy.tab", sep="\t")
head(MG_10)
MG_10 <- cSplit(MG_10, "V1", ":")
MG_10 <- rename(MG_10, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_10 <- subset(MG_10, select=-c(V6, V7, V8, V9, V10, V12))
MG_10 <- subset(MG_10, align_length >=30)
MG_10$type <- "metagenome"
MG_10$time <- 10
MG_10_GPM <- read.csv("metagenome_gpm/10_MG_GPM.csv")
head(MG_10_GPM)
MG_10 <- merge(MG_10, MG_10_GPM, all=FALSE)

MG_11 <- read.table("CAZy_results/11_MG_CAZy.tab", sep="\t")
head(MG_11)
MG_11 <- cSplit(MG_11, "V1", ":")
MG_11 <- rename(MG_11, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_11 <- subset(MG_11, select=-c(V6, V7, V8, V9, V10, V12))
MG_11 <- subset(MG_11, align_length >=30)
MG_11$type <- "metagenome"
MG_11$time <- 11
MG_11_GPM <- read.csv("metagenome_gpm/11_MG_GPM.csv")
head(MG_11_GPM)
MG_11 <- merge(MG_11, MG_11_GPM, all=FALSE)

MG_12A <- read.table("CAZy_results/12A_MG_CAZy.tab", sep="\t")
head(MG_12A)
MG_12A <- cSplit(MG_12A, "V1", ":")
MG_12A <- rename(MG_12A, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_12A <- subset(MG_12A, select=-c(V6, V7, V8, V9, V10, V12))
MG_12A <- subset(MG_12A, align_length >=30)
MG_12A$type <- "metagenome"
MG_12A$time <- 12
MG_12A_GPM <- read.csv("metagenome_gpm/12A_MG_GPM.csv")
head(MG_12A_GPM)
MG_12A <- merge(MG_12A, MG_12A_GPM, all=FALSE)

MG_12B <- read.table("CAZy_results/12B_MG_CAZy.tab", sep="\t")
head(MG_12B)
MG_12B <- cSplit(MG_12B, "V1", ":")
MG_12B <- rename(MG_12B, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_12B <- subset(MG_12B, select=-c(V6, V7, V8, V9, V10, V12))
MG_12B <- subset(MG_12B, align_length >=30)
MG_12B$type <- "metagenome"
MG_12B$time <- 12.5
MG_12B_GPM <- read.csv("metagenome_gpm/12B_MG_GPM.csv")
head(MG_12B_GPM)
MG_12B <- merge(MG_12B, MG_12B_GPM, all=FALSE)

MG_13A <- read.table("CAZy_results/13A_MG_CAZy.tab", sep="\t")
head(MG_13A)
MG_13A <- cSplit(MG_13A, "V1", ":")
MG_13A <- rename(MG_13A, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_13A <- subset(MG_13A, select=-c(V6, V7, V8, V9, V10, V12))
MG_13A <- subset(MG_13A, align_length >=30)
MG_13A$type <- "metagenome"
MG_13A$time <- 13
MG_13A_GPM <- read.csv("metagenome_gpm/13A_MG_GPM.csv")
head(MG_13A_GPM)
MG_13A <- merge(MG_13A, MG_13A_GPM, all=FALSE)

MG_13B <- read.table("CAZy_results/13B_MG_CAZy.tab", sep="\t")
head(MG_13B)
MG_13B <- cSplit(MG_13B, "V1", ":")
MG_13B <- rename(MG_13B, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MG_13B <- subset(MG_13B, select=-c(V6, V7, V8, V9, V10, V12))
MG_13B <- subset(MG_13B, align_length >=30)
MG_13B$type <- "metagenome"
MG_13B$time <- 13.5
MG_13B_GPM <- read.csv("metagenome_gpm/13B_MG_GPM.csv")
head(MG_13B_GPM)
MG_13B <- merge(MG_13B, MG_13B_GPM, all=FALSE)

MT_01 <- read.table("CAZy_results/01_MT_CAZy.m8", sep="\t")
head(MT_01)
MT_01 <- cSplit(MT_01, "V1", ":")
MT_01 <- rename(MT_01, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_01 <- subset(MT_01, select=-c(V6, V7, V8, V9, V10, V12))
MT_01 <- subset(MT_01, align_length >=30)
MT_01$type <- "metatranscriptome"
MT_01$time <- 1
MT_01_TPM <- read.csv("metatranscriptome_tpm/01_MT_TPM.csv")
head(MT_01_TPM)
MT_01_TPM <- rename(MT_01_TPM, GPM=TPM)
MT_01 <- merge(MT_01, MT_01_TPM, all=FALSE)

MT_02 <- read.table("CAZy_results/02_MT_CAZy.tab", sep="\t")
head(MT_02)
MT_02 <- cSplit(MT_02, "V1", ":")
MT_02 <- rename(MT_02, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_02 <- subset(MT_02, select=-c(V6, V7, V8, V9, V10, V12))
MT_02 <- subset(MT_02, align_length >=30)
MT_02$type <- "metatranscriptome"
MT_02$time <- 2
MT_02_TPM <- read.csv("metatranscriptome_tpm/02_MT_TPM.csv")
head(MT_02_TPM)
MT_02_TPM <- rename(MT_02_TPM, GPM=TPM)
MT_02 <- merge(MT_02, MT_02_TPM, all=FALSE)

MT_03 <- read.table("CAZy_results/03_MT_CAZy.tab", sep="\t")
head(MT_03)
MT_03 <- cSplit(MT_03, "V1", ":")
MT_03 <- rename(MT_03, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_03 <- subset(MT_03, select=-c(V6, V7, V8, V9, V10, V12))
MT_03 <- subset(MT_03, align_length >=30)
MT_03$type <- "metatranscriptome"
MT_03$time <- 3
MT_03_TPM <- read.csv("metatranscriptome_tpm/03_MT_TPM.csv")
head(MT_03_TPM)
MT_03_TPM <- rename(MT_03_TPM, GPM=TPM)
MT_03 <- merge(MT_03, MT_03_TPM, all=FALSE)

MT_05B <- read.table("CAZy_results/05B_MT_CAZy.m8", sep="\t")
head(MT_05B)
MT_05B <- cSplit(MT_05B, "V1", ":")
MT_05B <- rename(MT_05B, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_05B <- subset(MT_05B, select=-c(V6, V7, V8, V9, V10, V12))
MT_05B <- subset(MT_05B, align_length >=30)
MT_05B$type <- "metatranscriptome"
MT_05B$time <- 5
MT_05B_TPM <- read.csv("metatranscriptome_tpm/05B_MT_TPM.csv")
head(MT_05B_TPM)
MT_05B_TPM <- rename(MT_05B_TPM, GPM=TPM)
MT_05B <- merge(MT_05B, MT_05B_TPM, all=FALSE)

MT_05C <- read.table("CAZy_results/05C_MT_CAZy.tab", sep="\t")
head(MT_05C)
MT_05C <- cSplit(MT_05C, "V1", ":")
MT_05C <- rename(MT_05C, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_05C <- subset(MT_05C, select=-c(V6, V7, V8, V9, V10, V12))
MT_05C <- subset(MT_05C, align_length >=30)
MT_05C$type <- "metatranscriptome"
MT_05C$time <- 5.5
MT_05C_TPM <- read.csv("metatranscriptome_tpm/05C_MT_TPM.csv")
head(MT_05C_TPM)
MT_05C_TPM <- rename(MT_05C_TPM, GPM=TPM)
MT_05C <- merge(MT_05C, MT_05C_TPM, all=FALSE)

MT_06 <- read.table("CAZy_results/06_MT_CAZy.tab", sep="\t")
head(MT_06)
MT_06 <- cSplit(MT_06, "V1", ":")
MT_06 <- rename(MT_06, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_06 <- subset(MT_06, select=-c(V6, V7, V8, V9, V10, V12))
MT_06 <- subset(MT_06, align_length >=30)
MT_06$type <- "metatranscriptome"
MT_06$time <- 6
MT_06_TPM <- read.csv("metatranscriptome_tpm/06_MT_TPM.csv")
head(MT_06_TPM)
MT_06_TPM <- rename(MT_06_TPM, GPM=TPM)
MT_06 <- merge(MT_06, MT_06_TPM, all=FALSE)

MT_07 <- read.table("CAZy_results/07_MT_CAZy.tab", sep="\t")
head(MT_07)
MT_07 <- cSplit(MT_07, "V1", ":")
MT_07 <- rename(MT_07, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_07 <- subset(MT_07, select=-c(V6, V7, V8, V9, V10, V12))
MT_07 <- subset(MT_07, align_length >=30)
MT_07$type <- "metatranscriptome"
MT_07$time <- 7
MT_07_TPM <- read.csv("metatranscriptome_tpm/07_MT_TPM.csv")
head(MT_07_TPM)
MT_07_TPM <- rename(MT_07_TPM, GPM=TPM)
MT_07 <- merge(MT_07, MT_07_TPM, all=FALSE)

MT_08 <- read.table("CAZy_results/08_MT_CAZy.tab", sep="\t")
head(MT_08)
MT_08 <- cSplit(MT_08, "V1", ":")
MT_08 <- rename(MT_08, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_08 <- subset(MT_08, select=-c(V6, V7, V8, V9, V10, V12))
MT_08 <- subset(MT_08, align_length >=30)
MT_08$type <- "metatranscriptome"
MT_08$time <- 8
MT_08_TPM <- read.csv("metatranscriptome_tpm/08_MT_TPM.csv")
head(MT_08_TPM)
MT_08_TPM <- rename(MT_08_TPM, GPM=TPM)
MT_08 <- merge(MT_08, MT_08_TPM, all=FALSE)

MT_09A <- read.table("CAZy_results/09A_MT_CAZy.tab", sep="\t")
head(MT_09A)
MT_09A <- cSplit(MT_09A, "V1", ":")
MT_09A <- rename(MT_09A, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_09A <- subset(MT_09A, select=-c(V6, V7, V8, V9, V10, V12))
MT_09A <- subset(MT_09A, align_length >=30)
MT_09A$type <- "metatranscriptome"
MT_09A$time <- 9
MT_09A_TPM <- read.csv("metatranscriptome_tpm/09A_MT_TPM.csv")
head(MT_09A_TPM)
MT_09A_TPM <- rename(MT_09A_TPM, GPM=TPM)
MT_09A <- merge(MT_09A, MT_09A_TPM, all=FALSE)

MT_10 <- read.table("CAZy_results/10_MT_CAZy.tab", sep="\t")
head(MT_10)
MT_10 <- cSplit(MT_10, "V1", ":")
MT_10 <- rename(MT_10, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_10 <- subset(MT_10, select=-c(V6, V7, V8, V9, V10, V12))
MT_10 <- subset(MT_10, align_length >=30)
MT_10$type <- "metatranscriptome"
MT_10$time <- 10
MT_10_TPM <- read.csv("metatranscriptome_tpm/10_MT_TPM.csv")
head(MT_10_TPM)
MT_10_TPM <- rename(MT_10_TPM, GPM=TPM)
MT_10 <- merge(MT_10, MT_10_TPM, all=FALSE)

MT_11 <- read.table("CAZy_results/11_MT_CAZy.tab", sep="\t")
head(MT_11)
MT_11 <- cSplit(MT_11, "V1", ":")
MT_11 <- rename(MT_11, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_11 <- subset(MT_11, select=-c(V6, V7, V8, V9, V10, V12))
MT_11 <- subset(MT_11, align_length >=30)
MT_11$type <- "metatranscriptome"
MT_11$time <-11
MT_11_TPM <- read.csv("metatranscriptome_tpm/11_MT_TPM.csv")
head(MT_11_TPM)
MT_11_TPM <- rename(MT_11_TPM, GPM=TPM)
MT_11 <- merge(MT_11, MT_11_TPM, all=FALSE)

MT_12A <- read.table("CAZy_results/12A_MT_CAZy.tab", sep="\t")
head(MT_12A)
MT_12A <- cSplit(MT_12A, "V1", ":")
MT_12A <- rename(MT_12A, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_12A <- subset(MT_12A, select=-c(V6, V7, V8, V9, V10, V12))
MT_12A <- subset(MT_12A, align_length >=30)
MT_12A$type <- "metatranscriptome"
MT_12A$time <- 12
MT_12A_TPM <- read.csv("metatranscriptome_tpm/12A_MT_TPM.csv")
head(MT_12A_TPM)
MT_12A_TPM <- rename(MT_12A_TPM, GPM=TPM)
MT_12A <- merge(MT_12A, MT_12A_TPM, all=FALSE)

MT_12B <- read.table("CAZy_results/12B_MT_CAZy.tab", sep="\t")
head(MT_12B)
MT_12B <- cSplit(MT_12B, "V1", ":")
MT_12B <- rename(MT_12B, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_12B <- subset(MT_12B, select=-c(V6, V7, V8, V9, V10, V12))
MT_12B <- subset(MT_12B, align_length >=30)
MT_12B$type <- "metatranscriptome"
MT_12B$time <- 12.5
MT_12B_TPM <- read.csv("metatranscriptome_tpm/12B_MT_TPM.csv")
head(MT_12B_TPM)
MT_12B_TPM <- rename(MT_12B_TPM, GPM=TPM)
MT_12B <- merge(MT_12B, MT_12B_TPM, all=FALSE)

MT_13A <- read.table("CAZy_results/13A_MT_CAZy.tab", sep="\t")
head(MT_13A)
MT_13A <- cSplit(MT_13A, "V1", ":")
MT_13A <- rename(MT_13A, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_13A <- subset(MT_13A, select=-c(V6, V7, V8, V9, V10, V12))
MT_13A <- subset(MT_13A, align_length >=30)
MT_13A$type <- "metatranscriptome"
MT_13A$time <- 13
MT_13A_TPM <- read.csv("metatranscriptome_tpm/13A_MT_TPM.csv")
head(MT_13A_TPM)
MT_13A_TPM <- rename(MT_13A_TPM, GPM=TPM)
MT_13A <- merge(MT_13A, MT_13A_TPM, all=FALSE)

MT_13C <- read.table("CAZy_results/13C_MT_CAZy.tab", sep="\t")
head(MT_13C)
MT_13C <- cSplit(MT_13C, "V1", ":")
MT_13C <- rename(MT_13C, gene_id=V1_1, position=V1_2, CAZy=V2, percent_match=V3, align_length=V4, mismatch=V5, eval=V11)
MT_13C <- subset(MT_13C, select=-c(V6, V7, V8, V9, V10, V12))
MT_13C <- subset(MT_13C, align_length >=30)
MT_13C$type <- "metatranscriptome"
MT_13C$time <- 13.5
MT_13C_TPM <- read.csv("metatranscriptome_tpm/13C_MT_TPM.csv")
head(MT_13C_TPM)
MT_13C_TPM <- rename(MT_13C_TPM, GPM=TPM)
MT_13C <- merge(MT_13C, MT_13C_TPM, all=FALSE)

CAZy <- rbind(MG_01, MG_02, MG_03, MG_04, MG_05A, MG_05C, MG_06, MG_07, MG_09A, MG_09B, MG_10, MG_11, MG_12A, MG_12B, MG_13A, MG_13B, MT_01, MT_02, MT_03, MT_05B,
              MT_05C, MT_06, MT_07, MT_08, MT_09A, MT_10, MT_11, MT_12A, MT_12B, MT_13A, MT_13C)
head(CAZy)

CAZy <- read.csv("CAZy_GPM.csv")
head(CAZy)
CAZy <- rbind(CAZy, MT_13C, fill=TRUE)
write.csv(CAZy, "CAZy_GPM.csv", row.names=FALSE)
