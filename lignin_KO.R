require(tidyr)
require(dplyr)
require(ggplot2)
####load data for all metagenomes and subset to get aromatic genes####
MG_01 <- read.csv("metagenome_gpm/01_MG_GPM.csv")

MG_01_lig <- subset(MG_01, KO_num=="KO:K15733")

MG_01_lig$timepoint <- 1
head(MG_01_lig)


MG_02 <- read.csv("metagenome_gpm/02_MG_GPM.csv")

MG_02_lig <- subset(MG_02, KO_num=="KO:K15733")

MG_02_lig$timepoint <- 2
head(MG_02_lig)

MG_03 <- read.csv("metagenome_gpm/03_MG_GPM.csv")

MG_03_lig <- subset(MG_03, KO_num=="KO:K15733")

MG_03_lig$timepoint <- 3
head(MG_03_lig)

MG_04 <- read.csv("metagenome_gpm/04_MG_GPM.csv")

MG_04_lig <- subset(MG_04, KO_num=="KO:K15733")

MG_04_lig$timepoint <- 4
head(MG_04_lig)

MG_05A <- read.csv("metagenome_gpm/05A_MG_GPM.csv")

MG_05A_lig <- subset(MG_05A, KO_num=="KO:K15733")

MG_05A_lig$timepoint <- 5
head(MG_05A_lig)

MG_05C <- read.csv("metagenome_gpm/05C_MG_GPM.csv")

MG_05C_lig <- subset(MG_05C, KO_num=="KO:K15733")

MG_05C_lig$timepoint <- 5.5
head(MG_05C_lig)

MG_06 <- read.csv("metagenome_gpm/06_MG_GPM.csv")

MG_06_lig <- subset(MG_06, KO_num=="KO:K15733")

MG_06_lig$timepoint <- 6
head(MG_06_lig)

MG_07 <- read.csv("metagenome_gpm/07_MG_GPM.csv")

MG_07_lig <- subset(MG_07, KO_num=="KO:K15733")
MG_07_lig$timepoint <- 7
head(MG_07_lig)

MG_09A <- read.csv("metagenome_gpm/09A_MG_GPM.csv")

MG_09A_lig <- subset(MG_09A, KO_num=="KO:K15733")
MG_09A_lig$timepoint <- 9
head(MG_09A_lig)

MG_09B <- read.csv("metagenome_gpm/09B_MG_GPM.csv")

MG_09B_lig <- subset(MG_09B, KO_num=="KO:K15733")
MG_09B_lig$timepoint <- 9.5
head(MG_09B_lig)

MG_10 <- read.csv("metagenome_gpm/10_MG_GPM.csv")

MG_10_lig <- subset(MG_10, KO_num=="KO:K15733")
MG_10_lig$timepoint <- 10
head(MG_10_lig)

MG_11 <- read.csv("metagenome_gpm/11_MG_GPM.csv")

MG_11_lig <- subset(MG_11, KO_num=="KO:K15733") 
MG_11_lig$timepoint <- 11
head(MG_11_lig)

MG_12A <- read.csv("metagenome_gpm/12A_MG_GPM.csv")

MG_12A_lig <- subset(MG_12A, KO_num=="KO:K15733") 
MG_12A_lig$timepoint <- 12
head(MG_12A_lig)

MG_12B <- read.csv("metagenome_gpm/12B_MG_GPM.csv")

MG_12B_lig <- subset(MG_12B, KO_num=="KO:K15733") 
MG_12B_lig$timepoint <- 12.5
head(MG_12B_lig)

MG_13A <- read.csv("metagenome_gpm/13A_MG_GPM.csv")

MG_13A_lig <- subset(MG_13A, KO_num=="KO:K15733")
MG_13A_lig$timepoint <- 13
head(MG_13A_lig)

MG_13B <- read.csv("metagenome_gpm/13B_MG_GPM.csv")

MG_13B_lig <- subset(MG_13B, KO_num=="KO:K15733")
MG_13B_lig$timepoint <- 13.5
head(MG_13B_lig)

#####combine all aromatic data frames into one big one####
lignin <- rbind(MG_01_lig, MG_02_lig, MG_03_lig, MG_04_lig, MG_05A_lig, MG_05C_lig, MG_06_lig, MG_07_lig, MG_09A_lig,
                  MG_09B_lig, MG_10_lig, MG_11_lig, MG_12A_lig, MG_12B_lig, MG_13A_lig, MG_13B_lig)
head(lignin)
tail(lignin)
?facet_grid

lignin$type <- "metagenome"

MT_01 <- read.csv("metatranscriptome_tpm/01_MT_TPM.csv")

MT_01_lig <- subset(MT_01, KO_num=="KO:K15733")
MT_01_lig$timepoint <- 1
head(MT_01_lig)


MT_02 <- read.csv("metatranscriptome_tpm/02_MT_TPM.csv")

MT_02_lig <- subset(MT_02, KO_num=="KO:K15733")
MT_02_lig$timepoint <- 2
head(MT_02_lig)

MT_03 <- read.csv("metatranscriptome_tpm/03_MT_TPM.csv")

MT_03_lig <- subset(MT_03, KO_num=="KO:K15733")
MT_03_lig$timepoint <- 3
head(MT_03_lig)


MT_05B <- read.csv("metatranscriptome_tpm/05B_MT_TPM.csv")

MT_05B_lig <- subset(MT_05B, KO_num=="KO:K15733")
MT_05B_lig$timepoint <- 5
head(MT_05B_lig)

MT_05C <- read.csv("metatranscriptome_tpm/05C_MT_TPM.csv")

MT_05C_lig <- subset(MT_05C, KO_num=="KO:K15733")
MT_05C_lig$timepoint <- 5.5
head(MT_05C_lig)

MT_06 <- read.csv("metatranscriptome_tpm/06_MT_TPM.csv")

MT_06_lig <- subset(MT_06, KO_num=="KO:K15733")
MT_06_lig$timepoint <- 6
head(MT_06_lig)

MT_07 <- read.csv("metatranscriptome_tpm/07_MT_TPM.csv")

MT_07_lig <- subset(MT_07, KO_num=="KO:K15733")
MT_07_lig$timepoint <- 7
head(MT_07_lig)

MT_08 <- read.csv("metatranscriptome_tpm/08_MT_TPM.csv")

MT_08_lig <- subset(MT_08, KO_num=="KO:K15733")
MT_08_lig$timepoint <- 8
head(MT_08_lig)

MT_09A <- read.csv("metatranscriptome_tpm/09A_MT_TPM.csv")

MT_09A_lig <- subset(MT_09A, KO_num=="KO:K15733")
MT_09A_lig$timepoint <- 9
head(MT_09A_lig)

MT_10 <- read.csv("metatranscriptome_tpm/10_MT_TPM.csv")

MT_10_lig <- subset(MT_10, KO_num=="KO:K15733")
MT_10_lig$timepoint <- 10
head(MT_10_lig)

MT_11 <- read.csv("metatranscriptome_tpm/11_MT_TPM.csv")

MT_11_lig <- subset(MT_11, KO_num=="KO:K15733")
MT_11_lig$timepoint <- 11
head(MT_11_lig)

MT_12A <- read.csv("metatranscriptome_tpm/12A_MT_TPM.csv")

MT_12A_lig <- subset(MT_12A, KO_num=="KO:K15733")
MT_12A_lig$timepoint <- 12
head(MT_12A_lig)

MT_12B <- read.csv("metatranscriptome_tpm/12B_MT_TPM.csv")

MT_12B_lig <- subset(MT_12B, KO_num=="KO:K15733") 
MT_12B_lig$timepoint <- 12.5
head(MT_12B_lig)

MT_13A <- read.csv("metatranscriptome_tpm/13A_MT_TPM.csv")

MT_13A_lig <- subset(MT_13A, KO_num=="KO:K15733")
MT_13A_lig$timepoint <- 13
head(MT_13A_lig)

MT_13C <- read.csv("metatranscriptome_tpm/13C_MT_TPM.csv")

MT_13C_lig <- subset(MT_13C, KO_num=="KO:K15733")
MT_13C_lig$timepoint <- 13.5
head(MT_13C_lig)

#####combine all aromatic data frames into one big one####
lignin_mt <- rbind(MT_01_lig, MT_02_lig, MT_03_lig, MT_05B_lig, MT_05C_lig, MT_06_lig, MT_07_lig, MT_08_lig,
                     MT_09A_lig, MT_10_lig, MT_11_lig, MT_12A_lig, MT_12B_lig, MT_13A_lig, MT_13C_lig)
head(lignin_mt)

lignin_mt$type <- "metatranscriptome"

lignin_mt <- rename(lignin_mt, GPM=TPM)

lignin <- rbind(lignin, lignin_mt)
head(lignin)
tail(lignin)

lignin_agg <- aggregate(aromatic$GPM, by=list(aromatic$timepoint, aromatic$type), FUN=sum)
head(lignin_agg)

lignin_agg$Group.1 <- as.numeric(revalue(as.character(lignin_agg$Group.1), c("5.5" = "5", "9.5" = "9", "12.5"="12", "13.5"="13")))

lignin_agg <- rename(lignin_agg, Time='Group.1', Type='Group.2', Count='x')
head(lignin_agg)

lignin_agg$Type <- revalue(lignin_agg$Type, c('metagenome'='Metagenome', 'metatranscriptome'='Metatranscriptome'))

lignin_p <- ggplot(lignin_agg, aes(Time, Count, colour=Type)) +
  geom_point(aes(shape=Type), size=8) +
  xlab("Time Point") +
  ylab("Genes or Transcripts per Million (GPM or TPM)") +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_color_manual(values=c("springgreen", "darkgreen")) +
  theme_classic() +
  theme(text=element_text(size=20)) +
  theme(legend.position="bottom")
print(lignin_p)
