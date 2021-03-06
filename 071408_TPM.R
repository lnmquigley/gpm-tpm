require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MT_08 <- read.table("KO/071408_MT_assembled_derep.ko.txt")
reads_MT_08 <- read.table("covstats/08_mt_covstats.txt")
rl_MT_08 <- read.table("rl/08_MT_rl_sort.txt")
head(KO_MT_08)
head(reads_MT_08)
head(rl_MT_08)

####format rl data###
rl_MT_08 <- rename(rl_MT_08, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MT_08)
####split gene_id into gene_id and position####
rl_MT_08 <- cSplit(rl_MT_08, "gene_id", ":")
rl_MT_08 = rename(rl_MT_08, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MT_08 <- subset(rl_MT_08, gene_id!="*")
head(rl_MT_08)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MT_08 <-cSplit(reads_MT_08, "V1", ":")
head(reads_MT_08)
reads_MT_08$rg <-(reads_MT_08$V7+reads_MT_08$V8)
reads_MT_08 <-subset(reads_MT_08, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MT_08 = rename(reads_MT_08, flg=V3, gene_id=V1_1, position=V1_2)
reads_MT_08 <- reads_MT_08 %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MT_08, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MT_08)
rl_MT_08 <- as.data.frame(rl_MT_08)
reads_MT_08 <- merge(reads_MT_08, rl_MT_08, by=c("gene_id", "position"), all=TRUE)
head(reads_MT_08, n=15)
tail(reads_MT_08, n=15)
reads_MT_08 <- subset(reads_MT_08, select= -c(rg.x))
reads_MT_08 = rename(reads_MT_08, rg=rg.y)
tail(reads_MT_08, n=15)
reads_MT_08[is.na(reads_MT_08)] <- 0
###calculate T###
G <- ((sum(reads_MT_08$rg)*sum(reads_MT_08$rl))/sum(reads_MT_08$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
TPM <- transform(reads_MT_08, TPM=(rg*rl*1e6)/(flg*G))
head(TPM)
tail(TPM)



write.csv(TPM, file = "071408_MT_TPM.csv",row.names=FALSE)
TPM <- read.csv("metatranscriptome_tpm/071408_MT_TPM.csv")
dim(TPM)
head(TPM)
tail(TPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MT_08)
tail(KO_MT_08)
str(KO_MT_08)
KO_MT_08 <- rename(KO_MT_08, gene_id=V1, KO_num=V3, e_val=V9)
KO_MT_08 <- subset(KO_MT_08, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MT_08)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MT_08$CDS <- str_sub(KO_MT_08$gene_id, 18, 25)
KO_MT_08$gene_id <- substr(KO_MT_08$gene_id, 0, 17)
head(KO_MT_08)
tail(KO_MT_08)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MT_O8_KO <- merge(TPM, KO_MT_08, by=c("gene_id", "CDS"), all=T)
head(MT_O8_KO, n=100)
tail(MT_O8_KO, n=100)

phylo <- read.table("phylodist/08_MT_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MT_08 <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MT_08 <- subset(phylo_MT_08, select= -c(V2, V3))
phylo_MT_08$CDS <- str_sub(phylo_MT_08$gene_id, 18, 25)
phylo_MT_08$gene_id <- substr(phylo_MT_08$gene_id, 0, 17)

MT_08_TPM <- merge(phylo_MT_08, MT_O8_KO, by=c("gene_id", "CDS"), all=T)
head(MT_08_TPM)
tail(MT_08_TPM)

write.csv(MT_08_TPM, "metatranscriptome_tpm/08_MT_TPM.csv", row.names = FALSE)

TPM <- read.csv("metatranscriptome_tpm/08_MT_TPM.csv")
ec <- read.table("ec/08_MT_assembled.ec.txt")
head(TPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 18, 25)
ec$gene_id <- substr(ec$gene_id, 0, 17)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
TPM <- merge(MT_08_TPM, ec, by=c("gene_id", "CDS"), all=T)
head(TPM)
write.csv(TPM, "metatranscriptome_tpm/08_MT_TPM.csv", row.names=FALSE)

TPM_08 <- read.csv("metatranscriptome_tpm/V1/08_MT_TPM.csv")
head(TPM_08)
TPM_08_ko <- subset(TPM_08, KO_num !="<NA>")
TPM_08_ko <- subset(TPM_08_ko, select=-c(TPM, EC, percent_id_ec))
T8 <- ((sum(TPM_08_ko$rg)*sum(TPM_08_ko$rl))/sum(TPM_08_ko$flg))
T8
TPM_08_ko <- transform(TPM_08_ko, TPM=(rg*rl*1e6)/(flg*T8))
TPM_08_ko_75 <- subset(TPM_08_ko, percent_id_tax >=75)
write.csv(TPM_08_ko_75, "metatranscriptome_tpm/08_MT_TPM.csv", row.names=FALSE)
