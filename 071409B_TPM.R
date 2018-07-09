require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MT_09B <- read.table("KO/071409B_MT_assembled_derep.ko.txt")
reads_MT_09B <- read.table("covstats/09B_mt_covstats.txt")
rl_MT_09B <- read.table("rl/09B_MT_rl_sort.txt")
head(KO_MT_09B)
head(reads_MT_09B)
head(rl_MT_09B)

####format rl data###
rl_MT_09B <- rename(rl_MT_09B, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MT_09B)
####split gene_id into gene_id and position####
rl_MT_09B <- cSplit(rl_MT_09B, "gene_id", ":")
rl_MT_09B = rename(rl_MT_09B, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MT_09B <- subset(rl_MT_09B, gene_id!="*")
head(rl_MT_09B)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MT_09B <-cSplit(reads_MT_09B, "V1", ":")
head(reads_MT_09B)
reads_MT_09B$rg <-(reads_MT_09B$V7+reads_MT_09B$V8)
reads_MT_09B <-subset(reads_MT_09B, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MT_09B = rename(reads_MT_09B, flg=V3, gene_id=V1_1, position=V1_2)
reads_MT_09B <- reads_MT_09B %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MT_09B, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MT_09B)
rl_MT_09B <- as.data.frame(rl_MT_09B)
reads_MT_09B <- merge(reads_MT_09B, rl_MT_09B, by=c("gene_id", "position"), all=TRUE)
head(reads_MT_09B, n=15)
tail(reads_MT_09B, n=15)
reads_MT_09B <- subset(reads_MT_09B, select= -c(rg.x))
reads_MT_09B = rename(reads_MT_09B, rg=rg.y)
tail(reads_MT_09B, n=15)
reads_MT_09B[is.na(reads_MT_09B)] <- 0
###calculate T###
G <- ((sum(reads_MT_09B$rg)*sum(reads_MT_09B$rl))/sum(reads_MT_09B$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
TPM <- transform(reads_MT_09B, TPM=(rg*rl*1e6)/(flg*G))
head(TPM)
tail(TPM)

write.csv(TPM, file = "071409B_MT_TPM.csv",row.names=FALSE)
TPM <- read.csv("metatranscriptome_tpm/071409B_MT_TPM.csv")
dim(TPM)
head(TPM)
tail(TPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MT_09B)
tail(KO_MT_09B)
str(KO_MT_09B)
KO_MT_09B <- rename(KO_MT_09B, gene_id=V1, KO_num=V3, e_val=V9)
KO_MT_09B <- subset(KO_MT_09B, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MT_09B)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MT_09B$CDS <- str_sub(KO_MT_09B$gene_id, 19, 25)
KO_MT_09B$gene_id <- substr(KO_MT_09B$gene_id, 0, 18)
head(KO_MT_09B)
tail(KO_MT_09B)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MT_O9B_KO <- merge(TPM, KO_MT_09B, by=c("gene_id", "CDS"), all=T)
head(MT_O9B_KO, n=100)
tail(MT_O9B_KO, n=100)

phylo <- read.table("phylodist/09B_MT_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MT_09B <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MT_09B <- subset(phylo_MT_09B, select= -c(V2, V3))
phylo_MT_09B$CDS <- str_sub(phylo_MT_09B$gene_id, 19, 25)
phylo_MT_09B$gene_id <- substr(phylo_MT_09B$gene_id, 0, 18)

MT_09B_TPM <- merge(phylo_MT_09B, MT_O9B_KO, by=c("gene_id", "CDS"), all=T)
head(MT_09B_TPM)
tail(MT_09B_TPM)

write.csv(MT_09B_TPM, "metatranscriptome_tpm/09B_MT_TPM.csv", row.names = FALSE)

TPM <- read.csv("metatranscriptome_tpm/09B_MT_TPM.csv")
ec <- read.table("ec/09B_MT_assembled.ec.txt")
head(TPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
TPM <- merge(MT_09B_TPM, ec, by=c("gene_id", "CDS"), all=T)
head(TPM)
write.csv(TPM, "metatranscriptome_tpm/09B_MT_TPM.csv", row.names=FALSE)

