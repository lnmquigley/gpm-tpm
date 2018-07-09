require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MT_09A <- read.table("KO/071409A_MT_assembled_derep.ko.txt")
reads_MT_09A <- read.table("covstats/09A_mt_covstats.txt")
rl_MT_09A <- read.table("rl/09A_MT_rl_sort.txt")
head(KO_MT_09A)
head(reads_MT_09A)
head(rl_MT_09A)

####format rl data###
rl_MT_09A <- rename(rl_MT_09A, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MT_09A)
####split gene_id into gene_id and position####
rl_MT_09A <- cSplit(rl_MT_09A, "gene_id", ":")
rl_MT_09A = rename(rl_MT_09A, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MT_09A <- subset(rl_MT_09A, gene_id!="*")
head(rl_MT_09A)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MT_09A <-cSplit(reads_MT_09A, "V1", ":")
head(reads_MT_09A)
reads_MT_09A$rg <-(reads_MT_09A$V7+reads_MT_09A$V8)
reads_MT_09A <-subset(reads_MT_09A, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MT_09A = rename(reads_MT_09A, flg=V3, gene_id=V1_1, position=V1_2)
reads_MT_09A <- reads_MT_09A %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MT_09A, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MT_09A)
rl_MT_09A <- as.data.frame(rl_MT_09A)
reads_MT_09A <- merge(reads_MT_09A, rl_MT_09A, by=c("gene_id", "position"), all=TRUE)
head(reads_MT_09A, n=15)
tail(reads_MT_09A, n=15)
reads_MT_09A <- subset(reads_MT_09A, select= -c(rg.x))
reads_MT_09A = rename(reads_MT_09A, rg=rg.y)
tail(reads_MT_09A, n=15)
reads_MT_09A[is.na(reads_MT_09A)] <- 0
###calculate T###
G <- ((sum(reads_MT_09A$rg)*sum(reads_MT_09A$rl))/sum(reads_MT_09A$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
TPM <- transform(reads_MT_09A, TPM=(rg*rl*1e6)/(flg*G))
head(TPM)
tail(TPM)

write.csv(TPM, file = "071409A_MT_TPM.csv",row.names=FALSE)
TPM <- read.csv("metatranscriptome_tpm/071409A_MT_TPM.csv")
dim(TPM)
head(TPM)
tail(TPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MT_09A)
tail(KO_MT_09A)
str(KO_MT_09A)
KO_MT_09A <- rename(KO_MT_09A, gene_id=V1, KO_num=V3, e_val=V9)
KO_MT_09A <- subset(KO_MT_09A, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MT_09A)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MT_09A$CDS <- str_sub(KO_MT_09A$gene_id, 18, 25)
KO_MT_09A$gene_id <- substr(KO_MT_09A$gene_id, 0, 17)
head(KO_MT_09A)
tail(KO_MT_09A)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MT_O9A_KO <- merge(TPM, KO_MT_09A, by=c("gene_id", "CDS"), all=T)
head(MT_O9A_KO, n=100)
tail(MT_O9A_KO, n=100)

phylo <- read.table("phylodist/09A_MT_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MT_09A <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MT_09A <- subset(phylo_MT_09A, select= -c(V2, V3))
phylo_MT_09A$CDS <- str_sub(phylo_MT_09A$gene_id, 18, 25)
phylo_MT_09A$gene_id <- substr(phylo_MT_09A$gene_id, 0, 17)

MT_09A_TPM <- merge(phylo_MT_09A, MT_O9A_KO, by=c("gene_id", "CDS"), all=T)
head(MT_09A_TPM)
tail(MT_09A_TPM)

write.csv(MT_09A_TPM, "metatranscriptome_tpm/09A_MT_TPM.csv", row.names = FALSE)

TPM <- read.csv("metatranscriptome_tpm/09A_MT_TPM.csv")
ec <- read.table("ec/09A_MT_assembled.ec.txt")
head(TPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 18, 25)
ec$gene_id <- substr(ec$gene_id, 0, 17)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
TPM <- merge(MT_09A_TPM, ec, by=c("gene_id", "CDS"), all=T)
head(TPM)
write.csv(TPM, "metatranscriptome_tpm/09A_MT_TPM.csv", row.names=FALSE)

TPM_09A <- read.csv("metatranscriptome_tpm/V1/09A_MT_TPM.csv")
head(TPM_09A)
TPM_09A_ko <- subset(TPM_09A, KO_num !="<NA>")
TPM_09A_ko <- subset(TPM_09A_ko, select=-c(TPM, EC, percent_id_ec))
T9A <- ((sum(TPM_09A_ko$rg)*sum(TPM_09A_ko$rl))/sum(TPM_09A_ko$flg))
T9A
TPM_09A_ko <- transform(TPM_09A_ko, TPM=(rg*rl*1e6)/(flg*T9A))
TPM_09A_ko_75 <- subset(TPM_09A_ko, percent_id_tax >=75)
write.csv(TPM_09A_ko_75, "metatranscriptome_tpm/09A_MT_TPM.csv", row.names=FALSE)
