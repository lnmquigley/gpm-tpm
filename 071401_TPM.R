
require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MT_01 <- read.table("KO/071401_MT_assembled_derep_ko.txt")
reads_MT_01 <- read.table("covstats/01_MT_covstats.txt")
rl_MT_01 <- read.table("rl/01_MT_rl.txt")
head(KO_MT_01)
head(reads_MT_01)
head(rl_MT_01)

####format rl data###
rl_MT_01 <- rename(rl_MT_01, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MT_01)
####split gene_id into gene_id and position####
rl_MT_01 <- cSplit(rl_MT_01, "gene_id", ":")
rl_MT_01 = rename(rl_MT_01, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MT_01 <- subset(rl_MT_01, gene_id!="*")
head(rl_MT_01)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MT_01<-cSplit(reads_MT_01, "V1", ":")
head(reads_MT_01)
reads_MT_01$rg <-(reads_MT_01$V7+reads_MT_01$V8)
reads_MT_01 <-subset(reads_MT_01, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MT_01 = rename(reads_MT_01, flg=V3, gene_id=V1_1, position=V1_2)
reads_MT_01 <- reads_MT_01 %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MT_01, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MT_01)
rl_MT_01 <- as.data.frame(rl_MT_01)
reads_MT_01 <- merge(reads_MT_01, rl_MT_01, by=c("gene_id", "position"), all=TRUE)
head(reads_MT_01, n=15)
tail(reads_MT_01, n=15)
reads_MT_01 <- subset(reads_MT_01, select= -c(rg.x))
reads_MT_01 = rename(reads_MT_01, rg=rg.y)
tail(reads_MT_01, n=15)
reads_MT_01[is.na(reads_MT_01)] <- 0
###calculate T###
T <- ((sum(reads_MT_01$rg)*sum(reads_MT_01$rl))/sum(reads_MT_01$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
TPM <- transform(reads_MT_01, TPM=(rg*rl*1e6)/(flg*T))
head(TPM)
tail(TPM)



write.csv(TPM, file = "metatranscriptome_tpm/071401_MT_TPM.csv",row.names=FALSE)
TPM <- read.csv("071401_MT_TPM.csv")
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MT_01)
tail(KO_MT_01)
str(KO_MT_01)
KO_MT_01 <- rename(KO_MT_01, gene_id=V1, KO_num=V3, e_val=V9)
KO_MT_01 <- subset(KO_MT_01, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MT_01)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MT_01$CDS <- str_sub(KO_MT_01$gene_id,18, 25)
KO_MT_01$gene_id <- substr(KO_MT_01$gene_id, 0, 17)
head(KO_MT_01)
tail(KO_MT_01)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####
####I really only care about TPM for scaffolds that have a KO function assigned to them###
####After an hour of trying to merge the KO df and the TPM df I finally figured out the issue###
####Prepare yourself####
###HONESTLY JUST FUCK JGI WHY THE FUCK WOULD YOU GIVE SOMEONE WITH GENE_ID NUMBERS THAT###
###DON'T MATCH UP FOR THE SAME GODDAMN SEQUENCING PROJECT AND WHY THE FUCK WOULD YOU ####
####USE TWO DIFFERENT ANALYSIS PIPELINES ON FUCKING DATASETS BELONGING TO THE SAME####
####FUCKING PROJECT LIKE HONESTLY WHAT THE FUCK IS WRONG WITH THESE PEOPLE####

MT_01_KO <- merge(KO_MT_01, TPM, by=c("gene_id", "CDS"), all.x=T)
head(MT_01_KO, n=100)
tail(MT_01_KO, n=100)
write.csv(MT_01_KO, file = "071401_MT_TPM_KO.csv",row.names=FALSE)

####want to combine phylodist with KO and gene_id to create master spreadsheet for sample##
phylo_MT_01 <- read.table("phylodist/01_MT_assembled.phylodist.txt", sep="\t")
head(phylo_MT_01)
phylo_MT_01 <- rename(phylo_MT_01, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MT_01 <- subset(phylo_MT_01, select = -c(V2, V3))
phylo_MT_01$CDS <- str_sub(phylo_MT_01$gene_id, 18, 25)
phylo_MT_01$gene_id <- substr(phylo_MT_01$gene_id, 0, 17)
MT_01 <- merge(MT_01_KO, phylo_MT_01, by=c("gene_id", "CDS"), all.x=T)
head(MT_01)

write.csv(MT_01, file="metatranscriptome_tpm/01_MT_TPM.csv", row.names=FALSE)

