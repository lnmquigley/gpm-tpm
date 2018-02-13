require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_01 <- read.table("KO/071401_MG_assembled_derep.ko.txt")
reads_MG_01 <- read.table("covstats/01_mg_covstats.txt")
rl_MG_01 <- read.table("rl/01_MG_rl_sort.txt")
head(KO_MG_01)
head(reads_MG_01)
head(rl_MG_01)

####format rl data###
rl_MG_01 <- rename(rl_MG_01, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_01)
####split gene_id into gene_id and position####
rl_MG_01 <- cSplit(rl_MG_01, "gene_id", ":")
rl_MG_01 = rename(rl_MG_01, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_01 <- subset(rl_MG_01, gene_id!="*")
head(rl_MG_01)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_01<-cSplit(reads_MG_01, "V1", ":")
head(reads_MG_01)
reads_MG_01$rg <-(reads_MG_01$V7+reads_MG_01$V8)
reads_MG_01 <-subset(reads_MG_01, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_01 = rename(reads_MG_01, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_01 <- reads_MG_01 %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_01, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_01)
rl_MG_01 <- as.data.frame(rl_MG_01)
reads_MG_01 <- merge(reads_MG_01, rl_MG_01, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_01, n=15)
tail(reads_MG_01, n=15)
reads_MG_01 <- subset(reads_MG_01, select= -c(rg.x))
reads_MG_01 = rename(reads_MG_01, rg=rg.y)
tail(reads_MG_01, n=15)
reads_MG_01[is.na(reads_MG_01)] <- 0
###calculate T###
G <- ((sum(reads_MG_01$rg)*sum(reads_MG_01$rl))/sum(reads_MG_01$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_01, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071401_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071401_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_01)
tail(KO_MG_01)
str(KO_MG_01)
KO_MG_01 <- rename(KO_MG_01, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_01 <- subset(KO_MG_01, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_01)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_01$CDS <- str_sub(KO_MG_01$gene_id, 19, 25)
KO_MG_01$gene_id <- substr(KO_MG_01$gene_id, 0, 18)
head(KO_MG_01)
tail(KO_MG_01)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_O1_KO <- merge(KO_MG_01, GPM, by=c("gene_id", "CDS"), all.x=T)
head(MG_O1_KO, n=100)
tail(MG_O1_KO, n=100)

phylo <- read.table("phylodist/01_MG_assembled.phylodist.txt", sep="\t")
head(phylo)
head(phylo_MG_01)
phylo_MG_01 <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_01 <- subset(phylo_MG_01, select= -c(V2, V3))
phylo_MG_01$CDS <- str_sub(phylo_MG_01$gene_id, 19, 25)
phylo_MG_01$gene_id <- substr(phylo_MG_01$gene_id, 0, 18)

MG_01_GPM <- merge(GPM, phylo_MG_01, by=c("gene_id", "CDS"))
MG_01_GPM <- merge(phylo_MG_01, MG_O1_KO, by=c("gene_id", "CDS"))
head(MG_01_GPM)
tail(MG_01_GPM)

write.csv(MG_01_GPM, "metagenome_gpm/01_MG_GPM.csv", row.names = FALSE)
