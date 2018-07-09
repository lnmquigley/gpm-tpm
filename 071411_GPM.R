require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_11 <- read.table("KO/071411_MG_assembled_derep.ko.txt")
reads_MG_11 <- read.table("covstats/11_mg_covstats.txt")
rl_MG_11 <- read.table("rl/11_MG_rl_sort.txt")
head(KO_MG_11)
head(reads_MG_11)
head(rl_MG_11)

####format rl data###
rl_MG_11 <- rename(rl_MG_11, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_11)
####split gene_id into gene_id and position####
rl_MG_11 <- cSplit(rl_MG_11, "gene_id", ":")
rl_MG_11 = rename(rl_MG_11, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_11 <- subset(rl_MG_11, gene_id!="*")
head(rl_MG_11)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_11 <-cSplit(reads_MG_11, "V1", ":")
head(reads_MG_11)
reads_MG_11$rg <-(reads_MG_11$V7+reads_MG_11$V8)
reads_MG_11 <-subset(reads_MG_11, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_11 = rename(reads_MG_11, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_11 <- reads_MG_11 %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_11, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_11)
rl_MG_11 <- as.data.frame(rl_MG_11)
reads_MG_11 <- merge(reads_MG_11, rl_MG_11, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_11, n=15)
tail(reads_MG_11, n=15)
reads_MG_11 <- subset(reads_MG_11, select= -c(rg.x))
reads_MG_11 = rename(reads_MG_11, rg=rg.y)
tail(reads_MG_11, n=15)
reads_MG_11[is.na(reads_MG_11)] <- 0
###calculate T###
G <- ((sum(reads_MG_11$rg)*sum(reads_MG_11$rl))/sum(reads_MG_11$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_11, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071411_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071411_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_11)
tail(KO_MG_11)
str(KO_MG_11)
KO_MG_11 <- rename(KO_MG_11, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_11 <- subset(KO_MG_11, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_11)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_11$CDS <- str_sub(KO_MG_11$gene_id, 19, 25)
KO_MG_11$gene_id <- substr(KO_MG_11$gene_id, 0, 18)
head(KO_MG_11)
tail(KO_MG_11)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_11_KO <- merge(GPM, KO_MG_11, by=c("gene_id", "CDS"), all=T)
head(MG_11_KO, n=100)
tail(MG_11_KO, n=100)

phylo <- read.table("phylodist/11_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_11 <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_11 <- subset(phylo_MG_11, select= -c(V2, V3))
phylo_MG_11$CDS <- str_sub(phylo_MG_11$gene_id, 19, 25)
phylo_MG_11$gene_id <- substr(phylo_MG_11$gene_id, 0, 18)

MG_11_GPM <- merge(phylo_MG_11, MG_11_KO, by=c("gene_id", "CDS"), all=T)
head(MG_11_GPM)
tail(MG_11_GPM)

write.csv(MG_11_GPM, "metagenome_gpm/11_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/11_MG_GPM.csv")
ec <- read.table("ec/11_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/11_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_11 <- read.csv("metagenome_gpm/11_MG_GPM.csv")
GPM_11_sep <- separate(GPM_11, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                       sep=";", remove=TRUE)
head(GPM_11_sep)
Tthermo_11 <- subset(GPM_11_sep, species=="Thermus thermophilus")
GPM_11_sep <- subset(GPM_11_sep, species!="Thermus thermophilus" | is.na(species))

GPM_11_sep$taxonomy <- paste(GPM_11_sep$domain, GPM_11_sep$phylum, GPM_11_sep$class, GPM_11_sep$order, 
                             GPM_11_sep$family, GPM_11_sep$genus, GPM_11_sep$species, GPM_11_sep$species, sep=";")
GPM_11_sep <- subset(GPM_11_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_11, "fucking_thermo/11_thermo.csv", row.names = FALSE)

GPM_11_sep <- subset(GPM_11_sep, select=-c(GPM))
G11 <- ((sum(GPM_11_sep$rg)*sum(GPM_11_sep$rl))/sum(GPM_11_sep$flg))
GPM_11_sep <- transform(GPM_11_sep, GPM=(rg*rl*1e6)/(flg*G11))
write.csv(GPM_11_sep, "metagenome_gpm/11_MG_GPM_mod.csv", row.names=FALSE)

GPM_11 <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/11_MG_GPM_mod.csv")
GPM_11_ko <- subset(GPM_11, KO_num !="<NA>")
GPM_11_ko <- subset(GPM_11_ko, select=-c(GPM, EC, percent_id_ec))
G11 <- ((sum(GPM_11_ko$rg)*sum(GPM_11_ko$rl))/sum(GPM_11_ko$flg))
G11
GPM_11_ko <- transform(GPM_11_ko, GPM=(rg*rl*1e6)/(flg*G11))
GPM_11_ko_75 <- subset(GPM_11_ko, percent_id_tax >=75)
write.csv(GPM_11_ko_75, "metagenome_gpm/11_MG_GPM.csv", row.names=FALSE)
