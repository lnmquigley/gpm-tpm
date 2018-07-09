require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_07 <- read.table("KO/071407_MG_assembled_derep.ko.txt")
reads_MG_07 <- read.table("covstats/07_mg_covstats.txt")
rl_MG_07 <- read.table("rl/07_MG_rl_sort.txt")
head(KO_MG_07)
head(reads_MG_07)
head(rl_MG_07)

####format rl data###
rl_MG_07 <- rename(rl_MG_07, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_07)
####split gene_id into gene_id and position####
rl_MG_07 <- cSplit(rl_MG_07, "gene_id", ":")
rl_MG_07 = rename(rl_MG_07, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_07 <- subset(rl_MG_07, gene_id!="*")
head(rl_MG_07)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_07 <-cSplit(reads_MG_07, "V1", ":")
head(reads_MG_07)
reads_MG_07$rg <-(reads_MG_07$V7+reads_MG_07$V8)
reads_MG_07 <-subset(reads_MG_07, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_07 = rename(reads_MG_07, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_07 <- reads_MG_07 %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_07, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_07)
rl_MG_07 <- as.data.frame(rl_MG_07)
reads_MG_07 <- merge(reads_MG_07, rl_MG_07, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_07, n=15)
tail(reads_MG_07, n=15)
reads_MG_07 <- subset(reads_MG_07, select= -c(rg.x))
reads_MG_07 = rename(reads_MG_07, rg=rg.y)
tail(reads_MG_07, n=15)
reads_MG_07[is.na(reads_MG_07)] <- 0
###calculate T###
G <- ((sum(reads_MG_07$rg)*sum(reads_MG_07$rl))/sum(reads_MG_07$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_07, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071407_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071407_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_07)
tail(KO_MG_07)
str(KO_MG_07)
KO_MG_07 <- rename(KO_MG_07, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_07 <- subset(KO_MG_07, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_07)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_07$CDS <- str_sub(KO_MG_07$gene_id, 19, 25)
KO_MG_07$gene_id <- substr(KO_MG_07$gene_id, 0, 18)
head(KO_MG_07)
tail(KO_MG_07)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_O7_KO <- merge(GPM, KO_MG_07, by=c("gene_id", "CDS"), all=T)
head(MG_O7_KO, n=100)
tail(MG_O7_KO, n=100)

phylo <- read.table("phylodist/07_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_07 <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_07 <- subset(phylo_MG_07, select= -c(V2, V3))
phylo_MG_07$CDS <- str_sub(phylo_MG_07$gene_id, 19, 25)
phylo_MG_07$gene_id <- substr(phylo_MG_07$gene_id, 0, 18)

MG_07_GPM <- merge(phylo_MG_07, MG_O7_KO, by=c("gene_id", "CDS"), all=T)
head(MG_07_GPM)
tail(MG_07_GPM)

write.csv(MG_07_GPM, "metagenome_gpm/07_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/07_MG_GPM.csv")
ec <- read.table("ec/07_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/07_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_07 <- read.csv("metagenome_gpm/07_MG_GPM.csv")
GPM_07_sep <- separate(GPM_07, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                       sep=";", remove=TRUE)
head(GPM_07_sep)
Tthermo_07 <- subset(GPM_07_sep, species=="Thermus thermophilus")
GPM_07_sep <- subset(GPM_07_sep, species!="Thermus thermophilus" | is.na(species))

GPM_07_sep$taxonomy <- paste(GPM_07_sep$domain, GPM_07_sep$phylum, GPM_07_sep$class, GPM_07_sep$order, 
                             GPM_07_sep$family, GPM_07_sep$genus, GPM_07_sep$species, GPM_07_sep$species, sep=";")
GPM_07_sep <- subset(GPM_07_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_07, "fucking_thermo/07_thermo.csv", row.names = FALSE)

GPM_07_sep <- subset(GPM_07_sep, select=-c(GPM))
G7 <- ((sum(GPM_07_sep$rg)*sum(GPM_07_sep$rl))/sum(GPM_07_sep$flg))
GPM_07_sep <- transform(GPM_07_sep, GPM=(rg*rl*1e6)/(flg*G7))
write.csv(GPM_07_sep, "metagenome_gpm/07_MG_GPM_mod.csv", row.names=FALSE)

GPM_07 <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/07_MG_GPM_mod.csv")
GPM_07_ko <- subset(GPM_07, KO_num !="<NA>")
GPM_07_ko <- subset(GPM_07_ko, select=-c(GPM, EC, percent_id_ec))
G7 <- ((sum(GPM_07_ko$rg)*sum(GPM_07_ko$rl))/sum(GPM_07_ko$flg))
G7
GPM_07_ko <- transform(GPM_07_ko, GPM=(rg*rl*1e6)/(flg*G7))
GPM_07_ko_75 <- subset(GPM_07_ko, percent_id_tax >=75)
write.csv(GPM_07_ko_75, "metagenome_gpm/07_MG_GPM.csv", row.names=FALSE)
