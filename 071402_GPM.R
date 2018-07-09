require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_02 <- read.table("KO/071402_MG_assembled_derep.ko.txt")
reads_MG_02 <- read.table("covstats/02_mg_covstats.txt")
rl_MG_02 <- read.table("rl/02_MG_rl_sort.txt")
head(KO_MG_02)
head(reads_MG_02)
head(rl_MG_02)

####format rl data###
rl_MG_02 <- rename(rl_MG_02, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_02)
####split gene_id into gene_id and position####
rl_MG_02 <- cSplit(rl_MG_02, "gene_id", ":")
rl_MG_02 = rename(rl_MG_02, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_02 <- subset(rl_MG_02, gene_id!="*")
head(rl_MG_02)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_02<-cSplit(reads_MG_02, "V1", ":")
head(reads_MG_02)
reads_MG_02$rg <-(reads_MG_02$V7+reads_MG_02$V8)
reads_MG_02 <-subset(reads_MG_02, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_02 = rename(reads_MG_02, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_02 <- reads_MG_02 %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_02, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_02)
rl_MG_02 <- as.data.frame(rl_MG_02)
reads_MG_02 <- merge(reads_MG_02, rl_MG_02, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_02, n=15)
tail(reads_MG_02, n=15)
reads_MG_02 <- subset(reads_MG_02, select= -c(rg.x))
reads_MG_02 = rename(reads_MG_02, rg=rg.y)
tail(reads_MG_02, n=15)
reads_MG_02[is.na(reads_MG_02)] <- 0
###calculate T###
G <- ((sum(reads_MG_02$rg)*sum(reads_MG_02$rl))/sum(reads_MG_02$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_02, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071402_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071402_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_02)
tail(KO_MG_02)
str(KO_MG_02)
KO_MG_02 <- rename(KO_MG_02, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_02 <- subset(KO_MG_02, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_02)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_02$CDS <- str_sub(KO_MG_02$gene_id, 18, 25)
KO_MG_02$gene_id <- substr(KO_MG_02$gene_id, 0, 17)
head(KO_MG_02)
tail(KO_MG_02)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_O2_KO <- merge(GPM, KO_MG_02, by=c("gene_id", "CDS"), all=T)
head(MG_O2_KO, n=100)
tail(MG_O2_KO, n=100)

phylo <- read.table("phylodist/02_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_02 <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_02 <- subset(phylo_MG_02, select= -c(V2, V3))
phylo_MG_02$CDS <- str_sub(phylo_MG_02$gene_id, 18, 25)
phylo_MG_02$gene_id <- substr(phylo_MG_02$gene_id, 0, 17)

MG_02_GPM <- merge(phylo_MG_02, MG_O2_KO, by=c("gene_id", "CDS"), all=T)
head(MG_02_GPM)
tail(MG_02_GPM)

write.csv(MG_02_GPM, "metagenome_gpm/02_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/02_MG_GPM.csv")
ec <- read.table("ec/02_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 18, 25)
ec$gene_id <- substr(ec$gene_id, 0, 17)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/02_MG_GPM.csv", row.names=FALSE)

###remove leftover tthermo####
require(tidyr)
GPM_02 <- read.csv("metagenome_gpm/02_MG_GPM.csv")
GPM_02_sep <- separate(GPM_02, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                       sep=";", remove=TRUE)
head(GPM_02_sep)
Tthermo_02 <- subset(GPM_02_sep, species=="Thermus thermophilus")
GPM_02_sep <- subset(GPM_02_sep, species!="Thermus thermophilus" | is.na(species))

GPM_02_sep$taxonomy <- paste(GPM_02_sep$domain, GPM_02_sep$phylum, GPM_02_sep$class, GPM_02_sep$order, 
                             GPM_02_sep$family, GPM_02_sep$genus, GPM_02_sep$species, GPM_02_sep$species, sep=";")
GPM_02_sep <- subset(GPM_02_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(GPM_02_sep, "metagenome_gpm/02_MG_GPM_mod.csv", row.names=FALSE)
write.csv(Tthermo_02, "fucking_thermo/02_thermo.csv", row.names = FALSE)

GPM_02_sep <- subset(GPM_02_sep, select=-c(GPM))
G2 <- ((sum(GPM_02_sep$rg)*sum(GPM_02_sep$rl))/sum(GPM_02_sep$flg))
GPM_02_sep <- transform(GPM_02_sep, GPM=(rg*rl*1e6)/(flg*G2))
write.csv(GPM_02_sep, "metagenome_gpm/02_MG_GPM_mod.csv", row.names=FALSE)


GPM_02 <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/02_MG_GPM_mod.csv")
GPM_02_ko <- subset(GPM_02, KO_num !="<NA>")
GPM_02_ko <- subset(GPM_02_ko, select=-c(GPM, EC, percent_id_ec))
G2 <- ((sum(GPM_02_ko$rg)*sum(GPM_02_ko$rl))/sum(GPM_02_ko$flg))
G2
GPM_02_ko <- transform(GPM_02_ko, GPM=(rg*rl*1e6)/(flg*G2))
GPM_02_ko_75 <- subset(GPM_02_ko, percent_id_tax >=75)
write.csv(GPM_02_ko_75, "metagenome_gpm/02_MG_GPM.csv", row.names=FALSE)
