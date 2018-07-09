require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_03 <- read.table("KO/071403_MG_assembled_derep.ko.txt")
reads_MG_03 <- read.table("covstats/03_mg_covstats.txt")
rl_MG_03 <- read.table("rl/03_MG_rl_sort.txt")
head(KO_MG_03)
head(reads_MG_03)
head(rl_MG_03)

####format rl data###
rl_MG_03 <- rename(rl_MG_03, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_03)
####split gene_id into gene_id and position####
rl_MG_03 <- cSplit(rl_MG_03, "gene_id", ":")
rl_MG_03 = rename(rl_MG_03, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_03 <- subset(rl_MG_03, gene_id!="*")
head(rl_MG_03)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_03<-cSplit(reads_MG_03, "V1", ":")
head(reads_MG_03)
reads_MG_03$rg <-(reads_MG_03$V7+reads_MG_03$V8)
reads_MG_03 <-subset(reads_MG_03, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_03 = rename(reads_MG_03, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_03 <- reads_MG_03 %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_03, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_03)
rl_MG_03 <- as.data.frame(rl_MG_03)
reads_MG_03 <- merge(reads_MG_03, rl_MG_03, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_03, n=15)
tail(reads_MG_03, n=15)
reads_MG_03 <- subset(reads_MG_03, select= -c(rg.x))
reads_MG_03 = rename(reads_MG_03, rg=rg.y)
tail(reads_MG_03, n=15)
reads_MG_03[is.na(reads_MG_03)] <- 0
###calculate T###
G <- ((sum(reads_MG_03$rg)*sum(reads_MG_03$rl))/sum(reads_MG_03$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_03, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071403_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071403_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_03)
tail(KO_MG_03)
str(KO_MG_03)
KO_MG_03 <- rename(KO_MG_03, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_03 <- subset(KO_MG_03, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_03)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_03$CDS <- str_sub(KO_MG_03$gene_id, 19, 25)
KO_MG_03$gene_id <- substr(KO_MG_03$gene_id, 0, 18)
head(KO_MG_03)
tail(KO_MG_03)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_O3_KO <- merge(GPM, KO_MG_03, by=c("gene_id", "CDS"), all=T)
head(MG_O3_KO, n=100)
tail(MG_O3_KO, n=100)

phylo <- read.table("phylodist/03_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_03 <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_03 <- subset(phylo_MG_03, select= -c(V2, V3))
phylo_MG_03$CDS <- str_sub(phylo_MG_03$gene_id, 19, 25)
phylo_MG_03$gene_id <- substr(phylo_MG_03$gene_id, 0, 18)

MG_03_GPM <- merge(phylo_MG_03, MG_O3_KO, by=c("gene_id", "CDS"), all=T)
head(MG_03_GPM)
tail(MG_03_GPM)

write.csv(MG_03_GPM, "metagenome_gpm/03_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/03_MG_GPM.csv")
ec <- read.table("ec/03_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/03_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_03 <- read.csv("metagenome_gpm/03_MG_GPM.csv")
GPM_03_sep <- separate(GPM_03, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                       sep=";", remove=TRUE)
head(GPM_03_sep)
Tthermo_03 <- subset(GPM_03_sep, species=="Thermus thermophilus")
GPM_03_sep <- subset(GPM_03_sep, species!="Thermus thermophilus" | is.na(species))

GPM_03_sep$taxonomy <- paste(GPM_03_sep$domain, GPM_03_sep$phylum, GPM_03_sep$class, GPM_03_sep$order, 
                             GPM_03_sep$family, GPM_03_sep$genus, GPM_03_sep$species, GPM_03_sep$species, sep=";")
GPM_03_sep <- subset(GPM_03_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_03, "fucking_thermo/03_thermo.csv", row.names = FALSE)

GPM_03_sep <- subset(GPM_03_sep, select=-c(GPM))
G3 <- ((sum(GPM_03_sep$rg)*sum(GPM_03_sep$rl))/sum(GPM_03_sep$flg))
GPM_03_sep <- transform(GPM_03_sep, GPM=(rg*rl*1e6)/(flg*G3))
write.csv(GPM_03_sep, "metagenome_gpm/03_MG_GPM_mod.csv", row.names=FALSE)

GPM_03 <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/03_MG_GPM_mod.csv")
GPM_03_ko <- subset(GPM_03, KO_num !="<NA>")
GPM_03_ko <- subset(GPM_03_ko, select=-c(GPM, EC, percent_id_ec))
G3 <- ((sum(GPM_03_ko$rg)*sum(GPM_03_ko$rl))/sum(GPM_03_ko$flg))
G3
GPM_03_ko <- transform(GPM_03_ko, GPM=(rg*rl*1e6)/(flg*G3))
GPM_03_ko_75 <- subset(GPM_03_ko, percent_id_tax >=75)
write.csv(GPM_03_ko_75, "metagenome_gpm/03_MG_GPM.csv", row.names=FALSE)
