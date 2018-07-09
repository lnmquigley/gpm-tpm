require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_05C <- read.table("KO/071405C_MG_assembled_derep.ko.txt")
reads_MG_05C <- read.table("covstats/05C_mg_covstats.txt")
rl_MG_05C <- read.table("rl/05C_MG_rl_sort.txt")
head(KO_MG_05C)
head(reads_MG_05C)
head(rl_MG_05C)

####format rl data###
rl_MG_05C <- rename(rl_MG_05C, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_05C)
####split gene_id into gene_id and position####
rl_MG_05C <- cSplit(rl_MG_05C, "gene_id", ":")
rl_MG_05C = rename(rl_MG_05C, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_05C <- subset(rl_MG_05C, gene_id!="*")
head(rl_MG_05C)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_05C<-cSplit(reads_MG_05C, "V1", ":")
head(reads_MG_05C)
reads_MG_05C$rg <-(reads_MG_05C$V7+reads_MG_05C$V8)
reads_MG_05C <-subset(reads_MG_05C, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_05C = rename(reads_MG_05C, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_05C <- reads_MG_05C %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_05C, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_05C)
rl_MG_05C <- as.data.frame(rl_MG_05C)
reads_MG_05C <- merge(reads_MG_05C, rl_MG_05C, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_05C, n=15)
tail(reads_MG_05C, n=15)
reads_MG_05C <- subset(reads_MG_05C, select= -c(rg.x))
reads_MG_05C = rename(reads_MG_05C, rg=rg.y)
tail(reads_MG_05C, n=15)
reads_MG_05C[is.na(reads_MG_05C)] <- 0
###calculate T###
G <- ((sum(reads_MG_05C$rg)*sum(reads_MG_05C$rl))/sum(reads_MG_05C$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_05C, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071405C_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071405C_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_05C)
tail(KO_MG_05C)
str(KO_MG_05C)
KO_MG_05C <- rename(KO_MG_05C, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_05C <- subset(KO_MG_05C, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_05C)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_05C$CDS <- str_sub(KO_MG_05C$gene_id, 19, 25)
KO_MG_05C$gene_id <- substr(KO_MG_05C$gene_id, 0, 18)
head(KO_MG_05C)
tail(KO_MG_05C)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_O5C_KO <- merge(GPM, KO_MG_05C, by=c("gene_id", "CDS"), all=T)
head(MG_O5C_KO, n=100)
tail(MG_O5C_KO, n=100)

phylo <- read.table("phylodist/05C_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_05C <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_05C <- subset(phylo_MG_05C, select= -c(V2, V3))
phylo_MG_05C$CDS <- str_sub(phylo_MG_05C$gene_id, 19, 25)
phylo_MG_05C$gene_id <- substr(phylo_MG_05C$gene_id, 0, 18)

MG_05C_GPM <- merge(phylo_MG_05C, MG_O5C_KO, by=c("gene_id", "CDS"), all=T)
head(MG_05C_GPM)
tail(MG_05C_GPM)

write.csv(MG_05C_GPM, "metagenome_gpm/05C_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/05C_MG_GPM.csv")
ec <- read.table("ec/05C_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/05C_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_05C <- read.csv("metagenome_gpm/05C_MG_GPM.csv")
GPM_05C_sep <- separate(GPM_05C, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                        sep=";", remove=TRUE)
head(GPM_05C_sep)
Tthermo_05C <- subset(GPM_05C_sep, species=="Thermus thermophilus")
GPM_05C_sep <- subset(GPM_05C_sep, species!="Thermus thermophilus" | is.na(species))

GPM_05C_sep$taxonomy <- paste(GPM_05C_sep$domain, GPM_05C_sep$phylum, GPM_05C_sep$class, GPM_05C_sep$order, 
                              GPM_05C_sep$family, GPM_05C_sep$genus, GPM_05C_sep$species, GPM_05C_sep$species, sep=";")
GPM_05C_sep <- subset(GPM_05C_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_05C, "fucking_thermo/05C_thermo.csv", row.names = FALSE)

GPM_05C_sep <- subset(GPM_05C_sep, select=-c(GPM))
G5C <- ((sum(GPM_05C_sep$rg)*sum(GPM_05C_sep$rl))/sum(GPM_05C_sep$flg))
GPM_05C_sep <- transform(GPM_05C_sep, GPM=(rg*rl*1e6)/(flg*G5C))
write.csv(GPM_05C_sep, "metagenome_gpm/05C_MG_GPM_mod.csv", row.names=FALSE)

GPM_05C <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/05C_MG_GPM_mod.csv")
GPM_05C_ko <- subset(GPM_05C, KO_num !="<NA>")
GPM_05C_ko <- subset(GPM_05C_ko, select=-c(GPM, EC, percent_id_ec))
G5C <- ((sum(GPM_05C_ko$rg)*sum(GPM_05C_ko$rl))/sum(GPM_05C_ko$flg))
G5C
GPM_05C_ko <- transform(GPM_05C_ko, GPM=(rg*rl*1e6)/(flg*G5C))
GPM_05C_ko_75 <- subset(GPM_05C_ko, percent_id_tax >=75)
write.csv(GPM_05C_ko_75, "metagenome_gpm/05C_MG_GPM.csv", row.names=FALSE)
