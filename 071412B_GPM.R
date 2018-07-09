require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_12B <- read.table("KO/071412B_MG_assembled_derep.ko.txt")
reads_MG_12B <- read.table("covstats/12B_mg_covstats.txt")
rl_MG_12B <- read.table("rl/12B_MG_rl_sort.txt")
head(KO_MG_12B)
head(reads_MG_12B)
head(rl_MG_12B)

####format rl data###
rl_MG_12B <- rename(rl_MG_12B, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_12B)
####split gene_id into gene_id and position####
rl_MG_12B <- cSplit(rl_MG_12B, "gene_id", ":")
rl_MG_12B = rename(rl_MG_12B, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_12B <- subset(rl_MG_12B, gene_id!="*")
head(rl_MG_12B)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_12B <-cSplit(reads_MG_12B, "V1", ":")
head(reads_MG_12B)
reads_MG_12B$rg <-(reads_MG_12B$V7+reads_MG_12B$V8)
reads_MG_12B <-subset(reads_MG_12B, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_12B = rename(reads_MG_12B, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_12B <- reads_MG_12B %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_12B, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_12B)
rl_MG_12B <- as.data.frame(rl_MG_12B)
reads_MG_12B <- merge(reads_MG_12B, rl_MG_12B, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_12B, n=15)
tail(reads_MG_12B, n=15)
reads_MG_12B <- subset(reads_MG_12B, select= -c(rg.x))
reads_MG_12B = rename(reads_MG_12B, rg=rg.y)
tail(reads_MG_12B, n=15)
reads_MG_12B[is.na(reads_MG_12B)] <- 0
###calculate T###
G <- ((sum(reads_MG_12B$rg)*sum(reads_MG_12B$rl))/sum(reads_MG_12B$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_12B, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071412B_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071412B_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_12B)
tail(KO_MG_12B)
str(KO_MG_12B)
KO_MG_12B <- rename(KO_MG_12B, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_12B <- subset(KO_MG_12B, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_12B)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_12B$CDS <- str_sub(KO_MG_12B$gene_id, 19, 25)
KO_MG_12B$gene_id <- substr(KO_MG_12B$gene_id, 0, 18)
head(KO_MG_12B)
tail(KO_MG_12B)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_12B_KO <- merge(GPM, KO_MG_12B, by=c("gene_id", "CDS"), all=T)
head(MG_12B_KO, n=100)
tail(MG_12B_KO, n=100)

phylo <- read.table("phylodist/12B_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_12B <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_12B <- subset(phylo_MG_12B, select= -c(V2, V3))
phylo_MG_12B$CDS <- str_sub(phylo_MG_12B$gene_id, 19, 25)
phylo_MG_12B$gene_id <- substr(phylo_MG_12B$gene_id, 0, 18)

MG_12B_GPM <- merge(phylo_MG_12B, MG_12B_KO, by=c("gene_id", "CDS"), all=T)
head(MG_12B_GPM)
tail(MG_12B_GPM)

write.csv(MG_12B_GPM, "metagenome_gpm/12B_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/12B_MG_GPM.csv")
ec <- read.table("ec/12B_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/12B_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_12B <- read.csv("metagenome_gpm/12B_MG_GPM.csv")
GPM_12B_sep <- separate(GPM_12B, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                       sep=";", remove=TRUE)
head(GPM_12B_sep)
Tthermo_12B <- subset(GPM_12B_sep, species=="Thermus thermophilus")
GPM_12B_sep <- subset(GPM_12B_sep, species!="Thermus thermophilus" | is.na(species))

GPM_12B_sep$taxonomy <- paste(GPM_12B_sep$domain, GPM_12B_sep$phylum, GPM_12B_sep$class, GPM_12B_sep$order, 
                             GPM_12B_sep$family, GPM_12B_sep$genus, GPM_12B_sep$species, GPM_12B_sep$species, sep=";")
GPM_12B_sep <- subset(GPM_12B_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_12B, "fucking_thermo/12B_thermo.csv", row.names = FALSE)

GPM_12B_sep <- subset(GPM_12B_sep, select=-c(GPM))
G12B <- ((sum(GPM_12B_sep$rg)*sum(GPM_12B_sep$rl))/sum(GPM_12B_sep$flg))
GPM_12B_sep <- transform(GPM_12B_sep, GPM=(rg*rl*1e6)/(flg*G12B))
write.csv(GPM_12B_sep, "metagenome_gpm/12B_MG_GPM_mod.csv", row.names=FALSE)

GPM_12B <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/12B_MG_GPM_mod.csv")
GPM_12B_ko <- subset(GPM_12B, KO_num !="<NA>")
GPM_12B_ko <- subset(GPM_12B_ko, select=-c(GPM, EC, percent_id_ec))
G12B <- ((sum(GPM_12B_ko$rg)*sum(GPM_12B_ko$rl))/sum(GPM_12B_ko$flg))
G12B
GPM_12B_ko <- transform(GPM_12B_ko, GPM=(rg*rl*1e6)/(flg*G12B))
GPM_12B_ko_75 <- subset(GPM_12B_ko, percent_id_tax >=75)
write.csv(GPM_12B_ko_75, "metagenome_gpm/12B_MG_GPM.csv", row.names=FALSE)

GPM_12A <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/12A_MG_GPM_mod.csv")
GPM_12A_ko <- subset(GPM_12A, KO_num !="<NA>")
GPM_12A_ko <- subset(GPM_12A_ko, select=-c(GPM, EC, percent_id_ec))
G12A <- ((sum(GPM_12A_ko$rg)*sum(GPM_12A_ko$rl))/sum(GPM_12A_ko$flg))
G12A
GPM_12A_ko <- transform(GPM_12A_ko, GPM=(rg*rl*1e6)/(flg*G12A))
GPM_12A_ko_75 <- subset(GPM_12A_ko, percent_id_tax >=75)
write.csv(GPM_12A_ko_75, "metagenome_gpm/12A_MG_GPM.csv", row.names=FALSE)
