require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_09A <- read.table("KO/071409A_MG_assembled_derep.ko.txt")
reads_MG_09A <- read.table("covstats/09A_mg_covstats.txt")
rl_MG_09A <- read.table("rl/09A_MG_rl_sort.txt")
head(KO_MG_09A)
head(reads_MG_09A)
head(rl_MG_09A)

####format rl data###
rl_MG_09A <- rename(rl_MG_09A, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_09A)
####split gene_id into gene_id and position####
rl_MG_09A <- cSplit(rl_MG_09A, "gene_id", ":")
rl_MG_09A = rename(rl_MG_09A, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_09A <- subset(rl_MG_09A, gene_id!="*")
head(rl_MG_09A)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_09A <-cSplit(reads_MG_09A, "V1", ":")
head(reads_MG_09A)
reads_MG_09A$rg <-(reads_MG_09A$V7+reads_MG_09A$V8)
reads_MG_09A <-subset(reads_MG_09A, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_09A = rename(reads_MG_09A, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_09A <- reads_MG_09A %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_09A, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_09A)
rl_MG_09A <- as.data.frame(rl_MG_09A)
reads_MG_09A <- merge(reads_MG_09A, rl_MG_09A, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_09A, n=15)
tail(reads_MG_09A, n=15)
reads_MG_09A <- subset(reads_MG_09A, select= -c(rg.x))
reads_MG_09A = rename(reads_MG_09A, rg=rg.y)
tail(reads_MG_09A, n=15)
reads_MG_09A[is.na(reads_MG_09A)] <- 0
###calculate T###
G <- ((sum(reads_MG_09A$rg)*sum(reads_MG_09A$rl))/sum(reads_MG_09A$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_09A, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071409A_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071409A_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_09A)
tail(KO_MG_09A)
str(KO_MG_09A)
KO_MG_09A <- rename(KO_MG_09A, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_09A <- subset(KO_MG_09A, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_09A)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_09A$CDS <- str_sub(KO_MG_09A$gene_id, 19, 25)
KO_MG_09A$gene_id <- substr(KO_MG_09A$gene_id, 0, 18)
head(KO_MG_09A)
tail(KO_MG_09A)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_O9A_KO <- merge(GPM, KO_MG_09A, by=c("gene_id", "CDS"), all=T)
head(MG_O9A_KO, n=100)
tail(MG_O9A_KO, n=100)

phylo <- read.table("phylodist/09A_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_09A <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_09A <- subset(phylo_MG_09A, select= -c(V2, V3))
phylo_MG_09A$CDS <- str_sub(phylo_MG_09A$gene_id, 19, 25)
phylo_MG_09A$gene_id <- substr(phylo_MG_09A$gene_id, 0, 18)

MG_09A_GPM <- merge(phylo_MG_09A, MG_O9A_KO, by=c("gene_id", "CDS"), all=T)
head(MG_09A_GPM)
tail(MG_09A_GPM)

write.csv(MG_09A_GPM, "metagenome_gpm/09A_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/09A_MG_GPM.csv")
ec <- read.table("ec/09A_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/09A_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_09A <- read.csv("metagenome_gpm/09A_MG_GPM.csv")
GPM_09A_sep <- separate(GPM_09A, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                       sep=";", remove=TRUE)
head(GPM_09A_sep)
Tthermo_09A <- subset(GPM_09A_sep, species=="Thermus thermophilus")
GPM_09A_sep <- subset(GPM_09A_sep, species!="Thermus thermophilus" | is.na(species))

GPM_09A_sep$taxonomy <- paste(GPM_09A_sep$domain, GPM_09A_sep$phylum, GPM_09A_sep$class, GPM_09A_sep$order, 
                             GPM_09A_sep$family, GPM_09A_sep$genus, GPM_09A_sep$species, GPM_09A_sep$species, sep=";")
GPM_09A_sep <- subset(GPM_09A_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_09A, "fucking_thermo/09A_thermo.csv", row.names = FALSE)

GPM_09A_sep <- subset(GPM_09A_sep, select=-c(GPM))
G9A <- ((sum(GPM_09A_sep$rg)*sum(GPM_09A_sep$rl))/sum(GPM_09A_sep$flg))
GPM_09A_sep <- transform(GPM_09A_sep, GPM=(rg*rl*1e6)/(flg*G9A))
write.csv(GPM_09A_sep, "metagenome_gpm/09A_MG_GPM_mod.csv", row.names=FALSE)

GPM_09A <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/09A_MG_GPM_mod.csv")
GPM_09A_ko <- subset(GPM_09A, KO_num !="<NA>")
GPM_09A_ko <- subset(GPM_09A_ko, select=-c(GPM, EC, percent_id_ec))
G9A <- ((sum(GPM_09A_ko$rg)*sum(GPM_09A_ko$rl))/sum(GPM_09A_ko$flg))
G9A
GPM_09A_ko <- transform(GPM_09A_ko, GPM=(rg*rl*1e6)/(flg*G9A))
GPM_09A_ko_75 <- subset(GPM_09A_ko, percent_id_tax >=75)
write.csv(GPM_09A_ko_75, "metagenome_gpm/09A_MG_GPM.csv", row.names=FALSE)

