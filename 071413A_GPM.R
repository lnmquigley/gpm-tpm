require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)

###import data###
KO_MG_13A <- read.table("KO/071413A_MG_assembled_derep.ko.txt")
reads_MG_13A <- read.table("covstats/13A_mg_covstats.txt")
rl_MG_13A <- read.table("rl/13A_MG_rl_sort.txt")
head(KO_MG_13A)
head(reads_MG_13A)
head(rl_MG_13A)

####format rl data###
rl_MG_13A <- rename(rl_MG_13A, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_13A)
####split gene_id into gene_id and position####
rl_MG_13A <- cSplit(rl_MG_13A, "gene_id", ":")
rl_MG_13A = rename(rl_MG_13A, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_13A <- subset(rl_MG_13A, gene_id!="*")
head(rl_MG_13A)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_13A <-cSplit(reads_MG_13A, "V1", ":")
head(reads_MG_13A)
reads_MG_13A$rg <-(reads_MG_13A$V7+reads_MG_13A$V8)
reads_MG_13A <-subset(reads_MG_13A, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_13A = rename(reads_MG_13A, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_13A <- reads_MG_13A %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_13A, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_13A)
rl_MG_13A <- as.data.frame(rl_MG_13A)
reads_MG_13A <- merge(reads_MG_13A, rl_MG_13A, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_13A, n=15)
tail(reads_MG_13A, n=15)
reads_MG_13A <- subset(reads_MG_13A, select= -c(rg.x))
reads_MG_13A = rename(reads_MG_13A, rg=rg.y)
tail(reads_MG_13A, n=15)
reads_MG_13A[is.na(reads_MG_13A)] <- 0
###calculate T###
G <- ((sum(reads_MG_13A$rg)*sum(reads_MG_13A$rl))/sum(reads_MG_13A$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_13A, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071413A_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071413A_MG_GPM.csv")
dim(GPM)
head(GPM)
tail(GPM)
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_13A)
tail(KO_MG_13A)
str(KO_MG_13A)
KO_MG_13A <- rename(KO_MG_13A, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_13A <- subset(KO_MG_13A, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_13A)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_13A$CDS <- str_sub(KO_MG_13A$gene_id, 19, 25)
KO_MG_13A$gene_id <- substr(KO_MG_13A$gene_id, 0, 18)
head(KO_MG_13A)
tail(KO_MG_13A)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_13A_KO <- merge(GPM, KO_MG_13A, by=c("gene_id", "CDS"), all=T)
head(MG_13A_KO, n=100)
tail(MG_13A_KO, n=100)

phylo <- read.table("phylodist/13A_MG_assembled.phylodist.txt", sep="\t", quote= "")
head(phylo)
phylo_MG_13A <- rename(phylo, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_13A <- subset(phylo_MG_13A, select= -c(V2, V3))
phylo_MG_13A$CDS <- str_sub(phylo_MG_13A$gene_id, 19, 25)
phylo_MG_13A$gene_id <- substr(phylo_MG_13A$gene_id, 0, 18)

MG_13A_GPM <- merge(phylo_MG_13A, MG_13A_KO, by=c("gene_id", "CDS"), all=T)
head(MG_13A_GPM)
tail(MG_13A_GPM)

write.csv(MG_13A_GPM, "metagenome_gpm/13A_MG_GPM.csv", row.names = FALSE)

GPM <- read.csv("metagenome_gpm/13A_MG_GPM.csv")
ec <- read.table("ec/13A_MG_assembled.ec.txt")
head(GPM)
head(ec)
ec <- rename(ec, gene_id=V1, EC=V3, percent_id_ec=V4)
ec$CDS <- str_sub(ec$gene_id, 19, 25)
ec$gene_id <- substr(ec$gene_id, 0, 18)
head(ec)
ec <- subset(ec, select=-c(V2, V5, V6, V7, V8, V9, V10, V11))
GPM <- merge(GPM, ec, by=c("gene_id", "CDS"), all=T)
head(GPM)
write.csv(GPM, "metagenome_gpm/13A_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_13A <- read.csv("metagenome_gpm/13A_MG_GPM.csv")
GPM_13A_sep <- separate(GPM_13A, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                        sep=";", remove=TRUE)
head(GPM_13A_sep)
Tthermo_13A <- subset(GPM_13A_sep, species=="Thermus thermophilus")
GPM_13A_sep <- subset(GPM_13A_sep, species!="Thermus thermophilus" | is.na(species))

GPM_13A_sep$taxonomy <- paste(GPM_13A_sep$domain, GPM_13A_sep$phylum, GPM_13A_sep$class, GPM_13A_sep$order, 
                              GPM_13A_sep$family, GPM_13A_sep$genus, GPM_13A_sep$species, GPM_13A_sep$species, sep=";")
GPM_13A_sep <- subset(GPM_13A_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_13A, "fucking_thermo/13A_thermo.csv", row.names = FALSE)

GPM_13A_sep <- subset(GPM_13A_sep, select=-c(GPM))
G13A <- ((sum(GPM_13A_sep$rg)*sum(GPM_13A_sep$rl))/sum(GPM_13A_sep$flg))
GPM_13A_sep <- transform(GPM_13A_sep, GPM=(rg*rl*1e6)/(flg*G13A))
write.csv(GPM_13A_sep, "metagenome_gpm/13A_MG_GPM_mod.csv", row.names=FALSE)


GPM_13A <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/13A_MG_GPM_mod.csv")
GPM_13A_ko <- subset(GPM_13A, KO_num !="<NA>")
GPM_13A_ko <- subset(GPM_13A_ko, select=-c(GPM, EC, percent_id_ec))
G13A <- ((sum(GPM_13A_ko$rg)*sum(GPM_13A_ko$rl))/sum(GPM_13A_ko$flg))
G13A
GPM_13A_ko <- transform(GPM_13A_ko, GPM=(rg*rl*1e6)/(flg*G13A))
GPM_13A_ko_75 <- subset(GPM_13A_ko, percent_id_tax >=75)
write.csv(GPM_13A_ko_75, "metagenome_gpm/13A_MG_GPM.csv", row.names=FALSE)

