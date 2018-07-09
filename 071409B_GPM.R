
require(stringr)
require(splitstackshape)
require(dplyr)
require(data.table)


###import data###
KO_MG_09B <- read.table("KO/071409B_MG_assembled_derep.ko.txt")
reads_MG_09B <- read.table("covstats/09B_mg_covstats.txt")
rl_MG_09B <- read.table("rl/09B_MG_rl_sort.txt")
head(KO_MG_09B)
head(reads_MG_09B)
head(rl_MG_09B)

####format rl data###
rl_MG_09B <- rename(rl_MG_09B, gene_id=V1, rg=V2, total_length=V3, rl=V4)
head(rl_MG_09B)
####split gene_id into gene_id and position####
rl_MG_09B <- cSplit(rl_MG_09B, "gene_id", ":")
rl_MG_09B = rename(rl_MG_09B, gene_id = gene_id_1, position=gene_id_2)
####remove * which represents all unmapped reads####
rl_MG_09B <- subset(rl_MG_09B, gene_id!="*")
head(rl_MG_09B)

####generate column that is the occurence of gene_id which corresponds to###
####relative position of the CDS###
reads_MG_09B<-cSplit(reads_MG_09B, "V1", ":")
head(reads_MG_09B)
reads_MG_09B$rg <-(reads_MG_09B$V7+reads_MG_09B$V8)
reads_MG_09B <-subset(reads_MG_09B, select= -c(V2, V4, V5, V6, V7, V8, V9, V10, V11))
reads_MG_09B = rename(reads_MG_09B, flg=V3, gene_id=V1_1, position=V1_2)
reads_MG_09B <- reads_MG_09B %>% group_by(gene_id) %>% mutate (CDS=1:n())
tail(reads_MG_09B, n=12)

###need the covstats.txt and and rl.txt to be the same length###
###remove all scaffolds for which there was no coverage###

###combine covstats and rl to get a dataframe with all necessary components to calculate TPM###
head(rl_MG_09B)
rl_MG_09B <- as.data.frame(rl_MG_09B)
reads_MG_09B <- merge(reads_MG_09B, rl_MG_09B, by=c("gene_id", "position"), all=TRUE)
head(reads_MG_09B, n=15)
tail(reads_MG_09B, n=15)
reads_MG_09B <- subset(reads_MG_09B, select= -c(rg.x))
reads_MG_09B = rename(reads_MG_09B, rg=rg.y)
tail(reads_MG_09B, n=15)
reads_MG_09B[is.na(reads_MG_09B)] <- 0
###calculate T###
G <- ((sum(reads_MG_09B$rg)*sum(reads_MG_09B$rl))/sum(reads_MG_09B$flg))
###create a new variable with TPM for each row (scaffold/gene_id)####
GPM <- transform(reads_MG_09B, GPM=(rg*rl*1e6)/(flg*G))
head(GPM)
tail(GPM)



write.csv(GPM, file = "071409B_MG_GPM.csv",row.names=FALSE)
GPM <- read.csv("metagenome_gpm/071409B_MG_GPM.csv")
####there's a lot of extra information in the KO file that I don't necessary care about####
####like coordinates and GC content so lets get ride of that####
head(KO_MG_09B)
tail(KO_MG_09B)
str(KO_MG_09B)
KO_MG_09B <- rename(KO_MG_09B, gene_id=V1, KO_num=V3, e_val=V9)
KO_MG_09B <- subset(KO_MG_09B, select= -c(V2, V4, V5, V6, V7, V8, V10, V11))
head(KO_MG_09B)

###separate out CDS number from gene_ID and delete CDS number from gene_ID after)
KO_MG_09B$CDS <- str_sub(KO_MG_09B$gene_id, 19, 25)
KO_MG_09B$gene_id <- substr(KO_MG_09B$gene_id, 0, 18)
head(KO_MG_09B)
tail(KO_MG_09B)
####merge TPM df with KO df - this will be tricky since KO df is smaller than TPM df####

MG_09B_KO <- merge(KO_MG_09B, GPM, by=c("gene_id", "CDS"), all.x=T)
head(MG_09B_KO, n=100)
tail(MG_09B_KO, n=100)

####want to combine phylodist with KO and gene_id to create master spreadsheet for sample##
phylo_MG_09B <- read.table("phylodist/09B_MG_assembled.phylodist.txt", sep="\t")
head(phylo_MG_09B)
phylo_MG_09B <- rename(phylo_MG_09B, gene_id=V1, percent_id_tax=V4, taxonomy=V5)
phylo_MG_09B <- subset(phylo_MG_09B, select = -c(V2, V3))
phylo_MG_09B$CDS <- str_sub(phylo_MG_09B$gene_id, 19, 25)
phylo_MG_09B$gene_id <- substr(phylo_MG_09B$gene_id, 0, 18)
MG_09B <- merge(MG_09B_KO, phylo_MG_09B, by=c("gene_id", "CDS"), all.x=T)
head(MG_09B)

write.csv(MG_09B, file="metagenome_gpm/09B_MG_GPM.csv", row.names=FALSE)

###remove tthermo###
require(tidyr)
GPM_09B <- read.csv("metagenome_gpm/09B_MG_GPM.csv")
GPM_09B_sep <- separate(GPM_09B, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                        sep=";", remove=TRUE)
head(GPM_09B_sep)
Tthermo_09B <- subset(GPM_09B_sep, species=="Thermus thermophilus")
GPM_09B_sep <- subset(GPM_09B_sep, species!="Thermus thermophilus" | is.na(species))

GPM_09B_sep$taxonomy <- paste(GPM_09B_sep$domain, GPM_09B_sep$phylum, GPM_09B_sep$class, GPM_09B_sep$order, 
                              GPM_09B_sep$family, GPM_09B_sep$genus, GPM_09B_sep$species, GPM_09B_sep$species, sep=";")
GPM_09B_sep <- subset(GPM_09B_sep, select=-c(domain, phylum, class, order, family, genus, species, strain))
write.csv(Tthermo_09B, "fucking_thermo/09B_thermo.csv", row.names = FALSE)

GPM_09B_sep <- subset(GPM_09B_sep, select=-c(GPM))
G9B <- ((sum(GPM_09B_sep$rg)*sum(GPM_09B_sep$rl))/sum(GPM_09B_sep$flg))
GPM_09B_sep <- transform(GPM_09B_sep, GPM=(rg*rl*1e6)/(flg*G9A))
write.csv(GPM_09B_sep, "metagenome_gpm/09B_MG_GPM_mod.csv", row.names=FALSE)


GPM_09B <- read.csv("metagenome_gpm/V2 - GPM calculated based on all contigs (should only be based on annotated contigs/09B_MG_GPM_mod.csv")
GPM_09B_ko <- subset(GPM_09B, KO_num !="<NA>")
GPM_09B_ko <- subset(GPM_09B_ko, select=-c(GPM, EC, percent_id_ec))
G9B <- ((sum(GPM_09B_ko$rg)*sum(GPM_09B_ko$rl))/sum(GPM_09B_ko$flg))
G9B
GPM_09B_ko <- transform(GPM_09B_ko, GPM=(rg*rl*1e6)/(flg*G9B))
GPM_09B_ko_75 <- subset(GPM_09B_ko, percent_id_tax >=75)
write.csv(GPM_09B_ko_75, "metagenome_gpm/09B_MG_GPM.csv", row.names=FALSE)

