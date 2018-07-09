require(dplyr)
require(stringr)
require(splitstackshape)
####combine functional and taxonomic annotations for each scaffold###
ko <- read.table("raw data - acf:jgi files/041412BS_assembled_derep.ko", sep="\t")
phylo <- read.table("raw data - acf:jgi files/041412BS_assembled.phylodist", sep="\t", quote = "")
head(ko)
head(phylo)

phylo <- rename(phylo, gene=V1, percent_id=V4, taxonomy=V5)
phylo <- subset(phylo, select=-c(V2, V3))

ko <- rename(ko, gene=V1, KO_num=V3, evalue=V9)
ko <- subset(ko, select=-c(V2, V4, V5, V6, V7, V8, V10, V11))

genes <- merge(ko, phylo, all=TRUE)

####split out CDS from scaffolds###
genes$CDS <- str_sub(genes$gene, 18, 25)
genes$gene_id <- substr(genes$gene, 0, 17)
genes <- subset(genes, select=-c(gene))
head(genes)
tail(genes)

###combine covstats with average length###
cov <- read.table("raw data - acf:jgi files/12BS_covstats.txt")
cov$rg <- cov$V7+cov$V8
cov <- rename(cov, gene=V1, flg=V3)
cov <- subset(cov, select=c(gene, flg, rg))
head(cov)
rl <- read.table("raw data - acf:jgi files/12BS_rl_sort.txt")
rl <- subset(rl, select=-c(V2, V3))
rl <- rename(rl, gene=V1, rl=V4)
head(rl)
reads <- merge(cov, rl, all=TRUE)

###split out CDS from scaffolds in read information####
reads <- cSplit(reads, "gene", ":")
head(reads)
reads = rename(reads, gene_id = gene_1, position=gene_2)
reads <- reads %>% group_by(gene_id) %>% mutate (CDS=1:n())

###combine reads with annotations###
MT_12BS <- merge(genes, reads, all = TRUE)
MT_12BS[is.na(MT_12BS)] <- 0
###calculate T###
G <- ((sum(MT_12BS$rg)*sum(MT_12BS$rl))/sum(MT_12BS$flg))
G
###create a new variable with TPM for each row (scaffold/gene_id)####
TPM <- transform(MT_12BS, TPM=(rg*rl*1e6)/(flg*G))
TPM$KO_num <- as.character(TPM$KO_num)
TPM$KO_num[is.na(TPM$KO_num)] <- "CDS without KO annotation"

write.csv(TPM, "12BS_TPM_allCDS.csv", row.names = FALSE)
