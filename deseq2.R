require(DESeq2)

###load data####
ec_reads <- read.csv("ec_reads.csv", row.names=1)
head(ec_reads)

ec_meta <- read.csv("july_metadata.csv", row.names=1)
head(ec_meta)

phylo_reads <- read.csv("phylo_family_reads.csv", row.names=1)

ko_reads <- read.csv("MGMT_KOnum_reads.csv", row.names=1)
ko_reads_mt <- subset(ko_reads, select=c(MT_01, MT_02, MT_03, MT_05B, MT_05C, MT_06, MT_07, 
                                         MT_08, MT_09A, MT_10, MT_11, MT_12A, MT_12B, 
                                         MT_13A, MT_13C))


all(rownames(ec_meta)==colnames(ec_reads))
rownames(ec_meta)
colnames(ec_reads)
ec_reads <- subset(ec_reads, select=-c(name))

ec_meta$timepoint <- factor(ec_meta$timepoint)
####subset data so that we have metagenomes and metatranscriptomes####
####datasets for metatranscriptomes####
ec_reads_mt <- subset(ec_reads, select=c(MT_01, MT_02, MT_03, MT_05B, MT_05C, MT_06, MT_07, 
                                         MT_08, MT_09A, MT_09B, MT_10, MT_11, MT_12A, MT_12B, 
                                         MT_13A, MT_13C))

phylo_reads_mt <- subset(phylo_reads, select=c(MT_01, MT_02, MT_03, MT_05B, MT_05C, MT_06, MT_07, 
                                         MT_08, MT_09A, MT_09B, MT_10, MT_11, MT_12A, MT_12B, 
                                         MT_13A, MT_13C))
?rowsum

###remove single occurences which will throw off our analyses####
phylo_reads_mt$occur <- rowSums(phylo_reads_mt !=0)
ec_reads_mt$occur <- rowSums(ec_reads_mt !=0)
ko_reads_mt$occur <- rowSums(ko_reads_mt !=0)
head(ko_reads_mt)

ec_reads_mt <- subset(ec_reads_mt, occur > 1)
phylo_reads_mt <- subset(phylo_reads_mt, occur > 1)
ko_reads_mt <- subset(ko_reads_mt, occur > 1)

####get rid of occur column because metadata and samples need to be exactly the same####
ec_reads_mt <- subset(ec_reads_mt, select=-c(occur))
phylo_reads_mt <- subset(phylo_reads_mt, select=-c(occur))
ko_reads_mt <- subset(ko_reads_mt, select=-c(occur))


head(ec_reads_mt)
head(phylo_reads_mt)
meta_mt <- subset(ec_meta, type=="metatranscriptome")


###make datasets for just the metagenomes####
ec_reads_mg <- subset(ec_reads, select=c(MG_01, MG_09B, MG_13B, MG_02, MG_04, MG_05A, MG_03, MG_05C,
                                         MG_06, MG_07, MG_09A, MG_10, MG_11, MG_12A, MG_12B, MG_13A))
head(ec_reads_mg)


phylo_reads_mg <- subset(phylo_reads, select=c(MG_01, MG_09B, MG_13B, MG_02, MG_04, MG_05A, MG_03, MG_05C,
                                         MG_06, MG_07, MG_09A, MG_10, MG_11, MG_12A, MG_12B, MG_13A))
head(phylo_reads_mg)

####remove single occurences####
phylo_reads_mg$occur <- rowSums(phylo_reads_mg !=0)
ec_reads_mg$occur <- rowSums(ec_reads_mg !=0)

ec_reads_mg <- subset(ec_reads_mg, occur > 1)
phylo_reads_mg <- subset(phylo_reads_mg, occur > 1)

phylo_reads_mg <- subset(phylo_reads_mg, select=-c(occur))
ec_reads_mg <- subset(ec_reads_mg, select=-c(occur))

ec_meta_mg <- subset(ec_meta, type=="metagenome")
head(ec_meta_mg)

all(rownames(ec_meta_mg)==colnames(ec_reads_mg))
rownames(ec_meta_mg)
colnames(ec_reads_mg)
####deseq2 metagenome####
dds_mg <- DESeqDataSetFromMatrix(countData = ec_reads_mg,
                                 colData = ec_meta_mg,
                                 design = ~ tide)
dds_mg <- DESeq(dds_mg)
resultsNames(dds_mg)
res_mg <- results(dds_mg, name="tide_low_vs_high")
res_mg_LFC <- lfcShrink(dds_mg, coef="tide_low_vs_high")
res_mg
res_mg_ordered <- res_mg[order(res_mg$pvalue),]
summary(res_mg)
sum(res_mg$padj < 0.05, na.rm=TRUE)
res_mg_05 <- results(dds_mg, alpha=0.05)
summary(res_mg_05)
sum(res_mg_05$padj < 0.05, na.rm = TRUE)
plotMA(res_mg_LFC)
write.csv(as.data.frame(res_mg_ordered), file="tide_low_vs_high_ec_mg_qc.csv")


dds_mg_phylo <- DESeqDataSetFromMatrix(countData = phylo_reads_mg,
                                    colData = ec_meta_mg,
                                    design = ~ tide)
dds_mg_phylo <- DESeq(dds_mg_phylo)
resultsNames(dds_mg_phylo)
res_mg_phylo <- results(dds_mg_phylo, name="tide_low_vs_high")
res_mg_LFC_phylo <- lfcShrink(dds_mg_phylo, coef="tide_low_vs_high")
res_mg_phylo
res_mg_phylo_ordered <- res_mg_phylo[order(res_mg_phylo$log2FoldChange),]
summary(res_mg_phylo)
sum(res_mg_phylo$padj < 0.05, na.rm=TRUE)
res_mg_phylo_05 <- results(dds_mg_phylo, alpha=0.05)
summary(res_mg_phylo_05)
sum(res_mg_phylo_05$padj < 0.05, na.rm = TRUE)
plotMA(res_mg_LFC_phylo)
write.csv(as.data.frame(res_mg_phylo_ordered), file="tide_low_vs_high_phylo_mg_qc.csv")


###deseq2 metatranscriptome####
dds_mt <- DESeqDataSetFromMatrix(countData = ec_reads_mt,
                                 colData = ec_meta_mt,
                                 design = ~ tide)
dds_mt <- DESeq(dds_mt)
resultsNames(dds_mt)
res_mt <- results(dds_mt, name="tide_low_vs_high")
res_mt_LFC <- lfcShrink(dds_mt, coef="tide_low_vs_high")
res_mt
res_mt_ordered <- res_mt[order(res_mt$pvalue),]
summary(res_mt)
sum(res_mt$padj < 0.05, na.rm=TRUE)
res_mt_05 <- results(dds_mt, alpha=0.05)
summary(res_mt_05)
sum(res_mt_05$padj < 0.05, na.rm = TRUE)
plotMA(res_mt_LFC)
write.csv(as.data.frame(res_mt_ordered), file="tide_low_vs_high_ec_mt_qc.csv")


dds_mt_phylo <- DESeqDataSetFromMatrix(countData = phylo_reads_mt,
                                    colData = ec_meta_mt,
                                    design = ~ tide)
dds_mt_phylo <- DESeq(dds_mt_phylo)
resultsNames(dds_mt_phylo)
res_mt_phylo <- results(dds_mt_phylo, name="tide_low_vs_high")
res_mt_LFC_phylo <- lfcShrink(dds_mt_phylo, coef="tide_low_vs_high")
res_mt_phylo
res_mt_phylo_ordered <- res_mt_phylo[order(res_mt_phylo$log2FoldChange),]
summary(res_mt_phylo)
sum(res_mt_phylo$padj < 0.05, na.rm=TRUE)
res_mt_phylo_05 <- results(dds_mt_phylo, alpha=0.05)
summary(res_mt_phylo_05)
sum(res_mt_phylo_05$padj < 0.05, na.rm = TRUE)
plotMA(res_mt_LFC_phylo)
write.csv(as.data.frame(res_mt_phylo_ordered), file="tide_low_vs_high_phylo_mt_qc.csv")

dds_mt <- DESeqDataSetFromMatrix(countData = ko_reads_mt,
                                 colData = meta_mt,
                                 design = ~ tide)
dds_mt <- DESeq(dds_mt)
resultsNames(dds_mt)
res_mt <- results(dds_mt, name="tide_low_vs_high")
res_mt_LFC <- lfcShrink(dds_mt, coef="tide_low_vs_high")
res_mt
res_mt_ordered <- res_mt[order(res_mt$pvalue),]
summary(res_mt)
sum(res_mt$padj < 0.05, na.rm=TRUE)
res_mt_05 <- results(dds_mt, alpha=0.05)
summary(res_mt_05)
sum(res_mt_05$padj < 0.05, na.rm = TRUE)
plotMA(res_mt_LFC)



