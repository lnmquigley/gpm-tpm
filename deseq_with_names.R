require(dplyr)

ec <- read.csv("ec.csv")
ec <- subset(ec, select=c(EC, name))
ko <- read.csv("ko_metadata.csv")
head(ko)

ko$KO_number <- sub("^", "KO:", ko$KO_number)
head(ko)

ko <- subset(ko, select=c(KO_number, KO_name))
ko <- ko[!duplicated(ko$KO_number),]

deseq_ec_mg <- read.csv("tide_low_vs_high_ec_mg.csv")
deseq_ec_mt <- read.csv("tide_low_vs_high_ec_mt.csv")
deseq_ko_mg <- read.csv("tide_low_vs_high_ko_mg.csv")
deseq_ko_mt <- read.csv("tide_low_vs_high_ko_mt.csv")

deseq_ec_mg <- rename(deseq_ec_mg, EC=X)
head(deseq_ec_mg)
deseq_ec_mt <- rename(deseq_ec_mt, EC=X)
head(deseq_ec_mt)
deseq_ko_mg <- rename(deseq_ko_mg, KO_number=X)
head(deseq_ko_mg)
deseq_ko_mt <- rename(deseq_ko_mt, KO_number=X)
head(deseq_ko_mt)

?merge
deseq_ec_mg <- merge(deseq_ec_mg, ec)
deseq_ec_mg <- deseq_ec_mg[order(deseq_ec_mg$padj),]
tail(deseq_ec_mg)
deseq_ec_mg <- subset(deseq_ec_mg, baseMean!=0)
deseq_ec_mg <- subset(deseq_ec_mg, padj!="NA")
write.csv(deseq_ec_mg, "tide_low_vs_high_ec_mg.csv", row.names = FALSE)

deseq_ec_mt <- merge(deseq_ec_mt, ec)
deseq_ec_mt <- deseq_ec_mt[order(deseq_ec_mt$padj),]
deseq_ec_mt <- subset(deseq_ec_mt, baseMean!=0)
deseq_ec_mt <- subset(deseq_ec_mt, padj!="NA")
write.csv(deseq_ec_mt, "tide_low_vs_high_ec_mt.csv", row.names=FALSE)

deseq_ko_mg <- merge(deseq_ko_mg, ko)
deseq_ko_mg <- deseq_ko_mg[order(deseq_ko_mg$padj),]
deseq_ko_mg <- subset(deseq_ko_mg, baseMean!=0)
deseq_ko_mg <- subset(deseq_ko_mg, padj!="NA")
tail(deseq_ko_mg)
write.csv(deseq_ko_mg, "tide_low_vs_high_ko_mg.csv", row.names = FALSE)

deseq_ko_mt <- merge(deseq_ko_mt, ko)
deseq_ko_mt <- deseq_ko_mt[order(deseq_ko_mt$padj),]
deseq_ko_mt <- subset(deseq_ko_mt, baseMean!=0)
deseq_ko_mt <- subset(deseq_ko_mt, padj!="NA")
tail(deseq_ko_mt)
write.csv(deseq_ko_mt, "tide_low_vs_high_ko_mt.csv", row.names = FALSE)
