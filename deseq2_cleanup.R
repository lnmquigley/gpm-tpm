require(dplyr)
ec_mg <- read.csv("DESeq2/qc'ed_nosingleoccur/tide_low_vs_high_ec_mg_qc.csv")
head(ec_mg)
ec_mt <- read.csv("DESeq2/qc'ed_nosingleoccur/tide_low_vs_high_ec_mt_qc.csv")
head(ec_mt)

ec_mt <- rename(ec_mt, EC=X)
ec_mg <- rename(ec_mg, EC=X)

ec <- read.csv("ec.csv")
head(ec)

ec <- subset(ec, select=c(EC, name))

ec_mt <- merge(ec_mt, ec)
ec_mg <- merge(ec_mg, ec)

write.csv(ec_mt, "DESeq2/qc'ed_nosingleoccur/tide_low_vs_high_ec_mt_qc.csv")
write.csv(ec_mg, "DESeq2/qc'ed_nosingleoccur/tide_low_vs_high_ec_mg_qc.csv")
