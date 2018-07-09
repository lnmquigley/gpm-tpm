require(dplyr)
ko_meta <- read.csv("ko_metadata.csv")
ko_meta$KO_number <- sub("^", "KO:", ko_meta$KO_number)
ko_TPM <- read.csv("KOnum_TPM.csv")

head(ko_meta)
head(ko_TPM)

kobygroup2 <- count(ko_meta, group_2)

ko_group <-subset(ko_meta, select=c(group_1, group_2))

ko_path <- merge(ko_meta, ko_TPM)

ko_path <- subset(ko_path, select=-c(pathway, pathway_number, group_1_number, group_2_number, KO_name, KO_number,group_1))
head(ko_path)


ko_path <- aggregate(. ~ group_2, ko_path, sum)
ko_path <- merge(ko_path, ko_group)
ko_path <- merge(ko_path, kobygroup2)
head(ko_path)
?aggregate
?unique

ko_path <- unique(ko_path)
ko_path_bac <- subset(ko_path, group_1!="Organismal Systems" & group_1!="Human Diseases")
head(ko_path_bac)

write.csv(ko_path_bac, "log_fold_change_long_early.csv")



ko_path_bac$log6 <- log10(ko_path_bac$MT_06SA/ko_path_bac$MT_06TA)
ko_path_bac$log12 <- log10(ko_path_bac$MT_12SA/ko_path_bac$MT_12TA)
ko_path_bac$log13 <- log10(ko_path_bac$MT_13SA/ko_path_bac$MT_13TA)
ko_path_bac <- subset(ko_path_bac, group_2!="RNA family")

ko_path_bac$replog5 <- log10(ko_path_bac$MT_05CJ/ko_path_bac$MT_05BJ)
ko_path_bac$replog12 <- log10(ko_path_bac$MT_12BJ/ko_path_bac$MT_12AJ)

write.csv(ko_path_bac, "log_fold_change_ko.csv", row.names=FALSE)
