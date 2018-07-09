require(tidyr)
phylo <- read.csv("phylo_75_reads.csv")
head(phylo)

phylo_sep <- separate(phylo, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                      sep=";", remove=TRUE)
head(phylo_sep)

phylo_sep <- subset(phylo_sep, genus !="Thermus")

phylo_sep <- subset(phylo_sep, select=-c(genus, species, strain))
head(phylo_sep)

phylo_sep$taxonomy <- paste(phylo_sep$domain, phylo_sep$phylum, phylo_sep$class, phylo_sep$order, 
                            phylo_sep$family, sep=";")
head(phylo_sep)
phylo_sep <- subset(phylo_sep, select=-c(domain, phylum, class, order, family))
phylo_by_family <- aggregate(. ~ taxonomy, phylo_sep, sum)
head(phylo_by_family)

write.csv(phylo_by_family, "phylo_family_reads.csv", row.names = FALSE)

