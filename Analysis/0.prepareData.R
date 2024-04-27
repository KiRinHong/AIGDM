rm(list = ls())
setwd("~/Documents/Project/LONGITUDINAL/code/Analysis/")
library(tidyverse)
library(vegan)
library(phyloseq)
source("0a.utility.R")

# ECAM data ===========
#### OTU table ####
otu <- read.table(file = "../Data/Rawdata/ECAM/feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")
colnames(otu) <- unlist(strsplit(colnames(otu), split = "X"))[seq(2, 2*ncol(otu), 2)]

subjectid <- sapply(colnames(otu), function(x) unlist(strsplit(x, "\\."))[2])
# for children duplicates, keep the one with largest sequencing depth
id.duplicate <- colnames(otu)[which(endsWith(colnames(otu), ".r"))]
seq.depth.dup <- colSums(otu[,id.duplicate])
id.original <- sub(".r", "", id.duplicate)
seq.depth.ori <- colSums(otu[,id.original])
id.rm <- c(id.duplicate[seq.depth.dup < seq.depth.ori],
           id.original[seq.depth.ori < seq.depth.dup],
           colnames(otu)[which(startsWith(subjectid, "M"))]) # remove Mother subject
otu.keep <- otu[, !(colnames(otu) %in% id.rm)]
colnames(otu.keep) <- sub(".r", "", colnames(otu.keep))

#### meta data ####
metadata <- read.table(file = "../Data/Rawdata/ECAM/metadata.tsv", sep = "\t", header = T, row.names = 1)
metadata <- metadata[!(rownames(metadata) %in% id.rm),]
rownames(metadata) <- sub(".r", "", rownames(metadata))
metadata <- metadata[,-which(sapply(metadata, function(x) length(unique(x)))==1)]

length(table(metadata$studyid)) # children id 
length(table(metadata$host_subject_id)) # mother id

metadata$rowname <- rownames(metadata)

metadata.rm <- metadata %>%
  mutate(diet = ifelse(diet_3 %in% c("eb", "bd"), 1, 0)) %>% # 1:breast-feed dominant 0:formula-feed dominant
  mutate(baby_sex = ifelse(baby_sex == "female", 1, 0)) %>% # 1:female 0:male
  mutate(delivery = ifelse(delivery == "Vaginal", 1, 0)) %>% # 1:Vaginal 0:Cesarean
  mutate(life_month = sapply(diet_2_month, function(x) unlist(strsplit(x, "_"))[2])) %>%
  # period 1: 0-1 mon; 2: 1-3 mon; 3: 3-6 mon; 4: 6-12 mon; 5: 12-24 mon
  # 1 subject with 1 period (missing period 1-4), 9 subjects with 4 periods (missing period 5), 33 subjects with 5 periods
  mutate(period = ifelse(life_month %in% c(0,1), 1,
                         ifelse(life_month %in% c(2,3), 2,
                                ifelse(life_month %in% 4:6, 3,
                                       ifelse(life_month %in% 7:12, 4,
                                              ifelse(life_month %in% 12:24, 5, NA)))))) %>%
  select(rowname, studyid, antiexposedall, baby_sex, delivery, diet, day_of_life, description, life_month, period)

tmp <- metadata.rm %>%
  group_by(studyid) %>%
  summarise_at(vars(baby_sex, delivery, diet, antiexposedall, period), n_distinct) %>%
  ungroup()
lapply(tmp[,-1], sum)

tmp <- metadata.rm %>%
  group_by(studyid) %>%
  summarise(n_reps = n(), sex = first(baby_sex), delivery = first(delivery), diet = first(diet), 
            antibiotics = ifelse(n_distinct(antiexposedall)==2, 1, 0)) %>% # 1:yes 2:no
  ungroup()
# the result is consistent with table in the paper
table(tmp$delivery)
table(tmp$sex)
table(tmp$diet, tmp$delivery, tmp$antibiotics)


#### Taxonomy table ####
taxonomy <- read.table(file = "../Data/Rawdata/ECAM/taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)
tax.clean[tax.clean==""] <- NA
# first iterate through each taxonomic rank starting from the 'genus' and moving upward. 
# For each rank, it will group by the names in that rank and then check if there's more than one unique name in the rank above. 
# If this is the case, it will output the rows where these discrepancies occur.
tmp <- tax.clean[!duplicated(tax.clean),]
# Use the function on your data
tmp_incons <- na.omit(check_inconsistencies(tmp))

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""
# https://forum.qiime2.org/t/the-same-three-genus-belongs-to-different-family/13755/2
tax.clean$Genus[tax.clean$Family == "Lachnospiraceae" & tax.clean$Genus == "Clostridium"] <- "Lach_Clostridium"
tax.clean$Genus[tax.clean$Family == "Ruminococcaceae" & tax.clean$Genus == "Clostridium"] <- "Rumi_Clostridium"
tax.clean$Genus[tax.clean$Family == "Clostridiaceae" & tax.clean$Genus == "Clostridium"] <- "Clos_Clostridium"
tax.clean$Genus[tax.clean$Family == "Peptostreptococcaceae" & tax.clean$Genus == "Clostridium"] <- "Pept_Clostridium"
tax.clean$Genus[tax.clean$Family == "Erysipelotrichaceae" & tax.clean$Genus == "Clostridium"] <- "Erys_Clostridium"
tax.clean$Genus[tax.clean$Family == "Lachnospiraceae" & tax.clean$Genus == "Ruminococcus"] <- "Lach_Ruminococcus"
tax.clean$Genus[tax.clean$Family == "Ruminococcaceae" & tax.clean$Genus == "Ruminococcus"] <- "Rumi_Ruminococcus"

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

#tmp <- tax.clean[!duplicated(tax.clean),]
TREE <- read_tree("../Data/Rawdata/ECAM/tree.nwk")
TAX <- tax_table(as.matrix(tax.clean))
OTU <- otu_table(as.matrix(otu.keep), taxa_are_rows = TRUE)
META <- sample_data(metadata.rm)
ps <- phyloseq(OTU, TAX, META, TREE)
ps # 3835 taxa and 897 samples

ps_agg <- tax_glom(ps, taxrank = "Genus")
ps_agg <- subset_samples(ps_agg, sample_sums(ps_agg) > 1000)

ps_agg_subset <- subset_samples(ps_agg, period == 4)
ps_agg_subset <- subset_taxa(ps_agg_subset, taxa_sums(ps_agg_subset) != 0)
tmp.full <- otu_table(ps_agg_subset)@.Data
tmp.full <- t(tmp.full)
tmp.full <- tmp.full/rowSums(tmp.full)
sum(colMeans(tmp.full>0) > 0.1) # 0.1: 57
tmp.id <- names(which(colMeans(tmp.full>0) > 0.1))
sum(colMeans(tmp.full)[tmp.id]) # 0.1: 99.98%

ps.nonrare.10 <- get_nonrare_taxa(ps_agg_subset, threshold = 0.1, discard_other = T)
ps.nonrare.10 # 262 samples x 57 taxa
meta <- sample_data(ps.nonrare.10)
M <- otu_table(ps.nonrare.10)
ecam.p4 <- list(phy = ps.nonrare.10,
                meta = meta, otu.tab = M@.Data, 
                tax.tab = tax_table(ps.nonrare.10))

ps_agg_subset <- subset_samples(ps_agg, period %in% c(4,5))
ps_agg_subset <- subset_taxa(ps_agg_subset, taxa_sums(ps_agg_subset) != 0)
tmp.full <- otu_table(ps_agg_subset)@.Data
tmp.full <- t(tmp.full)
tmp.full <- tmp.full/rowSums(tmp.full)
sum(colMeans(tmp.full>0) > 0.1) # 0.1: 66
tmp.id <- names(which(colMeans(tmp.full>0) > 0.1))
sum(colMeans(tmp.full)[tmp.id]) # 0.1: 99.98%

ps.nonrare.10 <- get_nonrare_taxa(ps_agg_subset, threshold = 0.1, discard_other = T)
ps.nonrare.10 # 490 samples x 66 taxa
meta <- sample_data(ps.nonrare.10)
M <- otu_table(ps.nonrare.10)
ecam.p45 <- list(phy = ps.nonrare.10,
                 meta = meta, otu.tab = M@.Data, 
                 tax.tab = tax_table(ps.nonrare.10))

ps_agg_subset <- subset_samples(ps_agg, period %in% 1:5)
ps_agg_subset <- subset_taxa(ps_agg_subset, taxa_sums(ps_agg_subset) != 0)
tmp.full <- otu_table(ps_agg_subset)@.Data
tmp.full <- t(tmp.full)
tmp.full <- tmp.full/rowSums(tmp.full)
sum(colMeans(tmp.full>0) > 0.1) # 0.1: 61
tmp.id <- names(which(colMeans(tmp.full>0) > 0.1))
sum(colMeans(tmp.full)[tmp.id]) # 0.1: 99.98%

ps.nonrare.10 <- get_nonrare_taxa(ps_agg_subset, threshold = 0.1, discard_other = T)
ps.nonrare.10 # 840 samples x 61 taxa
meta <- sample_data(ps.nonrare.10)
M <- otu_table(ps.nonrare.10)
ecam.pall <- list(phy = ps.nonrare.10,
                  meta = meta, otu.tab = M@.Data, 
                  tax.tab = tax_table(ps.nonrare.10))

ecam <- list(p4 = ecam.p4, p45 = ecam.p45, pall = ecam.pall)
save(ecam, file = "../Data/Deriveddata/ECAM.Genus.filter10.rda")
