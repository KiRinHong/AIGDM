setwd("~/Documents/Project/LONGITUDINAL/code/Analysis/")
library(phyloseq)
library(data.tree)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggtree)
library(dplyr)
library(magick)
library(patchwork)
library(MASS)
library(picante)
library(matrixStats)

rm(list = ls())

####### ECAM (delivery mode) Fig 1########
load("../Data/Deriveddata/ECAM.Genus.filter10.rda")
ecam <- ecam$pall
dim(ecam$meta) # 840 x 10
dim(ecam$otu.tab) # 840 x 61
apply(ecam$tax.tab, 2, function(x) length(unique(x)))
# Kingdom  Phylum   Class   Order  Family   Genus Species 
# 1       6      12      18      31      61       1 
length(unique(ecam$meta$studyid)) # 43 (0:19, 1:24)
colSums(table(ecam$meta$studyid, ecam$meta$delivery) > 0)

ps <- ecam$phy
otu <- as.data.frame(t(ps@otu_table))
meta <- as.data.frame(ps@sam_data)

df.pd <- pd(otu, phy_tree(ps)) 
meta$PD <- df.pd$PD
meta$delivery[meta$delivery == 0] <- "C-section born"
meta$delivery[meta$delivery == 1] <- "Vaginally delivered"
meta$baby_sex[meta$baby_sex == 0] <- "Male"
meta$baby_sex[meta$baby_sex == 1] <- "Female"
meta$life_month[meta$life_month == "24"] <- "23"
meta$life_month <- factor(meta$life_month, labels = as.character(1:24), levels = 0:23)

meta.rmdup <- data.frame(as.matrix(meta[,c("studyid", "life_month", "delivery")])) # remove reps for generating plot, same studyID and same lifeMonth 
meta.rmdup <- meta.rmdup[!duplicated(meta.rmdup),] # 653
meta.timepoint <- expand.grid(StudyID = unique(meta.rmdup$studyid), LifeMonth = 0:23, Value = 0, Delivery = 0)
meta.timepoint$LifeMonth <- factor(meta.timepoint$LifeMonth, labels = 1:24, levels = 0:23)
meta.timepoint$Value[paste(meta.timepoint$StudyID, meta.timepoint$LifeMonth, sep = ".") %in% 
                       paste(meta.rmdup$studyid, meta.rmdup$life_month, sep = ".")] <- 1
c_ids <- unique(meta.rmdup$studyid[meta.rmdup$delivery == "C-section born"])
v_ids <- unique(meta.rmdup$studyid[meta.rmdup$delivery == "Vaginally delivered"])
meta.timepoint$Delivery[meta.timepoint$StudyID %in% c_ids] <- "C-section born"
meta.timepoint$Delivery[meta.timepoint$StudyID %in% v_ids] <- "Vaginally delivered"
meta.timepoint$fill <- ifelse(meta.timepoint$Value == 0, "white", 
                              ifelse(meta.timepoint$Delivery == "C-section born", "orange", "#673147"))
meta.timepoint$Value <- factor(meta.timepoint$Value, levels = 0:1)
g1 <- ggplot(meta.timepoint, aes(LifeMonth, factor(StudyID, levels = c(c_ids, v_ids)))) +
  geom_tile(aes(fill = I(fill)), color = "#2b2b2b", size=0.125) +
  scale_x_discrete() +
  scale_y_discrete() +
  theme(axis.text.x = element_text(size=14, angle = 90)) + 
  labs(title = "ECAM Study", x = "Age (months)", y = "Infants") + 
  theme_bw(base_size=14) +
  theme(axis.title = element_text(face = "bold",size = rel(1.5)),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = rel(1.5), angle = 90),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        strip.text.x = element_text(size = rel(1.8)),
        plot.title = element_text(face='bold', size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position = c(0.85, 0.1),
        legend.background = element_rect(fill = "white", color = "black"))
g1

meta.group <- meta %>%
  group_by(life_month, delivery) %>%
  summarise(PD.mean = mean(PD), PD.se = sd(PD)) %>%
  ungroup()
pd <- position_dodge(0.2)

g2 <- ggplot(meta.group, aes(life_month, PD.mean, 
                             color = delivery, fill = delivery, group = delivery)) + 
  #facet_wrap(~baby_sex, nrow = 1) +
  geom_errorbar(aes(ymin = PD.mean-PD.se, ymax = PD.mean+PD.se),
                width = 0.5, linewidth = 0.4, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + 
  scale_fill_manual(values = c("#673147", "orange")) +
  scale_color_manual(values = c("#673147", "orange"), guide = "none") +
  labs(x = "Age (months)", y = "Phylogenetic diversity") + 
  theme_bw(base_size=14) +
  theme(axis.title = element_text(face = "bold",size = rel(1.5)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5), angle = 90),
        axis.ticks.x = element_blank(),
        axis.text = element_text(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        strip.text.x = element_text(size = rel(1.8)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position = c(0.85, 0.1),
        legend.background = element_rect(fill = "white", color = "black")) + 
  guides(fill = guide_legend(title = "Delivery Modes",
                             title.theme = element_text(face = "bold", size = 15, angle = 0),
                             label.theme = element_text(size = 15, angle = 0)))
g2

Tax <- ecam$tax.tab@.Data[,1:6]
Tax.df <- as.data.frame(Tax)
Tax.df[] <- lapply(Tax.df, function(x) gsub("\\[Mogibacteriaceae\\]", "Clos_brack_Mogibacteriaceae", x))
Tax.df[] <- lapply(Tax.df, function(x) gsub("\\[Ruminococcus\\]", "Lach_brack_Ruminococcus", x))
Tax.df[] <- lapply(Tax.df, function(x) gsub("\\[Eubacterium\\]", "Erys_brack_Eubacterium", x))
names(Tax.df) <- paste0("Rank", 1:6)
Tax.df$Rank1 <- paste0("Rank1.", Tax.df$Rank1)
Tax.df$Rank2 <- paste0("Rank2.", Tax.df$Rank2)
Tax.df$Rank3 <- paste0("Rank3.", Tax.df$Rank3)
Tax.df$Rank4 <- paste0("Rank4.", Tax.df$Rank4)
Tax.df$Rank5 <- paste0("Rank5.", Tax.df$Rank5)
Tax.df$Rank6 <- paste0("Rank6.", Tax.df$Rank6)
Tax.df$pathString <- paste(Tax.df$Rank1, Tax.df$Rank2,
                           Tax.df$Rank3, Tax.df$Rank4,
                           Tax.df$Rank5, Tax.df$Rank6, sep = "/")
Tax.tree <- as.Node(Tax.df, pathDelimiter = "/")
tree.vis <- as.phylo.Node(Tax.tree)
tree.vis$node.label <- gsub("Clos_brack_Mogibacteriaceae", "\\[Mogibacteriaceae\\]", tree.vis$node.label)
tree.vis$node.label <- gsub("Lach_brack_Ruminococcus", "\\[Ruminococcus\\]", tree.vis$node.label)
tree.vis$node.label <- gsub("Erys_brack_Eubacterium", "\\[Eubacterium\\]", tree.vis$node.label)

cladedat <- data.frame(id= nodeid(tree.vis, "Rank5.Veillonellaceae"), extendto=110)

g3 <- ggtree(tree.vis, layout="circular", ) +
  geom_point2(shape=21, size=2, fill = "#2b2b2b") +
  geom_hilight(data=cladedat, mapping=aes(node=id, extendto =extendto), alpha=0.3, fill="#808080", color="#808080")
g3 <- rotate_tree(g3, angle = 100)
# f <- tempfile(fileext=".png")
# ggsave(filename = f, plot = g3, width=7, height=7)
# g3.trim <- image_trim(image_read(f, density = 300))

id <- ecam$meta$studyid
tax <- ecam$tax.tab@.Data[,1:6]
colnames(tax) <- paste0("Rank", 1:6)
Y <- t(ecam$otu.tab)
colOrders <- order(colMeans(Y/rowSums(Y), na.rm = T), decreasing = T)
Y <- Y[,colOrders]
tax <- tax[colOrders,]

pos <- order(id, decreasing = F)
Y <- Y[pos,,drop=F]
meta <- ecam$meta[pos,,drop=F]

keep = which(colSums(Y)>0)
Y = Y[, keep, drop=FALSE]
tax = tax[keep, ,drop=FALSE]
W.data = data.table(data.frame(tax, t(Y)))
n.rank = ncol(tax)
otucols = names(W.data)[-(1:n.rank)]

k <- n.rank-5 # n.rank - 1
Rank.low <- paste0("Rank",n.rank-k)  # change the low rank
Rank.high <- paste0("Rank",n.rank-k+1) # change the high rank
tt <- W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
setnames(tt, 1:2, c(Rank.low, Rank.high))
W.tax <- as.vector(unlist(tt[, Rank.low, with=FALSE]))
W.count <- data.matrix(tt[, otucols, with=FALSE])
Y.tmp <- t(W.count[which(W.tax == "Veillonellaceae"), , drop=FALSE]) # Bacteria
id.keep <- which(rowSums(Y.tmp) != 0)
Y.tmp <- Y.tmp[id.keep,]
colOrder <- order(colMeans(Y.tmp/rowSums(Y.tmp)), decreasing = TRUE)
Y.tmp <- Y.tmp[,colOrder]

meta.tmp <- data.frame(expo = unlist(meta[id.keep, "delivery"]), 
                       month = factor(unlist(meta[id.keep, "life_month"]), 
                                      levels = c(0:22, 24), labels = 1:24),
                       Y.tmp/rowSums(Y.tmp))
colnames(meta.tmp) <- c("exposure", "life_month", unlist(tt[W.tax == "Veillonellaceae", ..Rank.high])[colOrder]) # Rank.high

meta.group <- meta.tmp %>%
  group_by(exposure, life_month) %>%
  summarise_each(mean) %>%
  ungroup()

meta.group.l <- tidyr::gather(meta.group, taxa_child, count, unlist(tt[W.tax == "Veillonellaceae", ..Rank.high])[colOrder], factor_key=FALSE)
stackOrder <- unique(meta.group.l$taxa_child)
colorDef <- c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf")
id.unc <- grep("Unclassified", stackOrder)
newColor <- numeric(length(stackOrder))
if (length(id.unc) > 0) {
  newColor[id.unc] <- "#999999"
  newColor[-id.unc] <- colorDef[1:(length(stackOrder)-1)]
}else{
  newColor <- colorDef[1:length(stackOrder)]
}
meta.group.l$taxa_child <- factor(meta.group.l$taxa_child, levels=stackOrder, labels = c("g_Veillonella", "g_Dialister", "g_Phascolarctobacterium", "g_Megasphaera"))
meta.group.l$exposure[meta.group.l$exposure == 0] <- "C-section"
meta.group.l$exposure[meta.group.l$exposure == 1] <- "Vaginal birth"
g4 <- ggplot(meta.group.l, 
             aes(fill=taxa_child, y=count, x=life_month)) +
  geom_bar(position="fill", color="black", stat="identity", width=1, linewidth=0.3) +
  scale_fill_manual(values = newColor) + 
  scale_x_discrete() +
  facet_wrap( ~exposure, scale="free_x") +
  labs(title = "", x = "Age (months)", y = "Relative abundance") + # input[5] delivery mode
  theme_bw(base_size=20) +
  theme(axis.title = element_text(face = "bold",size = rel(1.5)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.1), angle = 90),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        strip.text.x = element_text(size = rel(1.8)),
        plot.title = element_text(face='bold', size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position = "bottom") + 
  guides(fill = guide_legend(title = "f_Veillonellaceae",
                             title.theme = element_text(face = "bold.italic", size = 25, angle = 0),
                             label.theme = element_text(face = "italic", size = 25, angle = 0),
                             byrow = TRUE, ncol = 2))
g4

# pdf("../Figs/Fig1c_tree.pdf", width = 10, height = 10)
# print(g3)
# dev.off()
# pdf("../Figs/Fig1c_subcomp.pdf", width = 15, height = 10)
# print(g4)
# dev.off()

# g1.w <- wrap_elements(g1 + plot_annotation(title = 'a') & theme(plot.title = element_text(face = "bold", size = 40)) )
# g2.w <- wrap_elements(g2 + plot_annotation(title = 'b') & theme(plot.title = element_text(face = "bold", size = 40)) )
# g6.w <- wrap_elements(image_ggplot(g3.trim)+g4  + plot_annotation(title = 'c') & theme(plot.title = element_text(face = "bold", size = 40)) )
# pdf("../Figs/Fig1_ECAM.pdf", width = 20, height = 20)
# print((g1.w+g2.w)/g6.w)
# dev.off()

g5 <- image_read("../Figs/Fig1c_treeNsubcomp.pdf", density = 1000) %>%
  image_trim() %>%
  image_ggplot()
layout <- "
AABB
AABB
CCCC
CCCC"
pdf("../Figs/Fig1_ECAM.pdf", width = 20, height = 20)
print(g1+g2+g5+plot_layout(design = layout)+plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = "bold", size = 40)))
dev.off()

load("../Data/Deriveddata/ECAM.full.rslt.delivery.woconf.omni.rda")
rslt <- rslt.Ecam.full.delivery.woconf$AIGDMC
mean(unlist(lapply(rslt$lineage.zi, function(x) sum(x)>0))) # 54.55%
length(unlist(lapply(rslt$lineage.zi, function(x) sum(x)>0))) # 22
load("../Data/Deriveddata/ECAM.Genus.filter10.rda")
ecam <- ecam$pall
Tax <- ecam$tax.tab@.Data[,1:6]
Tax.df <- as.data.frame(Tax)

Tax.df[] <- lapply(Tax.df, function(x) gsub("\\[Mogibacteriaceae\\]", "Clos_brack_Mogibacteriaceae", x))
Tax.df[] <- lapply(Tax.df, function(x) gsub("\\[Ruminococcus\\]", "Lach_brack_Ruminococcus", x))
Tax.df[] <- lapply(Tax.df, function(x) gsub("\\[Eubacterium\\]", "Erys_brack_Eubacterium", x))

names(Tax.df) <- paste0("Rank", 1:6)
Tax.df$Rank1 <- paste0("Rank1.", Tax.df$Rank1)
Tax.df$Rank2 <- paste0("Rank2.", Tax.df$Rank2)
Tax.df$Rank3 <- paste0("Rank3.", Tax.df$Rank3)
Tax.df$Rank4 <- paste0("Rank4.", Tax.df$Rank4)
Tax.df$Rank5 <- paste0("Rank5.", Tax.df$Rank5)
Tax.df$Rank6 <- paste0("Rank6.", Tax.df$Rank6)
Tax.df$pathString <- paste(Tax.df$Rank1, Tax.df$Rank2,
                           Tax.df$Rank3, Tax.df$Rank4,
                           Tax.df$Rank5, Tax.df$Rank6, sep = "/")
Tax.tree <- as.Node(Tax.df, pathDelimiter = "/")
tree.vis <- as.phylo.Node(Tax.tree)
tree.vis$node.label <- gsub("Clos_brack_Mogibacteriaceae", "\\[Mogibacteriaceae\\]", tree.vis$node.label)
tree.vis$node.label <- gsub("Lach_brack_Ruminococcus", "\\[Ruminococcus\\]", tree.vis$node.label)
tree.vis$node.label <- gsub("Erys_brack_Eubacterium", "\\[Eubacterium\\]", tree.vis$node.label)

rawp <- rep(NA, length(tree.vis$node.label))
names(rawp) <- tree.vis$node.label
node.name <- names(rslt$`Omni-Cauchy-Resampling`$lineage.pval)
rawp[node.name] <- rslt$`Omni-Cauchy-Resampling`$lineage.pval
nodeids <- nodeid(tree.vis, node.name[p.adjust(rslt$`Omni-Cauchy-Resampling`$lineage.pval, method = "BH") < 0.1])
tree.vis$node.label <- rawp

# getMRCA(tree.vis, Tax.df$Rank7[which(Tax.df$Rank4 == "Staphylococcales")]) # 135

dat <- data.frame(ID = gsub(" ","_",Tax.df[,"Rank6"]),
                  Phylum = gsub("Rank2.", "", Tax.df[,"Rank2"]))
nodelab <- tree.vis$node.label[node.name[p.adjust(rslt$`Omni-Cauchy-Resampling`$lineage.pval, method = "BH")<0.1]]
nodelab.taxa <- sub("Rank[0-9].", "", names(nodelab))
cladedat <- data.frame(id=nodeids, class=nodelab.taxa, pos=rep(3, length(nodelab)),
                       extendto=c(105, 105, 105, 110, 115, 120, 125),
                       labangle=rep(0, length(nodelab)),
                       vjust=rep(-1, length(nodelab)),
                       hjust=rep(0.3, length(nodelab)))

p <- ggtree(tree.vis, layout="circular") +
  geom_point2(aes(subset=!isTip), shape=21, size=-log10(rawp)*3, fill = ifelse(p.adjust(rawp, "BH")<0.1, "#e41a1c", "#808080")) +
  geom_hilight(data=cladedat, mapping=aes(node=id, extendto =extendto), alpha=0.3, fill="#808080", color="#808080") 
  # geom_cladelab(data=cladedat, mapping=aes(node=id, label=class, offset.txt=pos,
  #                                          angle=labangle, vjust=vjust, hjust=hjust),
  #               barsize=NA, fontsize=5, horizontal=FALSE, fontface=4)
p <- p %<+% dat + 
  geom_tippoint(mapping=aes(fill=Phylum), shape=22, size=5,stroke=0.05, position="identity")+
  scale_fill_manual(values=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#808080", "#f781bf"))+
  guides(fill = guide_legend(title = "Phylum",
                             title.theme = element_text(face = "bold", size = 25, angle = 0),
                             label.theme = element_text(face = "italic", size = 25, angle = 0),
                             ncol = 1))
p

# tree.vis <- as.phylo.Node(Tax.tree)
# tree.vis$node.label <- gsub("Clos_brack_Mogibacteriaceae", "\\[Mogibacteriaceae\\]", tree.vis$node.label)
# tree.vis$node.label <- gsub("Lach_brack_Ruminococcus", "\\[Ruminococcus\\]", tree.vis$node.label)
# tree.vis$node.label <- gsub("Erys_brack_Eubacterium", "\\[Eubacterium\\]", tree.vis$node.label)
# rawp <- rep(NA, length(tree.vis$node.label))
# names(rawp) <- tree.vis$node.label
# node.name <- names(rslt$`Disp-Resampling`$lineage.pval)
# rawp[node.name] <- rslt$`Disp-Resampling`$lineage.pval
# nodeids <- nodeid(tree.vis, node.name[p.adjust(rslt$`Disp-Resampling`$lineage.pval, method = "BH")<0.1])
# tree.vis$node.label <- rawp
# 
# nodelab <- tree.vis$node.label[node.name[p.adjust(rslt$`Disp-Resampling`$lineage.pval, method = "BH")<0.1]]
# nodelab.taxa <- sub("Rank[0-9].", "", names(nodelab))
# cladedat <- data.frame(id=nodeids, class=nodelab.taxa, pos=rep(3, length(nodelab.taxa)),
#                        extendto=c(105, 105, 120, 120, 125))
# p <- ggtree(tree.vis, layout="circular") +
#   geom_point2(aes(subset=!isTip), shape=21, size=-log10(rawp)*3, fill = ifelse(p.adjust(rawp, "BH")<0.1, "#e41a1c", "#808080")) +
#   geom_hilight(data=cladedat, mapping=aes(node=id, extendto =extendto), alpha=0.3, fill="#808080", color="#808080")
#   # geom_cladelab(data=cladedat, mapping=aes(node=id, label=class, offset.txt=pos),
#   #               barsize=NA, fontsize=3, horizontal=FALSE, fontface=4)
# p
# p <- p %<+% dat + 
#   geom_tippoint(mapping=aes(fill=Phylum), shape=22, size=5,stroke=0.05, position="identity")+
#   scale_fill_manual(values=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#808080", "#f781bf"))
# p_disp <- p + 
#   guides(fill = guide_legend(title = "Phylum",
#                              title.theme = element_text(face = "bold", size = 25, angle = 0),
#                              label.theme = element_text(face = "italic", size = 25, angle = 0),
#                              ncol = 4))
# p4 <- ggpubr::ggarrange(p_mean, p_disp, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
pdf("../Figs/ECAM_Tree.pdf", width = 15, height = 15)
print(p)
dev.off()


####### ECAM subcomposition figures (delivery mode, genus as the lowest level) ######
# Rank4.Pseudomonadales and Rank5.Veillonellaceae, no ZI
rm(list = ls())
source("1a.AIGDMC_all_utility.R")
source("3a.utility.R")
load("../Data/Deriveddata/ECAM.Genus.filter10.rda")
load("../Data/Deriveddata/ECAM.full.rslt.delivery.woconf.omni.rda")
rslt <- rslt.Ecam.full.delivery.woconf
tmp_ps <- ecam$pall
id <- tmp_ps$meta$studyid
tax <- tmp_ps$tax.tab@.Data[,1:6]
colnames(tax) <- paste0("Rank", 1:6)
Y <- t(tmp_ps$otu.tab)
colOrders <- order(colMeans(Y/rowSums(Y), na.rm = T), decreasing = T)
Y <- Y[,colOrders]
tax <- tax[colOrders,]
exposure <- "delivery"
covariate <- NULL
X.m <- as.matrix(tmp_ps$meta[,c(covariate, exposure)])
X.index <- W.index <- length(covariate) + 1:length(exposure)

pos <- order(id, decreasing = F)
id <- id[pos]; Y <- Y[pos,,drop=F]; X.m = X.m[pos,,drop=F];  exposure = X.m[,X.index]
meta <- tmp_ps$meta[pos,,drop=F]

keep = which(colSums(Y)>0)
Y = Y[, keep, drop=FALSE]
X = cbind(1, X.m) # add the intercept term
X.index = X.index + 1
X.r = X[, -X.index, drop=FALSE]
W = cbind(1, X.m) # add the intercept term
W.index = W.index + 1
W.r = W[, -W.index, drop=FALSE]
tax = tax[keep, ,drop=FALSE]
W.data = data.table(data.frame(tax, t(Y)))
n.rank = ncol(tax)
otucols = names(W.data)[-(1:n.rank)]
rslt.Ecam.full.delivery.woconf$AIGDMC$`Omni-Cauchy-Resampling`$sig.lineage
conf_taxa <-c("Rank5.Clostridiaceae", "Rank5.Veillonellaceae", "Rank4.Pseudomonadales",
              "Rank3.Actinobacteria", "Rank2.Proteobacteria", "Rank1.Bacteria")
for (i in conf_taxa) {
  plot_subcomp(i)
}
# source("3a.utility.R")
plot_subcomp("Rank5.Ruminococcaceae")

####### ECAM (delivery mode) Tab ##########

rm(list = ls())

analyze_taxa_zero_mean_disp <- function(rslt, adj_method = "BH", fdr.alpha = 0.1) {
  # Identify significant taxa for each model and metric combination
  aigdm.mean <- names(which(p.adjust(rslt$AIGDMC$`Mean-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  gdm.mean <- names(which(p.adjust(rslt$GDMC$`Mean-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  zigdm.mean <- names(which(p.adjust(rslt$ZIGDMC$`Mean-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  
  aigdm.disp <- names(which(p.adjust(rslt$AIGDMC$`Disp-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  gdm.disp <- names(which(p.adjust(rslt$GDMC$`Disp-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  zigdm.disp <- names(which(p.adjust(rslt$ZIGDMC$`Disp-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  
  aigdm.zero <- names(which(p.adjust(rslt$AIGDMC$`Zero-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  zigdm.zero <- names(which(p.adjust(rslt$ZIGDMC$`Zero-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  
  aigdm <- names(which(p.adjust(rslt$AIGDMC$`Omni-Cauchy-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  gdm <- names(which(p.adjust(rslt$GDMC$`Omni-Cauchy-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  zigdm <- names(which(p.adjust(rslt$ZIGDMC$`Omni-Cauchy-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  # Combine significant taxa and additional taxa of interest
  taxa_int <- union(c(aigdm,zigdm), gdm)
  taxa_int <- taxa_int[order(taxa_int, decreasing = T)]
  
  # Compile adjusted p-values for selected taxa across models and metrics
  tmp <- cbind(
    rslt$AIGDMC$`Mean-Resampling`$lineage.pval[taxa_int],
    rslt$AIGDMC$`Disp-Resampling`$lineage.pval[taxa_int],
    rslt$AIGDMC$`Zero-Resampling`$lineage.pval[taxa_int],
    rslt$GDMC$`Mean-Resampling`$lineage.pval[taxa_int],
    rslt$GDMC$`Disp-Resampling`$lineage.pval[taxa_int],
    rslt$ZIGDMC$`Mean-Resampling`$lineage.pval[taxa_int],
    rslt$ZIGDMC$`Disp-Resampling`$lineage.pval[taxa_int],
    rslt$ZIGDMC$`Zero-Resampling`$lineage.pval[taxa_int]
  )
  
  colnames(tmp) <- c("AIGDM.M", "AIGDM.D", "AIGDM.Z", "GDM.M", "GDM.D", 
                     "ZIGDM.M", "ZIGDM.D", "ZIGDM.Z")
  print(signif(tmp[nrow(tmp),], 4))
  return(round(tmp, 4))
}

analyze_taxa_omni <- function(rslt, rslt.qcatc1, rslt.qcatc2, adj_method = "BH", fdr.alpha = 0.1, more_taxa = NULL) {
  # Identify significant taxa for each model and metric combination
  aigdm <- names(which(p.adjust(rslt$AIGDMC$`Omni-Cauchy-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  gdm <- names(which(p.adjust(rslt$GDMC$`Omni-Cauchy-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  zigdm <- names(which(p.adjust(rslt$ZIGDMC$`Omni-Cauchy-Resampling`$lineage.pval, method = adj_method) < fdr.alpha))
  dm <- names(which(p.adjust(rslt$DMC$lineage.pval, method = adj_method) < fdr.alpha))
  qcatc1 <- names(which(p.adjust(rslt.qcatc1$lineage.pval[2,], method = adj_method) < fdr.alpha))
  qcatc2 <- names(which(p.adjust(rslt.qcatc2$`Two-Part`[2,], method = adj_method) < fdr.alpha))

  # Combine significant taxa and additional taxa of interest
  taxa_int <- union(
    c(aigdm, gdm, zigdm),
    c(dm, more_taxa)
  )
  taxa_int <- taxa_int[order(taxa_int, decreasing = T)]
  
  print(rslt$AIGDMC$lineage.zi[taxa_int])
  
  # Compile adjusted p-values for selected taxa across models and metrics
  tmp <- cbind(
    rslt$AIGDMC$`Omni-Cauchy-Resampling`$lineage.pval[taxa_int],
    rslt$GDMC$`Omni-Cauchy-Resampling`$lineage.pval[taxa_int],
    rslt$ZIGDMC$`Omni-Cauchy-Resampling`$lineage.pval[taxa_int],
    rslt$DMC$lineage.pval[taxa_int],
    rslt.qcatc1$lineage.pval[2,][sapply(taxa_int, function(x) unlist(strsplit(x, "\\."))[2])],
    rslt.qcatc2$lineage.pval$`Two-Part`[2,][sapply(taxa_int, function(x) unlist(strsplit(x, "\\."))[2])]
  )

  colnames(tmp) <- c("AIGDM","GDM", "ZIGDM", "DM", "QCATC1", "QCATC2")
  print(signif(tmp[nrow(tmp),], 4))
  return(round(tmp, 4))
}

load("../Data/Deriveddata/ECAM.full.rslt.delivery.woconf.omni.rda")
rslt <- rslt.Ecam.full.delivery.woconf
rslt.qcatc1 <- rslt$QCAT1
rslt.qcatc2 <- rslt$QCAT2
more_taxa <- c(c("Rank5.Enterobacteriaceae", "Rank1.Bacteria"), 
               c("Rank4.Bacteroidales", "Rank4.Clostridiales", "Rank3.Actinobacteria", "Rank1.Bacteria"))
analyze_taxa_omni(rslt, rslt.qcatc1, rslt.qcatc2, adj_method = "BH", more_taxa = more_taxa)
analyze_taxa_zero_mean_disp(rslt, adj_method = "BH")