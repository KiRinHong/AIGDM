rm(list = ls())
library(phyloseq)
library(HMP)
library(miLineage)
library(parallel)
library(matrixStats)
library(MASS)

seed <- as.numeric(commandArgs(trailingOnly=TRUE)) 

source("AIGDM_GEE.R")
source("1b.QCATC_utility.R")
source("1c.DMC_utility.R")

n.boot <- 2e3
n.perm <- 1e4
fdr.alpha <- 0.1
load("ECAM.Genus.filter10.rda")
period <- "all"
exposure <- "delivery"
covariate <- NULL
input_data <- ecam[[paste0("p",period)]]
id <- input_data$meta$studyid
tax <- input_data$tax.tab@.Data[,1:6]
colnames(tax) <- paste0("Rank", 1:6)
Y.m <- t(input_data$otu.tab)
colOrders <- sample(1:ncol(Y.m)) # order(colMeans(Y.m/rowSums(Y.m), na.rm = T), decreasing = TRUE)
Y.m <- Y.m[,colOrders]
tax <- tax[colOrders,]

X.m <- as.matrix(input_data$meta[,c(covariate, exposure)])
X.index <- W.index <- length(covariate) + 1:length(exposure)

pos <- order(id, decreasing = F)
id <- id[pos]; Y.m <- Y.m[pos,,drop=F]; X.m = X.m[pos,,drop=F]; exposure = X.m[,X.index]

rslt.GDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                           Tax = tax, min.depth = 0, order = "random", ZI.pos = "no", zi.check.pval.lst = NULL,
                           n.boot = n.boot, n.cores = detectCores()/2, 
                           fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.GDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.GDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### GDMC test done ####")

# load("ECAM.full.rslt.delivery.woconf.LBFGS.genus.rda")
# zi.check.pval.lst <- rslt.Ecam.full.delivery.woconf.testzero$AIGDMC.check$lineage.zi.pval
rslt.AIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                             Tax = tax, min.depth = 0, order = "random", ZI.pos = "adaptive", zi.check.pval.lst = NULL,
                             n.boot = n.boot, n.cores = detectCores()/2, 
                             fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.AIGDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.AIGDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### AIGDMC test done ####")

rslt.ZIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                             Tax = tax, min.depth = 0, order = "random", ZI.pos = "all", zi.check.pval.lst = NULL,
                             n.boot = n.boot, n.cores = detectCores()/2, 
                             fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.ZIGDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.ZIGDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### ZIGDMC test done ####")

rslt.Ecam.full.delivery.woconf <- list(AIGDMC = rslt.AIGDMC, GDMC = rslt.GDMC, ZIGDMC = rslt.ZIGDMC)
save(rslt.Ecam.full.delivery.woconf, file = paste0("ECAM.full.rslt.delivery.woconf.omni.random", seed, ".rda"))
