rm(list = ls())
library(phyloseq)
library(HMP)
library(miLineage)
library(parallel)
library(matrixStats)
library(MASS)

source("AIGDM_GEE.R")
source("1b.QCATC_utility.R")
source("1c.DMC_utility.R")

seed <- 12345
n.boot <- 2e3
n.perm <- 1e4
fdr.alpha <- 0.1
load("ECAM.Genus.filter10.rda")
period <- "all"
exposure <- "delivery"
covariate <- "baby_sex"
input_data <- ecam[[paste0("p",period)]]
id <- input_data$meta$studyid
tax <- input_data$tax.tab@.Data[,1:6]
colnames(tax) <- paste0("Rank", 1:6)
Y.m <- t(input_data$otu.tab)
colOrders <- order(colMeans(Y.m/rowSums(Y.m), na.rm = T), decreasing = TRUE)
Y.m <- Y.m[,colOrders]
tax <- tax[colOrders,]

X.m <- as.matrix(input_data$meta[,c(covariate, exposure)])
X.index <- W.index <- length(covariate) + 1:length(exposure)

pos <- order(id, decreasing = F)
id <- id[pos]; Y.m <- Y.m[pos,,drop=F]; X.m = X.m[pos,,drop=F]; exposure = X.m[,X.index]

# set.seed(seed)
# rslt.QCATC1 <- QCAT.Cluster(ID = id, OTU = Y.m, X = X.m, Tax = tax,
#                             perm.type = "BTW", n.perm = n.perm, fdr.alpha = fdr.alpha, test = "chisq")
# print("#### QCATC one-part done ####")
# 
# set.seed(seed)
# Y.m.rrf <- .rarefy(Y.m)
# rslt.QCATC2 <- QCAT_GEE.Cluster(ID = id, OTU = Y.m, X = X.m, X.index = X.index, Z = X.m, Z.index = W.index,
#                                 Tax = tax, perm.type = "BTW", n.perm = n.perm, fdr.alpha = fdr.alpha, test = "chisq")
# print("#### QCATC two-part done ####")

rslt.GDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                           Tax = tax, min.depth = 0, ZI.pos = "no", zi.check.pval.lst = NULL,
                           n.boot = n.boot, n.cores = detectCores()/2, 
                           fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.GDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.GDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### GDMC test done ####")

load("ECAM.full.rslt.delivery.wconf.LBFGS.genus.rda")
zi.check.pval.lst <- rslt.Ecam.full.delivery.wconf.testzero$AIGDMC.check$lineage.zi.pval
rslt.AIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                             Tax = tax, min.depth = 0, ZI.pos = "adaptive", zi.check.pval.lst = zi.check.pval.lst,
                             n.boot = n.boot, n.cores = detectCores()/2, 
                             fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.AIGDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.AIGDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### AIGDMC test done ####")

rslt.ZIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                             Tax = tax, min.depth = 0, ZI.pos = "all", zi.check.pval.lst = NULL,
                             n.boot = n.boot, n.cores = detectCores()/2, 
                             fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.ZIGDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.ZIGDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### ZIGDMC test done ####")

rslt.DMC <- DM.Cluster(ID = id, Y = Y.m, case = exposure, Tax = tax, min.depth = 0,
                       fdr.alpha = fdr.alpha, perm.type = "BTW", n.perm = n.perm, seed = seed)
print("#### DMC done ####")
rslt.Ecam.full.delivery.wconf <- list(AIGDMC = rslt.AIGDMC, GDMC = rslt.GDMC, ZIGDMC = rslt.ZIGDMC, DMC = rslt.DMC)
save(rslt.Ecam.full.delivery.wconf, file = "ECAM.full.rslt.delivery.wconf.omni.rda")

seed <- 12345
n.boot <- 2e3
n.perm <- 1e4
fdr.alpha <- 0.1
load("ECAM.Genus.filter10.rda")
period <- "all"
exposure <- "delivery"
covariate <- "baby_sex"
input_data <- ecam[[paste0("p",period)]]
id <- input_data$meta$studyid
tax <- input_data$tax.tab@.Data[,1:6]
colnames(tax) <- paste0("Rank", 1:6)
Y.m <- t(input_data$otu.tab)
colOrders <- order(colMeans(Y.m/rowSums(Y.m), na.rm = T), decreasing = TRUE)
Y.m <- Y.m[,colOrders]
tax <- tax[colOrders,]

X.m <- as.matrix(input_data$meta[,c(covariate, exposure)])
X.index <- W.index <- length(covariate) + 1:length(exposure)

pos <- order(id, decreasing = F)
id <- id[pos]; Y.m <- Y.m[pos,,drop=F]; X.m = X.m[pos,,drop=F]; exposure = X.m[,X.index]

# set.seed(seed)
# rslt.QCATC1 <- QCAT.Cluster(ID = id, OTU = Y.m, X = X.m, Tax = tax,
#                             perm.type = "BTW", n.perm = n.perm, fdr.alpha = fdr.alpha, test = "chisq")
# print("#### QCATC one-part done ####")
# 
# set.seed(seed)
# Y.m.rrf <- .rarefy(Y.m)
# rslt.QCATC2 <- QCAT_GEE.Cluster(ID = id, OTU = Y.m, X = X.m, X.index = X.index, Z = X.m, Z.index = W.index,
#                                 Tax = tax, perm.type = "BTW", n.perm = n.perm, fdr.alpha = fdr.alpha, test = "chisq")
# print("#### QCATC two-part done ####")

rslt.GDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                           Tax = tax, min.depth = 0, ZI.pos = "no", zi.check.pval.lst = NULL,
                           n.boot = n.boot, n.cores = detectCores()/2, 
                           fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.GDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.GDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### GDMC test done ####")

load("ECAM.full.rslt.diet.wconf.LBFGS.genus.rda")
zi.check.pval.lst <- rslt.Ecam.full.diet.wconf.testzero$AIGDMC.check$lineage.zi.pval
rslt.AIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                             Tax = tax, min.depth = 0, ZI.pos = "adaptive", zi.check.pval.lst = zi.check.pval.lst,
                             n.boot = n.boot, n.cores = detectCores()/2, 
                             fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.AIGDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.AIGDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### AIGDMC test done ####")

rslt.ZIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.m, X.m, X.index,
                             Tax = tax, min.depth = 0, ZI.pos = "all", zi.check.pval.lst = NULL,
                             n.boot = n.boot, n.cores = detectCores()/2, 
                             fdr.alpha = fdr.alpha, n.perm = n.perm, seed = seed)
print(rslt.ZIGDMC$`Omni-Cauchy-Resampling`$sig.lineage)
print(rslt.ZIGDMC$`Omni-MV-Resampling`$sig.lineage)
print("#### ZIGDMC test done ####")

rslt.DMC <- DM.Cluster(ID = id, Y = Y.m, case = exposure, Tax = tax, min.depth = 0,
                       fdr.alpha = fdr.alpha, perm.type = "BTW", n.perm = n.perm, seed = seed)
print("#### DMC done ####")
rslt.Ecam.full.diet.wconf <- list(AIGDMC = rslt.AIGDMC, GDMC = rslt.GDMC, ZIGDMC = rslt.ZIGDMC, DMC = rslt.DMC)
save(rslt.Ecam.full.diet.wconf, file = "ECAM.full.rslt.diet.wconf.omni.rda")
