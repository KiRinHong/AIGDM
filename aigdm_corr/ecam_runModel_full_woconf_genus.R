## ===================================================================
## ECAM genus-level driver (refactored)
##
## Fits the ECAM infant-gut microbiome genus data with the *merged*
## AIGDM implementation, which supports four GEE working-correlation
## structures (independence, exchangeable, ar1, toeplitz) and selects
## the best-fitting structure per taxon/lineage by GIC (Pan 2001 QIC).
##
## Changes vs. the original two-file workflow:
##   * sources the single merged AIGDM.R (no separate _Toep variant)
##   * passes `period` (ordinal study period) and `corstr = "GIC"`
##     into every AIGDM.Cluster call, so each taxon's working
##     correlation is chosen by GIC rather than fixed a priori
##   * writes a per-lineage GIC selection table (CSV) and the set of
##     selected structures for each of the three AIGDM variants
##
## Correlation-lag convention: lag = |period_s - period_t| on the
## ORDINAL study-period index (1..5), matching the original Toeplitz
## code. Same-period replicate samples are lag 0 (Toeplitz off-diagonal
## 0; AR-1 diagonal-adjacent handled by the |Δperiod| exponent). Max
## Toeplitz lag is data-driven (max within-subject ordinal distance).
## ===================================================================

rm(list = ls())
suppressMessages({
  ## phyloseq is only needed if the tax.tab slot must be resolved as an
  ## S4 taxonomyTable; we extract it robustly below without requiring it.
  if (requireNamespace("phyloseq", quietly = TRUE)) library(phyloseq)
  library(HMP)
  library(miLineage)
  library(parallel)
  library(matrixStats)
  library(MASS)
  library(data.table)
})

## ---- paths (edit AIGDM_DIR / DATA_FILE / OUT_DIR as needed) --------
AIGDM_DIR <- "~/Downloads/AIGDM/LONGITUDINAL/code/Analysis/aigdm_corr/"
DATA_FILE <- "../Data/Deriveddata/ECAM.Genus.filter10.rda"
OUT_DIR   <- "../Data/Deriveddata"

source(file.path(AIGDM_DIR, "AIGDM.R"))       # merged implementation
source(file.path(AIGDM_DIR, "1b.QCATC_utility.R"))
source(file.path(AIGDM_DIR, "1c.DMC_utility.R"))

## ---- run configuration --------------------------------------------
seed      <- 12345
n.boot    <- 2e3
n.perm    <- 1e4
fdr.alpha <- 0.1
corstr    <- "GIC"          # per-taxon structure selection

## ---- load + shape ECAM genus data ---------------------------------
load(DATA_FILE)             # provides `ecam`
period.tag <- "all"
exposure   <- "delivery"
covariate  <- NULL
input_data <- ecam[[paste0("p", period.tag)]]

id  <- input_data$meta$studyid
## Extract the taxonomy matrix robustly: the tax.tab slot is a phyloseq
## S4 taxonomyTable that simply wraps a character matrix, so unclass()
## yields the underlying matrix whether or not phyloseq is attached.
tax.mat <- tryCatch(input_data$tax.tab@.Data,
                    error = function(e) unclass(input_data$tax.tab))
tax <- tax.mat[, 1:6]
colnames(tax) <- paste0("Rank", 1:6)
Y.m <- t(input_data$otu.tab)                 # samples x taxa
period <- input_data$meta$period             # ordinal study period (1..5)

## order taxa by mean relative abundance (as in the original driver)
colOrders <- order(colMeans(Y.m / rowSums(Y.m), na.rm = TRUE), decreasing = TRUE)
Y.m <- Y.m[, colOrders]
tax <- tax[colOrders, ]

X.m <- as.matrix(input_data$meta[, c(covariate, exposure)])
X.index <- W.index <- length(covariate) + 1:length(exposure)

## sort by (subject, period) so within-subject blocks are contiguous
pos    <- order(id, period, decreasing = FALSE)
id     <- id[pos]
Y.m    <- Y.m[pos, , drop = FALSE]
X.m    <- X.m[pos, , drop = FALSE]
period <- period[pos]
exposure.vec <- X.m[, X.index]


## ---- AIGDM variants with GIC structure selection ------------------
## GDMC   : no zero inflation (ZI.pos="no")
## AIGDMC : adaptive zero inflation (ZI.pos="adaptive")
## ZIGDMC : full zero inflation    (ZI.pos="all")
set.seed(seed)
rslt.GDMC <- AIGDM.Cluster(id, Y.m, X.m, X.index, period = period, corstr = corstr,
                           Tax = tax, min.depth = 0, test.type = "Omni", ZI.pos = "no",
                           n.boot = n.boot, n.cores = 10,
                           fdr.alpha = fdr.alpha, n.perm = n.perm)
message("#### GDMC (GIC) done ####")

set.seed(seed)
rslt.ZIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.index, period = period, corstr = corstr,
                             Tax = tax, min.depth = 0, test.type = "Omni", ZI.pos = "all",
                             n.boot = n.boot, n.cores = 10,
                             fdr.alpha = fdr.alpha, n.perm = n.perm)
message("#### ZIGDMC (GIC) done ####")

set.seed(seed)
rslt.AIGDMC <- AIGDM.Cluster(id, Y.m, X.m, X.index, period = period, corstr = corstr,
                             Tax = tax, min.depth = 0, test.type = "Omni", ZI.pos = "adaptive",
                             n.boot = n.boot, n.cores = 10,
                             fdr.alpha = fdr.alpha, n.perm = n.perm)
message("#### AIGDMC (GIC) done ####")

## ---- QCAT one-part & two-part (unchanged) --------------------------
set.seed(seed)
rslt.QCATC1 <- QCAT.Cluster(ID = id, OTU = Y.m, X = X.m, Tax = tax,
                            perm.type = "BTW", n.perm = n.perm,
                            fdr.alpha = fdr.alpha, test = "chisq")
message("#### QCATC one-part done ####")

set.seed(seed)
Y.m.rrf <- .rarefy(Y.m)
rslt.QCATC2 <- QCAT_GEE.Cluster(ID = id, OTU = Y.m.rrf, X = X.m, X.index = X.index,
                                Z = X.m, Z.index = W.index, Tax = tax,
                                perm.type = "BTW", n.perm = n.perm,
                                fdr.alpha = fdr.alpha, test = "chisq")
message("#### QCATC two-part done ####")

## ---- Dirichlet-Multinomial baseline (unchanged) -------------------
rslt.DMC <- DM.Cluster(ID = id, Y = Y.m, case = exposure.vec, Tax = tax, min.depth = 0,
                       fdr.alpha = fdr.alpha, perm.type = "BTW", n.perm = n.perm, seed = seed)
message("#### DMC done ####")

## ---- assemble + save results --------------------------------------
rslt.Ecam.full.delivery.woconf <- list(
  AIGDMC = rslt.AIGDMC, GDMC = rslt.GDMC, ZIGDMC = rslt.ZIGDMC, DMC = rslt.DMC,
  QCATC1 = rslt.QCATC1, QCATC2 = rslt.QCATC2)
save(rslt.Ecam.full.delivery.woconf,
     file = file.path(OUT_DIR, "ECAM.full.rslt.delivery.woconf.omni.avg_dec.corr_str.rda"))

## ---- GIC selection tables (CSV deliverables) ----------------------
## Combine the per-lineage GIC tables from the three AIGDM variants,
## tagging each with the model that produced it.
.collect_gic <- function(fit, model) {
  gt <- fit$gic.table
  if (is.null(gt) || nrow(gt) == 0) return(NULL)
  cbind(model = model, gt)
}
gic.all <- do.call(rbind, list(
  .collect_gic(rslt.GDMC,   "GDMC"),
  .collect_gic(rslt.AIGDMC, "AIGDMC"),
  .collect_gic(rslt.ZIGDMC, "ZIGDMC")))
if (!is.null(gic.all)) {
  write.csv(gic.all,
            file = file.path(OUT_DIR, "ECAM.genus.GIC.table.csv"),
            row.names = FALSE)
  ## selected structure per (model, lineage, taxon)
  sel <- gic.all[gic.all$selected, c("model", "lineage", "taxon", "structure", "gic")]
  write.csv(sel,
            file = file.path(OUT_DIR, "ECAM.genus.GIC.selected.csv"),
            row.names = FALSE)
  ## quick tally of chosen structures per model
  message("#### GIC structure tallies (selected) ####")
  print(table(sel$model, sel$structure))
}
message("#### ALL DONE ####")
