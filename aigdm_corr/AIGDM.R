########################################
#                                      #
#   AIGDM test for clustered data      #
#             add in v5.0              #
#                                      #
#                                      #
########################################
#' @importFrom stats cor ecdf na.omit p.adjust pcauchy qnbinom rbeta rexp rmultinom runif optim pchisq anova coef rbinom
#' @importFrom parallel detectCores parSapply
#' @importFrom utils globalVariables
#' @importFrom MASS ginv
#' @importFrom geepack geeglm
#' @importFrom abind abind
#' @import data.table
NULL

#' Adaptively Zero-Inflated Generalized Dirichlet Multinomial (AIGDM) Tests for Clustered Microbiome Data
#'
#' The AIGDM tests utilize a marginal model for analyzing clustered compositional counts, such as those from a longitudinal study.
#' Under the model, we link taxa mean, dispersion, and proportion of absence-presence in the composition to the covariate of interest.
#' We introduce an association test for each set of parameters to delineate what aspect of the microbial composition distribution is perturbed by the covariate.
#' Additionally, we offer an omnibus test that combines these individual tests.
#' This function allows users to (a) perform AIGDM tests on multivariate taxon counts; (b) perform AIGDM tests on the taxonomic tree to localize the covariate-associated lineages; and (c) assess the overall association of the microbial community with the covariate of interest.
#'
#' Within-cluster dependence is modeled with a GEE-style working correlation
#' matrix. This merged implementation supports four structures --
#' \code{"independence"}, \code{"exchangeable"}, \code{"ar1"}, and
#' \code{"toeplitz"} -- selected via the \code{corstr} argument. Setting
#' \code{corstr = "GIC"} (the default) fits all four structures for each
#' taxon/lineage and chooses the best-fitting one by a QIC-type information
#' criterion (Pan 2001), so the working correlation need not be fixed a priori.
#' The correlation lag is the absolute difference in the ordinal \code{period}
#' index between two samples from the same cluster; same-period replicate
#' samples are lag 0. The maximum Toeplitz lag is data-driven (the largest
#' within-cluster ordinal period distance observed).
#'
#' @param ID a vector contains cluster IDs.
#' @param OTU a numeric matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxon. Column name is mandatory. No missing values are allowed.
#' @param X a numeric matrix contains covariates that link to mean abundance, dispersion level, and presence-absence frequency in microbial compositions. Each column pertains to one variable. Samples in the OTU and X matrices should be identical and in the same order. No missing values are allowed.
#' @param X.index X.index is a numeric value/vector indicates the columns in X for the covariate(s) of interest. The remaining columns in X will be treated as confounders in testing the differential composition. If X.index is not specified, all the columns in X will be treated as the covariates of interest.
#' @param period an integer vector giving the ordinal study period of each sample (one per row of \code{OTU}), used to define the correlation lag for the \code{"ar1"}, \code{"toeplitz"}, and \code{"GIC"} structures. Required when \code{corstr} is one of those values; ignored (defaults to a single period) for \code{"independence"} and \code{"exchangeable"}. Must be integer-valued and the same length as \code{ID}.
#' @param corstr the GEE working correlation structure: one of \code{"independence"}, \code{"exchangeable"}, \code{"ar1"}, \code{"toeplitz"}, or \code{"GIC"} (Default). \code{"GIC"} fits all four structures per taxon/lineage and selects the best by the QIC-type criterion. The \code{"exchangeable"} and \code{"toeplitz"} paths reproduce the original \code{AIGDM.R} and \code{AIGDM_Toep.R} respectively.
#' @param Tax a matrix define the taxonomic ranks with each row corresponds to an OTU or a taxon and each column corresponds to a rank (start from the higher taxonomic level, e.g., from kingdom to genus). Row name is mandatory and should be consistent with the column name of the OTU table. Column name should be formatted as "Rank1", "Rank2", etc, representing the taxonomic levels from highest to lowest.
#' @param test.type If \code{test.type = "Omni"}, the function will construct an omnibus test for differential composition (Default), which combines the following three tests. If \code{test.type = "Mean"}, the function will test for differential mean. If \code{test.type = "Disp"}, the function will test for differential dispersion. If \code{test.type = "Zero"}, the function will test for differential proportion of absence-presence.
#' @param min.depth lower bound of sample read depth. Samples with read depth less than min.depth will be removed before the analysis.
#' @param ZI.pos If \code{ZI.pos = "adaptive"}, the function will perform a likelihood ratio test to determine whether modeling zero inflation for each taxon is necessary (Default). If \code{ZI.pos = "no"}, the function will not model zero inflation for any taxon. If \code{ZI.pos = "all"}, the function will model zero inflation for all taxa.
#' @param n.boot number of bootstrap iterations. By default (\code{ZI.pos = "adaptive"}), bootstrapped Monte Carlo will be performed to estimate the p-value of the likelihood ratio test. This test assesses the necessity of modeling zero inflation for each taxon. If \code{ZI.pos} is set to any value other than \code{"adaptive"}, \code{n.boot} would be \code{NULL}.
#' @param n.perm number of permutations. By default (\code{n.perm = NULL}), only the asymptotic test will be performed.
#' @param n.cores number of CPU cores to be used in parallel operation. By default (\code{n.cores = detectCores() - 1}), \code{n.cores} available CPU cores will be used.
#' @param fdr.alpha desired false discovery rate for multiple tests on the lineages.
#'
#' @return
#' If \code{Tax = NULL} (Default), a test is performed using all the OTUs/taxa.
#'
#' If Tax is provided, tests are performed for lineages derived from the taxonomic hierarchy. The output is a list that contains the following components:
#' \item{lineage.pval}{p-values for all lineages. If \code{test.type = "Omni"} (Default), output the omnibus test p-value as well as the p-values of the three tests combined into the omnibus test; If \code{test.type = "Mean"}, output the mean test p-value; If \code{test.type = "Disp"}, output the dispersion test p-value; If \code{test.type = "Zero"}, output the absence-presence proportion test p-value. When \code{n.perm = NULL} (Default), only the asymptotic test will be output; otherwise output both the asymptotic and permutation test p-values.}
#' \item{sig.lineage}{a vector of significant lineages. It is based on the result of the omnibus test if \code{test.type = "Omni"} (Default), mean test if \code{test.type = "Mean"}, dispersion test if \code{test.type = "Disp"}, absence-presence proportion test if \code{test.type = "Zero"}. When n.perm=NULL (Default), it is based on the asymptotic test result; otherwise it is based on the permutation test result.}
#' \item{global.pval}{p-value of the global test, which varies depending on the \code{test.type} and \code{n.perm} settings, similar to \code{sig.lineage}}
#' \item{lineage.zi}{zero inflation status for taxa within each lineage. This output is provided only when \code{ZI.pos = "adaptive"}.}
#' \item{lineage.zi.pval}{p-values for the zero inflation likelihood ratio test. This output is provided only when \code{ZI.pos = "adaptive"}.}
#' \item{corr.structure}{the working correlation structure used for each taxon/lineage. When \code{corstr} is a single structure, this is that structure for all taxa; when \code{corstr = "GIC"}, it is the GIC-selected structure per taxon (a named vector when \code{Tax = NULL}, or a per-lineage list when \code{Tax} is provided).}
#' \item{gic.table}{provided only when \code{corstr = "GIC"}. A data frame with one row per (lineage, taxon, structure) giving the quasi-log-likelihood (\code{qhat}), the penalty (\code{penalty}), the criterion value (\code{gic}), and a logical \code{selected} flag marking the minimum-GIC (chosen) structure for each taxon.}
#'
#' @keywords parametric association test for clustered data
#'
#' @importFrom stats cor ecdf na.omit p.adjust pcauchy qnbinom rbeta rexp rmultinom runif nlm pnorm qnorm var
#' @importFrom parallel detectCores parSapply
#' @import data.table
#'
#' @examples
#' \donttest{
#' data(data.toy.clustered)
#' OTU.toy = data.toy.clustered$OTU.toy
#' Tax.toy = data.toy.clustered$Tax.toy
#' meta = data.toy.clustered$meta.toy
#' case = meta[,1,drop=FALSE]
#' id = meta[,2]
#' # the OTUs should be consistent between the OTU table and the taxonomy table
#' OTU.toy.reorder = OTU.toy[,match(rownames(Tax.toy), colnames(OTU.toy))]
#' # perform AIGDM.Cluster test with the exchangeable working correlation
#' AIGDM.Cluster(ID=id, OTU=OTU.toy.reorder, X=case, X.index=1, Tax=Tax.toy,
#'               corstr = "exchangeable", n.cores = 1)
#' # or let GIC select the best structure per taxon (requires a period index)
#' # period = meta$period   # ordinal study period per sample
#' # AIGDM.Cluster(ID=id, OTU=OTU.toy.reorder, X=case, X.index=1, Tax=Tax.toy,
#' #               period = period, corstr = "GIC", n.cores = 1)
#' }
#'
#' @export
AIGDM.Cluster <- function(ID, OTU, X, X.index, period = NULL, corstr = "GIC", Tax = NULL, test.type = "Omni", min.depth = 0,
                          ZI.pos = "adaptive", n.boot = 499,
                          n.perm = NULL, n.cores = detectCores() - 1, fdr.alpha = 0.05) {
  # ID = id; OTU = Y.m; X = X.m; period = period
  # Tax = tax; min.depth = 0; test.type = "Omni"; ZI.pos = "no"
  # n.cores = 10
  
  
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      library(matrixStats)
      library(stats)
      library(MASS)
    })
    
    static_funcs <- c(
      ".BB.log.lik", ".simData.GDM", ".BBnZIBB.log.lik",
      ".loglik_bb", ".loglik_zibb", ".gpd_approx",
      ".gpd_params_est", ".gpd_pval", ".gpd_goft",
      ".gp_test", ".amle_method", ".combined_method",
      ".R1", ".R2", ".rgp",
      ".score_test_zero", ".score_test_mean", ".score_test_disp",
      ".dpg", ".var_logit", ".BlockDiag", ".RowbyRow",
      ".dma", ".dmb", ".var_ZStar", ".DVD", ".colwise.cbind",
      ".rowCumprods",
      ".toeplitz_cor", ".corr_toep1", ".corr_toep2",
      ".corr_rho_toep1", ".corr_rho_toep2",
      ".corr_zero1", ".corr_zero2", ".corr_rho1", ".corr_rho2",
      ".ar1_zero1", ".ar1_zero2", ".ar1_rho1", ".ar1_rho2",
      ".build_corr", ".n_corr_params", ".getSubDiag", ".var_mu",
      ".est_para", ".est_para_select", ".gic_aigdm",
      ".aigdm_quasi_loglik", ".aigdm_mean_A_B"
    )
    parallel::clusterExport(cl, varlist = static_funcs, envir = environment())
  } else {
    cl <- NULL
  }
  
  
  if (is.null(X.index))
    X.index = 1:ncol(X.index)
  if (any(X.index < 1) || any(X.index > ncol(X)))
    stop("X.index out of bounds")
  if (!(test.type %in% c("Mean", "Disp", "Zero", "Omni")))
    stop("test.type should be one of c('Mean', 'Disp', 'Zero', 'Omni')!")
  if (!(ZI.pos %in% c("no", "all", "adaptive")))
    stop("ZI.pos should be one of c('no', 'all', 'adaptive')!")
  if (!(corstr %in% c("independence", "exchangeable", "ar1", "toeplitz", "GIC")))
    stop("corstr should be one of c('independence','exchangeable','ar1','toeplitz','GIC')!")
  .need.period <- corstr %in% c("ar1", "toeplitz", "GIC")
  if (.need.period) {
    if (missing(period) || is.null(period))
      stop("'period' must be provided for corstr in {ar1, toeplitz, GIC}: an integer vector, one per sample row.")
    if (length(period) != nrow(OTU))
      stop("'period' must have the same length as nrow(OTU).")
    if (any(is.na(period)) || any(period != round(period)))
      stop("'period' must be an integer-valued vector (ordinal period index).")
  } else {
    if (missing(period) || is.null(period)) period <- rep(1L, nrow(OTU))
    if (length(period) != nrow(OTU))
      stop("'period', when supplied, must have the same length as nrow(OTU).")
  }
  period <- as.integer(period)
  
  if (ZI.pos == "adaptive") {
    if (is.null(n.boot))
      stop("n.boot cannot be NULL when ZI.pos = 'adaptive'!")
  } else {
    n.boot = NULL
  }
  if (!is.null(Tax)) {
    colnames(Tax) <- paste0("Rank", 1:ncol(Tax))
  }
  # Remove subjects with read depth less than min.depth
  remove.subject = which(rowSums(OTU) < min.depth)
  if (length(remove.subject) > 0) {
    print(paste("Remove", length(remove.subject), "samples with read depth less than", min.depth))
    ID = ID[-remove.subject]
    period = period[-remove.subject]
    Y = OTU[-remove.subject, , drop = FALSE]
    X = X[-remove.subject, , drop = FALSE]
  } else {
    Y = OTU
  }
  n = length(ID)
  # Add intercept terms to design matrices Create reduced design matrices without the covariate of interest (only
  # when not NULL)
  X = .move_col(X, X.index)
  X.index = (ncol(X) - length(X.index) + 1):ncol(X)
  X.index = X.index + 1
  X = cbind(1, X)
  X.r = X[, -X.index, drop = FALSE]
  # Sort the subjects based on ID, period
  pos = order(ID, period, decreasing = F)
  ID = ID[pos]
  period = period[pos]
  Y = Y[pos, , drop = F]
  X = X[pos, , drop = F]
  X.r = X.r[pos, , drop = F]
  # data-driven Toeplitz max lag = max within-subject ordinal period distance
  maxlag <- max(1L, max(tapply(period, ID, function(p) if (length(p) > 1) max(p) - min(p) else 0L)))
  maxlag <- as.integer(maxlag)
  
  # Remove taxa with zero total counts (make sure the order is same? )
  keep = which(colSums(Y) > 0)
  Y = Y[, keep, drop = FALSE]
  if (is.null(Tax)) {
    R.sel = .choose_r(fdr.alpha, 0.05)
    colOrder = order(colMeans(Y/rowSums(Y), na.rm = T), decreasing = TRUE)
    Y = Y[, colOrder]
    K = ncol(Y) - 1
    if (ZI.pos == "no") {
      zi.check.pval.wo = rep(1, K)
    } else if (ZI.pos == "all") {
      zi.check.pval.wo = rep(0, K)
    } else if (ZI.pos == "adaptive") {
      # Zero-inflation diagnostic test!  TRUE means add zero-inflation
      zi.check.pval = sapply(1:K, function(x) {
        y.tmp = cbind(Y[, x], rowSums(Y[, (x + 1):(K + 1), drop = FALSE]))
        id.tmp = ID[rowSums(y.tmp) != 0]
        y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
        if (sum(y.tmp[, 1] == 0) == 0) {
          return(1)  # no need to perform ZI test
        } else {
          return(.check_zeroinflation(id.tmp, y.tmp, n.boot = n.boot, n.cores = n.cores, R.sel, cl = cl))
        }
      })
      zi.check.pval.wo = zi.check.pval
    }
    zi.check.pval.w = p.adjust(zi.check.pval.wo, method = "BH")
    zi.check = zi.check.pval.w < 0.05
    zi.id = which(zi.check)
    names(zi.check.pval.wo) = names(zi.check) = colnames(Y)[1:K]
    .sel.lst = lapply(1:K, function(x) {
      .est_para_select(ID, X.r, X.r, X.r,
                       cbind(Y[, x], rowSums(Y[, (x + 1):(K + 1), drop = FALSE])),
                       zi.check[x], period, corstr, maxlag,
                       taxon = colnames(Y)[x])
    })
    est.para.lst = lapply(.sel.lst, function(z) z$est.para)
    corr.structure.vec = vapply(.sel.lst, function(z) z$corr.structure, character(1))
    names(corr.structure.vec) = colnames(Y)[1:K]
    gic.table = do.call(rbind, lapply(.sel.lst, function(z) z$gic.row))
    if (test.type %in% c("Mean", "Omni")) {
      asym.mean = .score_test_mean(ID, X, X.index, Y, est.para.lst, period)
      stat.mean.sum = asym.mean$stat.sum
      pval.mean.a = asym.mean$pval
    }
    if (test.type %in% c("Disp", "Omni")) {
      asym.disp = .score_test_disp(ID, X, X.index, Y, est.para.lst, period)
      stat.disp.sum = asym.disp$stat.sum
      pval.disp.a = asym.disp$pval
    }
    if (test.type %in% c("Zero", "Omni")) {
      if (length(zi.id) > 0) {
        asym.zero = .score_test_zero(ID, X, X.index, Y, est.para.lst, zi.id, period)
        stat.zero.sum = asym.zero$stat.sum
        pval.zero.a = asym.zero$pval
      } else {
        pval.zero.a = NA
      }
    }
    if (is.null(n.perm)) {
      if (test.type == "Mean") {
        pval = pval.mean.a
        names(pval) = "Mean-Asymptotic"
      } else if (test.type == "Disp") {
        pval = pval.disp.a
        names(pval) = "Disp-Asymptotic"
      } else if (test.type == "Zero") {
        pval = pval.zero.a
        names(pval) = "Zero-Asymptotic"
      } else if (test.type == "Omni") {
        if (length(zi.id) > 0) {
          pval.omni.c.a = .ACAT(c(pval.zero.a, pval.mean.a, pval.disp.a))
        } else {
          pval.omni.c.a = .ACAT(c(pval.mean.a, pval.disp.a))
        }
        pval = c(pval.zero.a, pval.mean.a, pval.disp.a, pval.omni.c.a)
        names(pval) = c("Zero-Asymptotic", "Mean-Asymptotic", "Disp-Asymptotic", "Omni-Asymptotic")
      }
    } else {
      # compute residual forming matrix Rconf (Smith's method)
      if (ncol(X.r) == 1) {
        Rconf = diag(nrow(X.r))
      } else {
        Rconf = diag(nrow(X.r)) - X.r[, -1, drop = F] %*% solve(t(X.r[, -1, drop = F]) %*% X.r[, -1, drop = F]) %*%
          t(X.r[, -1, drop = F])
      }
      if (test.type %in% c("Mean", "Omni")) {
        pval.mean.p = .score_test_perm(ID, X, X.index, Y, period, est.para.lst, stat.mean.sum, test = "Mean", Rconf, n.perm,
                                       n.cores, R.sel, NULL, cl = cl)
      }
      if (test.type %in% c("Disp", "Omni")) {
        pval.disp.p = .score_test_perm(ID, X, X.index, Y, period, est.para.lst, stat.disp.sum, test = "Disp", Rconf, n.perm,
                                       n.cores, R.sel, NULL, cl = cl)
      }
      if (test.type %in% c("Zero", "Omni")) {
        if (length(zi.id) > 0) {
          pval.zero.p = .score_test_perm(ID, X, X.index, Y, period, est.para.lst, stat.zero.sum, test = "Zero", Rconf, n.perm,
                                         n.cores, R.sel, zi.id, cl = cl)
        } else {
          pval.zero.p = NA
        }
      }
      if (test.type == "Mean") {
        pval = c(pval.mean.a, pval.zero.p)
        names(pval) = c("Mean-Asymptotic", "Mean-Resampling")
      } else if (test.type == "Disp") {
        pval = c(pval.disp.a, pval.disp.p)
        names(pval) = c("Disp-Asymptotic", "Disp-Resampling")
      } else if (test.type == "Zero") {
        pval = c(pval.zero.a, pval.zero.p)
        names(pval) = c("Zero-Asymptotic", "Zero-Resampling")
      } else if (test.type == "Omni") {
        
        # Calculate asymptotic omni p-value if n.perm is NOT NULL
        if (length(zi.id) > 0) {
          pval.omni.c.a = .ACAT(c(pval.zero.a, pval.mean.a, pval.disp.a))
        } else {
          pval.omni.c.a = .ACAT(c(pval.mean.a, pval.disp.a))
        }
        
        if (length(zi.id) > 0) {
          pval.omni.c.p = .ACAT(c(pval.zero.p, pval.mean.p, pval.disp.p))
        } else {
          pval.omni.c.p = .ACAT(c(pval.mean.p, pval.disp.p))
        }
        pval = c(pval.zero.a, pval.zero.p, pval.mean.a, pval.mean.p, pval.disp.a, pval.disp.p, pval.omni.c.a, pval.omni.c.p)
        names(pval) = c("Zero-Asymptotic", "Zero-Resampling", "Mean-Asymptotic", "Mean-Resampling", "Disp-Asymptotic",
                        "Disp-Resampling", "Omni-Asymptotic", "Omni-Resampling")
      }
    }
    if (ZI.pos == "adaptive") {
      rslt = list(pval = pval, zi = zi.check, zi.pval = zi.check.pval.wo)
    } else {
      rslt = list(pval = pval)
    }
    rslt$corr.structure = corr.structure.vec
    if (corstr == "GIC") rslt$gic.table = gic.table
  } else {
    # Perform hierarchical analysis using taxonomic information
    tax = Tax[keep, , drop = FALSE]
    if (sum(colnames(Y) != rownames(tax)) > 0) {
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in TAX table")
    }
    W.data = data.table(data.frame(tax, t(Y)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    n.level = n.rank - 1
    subtree = NULL
    for (k in 1:n.level) {
      # k = 1
      Rank.low = paste("Rank", n.rank - k, sep = "")
      Rank.high = paste("Rank", n.rank - k + 1, sep = "")
      tmp = table(tax[, n.rank - k])
      level.uni = sort(names(tmp)[which(tmp > 1)])
      m.level = length(level.uni)
      tt = W.data[, lapply(.SD, sum, na.rm = TRUE), .SDcols = otucols, by = list(get(Rank.low), get(Rank.high))]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with = FALSE]))
      W.count = data.matrix(tt[, otucols, with = FALSE])
      for (m in 1:m.level) {
        Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop = FALSE])
        remove.index = which(colSums(Y.tmp) == 0)
        if (length(remove.index) == ncol(Y.tmp)) {
          next
        } else {
          if (length(remove.index) > 0) {
            Y.tmp = Y.tmp[, -remove.index, drop = FALSE]
          }
          if (ncol(Y.tmp) == 1) {
            next
          } else {
            subtree = c(subtree, paste(Rank.low, level.uni[m], sep = "."))
          }
        }
      }
    }
    R.sel = .choose_r(fdr.alpha/length(subtree), 0.05)
    zi.check.pval.lst = list()
    curr.ind = 1
    for (k in 1:n.level) {
      Rank.low = paste("Rank", n.rank - k, sep = "")
      Rank.high = paste("Rank", n.rank - k + 1, sep = "")
      tmp = table(tax[, n.rank - k])
      level.uni = sort(names(tmp)[which(tmp > 1)])
      m.level = length(level.uni)
      tt = W.data[, lapply(.SD, sum, na.rm = TRUE), .SDcols = otucols, by = list(get(Rank.low), get(Rank.high))]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with = FALSE]))
      W.count = data.matrix(tt[, otucols, with = FALSE])
      for (m in 1:m.level) {
        Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop = FALSE])
        remove.index = which(colSums(Y.tmp) == 0)
        if (length(remove.index) == ncol(Y.tmp)) {
          next
        } else {
          if (length(remove.index) > 0) {
            Y.tmp = Y.tmp[, -remove.index, drop = FALSE]
          }
          if (ncol(Y.tmp) == 1) {
            next
          } else {
            colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
            Y.tmp = Y.tmp[, colOrder]
            K = ncol(Y.tmp) - 1
            if (ZI.pos == "no") {
              zi.check.pval.lst[[curr.ind]] = rep(1, K)
            } else if (ZI.pos == "all") {
              zi.check.pval.lst[[curr.ind]] = rep(0, K)
            } else if (ZI.pos == "adaptive") {
              print(paste0("### Diagnostic Test ### Now processing Rank ", n.rank - k, ".", m))
              # TRUE means add zero-inflation
              zi.check.pval = sapply(1:K, function(x) {
                y.tmp = cbind(Y.tmp[, x], rowSums(Y.tmp[, (x + 1):(K + 1), drop = FALSE]))
                id.tmp = ID[rowSums(y.tmp) != 0]
                y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
                if (sum(y.tmp[, 1] == 0) == 0) {
                  return(1)  # no need to perform ZI test
                } else {
                  return(.check_zeroinflation(id.tmp, y.tmp, n.boot = n.boot, n.cores = n.cores, R.sel, cl = cl))
                }
              })
              zi.check.pval.lst[[curr.ind]] = zi.check.pval
            }
            names(zi.check.pval.lst[[curr.ind]]) = unlist(tt[W.tax == level.uni[m], get(Rank.high)])[colOrder][1:K]
            curr.ind = curr.ind + 1
          }
        }
      }
    }
    zi.check.pval.all.wo = unlist(zi.check.pval.lst)
    zi.check.pval.all.w = p.adjust(zi.check.pval.all.wo, method = "BH")
    zi.check.lst = list()
    pval = NULL
    curr.ind = 0
    gic.rows.all = list()      # per-lineage GIC tables (corstr='GIC')
    corr.struct.all = list()   # per-lineage selected structures
    for (k in 1:n.level) {
      Rank.low = paste("Rank", n.rank - k, sep = "")
      Rank.high = paste("Rank", n.rank - k + 1, sep = "")
      tmp = table(tax[, n.rank - k])
      level.uni = sort(names(tmp)[which(tmp > 1)])
      m.level = length(level.uni)
      tt = W.data[, lapply(.SD, sum, na.rm = TRUE), .SDcols = otucols, by = list(get(Rank.low), get(Rank.high))]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with = FALSE]))
      W.count = data.matrix(tt[, otucols, with = FALSE])
      for (m in 1:m.level) {
        Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop = FALSE])
        remove.index = which(colSums(Y.tmp) == 0)
        if (length(remove.index) == ncol(Y.tmp)) {
          next
        } else {
          if (length(remove.index) > 0) {
            Y.tmp = Y.tmp[, -remove.index, drop = FALSE]
          }
          if (ncol(Y.tmp) == 1) {
            next
          } else {
            print(paste0("### Hypothesis Test ### Now processing Rank ", n.rank - k, ".", m))
            keep.ind = which(rowSums(Y.tmp) != 0)
            Y.tmp = Y.tmp[keep.ind, , drop = FALSE]
            ID.tmp = ID[keep.ind]
            period.tmp = period[keep.ind]
            X.r.tmp = X.r[keep.ind, , drop = FALSE]
            X.tmp = X[keep.ind, , drop = FALSE]
            colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
            Y.tmp = Y.tmp[, colOrder]
            K = ncol(Y.tmp) - 1
            zi.check = zi.check.pval.all.w[curr.ind + 1:K] < 0.05
            zi.id = which(zi.check)
            curr.ind = curr.ind + K
            .sel.lst.tax = lapply(1:K, function(x) {
              .est_para_select(ID.tmp, X.r.tmp, X.r.tmp, X.r.tmp,
                               cbind(Y.tmp[, x], rowSums(Y.tmp[, (x + 1):(K + 1), drop = FALSE])),
                               zi.check[x], period.tmp, corstr, maxlag,
                               taxon = colnames(Y.tmp)[x])
            })
            est.para.lst = lapply(.sel.lst.tax, function(z) z$est.para)
            .lin.label = paste(Rank.low, level.uni[m], sep = ".")
            .tax.names = colnames(Y.tmp)[1:K]
            if (is.null(.tax.names) || any(is.na(.tax.names)))
              .tax.names = paste0(.lin.label, ".taxon", 1:K)
            .gr = do.call(rbind, lapply(seq_len(K), function(x) {
              g = .sel.lst.tax[[x]]$gic.row; g$taxon = .tax.names[x]; g }))
            if (!is.null(.gr)) { .gr$lineage = .lin.label
              gic.rows.all[[length(gic.rows.all) + 1]] = .gr }
            .cs = vapply(.sel.lst.tax, function(z) z$corr.structure, character(1))
            names(.cs) = .tax.names
            corr.struct.all[[.lin.label]] = .cs
            if (test.type %in% c("Mean", "Omni")) {
              asym.mean = .score_test_mean(ID.tmp, X.tmp, X.index, Y.tmp, est.para.lst, period.tmp)
              stat.mean.sum = asym.mean$stat.sum
              pval.mean.a = asym.mean$pval
            }
            if (test.type %in% c("Disp", "Omni")) {
              asym.disp = .score_test_disp(ID.tmp, X.tmp, X.index, Y.tmp, est.para.lst, period.tmp)
              stat.disp.sum = asym.disp$stat.sum
              pval.disp.a = asym.disp$pval
            }
            if (test.type %in% c("Zero", "Omni")) {
              if (length(zi.id) > 0) {
                asym.zero = .score_test_zero(ID.tmp, X.tmp, X.index, Y.tmp, est.para.lst, zi.id, period.tmp)
                stat.zero.sum = asym.zero$stat.sum
                pval.zero.a = asym.zero$pval
              } else {
                pval.zero.a = NA
              }
            }
            if (test.type == "Mean") {
              pval = cbind(pval, c(pval.mean.a))
            } else if (test.type == "Disp") {
              pval = cbind(pval, c(pval.disp.a))
            } else if (test.type == "Zero") {
              pval = cbind(pval, c(pval.zero.a))
            } else if (test.type == "Omni") {
              if (length(zi.id) > 0) {
                pval.omni.c.a = .ACAT(c(pval.zero.a, pval.mean.a, pval.disp.a))
              } else {
                pval.omni.c.a = .ACAT(c(pval.mean.a, pval.disp.a))
              }
            }
            if (is.null(n.perm)) {
              
              pval = cbind(pval, c(pval.zero.a, pval.mean.a, pval.disp.a, pval.omni.c.a))
              
            } else {
              # compute residual forming matrix Rconf (Smith's method)
              if (ncol(X.r.tmp) == 1) {
                Rconf = diag(nrow(X.r.tmp))
              } else {
                Rconf = diag(nrow(X.r.tmp)) - X.r.tmp[, -1, drop = F] %*% solve(t(X.r.tmp[, -1, drop = F]) %*% X.r.tmp[, -1, drop = F]) %*% t(X.r.tmp[, -1, drop = F])
              }
              if (test.type %in% c("Mean", "Omni"))
                pval.mean.p = .score_test_perm(ID.tmp, X.tmp, X.index, Y.tmp, period.tmp, est.para.lst, stat.mean.sum, test = "Mean",
                                               Rconf, n.perm, n.cores, R.sel, NULL, cl = cl)
              if (test.type %in% c("Disp", "Omni"))
                pval.disp.p = .score_test_perm(ID.tmp, X.tmp, X.index, Y.tmp, period.tmp, est.para.lst, stat.disp.sum, test = "Disp",
                                               Rconf, n.perm, n.cores, R.sel, NULL, cl = cl)
              if (test.type %in% c("Zero", "Omni")) {
                if (length(zi.id) > 0) {
                  pval.zero.p = .score_test_perm(ID.tmp, X.tmp, X.index, Y.tmp, period.tmp, est.para.lst, stat.zero.sum, test = "Zero",
                                                 Rconf, n.perm, n.cores, R.sel, zi.id, cl = cl)
                } else {
                  pval.zero.p = NA
                }
              }
              if (test.type == "Mean") {
                pval = cbind(pval, c(pval.mean.a, pval.mean.p))
              } else if (test.type == "Disp") {
                pval = cbind(pval, c(pval.disp.a, pval.disp.p))
              } else if (test.type == "Zero") {
                pval = cbind(pval, c(pval.zero.a, pval.zero.p))
              } else if (test.type == "Omni") {
                if (length(zi.id) > 0) {
                  pval.omni.c.p = .ACAT(c(pval.zero.p, pval.mean.p, pval.disp.p))
                } else {
                  pval.omni.c.p = .ACAT(c(pval.mean.p, pval.disp.p))
                }
                pval = cbind(pval, c(pval.zero.a, pval.zero.p, pval.mean.a, pval.mean.p, pval.disp.a, pval.disp.p,
                                     pval.omni.c.a, pval.omni.c.p))
              }
            }
            zi.check.lst[[ncol(pval)]] = zi.check
          }
        }
      }  # lineage loop
    }  # level loop
    colnames(pval) = subtree
    names(zi.check.lst) = names(zi.check.pval.lst) = subtree
    # identify significant lineages
    rslt = apply(pval, 1, .identifySigLineagesNsimesTest, fdr.alpha)
    if (is.null(n.perm)) {
      if (test.type == "Mean") {
        names(rslt) = c("Mean-Asymptotic")
        lineage.pval = matrix(rslt$`Mean-Asymptotic`$lineage.pval, nrow = 1)
        rownames(lineage.pval) = c("Mean-Asymptotic")
        sig.lineage = rslt$`Mean-Asymptotic`$sig.lineage
        global.pval = rslt$`Mean-Asymptotic`$global.pval
      } else if (test.type == "Disp") {
        names(rslt) = c("Disp-Asymptotic")
        lineage.pval = matrix(rslt$`Disp-Asymptotic`$lineage.pval, nrow = 1)
        rownames(lineage.pval) = c("Disp-Asymptotic")
        sig.lineage = rslt$`Disp-Asymptotic`$sig.lineage
        global.pval = rslt$`Disp-Asymptotic`$global.pval
      } else if (test.type == "Zero") {
        names(rslt) = c("Zero-Asymptotic")
        lineage.pval = matrix(rslt$`Zero-Asymptotic`$lineage.pval, nrow = 1)
        rownames(lineage.pval) = c("Zero-Asymptotic")
        sig.lineage = rslt$`Zero-Asymptotic`$sig.lineage
        global.pval = rslt$`Zero-Asymptotic`$global.pval
      } else if (test.type == "Omni") {
        names(rslt) = c("Zero-Asymptotic", "Mean-Asymptotic", "Disp-Asymptotic", "Omni-Asymptotic")
        lineage.pval = rbind(rslt$`Omni-Asymptotic`$lineage.pval, rslt$`Mean-Asymptotic`$lineage.pval, rslt$`Disp-Asymptotic`$lineage.pval,
                             rslt$`Zero-Asymptotic`$lineage.pval)
        rownames(lineage.pval) = c("Omni-Asymptotic", "Mean-Asymptotic", "Disp-Asymptotic", "Zero-Asymptotic")
        sig.lineage = rslt$`Omni-Asymptotic`$sig.lineage
        global.pval = rslt$`Omni-Asymptotic`$global.pval
      }
    } else {
      if (test.type == "Mean") {
        names(rslt) = c("Mean-Asymptotic", "Mean-Resampling")
        lineage.pval = rbind(rslt$`Mean-Asymptotic`$lineage.pval, rslt$`Mean-Resampling`$lineage.pval)
        rownames(lineage.pval) = c("Mean-Asymptotic", "Mean-Resampling")
        sig.lineage = rslt$`Mean-Resampling`$sig.lineage
        global.pval = rslt$`Mean-Resampling`$global.pval
      } else if (test.type == "Disp") {
        names(rslt) = c("Disp-Asymptotic", "Disp-Resampling")
        lineage.pval = rbind(rslt$`Disp-Asymptotic`$lineage.pval, rslt$`Disp-Resampling`$lineage.pval)
        rownames(lineage.pval) = c("Disp-Asymptotic", "Disp-Resampling")
        sig.lineage = rslt$`Disp-Resampling`$sig.lineage
        global.pval = rslt$`Disp-Resampling`$global.pval
      } else if (test.type == "Zero") {
        names(rslt) = c("Zero-Asymptotic", "Zero-Resampling")
        lineage.pval = rbind(rslt$`Zero-Asymptotic`$lineage.pval, rslt$`Zero-Resampling`$lineage.pval)
        rownames(lineage.pval) = c("Zero-Asymptotic", "Zero-Resampling")
        sig.lineage = rslt$`Zero-Resampling`$sig.lineage
        global.pval = rslt$`Zero-Resampling`$global.pval
      } else if (test.type == "Omni") {
        names(rslt) = c("Zero-Asymptotic", "Zero-Resampling", "Mean-Asymptotic", "Mean-Resampling", "Disp-Asymptotic",
                        "Disp-Resampling", "Omni-Asymptotic", "Omni-Resampling")
        lineage.pval = rbind(rslt$`Omni-Asymptotic`$lineage.pval, rslt$`Omni-Resampling`$lineage.pval, rslt$`Mean-Asymptotic`$lineage.pval,
                             rslt$`Mean-Resampling`$lineage.pval, rslt$`Disp-Asymptotic`$lineage.pval, rslt$`Disp-Resampling`$lineage.pval,
                             rslt$`Zero-Asymptotic`$lineage.pval, rslt$`Zero-Resampling`$lineage.pval)
        rownames(lineage.pval) = c("Omni-Asymptotic", "Omni-Resampling", "Mean-Asymptotic", "Mean-Resampling", "Disp-Asymptotic",
                                   "Disp-Resampling", "Zero-Asymptotic", "Zero-Resampling")
        sig.lineage = rslt$`Omni-Resampling`$sig.lineage
        global.pval = rslt$`Omni-Resampling`$global.pval
      }
    }
    # rslt = c(rslt, lineage.zi = list(zi.check.lst), lineage.zi.pval = list(zi.check.pval.lst))
    if (ZI.pos == "adaptive") {
      rslt = list(lineage.pval = lineage.pval, sig.lineage = sig.lineage, global.pval = global.pval, lineage.zi = zi.check.lst,
                  lineage.zi.pval = zi.check.pval.lst)
    } else {
      rslt = list(lineage.pval = lineage.pval, sig.lineage = sig.lineage, global.pval = global.pval)
    }
    rslt$corr.structure = corr.struct.all
    if (corstr == "GIC" && length(gic.rows.all) > 0) {
      gic.table = do.call(rbind, gic.rows.all)
      gic.table = gic.table[, c("lineage", "taxon", "structure", "gic", "qhat", "penalty", "selected")]
      rownames(gic.table) = NULL
      rslt$gic.table = gic.table
    }
  }
  return(rslt)
}

.move_col <- function(mat, x) {
  indices <- 1:ncol(mat)
  # Remove the xth index and append it at the end
  new_indices <- c(indices[-x], indices[x])
  # Rearrange the columns based on new indices
  return(mat[, new_indices, drop = FALSE])
}

.identifySigLineagesNsimesTest <- function(pval, fdr.alpha) {
  subtree.tmp = names(pval)
  score.tmp = pval
  index.na = which(is.na(score.tmp))
  if (length(index.na) == length(pval)) {
    return(list(lineage.pval = rep(NA, length(pval)), sig.lineage = NULL, global.pval = NA))
  }
  if (length(index.na) > 0) {
    score.tmp = score.tmp[-index.na]
    subtree.tmp = subtree.tmp[-index.na]
  }
  # score.tmp[score.tmp==0] = 1e-4
  m.test = length(score.tmp)
  # Benjamini-Hochberg FDR control
  index.p = order(score.tmp)
  p.sort = sort(score.tmp)
  reject = rep(0, m.test)
  tmp = which(p.sort <= (1:m.test) * fdr.alpha/m.test)
  if (length(tmp) > 0) {
    index.reject = index.p[1:max(tmp)]
    reject[index.reject] = 1
  }
  sig.lineage = subtree.tmp[reject == 1]
  # perform global test
  global.pval = .simes.test(score.tmp)
  names(global.pval) = "Simes"
  return(list(lineage.pval = pval, sig.lineage = sig.lineage, global.pval = global.pval))
}

.simes.test <- function(x) {
  return(min(length(x) * x/rank(x)))
}

.ZIeZ_vec <- function(pZv, av, bv, Y) {
  K = length(Y) - 1
  N = sum(Y)
  # Vectorized operations for faster computations
  av.prim = av + Y[1:K]
  bv.prim = bv + (N - cumsum(Y[1:K]))
  pZv.post = rep(0, K)
  # Compute pZv.post values in a vectorized way
  beta_values = ifelse(beta(av, bv) == 0, 1, beta(av.prim, bv.prim)/beta(av, bv))
  tmp_values = pZv + (1 - pZv) * beta_values
  pZv.post[Y[1:K] == 0] = ifelse(tmp_values[Y[1:K] == 0] == 0, 1, pZv[Y[1:K] == 0]/tmp_values[Y[1:K] == 0])
  return(list(pv.post = pZv.post, av.post = av.prim, bv.post = bv.prim))
}

.LogitNegLoglik <- function(par, data){
  tmp = as.numeric( exp(data$W %*% par) )
  p = tmp/(1+tmp)
  tmp = data$Del * log(p) + (1-data$Del) * log(1-p)
  index = which(p==0 | p==1)
  if(length(index)>0){
    tmp[index] = 0
  }
  return( -sum(tmp) )
}

.rowCumsums <- function(x) {
  if (!is.matrix(x))
    stop("Must be matrix")
  n_rows <- nrow(x)
  n_cols <- ncol(x)
  res <- matrix(0, nrow = n_rows, ncol = n_cols)
  for (i in 1:n_rows) {
    current_sum <- 0
    for (j in 1:n_cols) {
      current_sum <- current_sum + x[i, j]
      res[i, j] <- current_sum
    }
  }
  return(res)
}

.ZIeZ_mat <- function(pZv, av, bv, Y) {
  n = dim(Y)[1]
  K = dim(Y)[2] - 1
  # Vectorized operations for faster computations
  av.prim = av + Y[, 1:K, drop = FALSE]
  bv.prim = bv + (rowSums(Y) - .rowCumsums(Y[, 1:K, drop = FALSE]))
  pZv.post = matrix(0, n, K)
  # Compute pZv.post values in a vectorized way
  beta_values = ifelse(beta(av, bv) == 0, 1, beta(av.prim, bv.prim)/beta(av, bv))
  tmp_values = pZv + (1 - pZv) * beta_values
  zero_indices = Y[, 1:K, drop = FALSE] == 0
  pZv.post[zero_indices] = ifelse(tmp_values[zero_indices] == 0, 1, pZv[zero_indices]/tmp_values[zero_indices])
  return(list(pv.post = pZv.post, av.post = av.prim, bv.post = bv.prim))
}
.LogitOptim <- function(Del, W, gamma.ini) {
  # Del=DelZ.R[,1]; gamma.ini=gamma.last[,1]
  Logit.par.ini = gamma.ini
  Logit.data = list(Del = Del, W = W)
  return(optim(par = Logit.par.ini, fn = .LogitNegLoglik, gr = .LogitNegScore, data = Logit.data, method = "L-BFGS-B")$par)
}

.LogitNegScore <- function(par, data) {
  tmp = as.numeric(exp(data$W %*% par))
  p = tmp/(1 + tmp)
  return(-colSums((data$Del - p) * data$W))
}
.AIBetaOptim <- function(Del, A, B, Xa, Xb, alpha.ini, beta.ini) {
  # Del=Del.R;A=A.R;B=B.R;alpha.ini=anew;beta.ini=c(beta0)
  Beta.par.ini = c(alpha.ini, beta.ini)
  Beta.data = list(Del = Del, A = A, B = B, Xa = Xa, Xb = Xb)
  return(optim(par = Beta.par.ini, fn = .AIBetaNegLoglik, gr = .AIBetaNegScore, data = Beta.data, method = "BFGS")$par)
}
.AIBetaNegLoglik <- function(par, data) {
  # par = Beta.par.ini; data = Beta.data
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  tmp = exp(data$Xa %*% alpha)
  mu.tmp = as.numeric(tmp/(1 + tmp))
  tmp = exp(data$Xb %*% beta)
  phi.tmp = as.numeric(tmp/(1 + tmp))
  a = pmax(0, (1/phi.tmp - 1) * mu.tmp)
  b = pmax(0, (1/phi.tmp - 1) * (1 - mu.tmp))
  return(-sum((1 - data$Del) * (-lbeta(a, b) + data$A * (a - 1) + data$B * (b - 1))))
}
.AIBetaNegScore <- function(par, data) {
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  tmp.a = c(exp(data$Xa %*% alpha))  # matrix to vector
  mu.tmp = as.numeric(tmp.a/(1 + tmp.a))
  tmp.b = c(exp(data$Xb %*% beta))
  phi.tmp = as.numeric(tmp.b/(1 + tmp.b))
  a = pmax(0, (1/phi.tmp - 1) * mu.tmp)
  b = pmax(0, (1/phi.tmp - 1) * (1 - mu.tmp))
  zstar = data$A - data$B
  zdagger = data$B
  mustar = digamma(a) - digamma(b)
  mudagger = digamma(b) - digamma(a + b)
  one_minus_del = 1 - data$Del
  return(-c(colSums(one_minus_del * (1/tmp.b) * mu.tmp * (1/(1 + tmp.a)) * (zstar - mustar) * data$Xa),
            colSums(one_minus_del * (-1/tmp.b) * (mu.tmp * (zstar - mustar) + (zdagger - mudagger)) * data$Xb)))
}

# maximum likelihood estimators of gamma, alpha and phi
.AIGDM_EM <- function(Y, W, Xa, Xb, gamma0, alpha0, beta0, tol = 1e-04, max.iter = 1000) {
  n = nrow(Y)
  K = ncol(Y) - 1
  da = ncol(Xa)
  db = ncol(Xb)
  gamma.last = gamma.now = gamma0
  alpha.last = alpha.now = alpha0
  beta.last = beta.now = beta0
  Del.R = A.R = B.R = matrix(0, n, K)
  for (l in 1:max.iter) {
    # E-step print(paste('====== ', l, 'th ======', sep=''))
    expXa_alpha = exp(Xa %*% alpha.last)
    expXb_beta = exp(Xb %*% beta.last)
    expW_gamma = exp(W %*% gamma.last)
    # Compute tmpMv, tmpSv using vectorized operations
    tmpMv = expXa_alpha/(1 + expXa_alpha)
    tmpSv = expXb_beta/(1 + expXb_beta)
    # Compute tmpA and tmpB using vectorized operations
    tmpA = tmpMv * (1/tmpSv - 1)
    tmpB = (1 - tmpMv) * (1/tmpSv - 1)
    # Compute tmpZ using vectorized operations
    tmpZ = expW_gamma/(1 + expW_gamma)
    # Fix infinite values
    tmpZ[is.na(tmpZ)] = 0
    tmpZ[is.infinite(tmpZ) & tmpZ > 0] = 1
    par.post = .ZIeZ_mat(tmpZ, tmpA, tmpB, Y)
    tmp = par.post$av.post + par.post$bv.post
    Del.R = par.post$pv.post
    A.R = digamma(par.post$av.post) - digamma(tmp)
    B.R = digamma(par.post$bv.post) - digamma(tmp)
    # M-step
    for (j in 1:K) {
      if (!is.infinite(gamma.last[1, j])) {
        gamma.now[, j] = .LogitOptim(Del.R[, j], W, gamma.last[, j])
      }
      tmp = .AIBetaOptim(Del.R[, j], A.R[, j], B.R[, j], Xa, Xb, alpha.last[, j], beta.last[, j])
      alpha.now[, j] = tmp[1:da]
      beta.now[, j] = tmp[-(1:da)]
    }
    diffs = abs(c(gamma.now - gamma.last, alpha.now - alpha.last, beta.now - beta.last))
    if (sum(diffs[!is.infinite(diffs)], na.rm = TRUE) < tol)
      break
    gamma.last = gamma.now
    alpha.last = alpha.now
    beta.last = beta.now
  }
  return(list(gamma.est = gamma.now, alpha.est = alpha.now, beta.est = beta.now))
}

# .dma: dmu*/dalpha, first derivative of alpha in the mu*=digamma(mu*phi)-digamma((1-mu)*phi)
.dma <- function(xmat, alpha, mu, phi) {
  K = ncol(alpha)
  tmp = c(exp(xmat %*% alpha))  # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  xmat_expand = do.call(rbind, replicate(K, xmat, simplify = FALSE))
  tmp1 = (1 + tmp)^2
  sigma = 1/phi - 1
  return(.RowbyRow(xmat_expand, tmp/tmp1 * sigma * (trigamma(mu * sigma) + trigamma((1 - mu) * sigma))))
}

# .dmb: dmu*/dbeta, first derivative of beta in the mu*=digamma(mu*phi)-digamma((1-mu)*phi)
.dmb <- function(xmat, beta, mu, phi) {
  K = ncol(beta)
  tmp = c(exp(xmat %*% beta))  # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  xmat_expand = do.call(rbind, replicate(K, xmat, simplify = FALSE))
  sigma = 1/phi - 1
  return(.RowbyRow(xmat_expand, (-1/tmp) * sigma * (trigamma(mu * sigma) * mu - trigamma((1 - mu) * sigma) * (1 - mu))))
}

# .ddb: ddisp/dbeta, first derivative of beta in the dispersion model xmat 2*2 alpha 2*5
.ddb <- function(xmat, beta) {
  K = ncol(beta)
  tmp = c(exp(xmat %*% beta))  # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  xmat_expand = do.call(rbind, replicate(K, xmat, simplify = FALSE))
  tmp1 = (1 + tmp)^2
  return(.RowbyRow(xmat_expand, tmp/tmp1))
}

.dpg <- function(wmat, gamma) {
  K = ncol(gamma)
  tmp = c(exp(wmat %*% gamma))  # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  tmp[is.nan(tmp)] = 0
  wmat_expand = do.call(rbind, replicate(K, wmat, simplify = FALSE))
  tmp1 = (1 + tmp)^2
  return(.RowbyRow(wmat_expand, tmp/tmp1))
}

# ----------------------------------------------------------------
# Build a Toeplitz correlation matrix from observed period indices.
# periods_i: integer vector of observed period labels (e.g. c(1,3,5))
# rho_vec:   length-4 vector, rho_vec[h] = correlation at period-distance h
# Entry [s,t] = rho_vec[ |periods_i[s] - periods_i[t]| ]
# ----------------------------------------------------------------

## ===================================================================
## Unified working-correlation dispatch layer  (merged AIGDM)
## Structures: "independence","exchangeable","ar1","toeplitz"
## Lag metric: ordinal period distance h = |period_s - period_t|.
## Toeplitz max lag is data-driven; see AIGDM.Cluster (.aigdm.maxlag).
## ===================================================================

## Build a subject working correlation matrix.
##  corstr    : structure tag
##  params    : exchangeable/ar1 -> scalar ; toeplitz -> length-maxlag vector ;
##              independence -> ignored
##  periods_i : integer ordinal periods for this subject
## Number of correlation parameters for a structure (for sizing + GIC bookkeeping).
##  max.lag: data-driven max ordinal period distance within subjects.
.n_corr_params <- function(corstr, max.lag = NULL) {
  switch(corstr,
         independence = 0L,
         exchangeable = 1L,
         ar1          = 1L,
         toeplitz     = as.integer(max.lag),
         GIC          = as.integer(max.lag),
         stop("Unknown corstr: ", corstr))
}

.build_corr <- function(corstr, params, periods_i) {
  n.T <- length(periods_i)
  if (n.T == 1L) return(matrix(1, 1, 1))
  if (corstr == "independence") {
    diag(n.T)
  } else if (corstr == "exchangeable") {
    ## identical to original: diag(rho,1) %x% J + diag(1 - rep(rho,each=n.T))
    diag(params[1], 1) %x% matrix(1, n.T, n.T) + diag(1 - rep(params[1], each = n.T), n.T)
  } else if (corstr == "ar1") {
    h <- abs(outer(periods_i, periods_i, "-"))
    params[1]^h
  } else if (corstr == "toeplitz") {
    .toeplitz_cor(params, periods_i)
  } else {
    stop("Unknown corstr: ", corstr)
  }
}

## AR-1 moment estimator, ZERO (logit) part.
## Pools standardized residual cross-products over lag-1-adjacent ordinal pairs
## (h == 1) to estimate a single decay rho; returns num & den scalars.
.ar1_zero1 <- function(eDel, p, periods_i) {
  n.T <- length(periods_i)
  if (n.T < 2L) return(0)
  pairs <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  s <- pairs[,1]; t <- pairs[,2]
  h <- abs(periods_i[s] - periods_i[t])
  val <- (eDel[s]-p[s])*(eDel[t]-p[t])/sqrt(p[s]*(1-p[s])*p[t]*(1-p[t]))
  val[!is.finite(val)] <- 0
  sum(val[h == 1L])
}
.ar1_zero2 <- function(eDel, p, periods_i) {
  n.T <- length(periods_i)
  if (n.T < 2L) return(0)
  pairs <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  h <- abs(periods_i[pairs[,1]] - periods_i[pairs[,2]])
  sum(h == 1L)
}
.ar1_rho1 <- function(eDel, ezstar, mustar, mu, phi, periods_i) {
  n.T <- length(periods_i)
  if (n.T < 2L) return(0)
  varzstar <- diag(.var_ZStar(mu, phi))
  pairs <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  s <- pairs[,1]; t <- pairs[,2]
  h <- abs(periods_i[s] - periods_i[t])
  val <- (1-eDel[s])*(1-eDel[t])*(ezstar[s]-mustar[s])*(ezstar[t]-mustar[t])/
    sqrt(varzstar[s]*varzstar[t])
  val[!is.finite(val)] <- 0
  sum(val[h == 1L])
}
.ar1_rho2 <- function(eDel, ezstar, mustar, mu, phi, periods_i) {
  n.T <- length(periods_i)
  if (n.T < 2L) return(0)
  pairs <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  h <- abs(periods_i[pairs[,1]] - periods_i[pairs[,2]])
  sum(h == 1L)
}

## ---- exchangeable moment helpers (verbatim from symmetric AIGDM.R) ----
.corr_zero1 <- function(eDel, p, n.T) {
  uist = NULL
  K = length(p)/n.T
  for (i in seq(1, n.T * K, n.T)) {
    for (s in i:(i + n.T - 2)) {
      for (t in (s + 1):(i + n.T - 1)) {
        uist = c(uist, (eDel[s] - p[s]) * (eDel[t] - p[t])/sqrt(p[s] * (1 - p[s]) * p[t] * (1 - p[t])))
      }
    }
  }
  uist[is.infinite(uist) | is.nan(uist)] = 0
  uist.mat = matrix(uist, nrow = choose(n.T, 2), ncol = K)
  return(colSums(uist.mat, na.rm = T))
}
.corr_zero2 <- function(eDel, p, n.T) {
  K = length(p)/n.T
  res = (eDel - p)^2/(p * (1 - p))
  res[is.infinite(res) | is.nan(res)] = 0
  res.mat = matrix(res, nrow = n.T, ncol = K)
  return(colSums(res.mat, na.rm = T))
}
.corr_rho1 <- function(eDel, ezstar, mustar, mu, phi, n.T) {
  K = length(mu)/n.T
  uist = NULL
  varzstar = diag(.var_ZStar(mu, phi))
  for (i in seq(1, n.T * K, n.T)) {
    for (s in i:(i + n.T - 2)) {
      for (t in (s + 1):(i + n.T - 1)) {
        uist = c(uist, (1 - eDel[s]) * (1 - eDel[t]) * (ezstar[s] - mustar[s]) * (ezstar[t] - mustar[t])/sqrt(varzstar[s] * varzstar[t]))
      }
    }
  }
  uist[is.infinite(uist) | is.nan(uist)] = 0
  uist.mat = matrix(uist, nrow = choose(n.T, 2), ncol = K)
  return(colSums(uist.mat, na.rm = T))
}
.corr_rho2 <- function(eDel, ezstar, mustar, mu, phi, n.T) {
  K = length(mu)/n.T
  res = (1 - eDel)^2 * (ezstar - mustar)^2/diag(.var_ZStar(mu, phi))
  res[is.infinite(res) | is.nan(res)] = 0
  res.mat = matrix(res, nrow = n.T, ncol = K)
  return(colSums(res.mat, na.rm = T))
}

## ---- helpers only present in symmetric AIGDM.R (needed by exch path) ----
.getSubDiag <- function(m) {
  diag(m) = 0
  return(m/2)
  # diag(m[-nrow(m),-1], nrow(m)-1)
}
.var_mu <- function(mu) {
  return(diag(c(mu * (1 - mu)), length(mu)))
}

.toeplitz_cor <- function(rho_vec, periods_i) {
  h_mat <- abs(outer(periods_i, periods_i, "-"))      # period-distance matrix
  h_mat <- pmin(h_mat, length(rho_vec))               # cap at max lag
  mat   <- matrix(c(0, rho_vec)[h_mat + 1L],          # index rho_vec by distance
                  nrow = length(periods_i))
  diag(mat) <- 1
  return(mat)
}

# ----------------------------------------------------------------
# Toeplitz moment estimator numerator for sig (zero model)
# Accumulates cross-products for each lag h = 1..4 separately
# Returns a vector of length 4 (one element per lag)
# ----------------------------------------------------------------
.corr_toep1 <- function(eDel, p, periods_i, n.lags = 4L) {
  n.T    <- length(periods_i)
  # All unique pairs (s < t)
  pairs  <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  s_idx  <- pairs[, 1]; t_idx <- pairs[, 2]
  h      <- abs(periods_i[s_idx] - periods_i[t_idx])
  valid  <- h <= n.lags
  resid_s <- (eDel[s_idx] - p[s_idx])
  resid_t <- (eDel[t_idx] - p[t_idx])
  var_s   <- p[s_idx] * (1 - p[s_idx])
  var_t   <- p[t_idx] * (1 - p[t_idx])
  val     <- resid_s * resid_t / sqrt(var_s * var_t)
  val[!is.finite(val)] <- 0
  num <- vapply(1:n.lags, function(hh) sum(val[valid & h == hh]), numeric(1))
  return(num)
}

# ----------------------------------------------------------------
# Toeplitz moment estimator denominator for sig (zero model)
# Returns a vector of length 4
# ----------------------------------------------------------------
.corr_toep2 <- function(eDel, p, periods_i, n.lags = 4L) {
  n.T    <- length(periods_i)
  pairs  <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  h      <- abs(periods_i[pairs[,1]] - periods_i[pairs[,2]])
  valid  <- h <= n.lags
  den <- vapply(1:n.lags, function(hh) sum(valid & h == hh), numeric(1))
  return(den)
}

# ----------------------------------------------------------------
# Toeplitz moment estimator numerator for rho (mean model)
# Returns a vector of length 4
# ----------------------------------------------------------------
.corr_rho_toep1 <- function(eDel, ezstar, mustar, mu, phi, periods_i, n.lags = 4L) {
  n.T      <- length(periods_i)
  varzstar <- diag(.var_ZStar(mu, phi))
  pairs    <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  s_idx    <- pairs[, 1]; t_idx <- pairs[, 2]
  h        <- abs(periods_i[s_idx] - periods_i[t_idx])
  valid    <- h <= n.lags
  val <- (1 - eDel[s_idx]) * (1 - eDel[t_idx]) *
    (ezstar[s_idx] - mustar[s_idx]) * (ezstar[t_idx] - mustar[t_idx]) /
    sqrt(varzstar[s_idx] * varzstar[t_idx])
  val[!is.finite(val)] <- 0
  num <- vapply(1:n.lags, function(hh) sum(val[valid & h == hh]), numeric(1))
  return(num)
}

# ----------------------------------------------------------------
# Toeplitz moment estimator denominator for rho (mean model)
# Returns a vector of length 4
# ----------------------------------------------------------------
.corr_rho_toep2 <- function(eDel, ezstar, mustar, mu, phi, periods_i, n.lags = 4L) {
  n.T    <- length(periods_i)
  pairs  <- which(upper.tri(matrix(0, n.T, n.T)), arr.ind = TRUE)
  h      <- abs(periods_i[pairs[,1]] - periods_i[pairs[,2]])
  valid  <- h <= n.lags
  den <- vapply(1:n.lags, function(hh) sum(valid & h == hh), numeric(1))
  return(den)
}


# .var_Z: variance of Z in mean model (matrix Bi=Diag(mu_i*(1-mu_i)/(phi_i+1)))
.var_ZStar <- function(mu, phi) {
  sigma = (1/phi - 1)
  return(diag(trigamma(mu * sigma) + trigamma((1 - mu) * sigma), length(mu)))
}

# .var_logit
.var_logit <- function(p) {
  return(diag(p * (1 - p), length(p)))
}

# .BlockDiag: change 10*1/2/3 matrix to 10*5/10/15 diagonal block matrix
.BlockDiag <- function(mat, K, n.T) {
  tmp_nx = ncol(mat)
  J = matrix(1, n.T, tmp_nx)  # n.T denotes # of time points
  mat_expand = do.call(cbind, replicate(K, mat, simplify = FALSE))
  return(diag(K) %x% J * mat_expand)
}

.colwise.cbind <- function(mat1, mat2, para.ind) {
  tmp = matrix(0, nrow = length(para.ind), ncol = ncol(mat1) + ncol(mat2))
  tmp[, -para.ind] = mat1
  tmp[, para.ind] = mat2
  return(tmp)
}

.DVD <- function(xamat, xbmat, alpha, beta, pZ, mu, phi, one, two, K, n.T) {
  tmp.a = c(exp(xamat %*% alpha))  # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  tmp1.a = (1 + tmp.a)^2
  tmp.b = c(exp(xbmat %*% beta))
  sigma = 1/phi - 1
  psi1 = trigamma(mu * sigma)
  psi2 = trigamma((1 - mu) * sigma)
  # auxiliary transformations
  a = psi1 + psi2
  b = psi1 * mu^2 + psi2 * (1 - mu)^2 - trigamma(sigma)
  b.tmp = 1/tmp.b
  a0.tmp = 1/(1 + tmp.a)
  a1.tmp = tmp.a/(1 + tmp.a)
  a2.tmp = tmp.a/(1 + tmp.a)^2
  a3.tmp = tmp.a/(1 + tmp.a)^3
  a.a2 = b.tmp * a3.tmp * (1 - tmp.a)
  a.ab = -b.tmp * a2.tmp
  a.b2 = b.tmp * a1.tmp
  b.a2 = -b.tmp * a3.tmp * (1 - tmp.a)
  b.ab = b.tmp * a2.tmp
  b.b2 = b.tmp * a0.tmp
  waa = -(tmp.a/tmp1.a)^2 * sigma^2 * a * (1 - pZ) + one * a.a2 + two * b.a2
  xamat_expand = do.call(rbind, replicate(K, xamat, simplify = FALSE))
  Kaa = t(.BlockDiag(xamat_expand, K, n.T)) %*% diag(waa, length(waa)) %*% (.BlockDiag(xamat_expand, K, n.T))
  wab = -tmp.a/tmp1.a * (-1/tmp.b) * sigma * (mu * a - psi2) * (1 - pZ) + one * a.ab + two * b.ab
  Kab = t(xamat) %*% diag(wab, length(wab)) %*% xbmat
  wbb = -(1/tmp.b)^2 * b * (1 - pZ) + one * a.b2 + two * b.b2
  Kbb = t(xbmat) %*% diag(wbb, length(wbb)) %*% xbmat
  return(list(Kab = Kab, Kbb = Kbb))
}


.BlockDiagBind <- function(...) {
  d = list(...)
  nrows = sum(sapply(d, nrow))
  ncols = sum(sapply(d, ncol))
  ans = matrix(0, nrows, ncols)
  i1 = 1
  j1 = 1
  for (m in d) {
    i2 = i1 + nrow(m) - 1
    j2 = j1 + ncol(m) - 1
    ans[i1:i2, j1:j2] = m
    i1 = i2 + 1
    j1 = j2 + 1
  }
  return(ans)
}

.BlockDiagBind_list <- function(d){
  nrows = sum(sapply(d, nrow))
  ncols = sum(sapply(d, ncol))
  ans = matrix(0, nrows, ncols)
  i1 = 1
  j1 = 1
  for (m in d){
    i2 = i1 + nrow(m) - 1
    j2 = j1 + ncol(m) - 1
    ans[i1:i2, j1:j2] = m
    i1 = i2 + 1
    j1 = j2 + 1
  }
  return(ans)
}

# .RowbyRow: define a calculation between matrix A and vector b
# (the length of b is the same as the number of rows in A)
# For any row i, the result matrix has the same size as A containing A[i,]*b[i]
.RowbyRow <- function(A, b) {
  temp <- A
  for (i in 1:length(A[, 1])) temp[i, ] <- A[i, ] * b[i]
  return(temp)
}

# 1. Log-Likelihood for Beta-Binomial (H0)
.loglik_bb <- function(params, Y, N) {
  # params: log(a), log(b)
  a <- exp(params[1])
  b <- exp(params[2])
  
  if (is.na(a) || is.na(b) || is.infinite(a) || is.infinite(b)) return(1e20) # Return high cost for nlm minimization
  
  ll <- lchoose(N, Y) + lbeta(Y + a, N - Y + b) - lbeta(a, b)
  sum_ll <- sum(ll)
  
  if (!is.finite(sum_ll)) return(1e20)
  return(-sum_ll) # Return negative LogLik for minimization
}

# 2. Log-Likelihood for Zero-Inflated Beta-Binomial (H1)
.loglik_zibb <- function(params, Y, N) {
  # params: probit(pi), log(a), log(b)
  pi <- pnorm(params[1])
  a <- exp(params[2])
  b <- exp(params[3])
  
  if (is.na(a) || is.na(b) || is.infinite(a) || is.infinite(b)) return(1e20)
  
  # Log-likelihood for BB part
  ll_bb <- lchoose(N, Y) + lbeta(Y + a, N - Y + b) - lbeta(a, b)
  
  idx_pos <- Y > 0
  # Non-zero observations: log((1-pi) * P_bb(y)) = log(1-pi) + LL_bb
  val1 <- sum(ll_bb[idx_pos] + log(1 - pi))
  
  # Zero observations: log(pi + (1-pi) * P_bb(0))
  # Use log-sum-exp trick or direct calculation if stable
  term_zero_bb <- ll_bb[!idx_pos] # This is log(P_bb(0))
  # val2 = sum(log(pi + (1 - pi) * exp(term_zero_bb)))
  
  # More numerically stable approach for zero part:
  # log(pi + (1-pi)*exp(L)) = log(pi) + log(1 + exp(log(1-pi) + L - log(pi)))
  p_bb_zero <- exp(term_zero_bb)
  val2 <- sum(log(pi + (1 - pi) * p_bb_zero))
  
  total_ll <- val1 + val2
  
  if (!is.finite(total_ll)) return(1e20)
  return(-total_ll) # Return negative LogLik
}

# Helper to fit BB and return negative loglik and params
.BB.log.lik <- function(Y) {
  N = rowSums(Y)
  Y = Y[, 1]
  safe_ratio <- Y / (N + 1)
  mu_mom <- mean(safe_ratio)
  v_mom <- var(safe_ratio)
  if (!is.na(v_mom) && v_mom > 0 && v_mom < mu_mom * (1 - mu_mom)) {
    factor <- (mu_mom * (1 - mu_mom) / v_mom) - 1
    factor <- max(factor, 0.01)
    a_init <- mu_mom * factor
    b_init <- (1 - mu_mom) * factor
  } else {
    a_init <- 1; b_init <- 1
  }
  start_params <- c(log(a_init), log(b_init))
  fit = nlm(f = .loglik_bb, p = start_params, Y = Y, N = N, iterlim = 1000)
  ans = -fit$minimum
  a.mat = matrix(exp(fit$estimate)[1], nrow = length(N), ncol = 1)
  b.mat = matrix(exp(fit$estimate)[2], nrow = length(N), ncol = 1)
  return(list(ans = ans, a.mat = a.mat, b.mat = b.mat))
}

# Helper to fit ZIBB and return negative loglik
.BBnZIBB.log.lik <- function(Y) {
  
  # fit bb
  fit_bb = .BB.log.lik(Y)
  bb.ans = fit_bb$ans
  
  # Initialize for zibb
  bb_params = log(c(fit_bb$a.mat[1], fit_bb$b.mat[1]))
  pi_start <- mean(Y == 0) / 2
  logit_pi_init <- qnorm(max(0.01, min(0.99, pi_start)))
  
  # fit zibb
  fit_zibb <- nlm(f = .loglik_zibb,
                  p = c(logit_pi_init, bb_params),
                  Y = Y[, 1], N = rowSums(Y), iterlim = 1000)
  
  zibb.ans = -fit_zibb$minimum
  return(list(bb.ans = bb.ans, zibb.ans = zibb.ans))
}

.check_zeroinflation <- function(ID, Y, n.boot = 9999, n.cores, R.sel = R.sel, cl = NULL) {
  Y.seqDepth = rowSums(Y)
  obs.ans = .BBnZIBB.log.lik(Y)
  stat = -2 * (obs.ans$bb.ans - obs.ans$zibb.ans)
  
  if (is.infinite(stat)) {
    cat("Use Pseudo-ECDF approximation p-value\n")
    return(1/(n.boot + 1))
  }
  
  id.unique = unique(ID)
  indexList = split(1:nrow(Y), ID)
  
  n.boot = ceiling(n.boot/n.cores) * n.cores
  m = 0
  Nexc = 0
  
  boot_worker <- function(iter, current_m) {
    set.seed(12345 + current_m + iter)
    
    ID.tmp = sort(sample(id.unique, replace = TRUE))
    Y.b.list = lapply(ID.tmp, function(id) Y[indexList[[as.character(id)]], , drop = FALSE])
    Y.b.tmp = do.call(rbind, Y.b.list)
    
    bb.ll = .BB.log.lik(Y.b.tmp)
    a.mat = matrix(colMeans(bb.ll$a.mat), nrow(Y), 1)
    b.mat = matrix(colMeans(bb.ll$b.mat), nrow(Y), 1)
    
    Y.b = .simData.GDM(a.mat, b.mat, Y.seqDepth)
    boot.obs = .BBnZIBB.log.lik(Y.b)
    
    return(-2 * (boot.obs$bb.ans - boot.obs$zibb.ans))
  }
  
  if (!is.null(cl)) {
    vars_to_export <- c("id.unique", "Y", "indexList", "Y.seqDepth")
    parallel::clusterExport(cl, varlist = vars_to_export, envir = environment())
  }
  
  if (!is.null(cl)) {
    stat.boot = parallel::parSapply(cl, 1:n.boot, boot_worker, current_m = 0)
  } else {
    stat.boot = sapply(1:n.boot, boot_worker, current_m = 0)
  }
  
  Nexc = sum(stat.boot >= stat)
  if (Nexc <= 10) {
    pval = tryCatch(.gpd_approx(stat.boot, 250, stat), error = function(err) NA)
    if (is.na(pval)) {
      pval = (Nexc + 1)/(n.boot + 1)
    }
  } else {
    pval = Nexc/n.boot
  }
  return(pval)
}

.GenerateScoreFun <- function(par, data) {
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  tmp.a = as.numeric(exp(data$Xa %*% alpha))
  mu.tmp = as.numeric(tmp.a/(1 + tmp.a))
  tmp.b = as.numeric(exp(data$Xb %*% beta))
  phi.tmp = as.numeric(tmp.b/(1 + tmp.b))
  a = (1/phi.tmp - 1) * mu.tmp
  b = (1/phi.tmp - 1) * (1 - mu.tmp)
  a[a < 0] = 0
  b[b < 0] = 0
  # https://github.com/cran/betareg/blob/master/R/betareg.R
  zstar = data$A - data$B
  zdagger = data$B
  mustar = digamma(a) - digamma(b)
  mudagger = digamma(b) - digamma(a + b)
  return((1 - data$Del) * (-1/tmp.b) * (mu.tmp * (zstar - mustar) + (zdagger - mudagger)))
}

.est_para <- function(ID, W, Xa, Xb, Y, zi.check, period, corstr, maxlag) {
  # ID = ID.tmp; W = Xa = Xb = X.r.tmp;
  # Y = cbind(Y.tmp[, 1], rowSums(Y.tmp[, (1 + 1):(K + 1), drop = FALSE]))
  # zi.check = F; period = period.tmp
  da = ncol(Xa)
  db = ncol(Xb)
  dw = ncol(W)
  K = ncol(Y) - 1
  id_index_lst = split(seq_along(ID), ID)   # pre-built lookup list
  id.unique = sort(unique(ID))
  N = length(id.unique)
  n.reps = as.vector(table(ID))
  nt = nrow(Y)  # the total number of observed units
  N1 = sum(n.reps * (n.reps - 1)/2)  # total #within-subject pairs (exch normalizer)
  
  r.step = 0.05  # set up the tuning parameter for updating alpha
  max = 1000  # maximal iterations
  if (zi.check) {
    gamma0 = matrix(0.001, nrow = dw, ncol = K)
  } else {
    gamma0 = matrix(-Inf, nrow = dw, ncol = K)
  }
  alpha0 = matrix(0.001, nrow = da, ncol = K)
  beta0 = matrix(0.001, nrow = db, ncol = K)
  # correlation of same subject, all taxa, at different times
  .np <- .n_corr_params(corstr, maxlag)
  rho0 = sig0 = if (.np > 0) rep(0.1, .np) else 0
  aigdm.reg = .AIGDM_EM(Y = Y, W = W, Xa = Xa, Xb = Xb, gamma0 = gamma0, alpha0 = alpha0, beta0 = beta0)
  gamma0 = aigdm.reg$gamma.est
  alpha0 = aigdm.reg$alpha.est
  beta0 = aigdm.reg$beta.est
  ng = no = dw
  na = da
  for (l in 1:max) {
    all.ddg = matrix(0, ng, ng)  # derivative of GEE
    all.geeg = matrix(0, ng, 1)  # GEE
    all.sig1 = all.sig2 = if (.np > 0) rep(0, .np) else 0
    all.dda = matrix(0, na, na)
    all.geea = matrix(0, na, 1)
    all.rho1 = all.rho2 = if (.np > 0) rep(0, .np) else 0  # 1:numerator 2:denominator
    N2 = ntot = 0  # exch normalizers (accumulated over subjects)
    Del.R.lst = A.R.lst = B.R.lst = A2.R.lst = B2.R.lst = AB.R.lst = vector("list", N)
    for (i in 1:N) {
      Del.R.lst[[i]] = A.R.lst[[i]] = B.R.lst[[i]] = A2.R.lst[[i]] = B2.R.lst[[i]] = AB.R.lst[[i]] = matrix(NA, n.reps[i],
                                                                                                            1)
    }
    for (i in 1:N) {
      ID.sel = id_index_lst[[as.character(id.unique[i])]]
      w.tmp = W[ID.sel, , drop = FALSE]
      xa.tmp = Xa[ID.sel, , drop = FALSE]
      xb.tmp = Xb[ID.sel, , drop = FALSE]
      y.tmp = Y[ID.sel, , drop = FALSE]
      periods_i = period[ID.sel]
      tmp = exp(w.tmp %*% gamma0)
      tmp[is.na(tmp)] = 0
      pZ0.mat = tmp/(1 + tmp)
      pZ0.mat[is.infinite(tmp) & tmp > 0] = 1
      pZ0 = c(pZ0.mat)
      tmp = exp(xa.tmp %*% alpha0)
      mu0.mat = tmp/(1 + tmp)
      tmp = exp(xb.tmp %*% beta0)
      phi0.mat = tmp/(1 + tmp)
      av = mu0.mat * (1/phi0.mat - 1)
      bv = (1 - mu0.mat) * (1/phi0.mat - 1)
      n.T = n.reps[i]
      eDelZ0.mat = eZ0.mat = matrix(NA, n.T, 1)
      eZStar.mat = matrix(NA, n.T, 1)
      for (t in 1:n.T) {
        par.post = .ZIeZ_vec(pZ0.mat[t, , drop = FALSE], av[t, , drop = FALSE], bv[t, , drop = FALSE], y.tmp[t, ])
        eDelZ0.mat[t, ] = par.post$pv.post
        tmp = par.post$av.post + par.post$bv.post
        # post.expectation of logZ conditional on delta==0
        A.R.lst[[i]][t, ] = digamma(par.post$av.post) - digamma(tmp)
        # post.expectation of log(1-Z) conditional on delta==0
        B.R.lst[[i]][t, ] = digamma(par.post$bv.post) - digamma(tmp)
        eZ0.mat[t, ] = par.post$av.post/tmp
        # post expectation of log(Z/(1-Z)) conditional on delta==0
        eZStar.mat[t, ] = digamma(par.post$av.post) - digamma(par.post$bv.post)
      }
      eDelZ0 = c(eDelZ0.mat)
      Del.R.lst[[i]] = eDelZ0.mat
      if (corstr == "exchangeable" && n.T != 1) {
        ntot = ntot + colSums((1 - eDelZ0.mat)^2)
        N2 = N2 + sum(.getSubDiag((1 - eDelZ0.mat) %*% t(1 - eDelZ0.mat)))
      }
      
      # Vdel
      cor.DelZ = .build_corr(corstr, sig0, periods_i)
      
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      VdelZ_inv = ginv(VdelZ)
      tmp.g = t(.BlockDiag(.dpg(w.tmp, gamma0), 1, n.T)) %*% VdelZ_inv
      geeg = as.matrix(apply(.RowbyRow(t(tmp.g), eDelZ0 - pZ0), 2, sum))
      temp.g1 = tmp.g %*% .BlockDiag(.dpg(w.tmp, gamma0), 1, n.T) + geeg %*% t(geeg)
      all.ddg = all.ddg + temp.g1  ## sum up to N for derivative of gee for gamma
      all.geeg = all.geeg + geeg  ## sum up to N for gee for gamma
      # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
      eZ0 = c(eZ0.mat)
      mu0 = c(mu0.mat)
      phi0 = c(phi0.mat)
      eZStar = c(eZStar.mat)
      muStar = digamma(av) - digamma(bv)
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = .build_corr(corstr, rho0, periods_i)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      Vz_inv = ginv(Vz)
      tmp.a = t(.BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)) %*% Vz_inv %*% diag(c(1 - eDelZ0), length(eDelZ0))
      geea = as.matrix(apply(.RowbyRow(t(tmp.a), eZStar - muStar), 2, sum))
      temp.a1 = tmp.a %*% .BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T) + geea %*% t(geea)
      all.dda = all.dda + temp.a1  # sum up to N for derivative of gee for alpha
      all.geea = all.geea + geea  # sum up to N for gee for alpha
      if (n.T != 1) {
        if (corstr == "exchangeable") {
          all.sig1 = all.sig1 + .corr_zero1(eDelZ0, pZ0, n.T)
          all.sig2 = all.sig2 + .corr_zero2(eDelZ0, pZ0, n.T)
          all.rho1 = all.rho1 + .corr_rho1(eDelZ0, eZStar, muStar, mu0, phi0, n.T)
          all.rho2 = all.rho2 + .corr_rho2(eDelZ0, eZStar, muStar, mu0, phi0, n.T)
        } else if (corstr == "toeplitz") {
          all.sig1 = all.sig1 + .corr_toep1(eDelZ0, pZ0, periods_i, maxlag)
          all.sig2 = all.sig2 + .corr_toep2(eDelZ0, pZ0, periods_i, maxlag)
          all.rho1 = all.rho1 + .corr_rho_toep1(eDelZ0, eZStar, muStar, mu0, phi0, periods_i, maxlag)
          all.rho2 = all.rho2 + .corr_rho_toep2(eDelZ0, eZStar, muStar, mu0, phi0, periods_i, maxlag)
        } else if (corstr == "ar1") {
          all.sig1 = all.sig1 + .ar1_zero1(eDelZ0, pZ0, periods_i)
          all.sig2 = all.sig2 + .ar1_zero2(eDelZ0, pZ0, periods_i)
          all.rho1 = all.rho1 + .ar1_rho1(eDelZ0, eZStar, muStar, mu0, phi0, periods_i)
          all.rho2 = all.rho2 + .ar1_rho2(eDelZ0, eZStar, muStar, mu0, phi0, periods_i)
        }
        # independence: no correlation parameters to update
      }
    }
    ########## updating scheme #############
    if (l == 1) {
      all.geeg1 <- all.geeg
    } else {
      all.geeg1 <- ((l - 1) * all.geeg1 + all.geeg)/l
    }
    if (l == 1) {
      all.ddg1 <- all.ddg
    } else {
      all.ddg1 <- ((l - 1) * all.ddg1 + all.ddg)/l
    }
    gnew0 = c(gamma0) + r.step * ginv(all.ddg1) %*% all.geeg
    gnew = (l * c(gamma0) + gnew0)/(l + 1)
    if (l == 1) {
      all.geea1 <- all.geea
    } else {
      all.geea1 <- ((l - 1) * all.geea1 + all.geea)/l
    }
    if (l == 1) {
      all.dda1 <- all.dda
    } else {
      all.dda1 <- ((l - 1) * all.dda1 + all.dda)/l
    }
    anew0 = c(alpha0) + r.step * ginv(all.dda1) %*% all.geea  # matrix to vector
    anew = (l * c(alpha0) + anew0)/(l + 1)
    Del.R = do.call(c, Del.R.lst)
    A.R = do.call(c, A.R.lst)
    B.R = do.call(c, B.R.lst)
    tmp = .AIBetaOptim(Del.R, A.R, B.R, Xa, Xb, anew, c(beta0))
    bnew = tmp[-(1:da)]
    # ---- structure-specific correlation-parameter finalization ----
    if (corstr == "exchangeable") {
      sig.new = (all.sig1/N1)/(all.sig2/nt)
      nan.id = which(is.nan(sig.new) | is.infinite(sig.new) | sig.new == 1)
      if (length(nan.id) != 0) sig.new[nan.id] = 0
      rho.new = (all.rho1/N2)/(all.rho2/ntot)
      nan.id = which(is.nan(rho.new) | is.infinite(rho.new))
      if (length(nan.id) != 0) rho.new[nan.id] = 0
    } else if (corstr == "toeplitz" || corstr == "ar1") {
      sig.new = ifelse(all.sig2 > 0, all.sig1 / all.sig2, 0)
      sig.new = pmax(pmin(sig.new, 0.99), -0.99)   # keep correlation valid
      rho.new = ifelse(all.rho2 > 0, all.rho1 / all.rho2, 0)
      rho.new = pmax(pmin(rho.new, 0.99), -0.99)
    } else {  # independence
      sig.new = 0
      rho.new = 0
    }
    diffs = abs(c(gnew - gamma0, anew - alpha0, bnew - beta0, sig.new - sig0, rho.new - rho0))
    if (max(diffs[!is.infinite(diffs)], na.rm = TRUE) < 1e-04 | l == max)
      break  # max or sum?
    gamma0 = matrix(gnew, ncol = K)
    sig0 = sig.new
    alpha0 = matrix(anew, ncol = K)
    beta0 = matrix(bnew, ncol = K)
    rho0 = rho.new
  }
  gamma0 = matrix(gnew, ncol = K)
  sig0 = sig.new
  alpha0 = matrix(anew, ncol = K)
  beta0 = matrix(bnew, ncol = K)
  rho0 = rho.new
  Del.R.lst = eDelZ0.lst.store = pZ0.lst.store = A.R.lst = B.R.lst = A2.R.lst = B2.R.lst = AB.R.lst = eZ0.lst.store = mu0.lst.store = phi0.lst.store = eZStar.lst.store = muStar.lst.store = pseudoll.score.beta.lst.store = vector("list",
                                                                                                                                                                                                                                    N)
  for (i in 1:N) {
    Del.R.lst[[i]] = A.R.lst[[i]] = B.R.lst[[i]] = A2.R.lst[[i]] = B2.R.lst[[i]] = AB.R.lst[[i]] = matrix(NA, n.reps[i], 1)
    eDelZ0.lst.store[[i]] = eZ0.lst.store[[i]] = mu0.lst.store[[i]] = phi0.lst.store[[i]] = eZStar.lst.store[[i]] = muStar.lst.store[[i]] = numeric(n.reps[i])
    pseudoll.score.beta.lst.store[[i]] = numeric(nrow(beta0))
  }
  for (i in 1:N) {
    ID.sel = which(ID == id.unique[i])
    w.tmp = W[ID.sel, , drop = FALSE]
    xa.tmp = Xa[ID.sel, , drop = FALSE]
    xb.tmp = Xb[ID.sel, , drop = FALSE]
    y.tmp = Y[ID.sel, , drop = FALSE]
    tmp = exp(w.tmp %*% gamma0)
    tmp[is.na(tmp) | is.infinite(tmp)] = 0
    pZ0.mat = tmp/(1 + tmp)
    pZ0.mat[is.infinite(tmp) & tmp > 0] = 1
    pZ0 = c(pZ0.mat)
    pZ0.lst.store[[i]] = pZ0
    tmp = exp(xa.tmp %*% alpha0)
    mu0.mat = tmp/(1 + tmp)
    tmp = exp(xb.tmp %*% beta0)
    phi0.mat = tmp/(1 + tmp)
    av = mu0.mat * (1/phi0.mat - 1)
    bv = (1 - mu0.mat) * (1/phi0.mat - 1)
    muStar.mat = digamma(av) - digamma(bv)
    n.T = n.reps[i]
    eDelZ0.mat = eZ0.mat = eZStar.mat = matrix(NA, n.T, 1)
    for (t in 1:n.T) {
      par.post = .ZIeZ_vec(pZ0.mat[t, , drop = FALSE], av[t, , drop = FALSE], bv[t, , drop = FALSE], y.tmp[t, ])
      eDelZ0.mat[t, ] = par.post$pv.post
      tmp = par.post$av.post + par.post$bv.post
      A.R.lst[[i]][t, ] = digamma(par.post$av.post) - digamma(tmp)
      B.R.lst[[i]][t, ] = digamma(par.post$bv.post) - digamma(tmp)
      # post.expectation of logZ*logZ conditional on delta==0
      A2.R.lst[[i]][t, ] = trigamma(par.post$av.post) - trigamma(tmp) + (digamma(par.post$av.post) - digamma(tmp))^2
      # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0
      B2.R.lst[[i]][t, ] = trigamma(par.post$bv.post) - trigamma(tmp) + (digamma(par.post$bv.post) - digamma(tmp))^2
      # post.expectation of logZ*log(1-Z) conditional on delta==0
      AB.R.lst[[i]][t, ] = -trigamma(tmp) + (digamma(par.post$av.post) - digamma(tmp)) * (digamma(par.post$bv.post) - digamma(tmp))
      eZ0.mat[t, ] = par.post$av.post/tmp
      eZStar.mat[t, ] = digamma(par.post$av.post) - digamma(par.post$bv.post)
    }
    Del.R.lst[[i]] = eDelZ0.mat
    eDelZ0.lst.store[[i]] = c(eDelZ0.mat)
    beta.par = c(alpha0, beta0)
    beta.data = list(Del = c(Del.R.lst[[i]]), A = c(A.R.lst[[i]]), B = c(B.R.lst[[i]]), Xa = xa.tmp, Xb = xb.tmp)
    pseudoll.score.beta.lst.store[[i]] = .GenerateScoreFun(beta.par, beta.data)

    eZ0 = c(eZ0.mat)
    mu0 = c(mu0.mat)
    phi0 = c(phi0.mat)
    eZStar = c(eZStar.mat)
    muStar = c(muStar.mat)
    eZStar.lst.store[[i]] = eZStar
    muStar.lst.store[[i]] = muStar
    eZ0.lst.store[[i]] = eZ0
    mu0.lst.store[[i]] = mu0
    phi0.lst.store[[i]] = phi0
  }
  # print(paste0('Converge?', ifelse(l == max, 0, 1)))
  return(list(CONV = ifelse(l == max, 0, 1), W = W, Xa = Xa, Xb = Xb, gamma0 = matrix(gnew, ncol = K),
              alpha0 = matrix(anew, ncol = K), beta0 = matrix(bnew, ncol = K), sig0 = sig.new, rho0 = rho.new,  period = period,
              corstr = corstr, maxlag = maxlag, eDelZ0.lst = eDelZ0.lst.store, A.R.lst = A.R.lst,
              B.R.lst = B.R.lst, A2.R.lst = A2.R.lst, B2.R.lst = B2.R.lst, AB.R.lst = AB.R.lst, pZ0.lst = pZ0.lst.store, eZ0.lst = eZ0.lst.store,
              mu0.lst = mu0.lst.store, phi0.lst = phi0.lst.store, eZStar.lst = eZStar.lst.store, muStar.lst = muStar.lst.store,
              pseudoll.score.beta.lst = pseudoll.score.beta.lst.store))
}

.score_test_zero <- function(ID, X, X.index, Y, est.para.lst, zi.id, period) {
  K = ncol(Y) - 1
  K.zi = length(zi.id)
  stat = df = numeric(K.zi)
  for (j.zi in 1:K.zi) {
    j = zi.id[j.zi]
    Y.tmp = cbind(Y[, j], rowSums(Y[, (j + 1):(K + 1), drop = FALSE]))
    est.para = est.para.lst[[j]]
    ID.tmp = ID
    W.tmp = X
    W.r.tmp = est.para$W
    Xa.r.tmp = est.para$Xa
    Xb.r.tmp = est.para$Xb
    dw = ncol(W.tmp)
    dw.r = ncol(W.r.tmp)
    da.r = ncol(Xa.r.tmp)
    db.r = ncol(Xb.r.tmp)
    id.unique = sort(unique(ID.tmp))
    N = length(id.unique)
    n.reps = as.vector(table(ID.tmp))
    ng = dw
    para.ind = sort(sapply(X.index, function(x) seq(x, ng, ng)))
    gamma0 = est.para$gamma0
    alpha0 = est.para$alpha0
    beta0 = est.para$beta0
    sig0 = est.para$sig0
    rho0 = est.para$rho0
    corstr = est.para$corstr
    eDelZ0.lst = est.para$eDelZ0.lst
    A.R.lst = est.para$A.R.lst
    B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst
    B2.R.lst = est.para$B2.R.lst
    AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst
    eZ0.lst = est.para$eZ0.lst
    mu0.lst = est.para$mu0.lst
    phi0.lst = est.para$phi0.lst
    eZStar.lst = est.para$eZStar.ls
    muStar.lst = est.para$muStar.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, ng, ng)
    ind.g = 1:ng
    for (i in 1:N) {
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel, , drop = FALSE]
      w.tmp = W.tmp[ID.sel, , drop = FALSE]
      n.T = n.reps[i]
      periods_i = period[ID.sel]
      
      pZ0 = pZ0.lst[[i]]  # unlist(pZ0.lst[i,])
      eDelZ0 = eDelZ0.lst[[i]]
      cor.DelZ = .build_corr(corstr, sig0, periods_i)
      
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      VdelZ_inv = ginv(VdelZ)
      gamma0.w = matrix(0, dw, 1)
      gamma0.w[1:dw.r, ] = gamma0
      tmp.g.w = t(.BlockDiag(.dpg(w.tmp, gamma0.w), 1, n.T)) %*% VdelZ_inv
      dvd.g.w = t(.BlockDiag(.dpg(w.tmp, gamma0.w), 1, n.T)) %*% VdelZ_inv %*% .BlockDiag(.dpg(w.r.tmp, gamma0),
                                                                                            1, n.T)
      geeg.w = as.matrix(apply(.RowbyRow(t(tmp.g.w), eDelZ0 - pZ0), 2, sum))
      ES = geeg.w
      U = U + ES
      I = I + (-dvd.g.w)
      M = M + geeg.w %*% t(geeg.w) * (diag(1) %x% matrix(1, ng, ng))
    }
    U.r = U
    I.r = I
    M.r = M
    H.r = .colwise.cbind(-I.r[para.ind, -para.ind] %*% ginv(I.r[-para.ind, -para.ind]), diag(length(para.ind)), para.ind)
    stat[j.zi] = c(t(U.r[para.ind, , drop = F]) %*% ginv(H.r %*% M.r %*% t(H.r)) %*% U.r[para.ind, , drop = F])
    df[j.zi] = length(para.ind)
  }
  stat.sum = sum(stat)
  df.sum = sum(df)
  pval = pchisq(stat.sum, df.sum, lower.tail = F)
  return(list(stat.sum = stat.sum, pval.sum = pval, stat = stat))
}

.score_test_mean <- function(ID, X, X.index, Y, est.para.lst, period) {
  K = ncol(Y) - 1
  stat = df = numeric(K)
  # ind.tmp = NULL
  for (j in 1:K) {
    # j = 1
    Y.tmp = cbind(Y[, j], rowSums(Y[, (j + 1):(K + 1), drop = FALSE]))
    est.para = est.para.lst[[j]]
    ID.tmp = ID
    W.r.tmp = est.para$W
    Xa.tmp = X
    Xa.r.tmp = est.para$Xa
    Xb.r.tmp = est.para$Xb
    dw.r = ncol(W.r.tmp)
    db.r = ncol(Xb.r.tmp)
    da = ncol(Xa.tmp)
    da.r = ncol(Xa.r.tmp)
    id.unique = sort(unique(ID.tmp))
    N = length(id.unique)
    n.reps = as.vector(table(ID.tmp))
    na = da
    nb = db.r
    para.ind = sort(sapply(X.index, function(x) seq(x, na, na)))
    gamma0 = est.para$gamma0
    alpha0 = est.para$alpha0
    beta0 = est.para$beta0
    sig0 = est.para$sig0
    rho0 = est.para$rho0
    corstr = est.para$corstr
    eDelZ0.lst = est.para$eDelZ0.lst
    pZ0.lst = est.para$pZ0.lst
    A.R.lst = est.para$A.R.lst
    B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst
    B2.R.lst = est.para$B2.R.lst
    AB.R.lst = est.para$AB.R.lst
    eZ0.lst = est.para$eZ0.lst
    mu0.lst = est.para$mu0.lst
    phi0.lst = est.para$phi0.lst
    eZStar.lst = est.para$eZStar.ls
    muStar.lst = est.para$muStar.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, na + nb, na + nb)
    ind.a = 1:na
    ind.b = (na + 1):(na + nb)
    for (i in 1:N) {
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel, , drop = FALSE]
      xa.r.tmp = Xa.r.tmp[ID.sel, , drop = FALSE]
      xa.tmp = Xa.tmp[ID.sel, , drop = FALSE]
      xb.r.tmp = Xb.r.tmp[ID.sel, , drop = FALSE]
      n.T = n.reps[i]
      periods_i = period[ID.sel]
      
      pZ0 = pZ0.lst[[i]]  # unlist(pZ0.lst[i,])
      mu0 = mu0.lst[[i]]
      phi0 = phi0.lst[[i]]
      eZ0 = eZ0.lst[[i]]
      eDelZ0 = eDelZ0.lst[[i]]
      eZStar = eZStar.lst[[i]]
      muStar = muStar.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z    = .build_corr(corstr, rho0, periods_i)
      
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      Vz_inv = ginv(Vz)
      alpha0.x = matrix(0, da, 1)
      alpha0.x[1:da.r, ] = alpha0
      tmp.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% Vz_inv %*% diag(c(1 - eDelZ0), length(eDelZ0))
      dvd.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% Vz_inv%*% diag(c(1 - eDelZ0), length(eDelZ0)) %*%
        .BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZStar - muStar), 2, sum))
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.r.tmp)
      # dvd.b.x = t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.r.tmp)) %*%
      # t(t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.r.tmp)))
      geeb.x = as.matrix(pseudoll.score.beta)
      dvd.ab.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% Vz_inv %*% diag(c(1 - eDelZ0), length(eDelZ0)) %*%
        .BlockDiag(.dmb(xb.r.tmp, beta0, mu0, phi0), 1, n.T)
      ES = rbind(geea.x, geeb.x)
      U = U + ES
      a = (1/phi0 - 1) * mu0
      a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0)
      b[b < 0] = 0
      A = digamma(a) - digamma(a + b)
      B = digamma(b) - digamma(a + b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1 - eDelZ0) * (-A + A.post)
      two.R = (1 - eDelZ0) * (-B + B.post)
      tmp.dvd = .DVD(xa.tmp, xb.r.tmp, alpha0.x, beta0, pZ0, mu0, phi0, one.R, two.R, 1, n.T)
      dvd.ab.all = rbind(cbind(dvd.a.x, dvd.ab.x), cbind(t(tmp.dvd$Kab), tmp.dvd$Kbb))
      I = I + (-dvd.ab.all)
      M[ind.a, ind.a] = M[ind.a, ind.a] + geea.x %*% t(geea.x) * (diag(1) %x% matrix(1, na, na))
      M[ind.b, ind.b] = M[ind.b, ind.b] + geeb.x %*% t(geeb.x) * (diag(1) %x% matrix(1, nb, nb))
      M[ind.a, ind.b] = M[ind.a, ind.b] + geea.x %*% t(geeb.x) * (diag(1) %x% matrix(1, na, nb))
      M[ind.b, ind.a] = M[ind.b, ind.a] + geeb.x %*% t(geea.x) * (diag(1) %x% matrix(1, nb, na))
    }
    U.r = U
    I.r = I
    M.r = M
    H.r = .colwise.cbind(-I.r[para.ind, -para.ind] %*% ginv(I.r[-para.ind, -para.ind]), diag(length(para.ind)), para.ind)
    stat[j] = c(t(U.r[para.ind, , drop = F]) %*% ginv(H.r %*% M.r %*% t(H.r)) %*% U.r[para.ind, , drop = F])
    df[j] = length(para.ind)
  }
  stat.sum = sum(stat)
  df.sum = sum(df)
  pval = pchisq(stat.sum, df.sum, lower.tail = F)
  return(list(stat.sum = stat.sum, pval.sum = pval, stat = stat))
}

.score_test_disp <- function(ID, X, X.index, Y, est.para.lst, period) {
  # X = Xb
  K = ncol(Y) - 1
  stat = df = numeric(K)
  # ind.tmp = NULL
  for (j in 1:K) {
    Y.tmp = cbind(Y[, j], rowSums(Y[, (j + 1):(K + 1), drop = FALSE]))
    est.para = est.para.lst[[j]]
    ID.tmp = ID
    W.r.tmp = est.para$W
    Xb.tmp = X
    Xb.r.tmp = est.para$Xb
    Xa.r.tmp = est.para$Xa
    dw.r = ncol(W.r.tmp)
    da.r = ncol(Xa.r.tmp)
    db = ncol(Xb.tmp)
    db.r = ncol(Xb.r.tmp)
    id.unique = sort(unique(ID.tmp))
    N = length(id.unique)
    n.reps = as.vector(table(ID.tmp))
    ng = no = dw.r
    na = da.r
    nb = db
    para.ind = sort(sapply(X.index, function(x) seq(x, nb, nb))) + na
    gamma0 = est.para$gamma0
    alpha0 = est.para$alpha0
    beta0 = est.para$beta0
    sig0 = est.para$sig0
    rho0 = est.para$rho0
    corstr = est.para$corstr
    eDelZ0.lst = est.para$eDelZ0.lst
    A.R.lst = est.para$A.R.lst
    B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst
    B2.R.lst = est.para$B2.R.lst
    AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst
    eZ0.lst = est.para$eZ0.lst
    mu0.lst = est.para$mu0.lst
    phi0.lst = est.para$phi0.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, na + nb, na + nb)
    ind.a = 1:na
    ind.b = (na + 1):(na + nb)
    for (i in 1:N) {
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel, , drop = FALSE]
      xa.r.tmp = Xa.r.tmp[ID.sel, , drop = FALSE]
      xb.tmp = Xb.tmp[ID.sel, , drop = FALSE]
      n.T = n.reps[i]
      periods_i = period[ID.sel]
      pZ0 = pZ0.lst[[i]]  # unlist(pZ0.lst[i,])
      mu0 = mu0.lst[[i]]
      phi0 = phi0.lst[[i]]
      eZ0 = eZ0.lst[[i]]
      eDelZ0 = eDelZ0.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z    = .build_corr(corstr, rho0, periods_i)
      
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      Vz_inv = ginv(Vz)
      tmp.a.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% Vz_inv %*% diag(c(1 - eDelZ0), length(eDelZ0))
      dvd.a.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% Vz_inv %*% diag(c(1 - eDelZ0), length(eDelZ0)) %*%
        .BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZ0 - mu0), 2, sum))
      beta0.x = matrix(0, db, 1)
      beta0.x[1:db.r, ] = beta0
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.tmp)
      geeb.x = as.matrix(pseudoll.score.beta)
      dvd.ab.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% Vz_inv %*% diag(c(1 - eDelZ0), length(eDelZ0)) %*%
        .BlockDiag(.dmb(xb.tmp, beta0.x, mu0, phi0), 1, n.T)
      ES = rbind(geea.x, geeb.x)
      U = U + ES
      a = (1/phi0 - 1) * mu0
      a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0)
      b[b < 0] = 0
      A = digamma(a) - digamma(a + b)
      B = digamma(b) - digamma(a + b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1 - eDelZ0) * (-A + A.post)
      two.R = (1 - eDelZ0) * (-B + B.post)
      tmp.dvd = .DVD(xa.r.tmp, xb.tmp, alpha0, beta0.x, pZ0, mu0, phi0, one.R, two.R, 1, n.T)
      dvd.ab.all = rbind(cbind(dvd.a.x, dvd.ab.x), cbind(t(tmp.dvd$Kab), tmp.dvd$Kbb))
      I = I + (-dvd.ab.all)
      M[ind.a, ind.a] = M[ind.a, ind.a] + geea.x %*% t(geea.x) * (diag(1) %x% matrix(1, na, na))
      M[ind.b, ind.b] = M[ind.b, ind.b] + geeb.x %*% t(geeb.x) * (diag(1) %x% matrix(1, nb, nb))
      M[ind.a, ind.b] = M[ind.a, ind.b] + geea.x %*% t(geeb.x) * (diag(1) %x% matrix(1, na, nb))
      M[ind.b, ind.a] = M[ind.b, ind.a] + geeb.x %*% t(geea.x) * (diag(1) %x% matrix(1, nb, na))
    }
    U.r = U
    I.r = I
    M.r = M
    H.r = .colwise.cbind(-I.r[para.ind, -para.ind] %*% ginv(I.r[-para.ind, -para.ind]), diag(length(para.ind)), para.ind)
    stat[j] = c(t(U.r[para.ind, , drop = F]) %*% ginv(H.r %*% M.r %*% t(H.r)) %*% U.r[para.ind, , drop = F])
    df[j] = length(para.ind)
  }
  stat.sum = sum(stat)
  df.sum = sum(df)
  pval = pchisq(stat.sum, df.sum, lower.tail = F)
  return(list(stat.sum = stat.sum, pval.sum = pval, stat = stat))
}

.score_test_perm <- function(ID, X, X.index, Y, period, est.para.lst, stat.sum, test, Rconf, n.perm, n.cores, R.sel, zi.id = NULL, cl = NULL) {
  
  n.perm = ceiling(n.perm/n.cores) * n.cores
  m = 0
  Nexc = 0
  stat.perm = numeric(n.perm)
  
  id.unique = sort(unique(ID))
  N = length(id.unique)
  n.reps = as.vector(table(ID))
  
  X.int.uniq = matrix(NA, N, length(X.index))
  for (i in 1:N) {
    tmpX = X[which(ID == id.unique[i]), X.index, drop = FALSE]
    if (dim(unique(tmpX))[1] != 1) {
      stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
    }
    X.int.uniq[i, ] = unique(tmpX)
  }
  
  perm_worker <- function(iter, current_m, test_type) {
    set.seed(12345 + current_m + iter)
    
    cluster.p = sample(1:N)
    
    X.p.int = NULL
    for (i in 1:N) {
      n.T = n.reps[i]
      X.p.int = rbind(X.p.int, matrix(rep(X.int.uniq[cluster.p[i], ], n.T), nrow = n.T, byrow = TRUE))
    }
    
    X.curr = X
    X.curr[, X.index] = Rconf %*% X.p.int
    
    tmp.stat = 0
    
    if (test_type == "Zero")
      tmp.stat = .score_test_zero(ID, X.curr, X.index, Y, est.para.lst, zi.id, period)$stat.sum
    
    if (test_type == "Mean")
      tmp.stat = .score_test_mean(ID, X.curr, X.index, Y, est.para.lst, period)$stat.sum
    
    if (test_type == "Disp")
      tmp.stat = .score_test_disp(ID, X.curr, X.index, Y, est.para.lst, period)$stat.sum
    
    if (test_type == "Omni") {
      s.zero = 0
      if (length(zi.id) > 0) {
        s.zero = .score_test_zero(ID, X.curr, X.index, Y, est.para.lst, zi.id, period)$stat.sum
      }
      s.mean = .score_test_mean(ID, X.curr, X.index, Y, est.para.lst, period)$stat.sum
      s.disp = .score_test_disp(ID, X.curr, X.index, Y, est.para.lst, period)$stat.sum
      tmp.stat = sum(s.zero, s.mean, s.disp)
    }
    
    return(tmp.stat)
  }
  
  if (!is.null(cl)) {
    vars_to_export <- c("X", "X.index", "Y", "ID", "est.para.lst",
                        "zi.id", "Rconf", "X.int.uniq", "n.reps", "N")
    parallel::clusterExport(cl, varlist = vars_to_export, envir = environment())
  }
  
  if (!is.null(cl)) {
    stat.perm = parallel::parSapply(cl, 1:n.perm, perm_worker, current_m = 0, test_type = test)
  } else {
    stat.perm = sapply(1:n.perm, perm_worker, current_m = 0, test_type = test)
  }
  
  Nexc = sum(stat.perm >= stat.sum)
  
  if (Nexc <= 10) {
    pval = tryCatch(.gpd_approx(stat.perm, 250, stat.sum), error = function(err) NA)
    if (is.na(pval)) {
      pval = (Nexc + 1)/(n.perm + 1)
    }
  } else {
    pval = (Nexc + 1)/(n.perm + 1)
  }
  
  return(pval)
}




.rowCumprods <- function(x) {
  if (!is.matrix(x))
    stop("Must be matrix")
  n_rows <- nrow(x)
  n_cols <- ncol(x)
  res <- matrix(0, nrow = n_rows, ncol = n_cols)
  for (i in 1:n_rows) {
    current_prod <- 1
    for (j in 1:n_cols) {
      current_prod <- current_prod * x[i, j]
      res[i, j] <- current_prod
    }
  }
  return(res)
}

.simData.GDM <- function(a, b, SeqDepth) {
  N = nrow(a)
  K = ncol(a)
  # Generate Z matrix
  Z = matrix(rbeta(N * K, a, b), N, K)
  # Calculate P matrix
  cumprod_1_minus_Z = .rowCumprods(cbind(1, 1 - Z))
  P = Z * cumprod_1_minus_Z[, -(K + 1), drop = FALSE]
  P = cbind(P, 1 - rowSums(P))
  P[P < 0] = 0
  # Generate Y matrix using apply for rowwise operation
  Y = t(vapply(1:N, function(i) rmultinom(1, SeqDepth[i], P[i, ]), numeric(K + 1)))
  return(Y)
}

# http://motsinger-reif-lab.org/files/documents/choose-sampling-parameters-for-adaptive-permutation-1.r
# determines the number of test statistics that should
# be sampled in adaptive permutation
.choose_r <- function(alpha, c) {
  error <- alpha * c
  R <- 0
  foundR <- FALSE
  while (!foundR) {
    R <- R + 1
    brange <- qnbinom(c(0.1586553, 0.8413447), R, alpha)
    pvalRange <- R/(R + brange)
    diff <- max(abs(pvalRange - alpha))
    if (diff < error) {
      foundR <- TRUE
    }
  }
  return(R)
}
# use Nexc most exterme test statistic to approximate GPD
.gpd_approx <- function(yperm, nexc, obsts) {
  yperm = na.omit(yperm)
  y = sort(yperm, decreasing = TRUE)
  tmp = .gpd_params_est(nexc, y)
  a_hat = tmp[1]
  k_hat = tmp[2]
  t = tmp[3]
  ### goodness-fit test
  nexc_re = .gpd_goft(nexc, y[1:nexc] - t)
  if (nexc_re[2] == nexc) {
    z0 = obsts - t
    p = .gpd_pval(a_hat, k_hat, z0, length(yperm), nexc)
    return(p)
  } else {
    tmp = .gpd_params_est(nexc_re[2], y)
    a_hat = tmp[1]
    k_hat = tmp[2]
    t = tmp[3]
    z0 = obsts - t
    p = .gpd_pval(a_hat, k_hat, z0, length(yperm), nexc_re[2])
    return(p)
  }
}
# Parameter estimators for the generalized Pareto distribution
.gpd_params_est <- function(nexc, y) {
  t = (y[nexc] + y[nexc + 1])/2
  z = y[1:nexc] - t
  z2 = z^2
  m = mean(z)
  m2 = mean(z2)
  a_hat = m * m2/(m2 - m^2)/2
  k_hat = (m^2/(m2 - m^2) - 1)/2
  return(c(a_hat, k_hat, t))
}

# p-value for the generalized Pareto distribution approximation
.gpd_pval <- function(a_hat, k_hat, z0, nperm, nexc) {
  p = nexc/nperm * ((1 - k_hat * z0/a_hat)^(1/k_hat))
  return(p)
}
# iteratively reduce Nexc by 10 until a goodness-of-fit satisfy
.gpd_goft <- function(nexc, y) {
  # y: the sorted test statistic - threshold t
  nexc = length(y)
  p = .gp_test(y)
  nexc.c = seq(0, nexc - 10, 10)
  z = y
  i = 0
  re = c()
  for (i in nexc.c) {
    z = y[1:(nexc - i)]
    p = .gp_test(z)
    re = rbind(re, c(i, p))
    if (!is.na(p) & p > 0.05)
      break
    i = i + 10
  }
  if (nrow(re) >= 2) {
    nexc.c2 = seq(re[(nrow(re) - 1), 1] + 1, re[nrow(re), 1], 1)
    re = c()
    for (i in nexc.c2) {
      z = y[1:(nexc - i)]
      p = .gp_test(z)
      re = rbind(re, c(i, p))
      if (!is.na(p) & p > 0.05)
        break
      i = i + 10
    }
  }
  p = re[nrow(re), 2]
  len = nexc - re[nrow(re), 1]
  return(c(p, len))
}
# Bootstrap goodness-of-fit test for the Generalized Pareto distribution test of fit for the GPD with unknown parameter
# refer to 'goft'
.gp_test <- function(x, B = 999) {
  x = x[!is.na(x)]
  x = as.vector(x)
  n = length(x)  #  sample size without NA values
  samplerange = max(x) - min(x)
  gammap = .amle_method(x, k = ceiling(0.2 * n))[1]
  gamman = .combined_method(x)[1]
  r1 = .R1(x)  # observed value of R^-
  r2 = .R2(x)  # observed value of R^+
  # use 'combined' for fitting GPD with negative shape parameter use 'amle' for fitting GPD with non-negative shape
  # parameter
  p.value1 = sum(replicate(B, .R1(.rgp(n, shape = gamman))) < r1)/B  # bootstrap p-value for H_0^-
  p.value2 = sum(replicate(B, .R2(.rgp(n, shape = gammap))) < r2)/B  # bootstrap p-value for H_0^+
  p.value = max(p.value1, p.value2)  # p-value of the intersection-union test
  return(p.value)
}

# Asymptotic maximum likelihood estimators
.amle_method <- function(x, k) {
  x = sort(x)
  n = length(x)
  nk = n - k
  x1 = x[(nk + 1):n]
  w = log(x1)
  g = -(w[1] - sum(w)/k)
  sigma = g * exp(w[1] + g * log(k/n))
  return(c(g, sigma))
}
# Combined estimators
.combined_method <- function(x) {
  m = mean(x)
  maxi = max(x)
  g = m/(m - maxi)
  sigma = -g * maxi
  return(c(g, sigma))
}
# Test statistic for H_0^-
.R1 <- function(x) {
  gamma_neg = .combined_method(x)[1]
  Fn = ecdf(x)
  x1 = x[x != max(x)]
  z1 = (1 - Fn(x1))^(-gamma_neg)
  return(abs(cor(x1, z1)))
}
# Test statistic for H_0^+
.R2 <- function(x) {
  n = length(x)
  Fn = ecdf(x)
  gamma_positive = .amle_method(x, ceiling(0.2 * n))[1]
  x1 = x[x != max(x)]
  y1 = (1 - Fn(x1))^(-gamma_positive)
  x.star = log(x1)
  y.star = log(y1 - 1)
  if (gamma_positive <= 0.5)
    return(cor(x1, y1))
  if (gamma_positive > 0.5)
    return((cor(x.star, y.star)))
}
# Simulation of random numbers from the gPd
.rgp <- function(n, shape) {
  if (shape != 0)
    return((1/shape) * (runif(n)^(-shape) - 1)) else return(rexp(n, 1))
}

# ACAT function from ACAT package
# https://github.com/yaowuliu/ACAT
# date: 4/16/2024
.ACAT <- function(Pvals, weights = NULL, is.check = TRUE) {
  Pvals <- as.matrix(Pvals)
  if (is.check) {
    #### check if there is NA
    if (sum(is.na(Pvals)) > 0) {
      stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals < 0) + sum(Pvals > 1)) > 0) {
      stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero <- (colSums(Pvals == 0) >= 1)
    is.one <- (colSums(Pvals == 1) >= 1)
    if (sum((is.zero + is.one) == 2) > 0) {
      stop("Cannot have both 0 and 1 p-values in the same column!")
    }
    if (sum(is.zero) > 0) {
      warning("There are p-values that are exactly 0!")
    }
    if (sum(is.one) > 0) {
      warning("There are p-values that are exactly 1!")
    }
  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights)) {
    is.weights.null <- TRUE
  } else {
    is.weights.null <- FALSE
    weights <- as.matrix(weights)
    if (sum(dim(weights) != dim(Pvals)) > 0) {
      stop("The dimensions of weights and Pvals must be the same!")
    } else if (is.check & (sum(weights < 0) > 0)) {
      stop("All the weights must be nonnegative!")
    } else {
      w.sum <- colSums(weights)
      if (sum(w.sum <= 0) > 0) {
        stop("At least one weight should be positive in each column!")
      } else {
        for (j in 1:ncol(weights)) {
          weights[, j] <- weights[, j]/w.sum[j]
        }
      }
    }
  }
  #### check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small <- (Pvals < 1e-15)
  if (is.weights.null) {
    Pvals[!is.small] <- tan((0.5 - Pvals[!is.small]) * pi)
    Pvals[is.small] <- 1/Pvals[is.small]/pi
    cct.stat <- colMeans(Pvals)
  } else {
    Pvals[!is.small] <- weights[!is.small] * tan((0.5 - Pvals[!is.small]) * pi)
    Pvals[is.small] <- (weights[is.small]/Pvals[is.small])/pi
    cct.stat <- colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval <- pcauchy(cct.stat, lower.tail = F)
  return(pval)
}


## ===================================================================
## GIC (QIC/CIC-style) model-selection layer for merged AIGDM
##
## Criterion (per taxon, per candidate working-correlation structure R):
##
##   GIC(R) = -2 * Qhat + 2 * trace( ginv(Omega_I) %*% V_R )
##
## following Pan (2001, Biometrics) QIC and Cui (2007, Stata Journal):
##   Qhat     : the two-part (zero-inflation + beta) INDEPENDENCE
##              quasi-(pseudo)log-likelihood evaluated at the fitted
##              estimates. The working correlation influences the fit's
##              variance, not this goodness-of-fit term, so Qhat is
##              essentially shared across structures for a given taxon.
##   Omega_I  : the model-based (naive) covariance of the GEE estimator
##              computed under the INDEPENDENCE working correlation,
##              i.e. Omega_I = ginv(A_indep) where A_indep is the GEE
##              sensitivity (bread) built with V_i = identity-correlation.
##   V_R      : the robust sandwich covariance under the candidate
##              structure R, V_R = ginv(A_R) %*% B_R %*% ginv(A_R),
##              A_R = sensitivity (bread), B_R = variability (meat).
##
## trace(ginv(Omega_I) %*% V_R) is the QIC/CIC generalized number of
## parameters (effective d.f. penalty). Because A_R is assembled as a
## Gram matrix (sum of D' V^{-1} W D), it is positive semi-definite by
## construction -- this is the fix for the earlier .DVD-based bread that
## produced negative eigenvalues and negative penalties.
##
## GEE ingredients for the MEAN model (this is where corstr enters),
## matching the estimating equations solved inside .est_para:
##   D_i    = .BlockDiag(.dma(Xa_i, alpha0, mu0, phi0), 1, n.T)   (dmu*/dalpha)
##   V_i    = sqrt(.var_ZStar) %*% cor.Z(R) %*% sqrt(.var_ZStar)
##   W_i    = diag(1 - eDelZ0)     (posterior non-zero weight)
##   r_i    = eZStar - muStar      (working residual of mu* GEE)
## so the per-subject score is  U_i = D_i' V_i^{-1} W_i r_i  (== geea),
## the sensitivity contribution is  D_i' V_i^{-1} W_i D_i,
## and the meat contribution is  U_i U_i'.
##
## The mean-model block is the identified/estimated correlation target
## for structure selection (dispersion/zero share the same corstr in
## the merged model). Restricting the criterion to the mean-model GEE
## keeps the penalty well-defined and interpretable.
## ===================================================================

## ---- fitted two-part quasi-log-likelihood (independence) ----------
## Sums, over all observations of one taxon fit, the working
## log-likelihood of the beta component (mean/dispersion) plus the
## zero-inflation Bernoulli component. Uses fitted mu0, phi0, pZ0 and
## posterior eDelZ0 / posterior E[logZ], E[log(1-Z)]. This is the
## independence quasi-loglik consistent with the EM objective; it does
## NOT depend on the working correlation, so it is (essentially) shared
## across structures -- only the penalty differs.
.aigdm_quasi_loglik <- function(est.para) {
  N <- length(est.para$mu0.lst)
  ll <- 0
  for (i in 1:N) {
    mu0  <- est.para$mu0.lst[[i]]
    phi0 <- est.para$phi0.lst[[i]]
    pZ0  <- est.para$pZ0.lst[[i]]
    eDel <- est.para$eDelZ0.lst[[i]]
    A.post <- c(est.para$A.R.lst[[i]])   # E[log Z | delta=0]
    B.post <- c(est.para$B.R.lst[[i]])   # E[log(1-Z) | delta=0]
    sigma <- (1/phi0 - 1)
    a <- pmax((sigma) * mu0, 0)
    b <- pmax((sigma) * (1 - mu0), 0)
    ## beta working loglik contribution for the non-zero part,
    ## E_delta[ (1-delta) * logBeta(Z; a,b) ] using posterior moments:
    ##   (a-1)E[logZ] + (b-1)E[log(1-Z)] - lbeta(a,b)
    beta_ll <- (a - 1) * A.post + (b - 1) * B.post - lbeta(a, b)
    beta_ll[!is.finite(beta_ll)] <- 0
    ## zero-inflation Bernoulli part (guard against pZ0 in {0,1})
    p <- pmin(pmax(pZ0, 1e-12), 1 - 1e-12)
    zi_ll <- eDel * log(p) + (1 - eDel) * log(1 - p)
    zi_ll[!is.finite(zi_ll)] <- 0
    ll <- ll + sum((1 - eDel) * beta_ll) + sum(zi_ll)
  }
  ll
}

## ---- mean-model GEE sensitivity (bread A) and meat (B) -------------
## For a fitted est.para, rebuild per-subject GEE quantities of the
## MEAN model and accumulate:
##   A = sum_i D_i' V_i^{-1} W_i D_i     (sensitivity / bread; PSD)
##   B = sum_i U_i U_i'                  (variability / meat)
##   U = D_i' V_i^{-1} W_i r_i           (score contribution)
## The working correlation used to form V_i is chosen by `corstr.use`
## (defaults to the fitted structure). Passing corstr.use="independence"
## yields the independence-working sensitivity A_indep used for Omega_I.
##
## D_i = .BlockDiag(.dma(...), 1, n.T) exactly matches the alpha GEE
## Jacobian block used inside .est_para (see the EM loop:
##   tmp.a = t(.BlockDiag(.dma(...),1,n.T)) %*% Vz_inv %*% diag(1-eDelZ0)
##   geea  = sum rows of t(tmp.a) * (eZStar - muStar)
##   sensitivity term = tmp.a %*% .BlockDiag(.dma(...),1,n.T)  ).
.aigdm_mean_A_B <- function(est.para, corstr.use = NULL) {
  corstr <- if (is.null(corstr.use)) est.para$corstr else corstr.use
  rho0   <- est.para$rho0
  alpha0 <- est.para$alpha0
  period <- est.para$period
  na <- ncol(alpha0) * nrow(alpha0)   # length of vec(alpha) == d.f. of mean GEE
  ## per-subject bookkeeping from stored list lengths
  N    <- length(est.para$mu0.lst)
  reps <- vapply(est.para$mu0.lst, length, integer(1))
  ends <- cumsum(reps); starts <- c(1L, head(ends, -1) + 1L)
  p_a  <- ncol(est.para$Xa)            # number of alpha covariates (per taxon col)
  ## dimension of the stacked alpha GEE == p_a * K, with K = ncol(alpha0)
  A <- matrix(0, na, na)
  B <- matrix(0, na, na)
  for (i in 1:N) {
    ID.sel    <- starts[i]:ends[i]
    xa.tmp    <- est.para$Xa[ID.sel, , drop = FALSE]
    n.T       <- reps[i]
    periods_i <- period[ID.sel]
    mu0    <- est.para$mu0.lst[[i]]
    phi0   <- est.para$phi0.lst[[i]]
    eDelZ0 <- est.para$eDelZ0.lst[[i]]
    eZStar <- est.para$eZStar.lst[[i]]
    muStar <- est.para$muStar.lst[[i]]
    cor.Z  <- .build_corr(corstr, rho0, periods_i)
    sdz    <- sqrt(.var_ZStar(mu0, phi0))
    Vz     <- sdz %*% cor.Z %*% sdz
    Vz_inv <- ginv(Vz)
    Dblk   <- .BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)
    tmp.a  <- t(Dblk) %*% Vz_inv %*% diag(c(1 - eDelZ0), length(eDelZ0))
    Ui     <- as.matrix(apply(.RowbyRow(t(tmp.a), eZStar - muStar), 2, sum))
    Ai     <- tmp.a %*% Dblk
    A <- A + Ai
    B <- B + Ui %*% t(Ui)
  }
  list(A = A, B = B)
}

## ---- GIC for a single fitted est.para (Pan QIC) --------------------
## Requires an independence-working sensitivity (Omega_I = ginv(A_indep))
## and the structure-R sandwich V_R = ginv(A_R) B_R ginv(A_R).
## A_indep can be supplied (shared across the 4 structures for a taxon)
## or is computed here from est.para with corstr.use="independence".
.gic_aigdm <- function(est.para, A.indep = NULL) {
  qhat <- .aigdm_quasi_loglik(est.para)
  AB <- tryCatch(.aigdm_mean_A_B(est.para), error = function(e) NULL)
  if (is.null(AB)) return(list(gic = NA_real_, qhat = qhat, penalty = NA_real_))
  if (is.null(A.indep)) {
    ABi <- tryCatch(.aigdm_mean_A_B(est.para, corstr.use = "independence"),
                    error = function(e) NULL)
    A.indep <- if (is.null(ABi)) NULL else ABi$A
  }
  penalty <- tryCatch({
    V_R    <- ginv(AB$A) %*% AB$B %*% ginv(AB$A)   # robust sandwich under R
    Omega_I <- ginv(A.indep)                        # naive cov under independence
    sum(diag(ginv(Omega_I) %*% V_R))                # trace(A_indep %*% V_R)
  }, error = function(e) NA_real_)
  gic <- -2 * qhat + 2 * penalty
  list(gic = gic, qhat = qhat, penalty = penalty)
}

## ---- per-taxon structure selection wrapper ------------------------
## Fits the requested structure(s) for ONE taxon and (for corstr="GIC")
## selects the structure minimizing GIC. Shares the independence-working
## sensitivity A_indep across the 4 candidate structures (Omega_I is,
## by definition, the same independence bread for a given taxon), so the
## Pan penalty is comparable across structures.
##
## Returns list(est.para = <chosen fit>, corr.structure = <name>,
##              gic.row = data.frame(taxon, structure, gic, qhat,
##                                   penalty, selected))
.est_para_select <- function(ID, W, Xa, Xb, Y, zi.check, period, corstr, maxlag,
                             taxon = NA_character_) {
  ## guard: taxon must be a length-1 label (colnames(Y)[x] can be NULL)
  if (length(taxon) == 0 || is.null(taxon)) taxon <- NA_character_
  taxon <- as.character(taxon)[1]
  structs_all <- c("independence", "exchangeable", "ar1", "toeplitz")
  if (corstr != "GIC") {
    ep <- .est_para(ID, W, Xa, Xb, Y, zi.check, period, corstr, maxlag)
    return(list(est.para = ep, corr.structure = corstr,
                gic.row = data.frame(taxon = taxon, structure = corstr,
                                     gic = NA_real_, qhat = NA_real_,
                                     penalty = NA_real_, selected = TRUE,
                                     stringsAsFactors = FALSE)))
  }
  ## corstr == "GIC": fit all four structures.
  fits <- vector("list", length(structs_all)); names(fits) <- structs_all
  gics <- setNames(rep(NA_real_, length(structs_all)), structs_all)
  qh   <- setNames(rep(NA_real_, length(structs_all)), structs_all)
  pen  <- setNames(rep(NA_real_, length(structs_all)), structs_all)
  for (s in structs_all) {
    fit <- tryCatch(.est_para(ID, W, Xa, Xb, Y, zi.check, period, s, maxlag),
                    error = function(e) NULL)
    fits[[s]] <- fit
  }
  ## independence-working sensitivity (Omega_I bread) shared across
  ## structures; take it from the independence fit when available, else
  ## from any available fit (A_indep depends only on fitted alpha/mu/phi,
  ## which are consistent regardless of working corr, per the paper).
  A.indep <- NULL
  ref.fit <- if (!is.null(fits[["independence"]])) fits[["independence"]] else
    fits[[which(!vapply(fits, is.null, logical(1)))[1]]]
  if (!is.null(ref.fit)) {
    ABi <- tryCatch(.aigdm_mean_A_B(ref.fit, corstr.use = "independence"),
                    error = function(e) NULL)
    A.indep <- if (is.null(ABi)) NULL else ABi$A
  }
  for (s in structs_all) {
    fit <- fits[[s]]
    if (!is.null(fit)) {
      g <- .gic_aigdm(fit, A.indep = A.indep)
      gics[s] <- g$gic; qh[s] <- g$qhat; pen[s] <- g$penalty
    }
  }
  ## choose min finite GIC; fall back to exchangeable if all NA
  if (all(!is.finite(gics))) {
    best <- "exchangeable"
  } else {
    best <- names(which.min(gics))
  }
  gic.row <- data.frame(taxon = taxon, structure = structs_all,
                        gic = as.numeric(gics[structs_all]),
                        qhat = as.numeric(qh[structs_all]),
                        penalty = as.numeric(pen[structs_all]),
                        selected = structs_all == best,
                        stringsAsFactors = FALSE)
  list(est.para = fits[[best]], corr.structure = best, gic.row = gic.row)
}
