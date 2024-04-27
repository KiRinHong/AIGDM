rm(list = ls())
########### input arguments
args = (commandArgs(trailingOnly=TRUE))
if(length(args) == 5){
  test = as.numeric(args[1]) # 1: mean 2: disp 3: both mean/disp equal to 0
  mod = args[2] # ZIGDM1 ZIGDM2 ZIGDM0 ZILN1 ZILN2 ZILN0
  sim = as.numeric(args[3]) # replicate
  r = as.numeric(args[4]) # correlation parameter between two time points # 0 0.2 0.4 0.6 0.8
  n.T = as.numeric(args[5])
} else {
  cat('usage: Rscript AIgetRslt_Omni.R <test> <mod> <sim> <r> <n.T>\n', file = stderr())
  stop()
}

# setwd("~/Documents/Project/LONGITUDINAL/code/Simulation/ProposedModel/shared/")
# rm(list = ls())
# test = 3; mod = "ZIGDM1"; sim = 1; r = 0; n.T = 5

library(MASS)
library(Matrix)
library(matrixStats)
library(HMP)
# library(miLineage)
source("sim_utility_Omni.R")
# source("ZIGDM.R")

set.seed(12+sim)
N <- 50 # number of clusters
# n.T <- 20 # number of time points
n.perm <- 1000

# Generates random samples from generalized Dirichlet Multinomial distribution n*(K+1)
# xmat: predictor for the mean model N * k.alpha
# alpha.mat: coef for the mean model k.alpha * K
xmat <- matrix(rep(c(0,1), each = N/2), ncol = 1)
xmat.m <- do.call(rbind, replicate(n.T, xmat, simplify=FALSE))
# if(type_perm == 1){ # WTH
#   xmat.1 <- matrix(rep(0,N), ncol = 1)
#   xmat.2 <- matrix(rep(1,N), ncol = 1)
#   xmat.m <- rbind(xmat.1, xmat.2)
# }else if(type_perm == 2){ # BTW
#   xmat.1 <- matrix(rep(c(0,1), each = N/2), ncol = 1)
#   xmat.2 <- xmat.1
#   xmat.m <- rbind(xmat.1, xmat.2)
# }
dx <- ncol(xmat.m)+1
ID <- rep(1:N, n.T)

load("para.Rdata")
K <- 5 # number of taxa - 1
taxa.id <- sample(1:K, 2)

if(mod %in% c("ZIGDM1", "ZILNM1")){
  zi.id <- 1
  pstr <- c(0.4, 0, 0, 0, 0)
  # pstr <- c(0.1, 0.15, 0.2, 0.25, 0.3)
  # pstr <- c(0.4, 0.5, 0.6, 0.7, 0.8)
  # pstr[-zi.id] <- rep(0, 3)
}else if(mod %in% c("ZIGDM2", "ZILNM2")){
  zi.id <- 2:K
  pstr <- c(0, 0.4, 0.4, 0.4, 0.4)
}else if(mod %in% c("ZIGDM0", "ZILNM0")){
  zi.id <- NULL
  pstr <- rep(0, K)
}

if(mod %in% c("ZIGDM0", "ZIGDM1", "ZIGDM2")){
  gamma.int <- -log(1/pstr-1)
  a.real <- para.list$gdm$a.real
  b.real <- para.list$gdm$b.real
  mu.real <- a.real/(a.real + b.real)
  sigma.real <- 1/(a.real + b.real + 1)
  
  alpha.int <- log(mu.real/(1-mu.real)) # intercept term
  beta.int <- log(sigma.real/(1-sigma.real))
  # if(mod == "ZIGDM0"){
  #   gamma.coef <- 0
  # }else{
  #   gamma.coef <- 0.4
  # }
  gamma.coef <- 0
  
  if(n.T == 5){
    if(test == 1){
      tmp.val <- 1
    }else if(test == 2){
      tmp.val <- 0.8
    }
  }else if(n.T == 20){
    if(test == 1){
      tmp.val <- 0.8
    }else if(test == 2){
      tmp.val <- 0.6
    }
  }
  
  if(test == 1){ # type1 disp; power mean
    alpha.coef <- tmp.val; beta.coef <- 0
  }else if(test == 2){ # type1 mean; power disp
    alpha.coef <- 0; beta.coef <- tmp.val
  }else if(test == 3){
    alpha.coef <- beta.coef <- 0
  }
  # gamma.trt <- alpha.trt <- beta.trt <- rep(0, K)
  # gamma.trt[zi.id] <- gamma.coef; alpha.trt[taxa.id] <- alpha.coef; beta.trt[taxa.id] <- beta.coef
  # 
  # gamma.mat <- matrix(c(gamma.int, gamma.trt), nrow = dx, ncol = K, byrow = TRUE)
  # alpha.mat <- matrix(c(alpha.int, alpha.trt), nrow = dx, ncol = K, byrow = TRUE) 
  # beta.mat <- matrix(c(beta.int, beta.trt), nrow = dx, ncol = K, byrow = TRUE) 
  
  alpha.trt <- beta.trt <- rep(0, K)
  alpha.trt[taxa.id] <- alpha.coef; beta.trt[taxa.id] <- beta.coef
  
  alpha.mat <- matrix(c(alpha.int, alpha.trt), nrow = dx, ncol = K, byrow = TRUE) 
  beta.mat <- matrix(c(beta.int, beta.trt), nrow = dx, ncol = K, byrow = TRUE) 
  
  # mu.mat: mean parameter N * K
  # phi.mat: dispersion parameter
  mu.mat <- 1/(1+exp(-(cbind(1, xmat.m) %*% alpha.mat)))
  sigma.mat <- 1/(1+exp(-(cbind(1, xmat.m) %*% beta.mat)))
  
  # a.mat, b.mat: mu = a/(a+b), sigma = 1/(a+b+1)
  # a, b refer to the parameters in a generalized Dirichlet distribution
  a.mat <- mu.mat * (1/sigma.mat - 1)
  b.mat <- (1 - mu.mat) * (1/sigma.mat - 1)
  
  # pstr.mat <- 1/(1+exp(-(cbind(1, xmat.m) %*% gamma.mat)))
  gamma.mat <- matrix(gamma.int, nrow = 1, ncol = K, byrow = TRUE) 
  pstr.mat <- 1/(1+exp(-(matrix(1, nrow(xmat.m), 1) %*% gamma.mat)))
  
  if(n.T == 5){
    Z.datalist <- .corrData.ZIBeta(n.T, r, a.mat, b.mat, pstr.mat) #ZIB
  }else if(n.T == 20){
    Z.datalist <- .corrData.ZIBeta.diab(n.T, r, a.mat, b.mat, pstr.mat) #ZIB
  }
  
  Y.datalist <- lapply(Z.datalist, .simCount.ZIGDM, SeqDepth = 1000)
  
  Y.m <- do.call(rbind, Y.datalist)
  
}else if(mod %in% c("ZILNM0", "ZILNM1", "ZILNM2")){
  mu.real <- para.list$ln.lt$mean
  sigma.real <- 1/(diag(para.list$ln.lt$varcov)+1)
  cor.real <- cov2cor(para.list$ln.lt$varcov)
  alpha.int <- mu.real # intercept term
  beta.int <- log(sigma.real/(1-sigma.real))
  if(n.T == 5){
    if(test == 1){
      tmp.val <- 1.6
    }else if(test == 2){
      tmp.val <- 1.2
    }
  }else if(n.T == 20){
    if(test == 1){
      tmp.val <- 1.2
    }else if(test == 2){
      tmp.val <- 0.8
    }
  }
  
  if(test == 1){ # type1 disp; power mean
    alpha.coef <- tmp.val; beta.coef <- 0
  }else if(test == 2){ # type1 mean; power disp
    alpha.coef <- 0; beta.coef <- tmp.val
  }else if(test == 3){
    alpha.coef <- beta.coef <- 0
  }
  # gamma.int <- -log(1/pstr-1)
  # if(mod == "ZILNM0"){
  #   gamma.coef <- 0
  # }else{
  #   gamma.coef <- 0.4
  # }
  # gamma.trt <- alpha.trt <- beta.trt <- rep(0, K)
  # gamma.trt[taxa.id] <- gamma.coef; alpha.trt[taxa.id] <- alpha.coef; beta.trt[taxa.id] <- beta.coef
  # gamma.mat <- matrix(c(gamma.int, gamma.trt), nrow = dx, ncol = K, byrow = TRUE)
  
  alpha.trt <- beta.trt <- rep(0, K)
  alpha.trt[taxa.id] <- alpha.coef; beta.trt[taxa.id] <- beta.coef
  
  alpha.mat <- matrix(c(alpha.int, alpha.trt), nrow = dx, ncol = K, byrow = TRUE) 
  beta.mat <- matrix(c(beta.int, beta.trt), nrow = dx, ncol = K, byrow = TRUE) 
  
  mu.mat <- cbind(1, xmat.m) %*% alpha.mat
  var.mat <- exp(-(cbind(1, xmat.m) %*% beta.mat))
  pstr.mat <- matrix(pstr, nrow = length(ID), ncol = K, byrow = T)
  
  # pstr.mat <- 1/(1+exp(-(cbind(1, xmat.m) %*% gamma.mat)))
  
  if(n.T == 5){
    Y.m <- .simCount.ZILN(mu.mat, var.mat, pstr.mat, cor.real, n.T, r, 1000)
  }else if(n.T == 20){
    Y.m <- .simCount.ZILN.diab(mu.mat, var.mat, pstr.mat, cor.real, n.T, r, 1000)
  }
  # Y.m <- Y.m[,order(colMeans(Y.m/rowSums(Y.m)), decreasing = TRUE)]
}

# remember to reorder the ID
id.order <- order(ID, decreasing = FALSE)
ID <- ID[id.order]
Y.m <- Y.m[id.order,]
xmat.m <- xmat.m[id.order,,drop=FALSE]

azigdm.est.mod <- AIGDM_GEE.Cluster(ID, Y.m, xmat.m, xmat.m, 1, xmat.m, "AZIGDM", zi.id,
                                    n.boot = 500, fdr.alpha = 0.05, perm.type = "BTW", n.perm = n.perm, seed = 12+sim)
zigdm.est.mod <- AIGDM_GEE.Cluster(ID, Y.m, xmat.m, xmat.m, 1, xmat.m, "ZIGDM", zi.id, 
                                   n.boot = NULL, fdr.alpha = 0.05, perm.type = "BTW", n.perm = n.perm, seed = 12+sim)
gdm.est.mod <- AIGDM_GEE.Cluster(ID, Y.m, xmat.m, xmat.m, 1, xmat.m, "GDM", zi.id, 
                                 n.boot = NULL, fdr.alpha = 0.05, perm.type = "BTW", n.perm = n.perm, seed = 12+sim)

pval.dm.resample <- DM_resample(ID, Y.m, c(xmat.m), perm.type = "BTW", fdr.alpha = 0.05, n.perm = n.perm)
mod.qcatc.1p <- QCAT.Cluster(ID = ID, OTU = Y.m, X = xmat.m, perm.type = "BTW", n.perm = n.perm, fdr.alpha=0.05, test="chisq")
Y.m.rff <- .rarefy(Y.m)
mod.qcatc.2p <- QCAT_GEE.Cluster(ID = ID, OTU = Y.m.rff, X = xmat.m, perm.type = "BTW", n.perm = n.perm, fdr.alpha=0.05, test="chisq")

check.azi1 <- azigdm.est.mod$zi == (pstr != 0)
check.azi2 <- pstr != 0

pval <- c(test, mod, sim, r, n.T, azigdm.est.mod$pval, zigdm.est.mod$pval, gdm.est.mod$pval, 
          pval.dm.resample$pval, mod.qcatc.1p$pval, mod.qcatc.2p$pval, check.azi1, check.azi2)

filename <- paste0("Test", test, "Mod", mod, "Sim", sim, "Corr", r, "T", n.T, ".txt")
write.table(t(pval), file = filename, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)