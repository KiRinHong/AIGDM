rm(list = ls())
setwd("~/Documents/Project/LONGITUDINAL/code/Simulation/")

library(phyloseq)
library(data.table)
library(matrixStats)
library(MGLM)
#### ECAM #####
load("../Data/Deriveddata/ECAM.Species.filter10.rda")
input_data <- ecam$pall
id <- input_data$meta$studyid
tax <- input_data$tax.tab@.Data# [,1:6]
colnames(tax) <- paste0("Rank", 1:7)
Y <- input_data$otu.tab # t(input_data$otu.tab)
colOrders <- order(colMeans(Y/rowSums(Y), na.rm = T), decreasing = TRUE)
Y <- Y[,colOrders]
tax <- tax[colOrders,]
Xa <- Xb <- W <- matrix(1, nrow(Y), 1)

remove.subject <- which(rowSums(Y) <= 0)
if(length(remove.subject) > 0){
  print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth))
  Xa <- Xa[-remove.subject, , drop=FALSE]
  Xb <- Xb[-remove.subject, , drop=FALSE]
  W <- W[-remove.subject, , drop=FALSE]
  Y <- Y[-remove.subject, , drop=FALSE]
}

keep <- which(colSums(Y)>0)
Y <- Y[, keep, drop=FALSE]
tax <- tax[keep, ,drop=FALSE]

W.data <- data.table(data.frame(tax, t(Y)))
n.rank <- ncol(tax)
otucols <- names(W.data)[-(1:n.rank)]

n.level <- n.rank-1
subtree <- subtree.K <- subtree.seqDepth <- NULL
for(k in 1:n.level){
  # k <- 1
  Rank.low <- paste("Rank", n.rank-k,sep="")
  Rank.high <- paste("Rank", n.rank-k+1,sep="")
  
  tmp <- table(tax[,n.rank-k])
  level.uni <- sort( names(tmp)[which(tmp>1)] )
  m.level <- length(level.uni)    
  
  tt <- W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
  setnames(tt, 1:2, c(Rank.low, Rank.high))
  W.tax <- as.vector(unlist(tt[, Rank.low, with=FALSE]))
  W.count <- data.matrix(tt[, otucols, with=FALSE])      
  
  for(m in 1:m.level){
    # m <- 1
    Y.tmp <- t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
    remove.index <- which(colSums(Y.tmp)==0)
    if(length(remove.index)==ncol(Y.tmp)){
      next
    }else{
      if(length(remove.index)>0){
        Y.tmp <- Y.tmp[, -remove.index, drop=FALSE] 
      }
      if(ncol(Y.tmp)==1){
        next
      }else{
        print(paste0("Now processing Rank ", n.rank-k, ".", m))
        keep.ind <- which(rowSums(Y.tmp) != 0)
        Y.tmp <- Y.tmp[keep.ind,,drop=FALSE]
        print(colMeans(Y.tmp/rowSums(Y.tmp)))
        subtree <- c(subtree, paste(Rank.low, level.uni[m], sep = "."))
        subtree.K <- c(subtree.K, ncol(Y.tmp)-1)
        subtree.seqDepth <- c(subtree.seqDepth, mean(rowSums(Y.tmp)))
      }
    }
  }
}
length(subtree) # 22
mean(subtree.K) # 2.73
mean(subtree.seqDepth) # 2405


k <- 1 # 2
Rank.low <- paste("Rank", n.rank-k,sep="")
Rank.high <- paste("Rank", n.rank-k+1,sep="")

tmp <- table(tax[,n.rank-k])
level.uni <- sort( names(tmp)[which(tmp>1)] )
m.level <- length(level.uni)    

tt <- W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
setnames(tt, 1:2, c(Rank.low, Rank.high))
W.tax <- as.vector(unlist(tt[, Rank.low, with=FALSE]))
W.count <- data.matrix(tt[, otucols, with=FALSE])      

m <- 14 # 3
Y.tmp <- t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
# Y.tmp <- cbind(Y.tmp[,1:4], rowSums(Y.tmp[,5:6]))
keep.ind <- which(rowSums(Y.tmp) != 0)

Y.tmp <- Y.tmp[keep.ind,,drop=FALSE]
W.r.tmp <- Xa.r.tmp <- Xb.r.tmp <- matrix(1, nrow(Y.tmp), 1)

colOrder <- order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
Y.tmp <- Y.tmp[,colOrder]
Y.tmp <- cbind(Y.tmp[,1:4], rowSums(Y.tmp[,6:7]), Y.tmp[,5])
K <- ncol(Y.tmp) - 1
source("../Analysis/1a.AIGDMC_all_utility.R")
aigdm.reg <- .AIGDM_EM(Y.tmp, W.r.tmp, Xa.r.tmp, Xb.r.tmp,
                      matrix(-Inf, 1, K), matrix(-Inf, 1, K),
                      matrix(0.001, 1, K), matrix(0.001, 1, K))
tmp.a <- exp(Xa.r.tmp %*% aigdm.reg$alpha.est)
mu.real <- tmp.a/(1+tmp.a)

tmp.b <- exp(Xb.r.tmp %*% aigdm.reg$beta.est)
phi.real <- 1/tmp.b

colMeans(mu.real*phi.real)
colMeans((1-mu.real)*phi.real)

# # fit GDM
# gdmFit <- MGLMfit(Y.tmp, dist = "GDM")
# ab.real <- round(gdmFit@estimate, 4)
# a.real <- ab.real[1:5]
# b.real <- ab.real[6:10]
para.gdm <- list(a.real = colMeans(mu.real*phi.real), 
                 b.real = colMeans((1-mu.real)*phi.real), 
                 seqDepth.real = rowSums(Y.tmp),
                 W.real = Y.tmp)


# fit LN
library(foreach)
library(doParallel)
source("realpara_utility.R")
Y.tmp.reloc <- Y.tmp[,c(2:6,1)]
n <- dim(Y.tmp.reloc)[1] # number of samples
K <- dim(Y.tmp.reloc)[2] - 1 # number of OTUs (reference OTU not counted)
M <- as.numeric(apply(Y.tmp.reloc, 1, sum))

num_cores_StARS <- 10
iter_sub <- 100
n_sub <- floor(7 * sqrt(n))
beta <- 0.05
length_rholist <- 30

# Perform log-ratio transformation
p.hat <- colMeans(Y.tmp.reloc/rowSums(Y.tmp.reloc))
offset <- K + 1
Y.tmp.adj <- t(t(Y.tmp.reloc) + p.hat * offset)
z.hat <- log(Y.tmp.adj[, -(K + 1)]/Y.tmp.adj[, (K + 1)])


# Compo-glasso with StARS selection
Sigma.2 <- cov(z.hat)
rho.list.1 <- exp(seq(log(max(Sigma.2) * 1), log(max(Sigma.2) / 300), length = length_rholist))
n1.rho <- length(rho.list.1)
Omegas.1 <- array(0, dim = c(K, K, n1.rho))

Omegas.1.orig <- Compo_glasso(Y.tmp.reloc, rho.list = rho.list.1)$Omegas.1
Omegas.1.prob <- array(0, dim = c(K, K, n1.rho))
Omegas.1.ksi <- array(0, dim = c(K, K, n1.rho))
D_b.1 <- array(0, n1.rho)
sup_D_b.1 <- array(0, n1.rho)

for(i.rho.1 in 1 : n1.rho)    {
  cat("i.rho.1 =", i.rho.1, "\n")
  rho <- rho.list.1[i.rho.1]
  
  cl<-makeCluster(num_cores_StARS)
  registerDoParallel(cl)
  Omegas.1.sub <- foreach (i_iter_sub = 1:iter_sub, combine = cbind, .export = "obj", .packages = c("huge", "propagate")) %dopar%
    {
      source("realpara_utility.R")
      cat("StARS selection subsampling number =", i_iter_sub, "\n")
      Omega.0 <- Omegas.1.orig[, , i.rho.1]
      x_sub <- Y.tmp.adj[sample(1:n, size = n_sub), ]
      Omega.1.sub <- Compo_glasso(x_sub, option = 0, rho.list = rho, O_ratio = 100, z_ratio = 100, max_iter = 100)$Omega.1
      Omega.1.sub
    }
  stopCluster(cl)
  Omegas.1.sub.tmp = array(dim = c(K, K, iter_sub))
  for (i in 1:iter_sub){
    Omegas.1.sub.tmp[, , i] = Omegas.1.sub[[i]]
  }
  Omegas.1.sub = Omegas.1.sub.tmp
  
  Omega.1.prob <- rowMeans((Omegas.1.sub > 0) * 1, dims = 2) 
  Omegas.1.prob[, , i.rho.1] <- Omega.1.prob
  Omega.1.ksi <- 2 * Omega.1.prob * (1 - Omega.1.prob)
  Omegas.1.ksi[, , i.rho.1] <- Omega.1.ksi
  D_b.1[i.rho.1] <- (sum(Omega.1.ksi) / 2) / (K * (K - 1) / 2)
  sup_D_b.1[i.rho.1] <- max(D_b.1)
  
  if (sup_D_b.1[i.rho.1] > beta) break
}

rho.1 <- rho.list.1[i.rho.1 - 1]
Omega.1 <- Omegas.1.orig[, , i.rho.1 - 1]
compGlasso.mod <- Compo_glasso(Y.tmp.reloc, rho.list = rho.1)
colMeans(compGlasso.mod$Z.mat)
compGlasso.mod$Sigma.1
z <- compGlasso.mod$Z.mat
p <- matrix(0, nrow = n, ncol = K + 1)
p[, -(K + 1)] <- exp(z)/((apply(exp(z), 1, sum) + 1) %*% matrix(1, ncol = K))
p[, (K + 1)] <- 1/(apply(exp(z), 1, sum) + 1)

para.ln <- list(p = colMeans(p),
                mean = colMeans(compGlasso.mod$Z.mat), 
                varcov = cov(compGlasso.mod$Z.mat))

# # fit DM
# dmFit <- MGLMfit(W.rarefy, dist = "DM")
# para.dm <- list(a.real = round(dmFit@estimate, 4))
tmp.p <- para.ln$p
tmp.mean <- para.ln$mean
tmp.cov <- para.ln$varcov
B <- matrix(c(0,0,0,0,-1,
              1,0,0,0,-1,
              0,1,0,0,-1,
              0,0,1,0,-1,
              0,0,0,1,-1), nrow = 5, ncol = 5, byrow = T)
para.ln.lt <- list(p = c(tmp.p[6], tmp.p[-6]),
                   mean = c(B %*% tmp.mean),
                   varcov = B %*% tmp.cov %*% t(B))
para.list <- list(gdm = para.gdm, ln = para.ln, ln.lt = para.ln.lt)

save(para.list, file = "ProposedModel/shared/para.Rdata")
