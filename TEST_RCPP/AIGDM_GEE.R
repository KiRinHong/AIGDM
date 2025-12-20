# library(extraDistr)
library(parallel)
library(matrixStats)
library(MASS)
library(data.table)

AIGDM.Cluster <- function(ID, OTU, X4mean, X4disp, X4zero, X.index, 
                      Tax = NULL, min.depth = 0, order = "random", ZI.pos = "adaptive", zi.check.pval.lst = NULL, n.boot = NULL, 
                      n.cores = 1, fdr.alpha = 0.05, n.perm = NULL, seed = 123){
  # ID=id; OTU=Y.m; X4mean=X4disp=X4zero=X.m;X.index=1;ZI.pos="adaptive";Tax=tax;min.depth=0;n.cores=10
  
  # Remove subjects with read depth less than min.depth #.seed = 123, verbose = FALSE
  remove.subject = which(rowSums(OTU) < min.depth)
  if(length(remove.subject) > 0){
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth))
    ID = ID[-remove.subject]
    Xa = X4mean[-remove.subject, , drop=FALSE]
    Xb = X4disp[-remove.subject, , drop=FALSE]
    W = X4zero[-remove.subject, , drop=FALSE]
    Y = OTU[-remove.subject, , drop=FALSE]
  }else{
    Xa = X4mean
    Xb = X4disp
    W = X4zero
    Y = OTU
  }
  
  # Remove taxa with zero total counts (make sure the order is same? )
  keep = which(colSums(Y)>0)
  Y = Y[, keep, drop=FALSE]
  
  # Add intercept terms to design matrices (check if the intercept exist? )
  Xa = cbind(1, Xa) # add the intercept term
  Xb = cbind(1, Xb) 
  W = cbind(1, W)
  
  # Adjust X.index due to intercept term
  X.index = X.index + 1
  # Create reduced design matrices without the covariate of interest (only when not NULL)
  Xa.r = Xa[, -X.index, drop=FALSE]
  Xb.r = Xb[, -X.index, drop=FALSE]
  W.r = W[, -X.index, drop=FALSE]
  
  # Perform hierarchical analysis using taxonomic information
  tax = Tax[keep, ,drop=FALSE]
  if( sum(colnames(Y)!=rownames(tax))>0 ){
    stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
  }
  W.data = data.table(data.frame(tax, t(Y)))
  n.rank = ncol(tax)
  otucols = names(W.data)[-(1:n.rank)]
  
  n.level = n.rank-1
  subtree = NULL
  for(k in 1:n.level){
    # k = 1
    Rank.low = paste("Rank", n.rank-k,sep="")
    Rank.high = paste("Rank", n.rank-k+1,sep="")
    
    tmp = table(tax[,n.rank-k])
    level.uni = sort( names(tmp)[which(tmp>1)] )
    m.level = length(level.uni)    
    
    tt = W.data[, lapply(.SD , sum, na.rm = TRUE), .SDcols = otucols, by = list( get(Rank.low), get(Rank.high) )]
    setnames(tt, 1:2, c(Rank.low, Rank.high))
    W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
    W.count = data.matrix(tt[, otucols, with=FALSE])      
    
    for(m in 1:m.level){
      # m = 1
      Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
      remove.index = which(colSums(Y.tmp)==0)
      if(length(remove.index)==ncol(Y.tmp)){
        next
      }else{
        if(length(remove.index)>0){
          Y.tmp = Y.tmp[, -remove.index, drop=FALSE] 
        }
        if(ncol(Y.tmp)==1){
          next
        }else{
          subtree = c(subtree, paste(Rank.low, level.uni[m], sep = "."))
        }
      }
    }
  } 
  R.sel = .choose_r(fdr.alpha/length(subtree), 0.05)
  
  if(is.null(zi.check.pval.lst)){
    zi.check.pval.lst = list()
    curr.ind = 1
    for(k in 1:n.level){
      # k = 4
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])
      
      for(m in 1:m.level){
        # m = 1
        Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
        remove.index = which(colSums(Y.tmp)==0)
        if(length(remove.index)==ncol(Y.tmp)){
          next
        }else{
          if(length(remove.index)>0){
            Y.tmp = Y.tmp[, -remove.index, drop=FALSE]
          }
          if(ncol(Y.tmp)==1){
            next
          }else{
            if(order == "avg_dec"){
              colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
            }else if(order == "avg_inc"){
              colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = FALSE)
            }else if(order == "prev_dec"){
              colOrder = order(colMeans(Y.tmp > 0, na.rm = T), decreasing = TRUE)
            }else if(order == "prev_inc"){
              colOrder = order(colMeans(Y.tmp > 0, na.rm = T), decreasing = FALSE)
            }else if(order == "random"){
              # colOrder = 1:ncol(Y.tmp)
              set.seed(seed)
              colOrder = sample(1:ncol(Y.tmp))
            }
            Y.tmp = Y.tmp[,colOrder]
            K = ncol(Y.tmp) - 1
            
            if(ZI.pos == "no"){
              zi.check.pval.lst[[curr.ind]] = rep(1, K) 
            }else if(ZI.pos == "all"){
              zi.check.pval.lst[[curr.ind]] = rep(0, K)
            }else if(ZI.pos == "adaptive"){
              print(paste0("### Diagnostic Test ### Now processing Rank ", n.rank-k, ".", m))
              print("Zero-inflation diagnostic test!")
              # TRUE means add zero-inflation
              zi.check.pval = sapply(1:K, function(x){
                y.tmp = cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE]))
                id.tmp = ID[rowSums(y.tmp) != 0]
                y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
                return(.check_zeroinflation(id.tmp, y.tmp, n.boot = n.boot, n.cores = n.cores, R.sel, seed))
              })
              zi.check.pval.lst[[curr.ind]] = zi.check.pval
              print(zi.check.pval)
            }else{
              stop("ZI.pos should be in c('no', 'all', 'adaptive')!")
            }
            names(zi.check.pval.lst[[curr.ind]]) = unlist(tt[W.tax == level.uni[m], ..Rank.high])[colOrder][1:K]
            curr.ind = curr.ind + 1
          }
        }
      }
    }
  }
  
  zi.check.pval.all.wo = unlist(zi.check.pval.lst)
  zi.check.pval.all.w = p.adjust(zi.check.pval.all.wo, method = "BH") 
  # zi.check.pval.all.w = unlist(lapply(zi.check.pval.lst, p.adjust, method = "BH")) 
  
  zi.check.lst = list() # stat.taxonwise.lst = list()
  pval = NULL
  curr.ind = 0
  for(k in 1:n.level){
    # k = 1
    Rank.low = paste("Rank", n.rank-k,sep="")
    Rank.high = paste("Rank", n.rank-k+1,sep="")
    
    tmp = table(tax[,n.rank-k])
    level.uni = sort( names(tmp)[which(tmp>1)] )
    m.level = length(level.uni)    
    
    tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
    setnames(tt, 1:2, c(Rank.low, Rank.high))
    W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
    W.count = data.matrix(tt[, otucols, with=FALSE])      
    
    for(m in 1:m.level){
      # m = 6
      Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
      remove.index = which(colSums(Y.tmp)==0)
      if(length(remove.index)==ncol(Y.tmp)){
        next
      }else{
        if(length(remove.index)>0){
          Y.tmp = Y.tmp[, -remove.index, drop=FALSE] 
        }
        if(ncol(Y.tmp)==1){
          next
        }else{
          print(paste0("### Hypothesis Test ### Now processing Rank ", n.rank-k, ".", m))
          keep.ind = which(rowSums(Y.tmp) != 0)
          
          Y.tmp = Y.tmp[keep.ind,,drop=FALSE]
          ID.tmp = ID[keep.ind]
          W.r.tmp = W.r[keep.ind,,drop=FALSE]; W.tmp = W[keep.ind,,drop=FALSE]
          Xa.r.tmp = Xa.r[keep.ind,,drop=FALSE]; Xa.tmp = Xa[keep.ind,,drop=FALSE]
          Xb.r.tmp = Xb.r[keep.ind,,drop=FALSE]; Xb.tmp = Xb[keep.ind,,drop=FALSE]
          
          if(order == "avg_dec"){
            colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
          }else if(order == "avg_inc"){
            colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = FALSE)
          }else if(order == "prev_dec"){
            colOrder = order(colMeans(Y.tmp > 0, na.rm = T), decreasing = TRUE)
          }else if(order == "prev_inc"){
            colOrder = order(colMeans(Y.tmp > 0, na.rm = T), decreasing = FALSE)
          }else if(order == "random"){
            # colOrder = 1:ncol(Y.tmp)
            set.seed(seed)
            colOrder = sample(1:ncol(Y.tmp))
          }
          Y.tmp = Y.tmp[,colOrder]
          K = ncol(Y.tmp) - 1
          
          zi.check = zi.check.pval.all.w[curr.ind+1:K] < 0.05
          # zi.check = zi.check.pval.all.w[which(startsWith(names(zi.check.pval.all.w), "Rank5.Pasteurellaceae"))] < 0.05
          zi.id = which(zi.check)
          curr.ind = curr.ind+K
          
          est.para.lst = lapply(1:K, function(x){
            .est_para(ID.tmp, W.r.tmp, Xa.r.tmp, Xb.r.tmp, cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])), zi.check[x])
          })
          
          # est.para.lst.tmp = lapply(1:K, function(x){
          #   .est_para(ID.tmp, W.r.tmp[,-1,drop=F], Xa.r.tmp[,-1,drop=F], Xb.r.tmp[,-1,drop=F], cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])), zi.check[x])
          # })
          # est.para.lst[[1]]$gamma0
          # est.para.lst.tmp[[1]]$gamma0

          asym.mean = .score_test_mean(ID.tmp, Xa.tmp, X.index, Y.tmp, est.para.lst)
          stat.mean.sum = asym.mean$stat.sum
          pval.mean.a = asym.mean$pval
          
          asym.disp = .score_test_disp(ID.tmp, Xb.tmp, X.index, Y.tmp, est.para.lst)
          stat.disp.sum = asym.disp$stat.sum
          pval.disp.a = asym.disp$pval
          
          if(length(zi.id) > 0){
            asym.zero = .score_test_zero(ID.tmp, W.tmp, X.index, Y.tmp, est.para.lst, zi.id)
            stat.zero.sum = asym.zero$stat.sum
            pval.zero.a = asym.zero$pval
            
            pval.omni.c.a = .ACAT(c(pval.zero.a, pval.mean.a, pval.disp.a))
            
            # stat.omni.sum = sum(stat.zero.sum, stat.mean.sum, stat.disp.sum)
          }else{
            pval.zero.a = NA
            
            pval.omni.c.a = .ACAT(c(pval.mean.a, pval.disp.a))
            
            # stat.omni.sum = sum(stat.mean.sum, stat.disp.sum)
          }
          # pval.omni.m.a = 1 - pchisq(stat.omni.sum, length(zi.id) + K*2)
          

          if(is.null(n.perm)){
            pval = cbind(pval, c(pval.zero.a, pval.mean.a, pval.disp.a, pval.omni.c.a)) # pval.omni.m.a
          }else{
            
            # compute residual forming matrix Rconf (Smith's method)
            if(ncol(W.r.tmp) == 1){
              Rconf = diag(nrow(W.r.tmp))
            }else{
              Rconf = diag(nrow(W.r.tmp)) - W.r.tmp[,-1,drop=F] %*% solve(t(W.r.tmp[,-1,drop=F]) %*% W.r.tmp[,-1,drop=F]) %*% t(W.r.tmp[,-1,drop=F])
            }
            pval.mean.p = .score_test_perm(ID.tmp, Xa.tmp, X.index, Y.tmp, est.para.lst, stat.mean.sum, test = "Mean", Rconf, n.perm, n.cores, R.sel, NULL, seed)
            pval.disp.p = .score_test_perm(ID.tmp, Xb.tmp, X.index, Y.tmp, est.para.lst, stat.disp.sum, test = "Disp", Rconf, n.perm, n.cores, R.sel, NULL, seed)
            if(length(zi.id) > 0){
              pval.zero.p = .score_test_perm(ID.tmp, W.tmp, X.index, Y.tmp, est.para.lst, stat.zero.sum, test = "Zero", Rconf, n.perm, n.cores, R.sel, zi.id, seed)
              pval.omni.c.p = .ACAT(c(pval.zero.p, pval.mean.p, pval.disp.p))
            }else{
              pval.zero.p = NA
              pval.omni.c.p = .ACAT(c(pval.mean.p, pval.disp.p))
            }
            
            # pval.omni.m.p = .score_test_perm(ID.tmp, W.tmp, X.index, Y.tmp, est.para.lst, stat.omni.sum, "Omni", Rconf, n.perm, n.cores, R.sel, zi.id, seed)
            
            pval = cbind(pval, c(pval.zero.a, pval.zero.p,
                                 pval.mean.a, pval.mean.p,
                                 pval.disp.a, pval.disp.p,
                                 pval.omni.c.a, pval.omni.c.p)) #pval.omni.m.a, pval.omni.m.p
        
          }
          
          zi.check.lst[[ncol(pval)]] = zi.check
          # stat.taxonwise.lst[[ncol(pval)]] = rbind(stat.zero.taxonwise, stat.mean.taxonwise, stat.disp.taxonwise)
          # In rbind(stat.zero.taxonwise, stat.mean.taxonwise, stat.disp.taxonwise) :
          #   number of columns of result is not a multiple of vector length (arg 1)
        }
      }
    }# lineage loop
  }# level loop
  colnames(pval) = subtree
  names(zi.check.lst) = names(zi.check.pval.lst) = subtree # names(stat.taxonwise.lst) = subtree
  # identify significant lineages
  # na.omit(pval) for .identifySigLineagesNsimesTest
  rslt = apply(pval, 1, .identifySigLineagesNsimesTest, fdr.alpha)
  
  if(is.null(n.perm)){
    names(rslt) = c("Zero-Asymptotic", "Mean-Asymptotic", "Disp-Asymptotic", "Omni-Cauchy-Asymptotic") # "Omni-MV-Asymptotic"
  }else{
    names(rslt) = c("Zero-Asymptotic", "Zero-Resampling",
                    "Mean-Asymptotic", "Mean-Resampling",
                    "Disp-Asymptotic", "Disp-Resampling",
                    "Omni-Cauchy-Asymptotic", "Omni-Cauchy-Resampling") # "Omni-MV-Asymptotic", "Omni-MV-Resampling"
  }
  rslt = c(rslt, lineage.zi = list(zi.check.lst), 
           lineage.zi.pval = list(zi.check.pval.lst)) 
           # stat.taxonwise = list(stat.taxonwise.lst)
  return(rslt)
}

.identifySigLineagesNsimesTest <- function(pval, fdr.alpha){
  subtree.tmp = names(pval)
  score.tmp = pval
  index.na = which(is.na(score.tmp))
  if(length(index.na)>0){
    score.tmp = score.tmp[-index.na]
    subtree.tmp = subtree.tmp[-index.na]
  }
  #score.tmp[score.tmp==0] = 1e-4
  m.test = length(score.tmp)
  # Benjamini-Hochberg FDR control
  index.p = order(score.tmp)
  p.sort = sort(score.tmp)
  
  reject = rep(0, m.test)
  tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
  if(length(tmp)>0){
    index.reject = index.p[1:max(tmp)]
    reject[index.reject] = 1      
  }
  
  sig.lineage = subtree.tmp[reject==1]
  
  # perform global test
  global.pval = .simes.test(score.tmp)
  names(global.pval) = "Simes"
  return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
}

.simes.test <- function(x){
  return( min(length(x) * x/rank(x)) )
}

.ZIeZ_vec <- function(pZv, av, bv, Y){
  K = length(Y)-1
  N = sum(Y)
  
  # Vectorized operations for faster computations
  av.prim = av + Y[1:K]
  bv.prim = bv + (N - cumsum(Y[1:K]))
  pZv.post = rep(0, K)
  
  # Compute pZv.post values in a vectorized way
  beta_values = ifelse(beta(av, bv) == 0, 
                       1,
                       beta(av.prim, bv.prim) / beta(av, bv))
  tmp_values = pZv + (1 - pZv) * beta_values
  pZv.post[Y[1:K] == 0] = ifelse(tmp_values[Y[1:K] == 0] == 0, 
                                 1, 
                                 pZv[Y[1:K] == 0] / tmp_values[Y[1:K] == 0])
  return(list(pv.post = pZv.post, av.post = av.prim, bv.post = bv.prim))
}
.ZIeZ_mat <- function(pZv, av, bv, Y){
  n = dim(Y)[1]
  K = dim(Y)[2] - 1
  
  # Vectorized operations for faster computations
  av.prim = av + Y[,1:K, drop=FALSE]
  bv.prim = bv + (rowSums(Y) - rowCumsums(Y[,1:K,drop=FALSE]))
  pZv.post = matrix(0, n, K)
  
  # Compute pZv.post values in a vectorized way
  beta_values = ifelse(beta(av, bv) == 0, 1, beta(av.prim, bv.prim) / beta(av, bv))
  
  tmp_values = pZv + (1 - pZv) * beta_values
  zero_indices = Y[, 1:K, drop=FALSE] == 0
  pZv.post[zero_indices] = ifelse(tmp_values[zero_indices] == 0, 1, pZv[zero_indices] / tmp_values[zero_indices])
  
  return(list(pv.post = pZv.post, av.post = av.prim, bv.post = bv.prim))
}

.LogitOptim <- function(Del, W, gamma.ini){
  # Del=DelZ.R[,1]; gamma.ini=gamma.last[,1]
  Logit.par.ini = gamma.ini
  Logit.data = list(Del=Del, W=W)
  
  return(optim(par=Logit.par.ini, fn=.LogitNegLoglik, gr=.LogitNegScore, data = Logit.data, method="L-BFGS-B")$par)
  # return(tryCatch(optim(par=Logit.par.ini, fn=.LogitNegLoglik, gr=.LogitNegScore, data = Logit.data, method="BFGS")$par,
  #                 error=function(err) optim(par=Logit.par.ini, fn=.LogitNegLoglik, gr=.LogitNegScore, data = Logit.data)$par))
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


.LogitNegScore <- function(par, data){
  tmp = as.numeric(exp(data$W %*% par))
  p = tmp/(1+tmp)
  return( -colSums((data$Del - p) * data$W))
}

.AIBetaOptim <- function(Del, A, B, Xa, Xb, alpha.ini, beta.ini){
  # Del=Del.R;A=A.R;B=B.R;alpha.ini=anew;beta.ini=c(beta0)
  Beta.par.ini = c(alpha.ini, beta.ini)
  Beta.data = list(Del = Del, A = A, B = B, Xa = Xa, Xb = Xb)
  return(tryCatch(optim(par=Beta.par.ini, fn = .AIBetaNegLoglik, gr = .AIBetaNegScore, data = Beta.data, method="BFGS")$par,
                  error=function(err) optim(par=Beta.par.ini, fn = .AIBetaNegLoglik, gr = .AIBetaNegScore, data = Beta.data)$par))
}

.AIBetaNegLoglik <- function(par, data){
  # par = Beta.par.ini; data = Beta.data
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  
  tmp = exp(data$Xa %*% alpha)
  mu.tmp = as.numeric( tmp/(1+tmp) )
  
  tmp = exp(data$Xb %*% beta)
  phi.tmp = as.numeric( tmp/(1+tmp) )
  
  a = pmax(0, (1/phi.tmp - 1) * mu.tmp)
  b = pmax(0, (1/phi.tmp - 1) * (1 - mu.tmp))
  return(-sum( (1-data$Del) * (-lbeta(a, b) + data$A * (a - 1) + data$B * (b - 1))))
}

.AIBetaNegScore <- function(par, data){
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  
  tmp.a = c(exp(data$Xa %*% alpha)) # matrix to vector
  mu.tmp = as.numeric( tmp.a/(1+tmp.a) )
  
  tmp.b = c(exp( data$Xb %*% beta))
  phi.tmp = as.numeric( tmp.b/(1+tmp.b) )
  
  a = pmax(0, (1/phi.tmp - 1) * mu.tmp)
  b = pmax(0, (1/phi.tmp - 1) * (1 - mu.tmp))
  # https://github.com/cran/betareg/blob/master/R/betareg.R
  
  zstar = data$A - data$B
  zdagger = data$B
  mustar = digamma(a) - digamma(b)
  mudagger = digamma(b) - digamma(a + b)
  one_minus_del = 1 - data$Del
  
  return(-c(colSums(one_minus_del * (1/tmp.b) * mu.tmp * (1/(1+tmp.a)) * (zstar - mustar) * data$Xa),
            colSums(one_minus_del * (-1/tmp.b) * (mu.tmp * (zstar - mustar) + (zdagger - mudagger)) * data$Xb)))
  
}

# maximum likelihood estimators of gamma, alpha and phi
.AIGDM_EM <- function(Y, W, Xa, Xb, gamma0, alpha0, beta0, tol = 1e-4, max.iter = 1000){
  n = nrow(Y)
  K = ncol(Y) - 1
  da = ncol(Xa); db = ncol(Xb)
  
  gamma.last = gamma.now = gamma0; alpha.last = alpha.now = alpha0; beta.last = beta.now = beta0
  Del.R = A.R = B.R = matrix(0, n, K)
  for (l in 1:max.iter) {
    # E-step
    # print(paste("====== ", l, "th ======", sep=""))
    expXa_alpha = exp(Xa %*% alpha.last)
    expXb_beta = exp(Xb %*% beta.last)
    expW_gamma = exp(W %*% gamma.last)
    # Compute tmpMv, tmpSv using vectorized operations
    tmpMv = expXa_alpha / (1 + expXa_alpha)
    tmpSv = expXb_beta / (1 + expXb_beta)
    # Compute tmpA and tmpB using vectorized operations
    tmpA = tmpMv * (1 / tmpSv - 1)
    tmpB = (1 - tmpMv) * (1 / tmpSv - 1)
    # Compute tmpZ using vectorized operations
    tmpZ = expW_gamma / (1 + expW_gamma)
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
      if( !is.infinite(gamma.last[1,j]) ){
        gamma.now[,j] = .LogitOptim(Del.R[,j], W, gamma.last[,j])
      }
      tmp = .AIBetaOptim(Del.R[,j], A.R[,j], B.R[,j], Xa, Xb, alpha.last[,j], beta.last[,j])
      alpha.now[,j] = tmp[1:da]; beta.now[,j] = tmp[-(1:da)]
    }
    
    diffs = abs(c(gamma.now - gamma.last, alpha.now - alpha.last, beta.now - beta.last))
    if (sum(diffs[!is.infinite(diffs)], na.rm=TRUE) < tol) break
    gamma.last = gamma.now; alpha.last = alpha.now; beta.last = beta.now
  }
  return(list(gamma.est = gamma.now, alpha.est = alpha.now, beta.est = beta.now))
}


# .dma: dmu*/dalpha, first derivative of alpha in the mu*=digamma(mu*phi)-digamma((1-mu)*phi) 
.dma <- function(xmat, alpha, mu, phi){
  K = ncol(alpha)
  tmp = c(exp(xmat %*% alpha)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  xmat_expand = do.call(rbind, replicate(K, xmat, simplify = FALSE))
  tmp1 = (1 + tmp)^2
  sigma = 1/phi-1
  return(.RowbyRow(xmat_expand, 
                   tmp/tmp1*sigma*(trigamma(mu*sigma)+trigamma((1-mu)*sigma))))
}
# .dmb: dmu*/dbeta, first derivative of beta in the mu*=digamma(mu*phi)-digamma((1-mu)*phi)
.dmb <- function(xmat, beta, mu, phi){
  K = ncol(beta)
  tmp = c(exp(xmat %*% beta)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  xmat_expand = do.call(rbind, replicate(K, xmat, simplify = FALSE))
  sigma = 1/phi-1
  return(.RowbyRow(xmat_expand, 
                   (-1/tmp)*sigma*(trigamma(mu*sigma)*mu-trigamma((1-mu)*sigma)*(1-mu))))
}
# .ddb: ddisp/dbeta, first derivative of beta in the dispersion model 
.ddb <- function(xmat, beta){ # xmat 2*2 alpha 2*5 
  K = ncol(beta)
  tmp = c(exp(xmat %*% beta)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  xmat_expand = do.call(rbind, replicate(K, xmat, simplify = FALSE))
  tmp1 = (1 + tmp)^2
  return(.RowbyRow(xmat_expand, tmp/tmp1))
}

.dpg <- function(wmat, gamma){
  K = ncol(gamma)
  tmp = c(exp(wmat %*% gamma)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  tmp[is.nan(tmp)] = 0
  wmat_expand = do.call(rbind, replicate(K, wmat, simplify = FALSE))
  tmp1 = (1 + tmp)^2
  return(.RowbyRow(wmat_expand, tmp/tmp1))
}
.corr_zero1 <- function(eDel, p, n.T){
  uist = NULL
  K = length(p)/n.T
  for(i in seq(1,n.T*K,n.T)){
    for (s in i:(i+n.T-2)) {
      for (t in (s+1):(i+n.T-1)) {
        uist = c(uist, (eDel[s]-p[s])*(eDel[t]-p[t])/sqrt(p[s]*(1-p[s])*p[t]*(1-p[t])))
      }
    }
  }
  uist[is.infinite(uist) | is.nan(uist)] = 0
  uist.mat = matrix(uist, nrow = choose(n.T, 2), ncol = K)
  return(colSums(uist.mat, na.rm = T))
}

.corr_zero2 <- function(eDel, p, n.T){
  K = length(p)/n.T
  res = (eDel-p)^2/(p*(1-p))
  res[is.infinite(res) | is.nan(res)] = 0
  res.mat = matrix(res, nrow = n.T, ncol = K)
  return(colSums(res.mat, na.rm = T))
}

.corr_rho1 <- function(eDel, ezstar, mustar, mu, phi, n.T){
  K = length(mu)/n.T
  uist = NULL
  varzstar = diag(.var_ZStar(mu, phi))
  for(i in seq(1,n.T*K,n.T)){
    for (s in i:(i+n.T-2)) {
      for (t in (s+1):(i+n.T-1)) {
        uist = c(uist, (1-eDel[s])*(1-eDel[t])*(ezstar[s]-mustar[s])*(ezstar[t]-mustar[t])/sqrt(varzstar[s]*varzstar[t]))
      }
    }
  }
  uist[is.infinite(uist) | is.nan(uist)] = 0
  uist.mat = matrix(uist, nrow = choose(n.T, 2), ncol = K)
  return(colSums(uist.mat, na.rm = T))
}
.corr_rho2 <- function(eDel, ezstar, mustar, mu, phi, n.T){
  K = length(mu)/n.T
  res = (1-eDel)^2*(ezstar-mustar)^2/diag(.var_ZStar(mu, phi))
  res[is.infinite(res) | is.nan(res)] = 0
  res.mat = matrix(res, nrow = n.T, ncol = K)
  return(colSums(res.mat, na.rm = T))
}

# .var_Z: variance of Z in mean model (matrix Bi=Diag(mu_i*(1-mu_i)/(phi_i+1)))
.var_ZStar <- function(mu, phi){
  sigma = (1/phi-1)
  return(diag(trigamma(mu*sigma) + trigamma((1-mu)*sigma), length(mu)))
}
.var_mu <- function(mu){
  return(diag(c(mu*(1-mu)), length(mu)))
}
# .var_logit
.var_logit <- function(p) {
  return(diag(p*(1-p), length(p)))
}
# .BlockDiag: change 10*1/2/3 matrix to 10*5/10/15 diagonal block matrix
.BlockDiag <- function(mat, K, n.T){
  tmp_nx = ncol(mat)
  J = matrix(1, n.T, tmp_nx) # n.T denotes # of time points
  mat_expand = do.call(cbind, replicate(K, mat, simplify = FALSE))
  return(diag(K) %x% J * mat_expand)
}
.colwise.cbind <- function(mat1, mat2, para.ind){
  tmp = matrix(0, nrow = length(para.ind), ncol = ncol(mat1)+ncol(mat2))
  tmp[,-para.ind] = mat1
  tmp[,para.ind] = mat2
  return(tmp)
}

.DVD <- function(xamat, xbmat, alpha, beta, pZ, mu, phi, one, two, K, n.T) {

  tmp.a = c(exp(xamat %*% alpha)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  tmp1.a = (1 + tmp.a)^2
  tmp.b = c(exp(xbmat %*% beta))
  sigma = 1/phi-1
  
  psi1 = trigamma(mu * sigma)
  psi2 = trigamma((1 - mu) * sigma)
  # auxiliary transformations
  a = psi1 + psi2
  b = psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(sigma)
  
  b.tmp = 1/tmp.b
  a0.tmp = 1/(1+tmp.a)
  a1.tmp = tmp.a/(1+tmp.a)
  a2.tmp = tmp.a/(1+tmp.a)^2
  a3.tmp = tmp.a/(1+tmp.a)^3
  a.a2 = b.tmp * a3.tmp * (1-tmp.a) 
  a.ab = - b.tmp * a2.tmp
  a.b2 = b.tmp * a1.tmp
  b.a2 = - b.tmp * a3.tmp * (1-tmp.a) 
  b.ab = b.tmp * a2.tmp
  b.b2 = b.tmp * a0.tmp
  
  waa = - (tmp.a/tmp1.a)^2 * sigma^2 * a * (1-pZ) + one * a.a2 + two * b.a2
  xamat_expand = do.call(rbind, replicate(K, xamat, simplify = FALSE))
  Kaa = t(.BlockDiag(xamat_expand, K, n.T)) %*% diag(waa, length(waa)) %*% (.BlockDiag(xamat_expand, K, n.T))
  
  wab = - tmp.a/tmp1.a * (-1/tmp.b) * sigma * (mu * a - psi2) * (1-pZ) + one * a.ab + two * b.ab 
  Kab = t(xamat) %*% diag(wab, length(wab)) %*% xbmat
  
  wbb = - (1/tmp.b)^2 * b * (1-pZ) + one * a.b2 + two * b.b2
  Kbb = t(xbmat) %*% diag(wbb, length(wbb)) %*% xbmat
  
  return(list(Kab = Kab, Kbb = Kbb))
}


.getSubDiag <- function(m){
  diag(m) = 0
  return(m/2)
  # diag(m[-nrow(m),-1], nrow(m)-1)
}

.BlockDiagBind <- function(...){
  d = list(...)
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
.RowbyRow<-function(A, b){  
  temp <- A
  for (i in 1:length(A[,1])) temp[i,] <- A[i,]*b[i]
  return(temp)
}


.BB.log.lik <- function(Y){ # fit beta-binomial with intercept only
  X = matrix(1, nrow(Y), 1)
  alpha0 = beta0 = matrix(0.001, 1, 1)
  y.ori = Y[,1]
  m.ori = rowSums(Y)
  
  # BB Model
  gamma0 = matrix(-Inf, 1, 1)
  gdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, alpha0 = alpha0, beta0 = beta0)
  mu = exp(X %*% gdm.reg$alpha.est) / (1 + exp(X %*% gdm.reg$alpha.est))
  phi = exp(X %*% gdm.reg$beta.est) / (1 + exp(X %*% gdm.reg$beta.est))
  a.mat = pmax(mu * (1/phi - 1), 0)  #  conditional assignments with pmax
  b.mat = pmax((1 - mu) * (1/phi - 1), 0)
  
  a = a.mat[,1]
  b = b.mat[,1]
  ans = sum( lgamma(m.ori + 1) + lgamma(a + b) - lgamma(a) - lgamma(b) - 
               lgamma(a + m.ori + b) + lgamma(y.ori + a) + lgamma(m.ori - y.ori + b) - 
               lgamma(y.ori + 1) - lgamma(m.ori - y.ori + 1) )
  return(list(ans = ans, a.mat = a.mat, b.mat = b.mat))
}

.BBnZIBB.log.lik <- function(Y){ # fit both beta-binomial and zero-inflated beta-binomial with intercept only
  X = matrix(1, nrow(Y), 1)
  alpha0 = beta0 = matrix(0.001, 1, 1)
  y.ori = Y[,1]
  m.ori = rowSums(Y)
  
  # BB Model
  gamma0 = matrix(-Inf, 1, 1)
  gdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, alpha0 = alpha0, beta0 = beta0)
  mu = exp(X %*% gdm.reg$alpha.est) / (1 + exp(X %*% gdm.reg$alpha.est))
  phi = exp(X %*% gdm.reg$beta.est) / (1 + exp(X %*% gdm.reg$beta.est))
  a.mat = pmax(mu * (1/phi - 1), 0)  #  conditional assignments with pmax
  b.mat = pmax((1 - mu) * (1/phi - 1), 0)
  
  a = a.mat[,1]
  b = b.mat[,1]
  bb.ans = sum( lgamma(m.ori + 1) + lgamma(a + b) - lgamma(a) - lgamma(b) - 
                  lgamma(a + m.ori + b) + lgamma(y.ori + a) + lgamma(m.ori - y.ori + b) - 
                  lgamma(y.ori + 1) - lgamma(m.ori - y.ori + 1) )
  
  # ZIBB Model
  gamma0 = matrix(0.001, 1, 1)
  zigdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, alpha0 = alpha0, beta0 = beta0)
  mu = exp(X %*% zigdm.reg$alpha.est) / (1 + exp(X %*% zigdm.reg$alpha.est))
  phi = exp(X %*% zigdm.reg$beta.est) / (1 + exp(X %*% zigdm.reg$beta.est))
  a.mat = pmax(mu * (1/phi - 1), 0)
  b.mat = pmax((1 - mu) * (1/phi - 1), 0)
  
  a.ori = a.mat[,1]
  b.ori = b.mat[,1]
  id.sel = which(y.ori > 0)
  y = y.ori[id.sel]
  a = a.ori[id.sel]
  b = b.ori[id.sel]
  m = m.ori[id.sel]
  
  logA = lgamma(a+m+b) + lgamma(b)
  logB = lgamma(a+b) + lgamma(m+b)
  zibb.ans = sum( - log(pmax(1 - exp(logB - logA), .Machine$double.eps)) - logA + lgamma(m+1) +  # Prevent potential numerical underflow
                    lgamma(y+a) + lgamma(m-y+b) + lgamma(a+b) - lgamma(y+1) -
                    lgamma(m-y+1) - lgamma(a) )
  
  return(list(bb.ans = bb.ans, zibb.ans = zibb.ans))
}

.check_zeroinflation <- function(ID, Y, n.boot = 9999, n.cores, R.sel = R.sel, seed){
  Y.seqDepth = rowSums(Y)
  obs.ans = .BBnZIBB.log.lik(Y)
  stat = -2 * (obs.ans$bb.ans - obs.ans$zibb.ans)
  if (is.infinite(stat)) {
    cat("Use Pseudo-ECDF approximation p-value\n")
    return(1/(n.boot+1))
  }
  
  id.unique = unique(ID)
  # Create an index list for each unique ID
  indexList = split(1:nrow(Y), ID)
  # Parallel calculation using mclapply
  n.boot = ceiling(n.boot/n.cores)*n.cores
  m = 0
  Nexc = 0
  stat.boot = numeric(n.boot)
  while (all(Nexc < R.sel, m < n.boot)) {
    results = mclapply(1:n.cores, function(iter){
      set.seed(seed+m+iter)
      ID.tmp = sort(sample(id.unique, replace = TRUE))
      Y.b.list = lapply(ID.tmp, function(id) Y[indexList[[as.character(id)]],, drop=FALSE])
      Y.b.tmp = do.call(rbind, Y.b.list)
      
      bb.ll = .BB.log.lik(Y.b.tmp)
      a.mat = matrix(colMeans(bb.ll$a.mat), nrow(Y), 1); b.mat = matrix(colMeans(bb.ll$b.mat), nrow(Y), 1)
      Y.b = .simData.GDM(a.mat, b.mat, Y.seqDepth)
      
      boot.obs = .BBnZIBB.log.lik(Y.b)
      return(-2 * (boot.obs$bb.ans - boot.obs$zibb.ans))
    }, mc.cores = n.cores)
    Nexc = Nexc + sum(unlist(results) >= stat)
    stat.boot[m + 1:n.cores] = unlist(results)
    m = m + n.cores
  }
  
  if (m < n.boot) {
    pval = Nexc/m
    cat(sprintf("# of bootstraps: %g\n", m))
    cat("Use ECDF approximation p-value\n")
  } else {
    if (Nexc <= 10) {
      pval = tryCatch(.gpd_approx(stat.boot, 250, stat), error=function(err) NA)
      cat(sprintf("# of bootstraps: %g\n", n.boot))
      if (is.na(pval)) {
        pval = (Nexc+1)/(n.boot+1)
        cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
      } else {
        cat("Use GPD approximation p-value\n")
      }
    } else {
      pval = Nexc/n.boot
      cat(sprintf("# of bootstraps: %g\n", n.boot))
      cat("Use ECDF approximation p-value\n")
    }
  }
  # print(stat)
  # print(stat.boot)
  return(pval)
}

.GenerateScoreFun <- function(par, data){
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  
  tmp.a = as.numeric(exp(data$Xa %*% alpha))
  mu.tmp = as.numeric( tmp.a/(1+tmp.a) )
  
  tmp.b = as.numeric( exp( data$Xb %*% beta) )
  phi.tmp = as.numeric( tmp.b/(1+tmp.b) )
  
  a = (1/phi.tmp - 1) * mu.tmp
  b = (1/phi.tmp - 1) * (1 - mu.tmp)
  a[a < 0] = 0
  b[b < 0] = 0
  # https://github.com/cran/betareg/blob/master/R/betareg.R
  
  zstar = data$A - data$B
  zdagger = data$B
  mustar = digamma(a) - digamma(b)
  mudagger = digamma(b) - digamma(a + b)
  
  return((1-data$Del) * (-1/tmp.b) * (mu.tmp * (zstar - mustar) + (zdagger - mudagger)))
}

.est_para <- function(ID, W, Xa, Xb, Y, zi.check){
  #ID=ID.tmp; W=W.r.tmp; Xa=Xa.r.tmp; Xb=Xb.r.tmp; Y=cbind(Y.tmp[,1], rowSums(Y.tmp[,(1+1):(K+1),drop=FALSE])); zi.check=T
  da = ncol(Xa); db = ncol(Xb); dw = ncol(W); K = ncol(Y)-1
  
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
  nt = nrow(Y) # the total number of observed units
  N1 = sum(n.reps*(n.reps-1)/2)
  
  r.step = 0.05 # set up the tuning parameter for updating alpha
  max = 1e3 # maximal iterations
  
  if(zi.check){
    gamma0 = matrix(0.001, nrow = dw, ncol = K)
  }else{
    gamma0 = matrix(-Inf, nrow = dw, ncol = K)
  }
  alpha0 = matrix(0.001, nrow = da, ncol = K)
  beta0 = matrix(0.001, nrow = db, ncol = K)
  # correlation of same subject, all taxa, at different times
  rho0 = sig0 = 0.1

  aigdm.reg = .AIGDM_EM(Y = Y, W = W, Xa = Xa, Xb = Xb,
                        gamma0 = gamma0, alpha0 = alpha0, beta0 = beta0)
  gamma0 = aigdm.reg$gamma.est
  alpha0 = aigdm.reg$alpha.est
  beta0 = aigdm.reg$beta.est

  ng = no = dw; na = da
  for (l in 1:max) {
    all.ddg = matrix(0, ng, ng) # derivative of GEE
    all.geeg = matrix(0, ng, 1) # GEE
    all.sig1 = all.sig2 = 0
    
    all.dda = matrix(0, na, na)
    all.geea = matrix(0, na, 1)
    all.rho1 = all.rho2 = 0 # 1:numerator 2:denominator
    
    Del.R.lst = A.R.lst = B.R.lst = A2.R.lst = B2.R.lst = AB.R.lst = vector("list", N)
    for (i in 1:N) {
      Del.R.lst[[i]] = A.R.lst[[i]] = B.R.lst[[i]] = 
        A2.R.lst[[i]] = B2.R.lst[[i]] = AB.R.lst[[i]] = matrix(NA, n.reps[i], 1)
    }
    
    N2 = ntot = 0
    for(i in 1:N){
      ID.sel = which(ID == id.unique[i])
      w.tmp = W[ID.sel,,drop=FALSE]
      xa.tmp = Xa[ID.sel,,drop=FALSE]
      xb.tmp = Xb[ID.sel,,drop=FALSE]
      y.tmp = Y[ID.sel,,drop=FALSE]
      
      tmp = exp(w.tmp %*% gamma0)
      tmp[is.na(tmp)] = 0
      pZ0.mat = tmp/(1+tmp)
      pZ0.mat[is.infinite(tmp) & tmp>0] = 1
      pZ0 = c(pZ0.mat)
      
      tmp = exp(xa.tmp %*% alpha0)
      mu0.mat = tmp/(1+tmp)
      tmp = exp(xb.tmp %*% beta0)
      phi0.mat = tmp/(1+tmp)
      
      av = mu0.mat * (1/phi0.mat - 1)
      bv = (1 - mu0.mat) * (1/phi0.mat - 1)
      
      n.T = n.reps[i]
      eDelZ0.mat = eZ0.mat = matrix(NA, n.T, 1)
      eZStar.mat = matrix(NA, n.T, 1)
      for (t in 1:n.T) {
        par.post = .ZIeZ_vec(pZ0.mat[t,,drop=FALSE], av[t,,drop=FALSE], bv[t,,drop=FALSE], y.tmp[t,])
        eDelZ0.mat[t,] = par.post$pv.post
        
        tmp = par.post$av.post + par.post$bv.post
        # post.expectation of logZ conditional on delta==0
        A.R.lst[[i]][t,] = digamma(par.post$av.post) - digamma(tmp)
        # post.expectation of log(1-Z) conditional on delta==0
        B.R.lst[[i]][t,] = digamma(par.post$bv.post) - digamma(tmp)
        eZ0.mat[t,] = par.post$av.post/tmp
        # post expectation of log(Z/(1-Z)) conditional on delta==0
        eZStar.mat[t,] = digamma(par.post$av.post) - digamma(par.post$bv.post)
      }
      
      eDelZ0 = c(eDelZ0.mat); Del.R.lst[[i]] = eDelZ0.mat
      
      ntot = ntot + colSums((1-eDelZ0.mat)^2)
      N2 = N2 + sum(.getSubDiag( (1-eDelZ0.mat)%*%t(1-eDelZ0.mat)))
      
      # Vdel
      cor.DelZ = diag(sig0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(sig0,each=n.T), n.T)
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      tmp.g = t(.BlockDiag(.dpg(w.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ)
      geeg = as.matrix(apply(.RowbyRow(t(tmp.g), eDelZ0 - pZ0), 2, sum))
      temp.g1 = tmp.g %*% .BlockDiag(.dpg(w.tmp, gamma0), 1, n.T) + geeg %*% t(geeg)
      all.ddg = all.ddg + temp.g1 ## sum up to N for derivative of gee for gamma
      all.geeg = all.geeg + geeg  ## sum up to N for gee for gamma
      
      # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
      eZ0 = c(eZ0.mat); mu0 = c(mu0.mat); phi0 = c(phi0.mat)
      eZStar = c(eZStar.mat); muStar = digamma(av) - digamma(bv)
      
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      tmp.a = t(.BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDelZ0), length(eDelZ0))
      geea = as.matrix(apply(.RowbyRow(t(tmp.a), eZStar - muStar), 2, sum))
      temp.a1 = tmp.a %*% .BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T) + geea %*% t(geea)
      all.dda = all.dda + temp.a1 # sum up to N for derivative of gee for alpha
      all.geea = all.geea + geea # sum up to N for gee for alpha
      if(n.T != 1){
        # update sig estimate
        all.sig1 = all.sig1 + .corr_zero1(eDelZ0, pZ0, n.T)
        all.sig2 = all.sig2 + .corr_zero2(eDelZ0, pZ0, n.T)
        
        # update rho estimate
        all.rho1 = all.rho1 + .corr_rho1(eDelZ0, eZStar, muStar, mu0, phi0, n.T)
        all.rho2 = all.rho2 + .corr_rho2(eDelZ0, eZStar, muStar, mu0, phi0, n.T)
      }
    }
    
    ########## updating scheme #############
    if(l==1){all.geeg1 <- all.geeg} else{ all.geeg1 <- ((l-1)*all.geeg1 + all.geeg)/l} 
    if(l==1){all.ddg1 <- all.ddg} else{ all.ddg1 <- ((l-1)*all.ddg1 + all.ddg)/l} 
    gnew0 = c(gamma0) + r.step * ginv(all.ddg1) %*% all.geeg
    gnew = (l * c(gamma0) + gnew0)/(l + 1)
    
    if(l==1){all.geea1 <- all.geea} else{ all.geea1 <- ((l-1)*all.geea1 + all.geea)/l} 
    if(l==1){all.dda1 <- all.dda} else{ all.dda1 <- ((l-1)*all.dda1 + all.dda)/l}
    anew0 = c(alpha0) + r.step * ginv(all.dda1) %*% all.geea # matrix to vector
    anew = (l * c(alpha0) + anew0)/(l + 1)
    
    Del.R = do.call(c, Del.R.lst)
    A.R = do.call(c, A.R.lst)
    B.R = do.call(c, B.R.lst)
    tmp = .AIBetaOptim(Del.R, A.R, B.R, Xa, Xb, anew, c(beta0))
    bnew = tmp[-(1:da)]
    
    sig.new = (all.sig1/N1)/(all.sig2/nt) #Note: p can be zero then sig.new = NaN
    nan.id = which(is.nan(sig.new) | is.infinite(sig.new) | sig.new == 1)
    if(length(nan.id) != 0)
      sig.new[nan.id] = 0
    
    rho.new = (all.rho1/N2)/(all.rho2/ntot)
    nan.id = which(is.nan(rho.new) | is.infinite(rho.new))
    if(length(nan.id) != 0)
      rho.new[nan.id] = 0
    
    diffs = abs(c(gnew - gamma0, anew - alpha0, bnew - beta0, sig.new - sig0, rho.new - rho0))
    if (max(diffs[!is.infinite(diffs)], na.rm=TRUE) < 1e-4 | l == max) break # max or sum?
    
    gamma0 = matrix(gnew, ncol = K); sig0 = sig.new
    alpha0 = matrix(anew, ncol = K); beta0 = matrix(bnew, ncol = K); rho0 = rho.new
  }
  
  gamma0 = matrix(gnew, ncol = K); sig0 = sig.new
  alpha0 = matrix(anew, ncol = K); beta0 = matrix(bnew, ncol = K); rho0 = rho.new
  
  Del.R.lst = eDelZ0.lst.store = pZ0.lst.store = 
    A.R.lst = B.R.lst = A2.R.lst = B2.R.lst = AB.R.lst = 
    eZ0.lst.store = mu0.lst.store = phi0.lst.store = 
    eZStar.lst.store = muStar.lst.store =
    pseudoll.score.beta.lst.store = vector("list", N)
  for (i in 1:N) {
    Del.R.lst[[i]]  = A.R.lst[[i]] = B.R.lst[[i]] = 
      A2.R.lst[[i]] = B2.R.lst[[i]] = AB.R.lst[[i]] = matrix(NA, n.reps[i], 1)
      eDelZ0.lst.store[[i]] = eZ0.lst.store[[i]] = mu0.lst.store[[i]] = phi0.lst.store[[i]] =
      eZStar.lst.store[[i]] = muStar.lst.store[[i]] = numeric(n.reps[i])
    pseudoll.score.beta.lst.store[[i]] = numeric(nrow(beta0))
  }
  
  N2 = ntot = 0
  for(i in 1:N){
    ID.sel = which(ID == id.unique[i])
    w.tmp = W[ID.sel,,drop=FALSE]
    xa.tmp = Xa[ID.sel,,drop=FALSE]
    xb.tmp = Xb[ID.sel,,drop=FALSE]
    y.tmp = Y[ID.sel,,drop=FALSE]
    
    tmp = exp(w.tmp %*% gamma0)
    tmp[is.na(tmp) | is.infinite(tmp)] = 0
    pZ0.mat = tmp/(1+tmp)
    pZ0.mat[is.infinite(tmp) & tmp>0] = 1
    pZ0 = c(pZ0.mat)
    pZ0.lst.store[[i]] = pZ0
    
    tmp = exp(xa.tmp %*% alpha0)
    mu0.mat = tmp/(1+tmp)
    tmp = exp(xb.tmp %*% beta0)
    phi0.mat = tmp/(1+tmp)
    av = mu0.mat * (1/phi0.mat - 1)
    bv = (1 - mu0.mat) * (1/phi0.mat - 1)
    muStar.mat = digamma(av) - digamma(bv)
    
    n.T = n.reps[i]
    eDelZ0.mat = eZ0.mat = eZStar.mat = matrix(NA, n.T, 1)
    for (t in 1:n.T) {
      par.post = .ZIeZ_vec(pZ0.mat[t,,drop=FALSE], av[t,,drop=FALSE], bv[t,,drop=FALSE], y.tmp[t,])
      eDelZ0.mat[t,] = par.post$pv.post
      
      tmp = par.post$av.post + par.post$bv.post
      A.R.lst[[i]][t,] = digamma(par.post$av.post) - digamma(tmp)
      B.R.lst[[i]][t,] = digamma(par.post$bv.post) - digamma(tmp)
      # post.expectation of logZ*logZ conditional on delta==0  
      A2.R.lst[[i]][t,] = trigamma(par.post$av.post) - trigamma(tmp) + ( digamma(par.post$av.post) - digamma(tmp) )^2
      # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
      B2.R.lst[[i]][t,] = trigamma(par.post$bv.post) - trigamma(tmp) + ( digamma(par.post$bv.post) - digamma(tmp) )^2
      # post.expectation of logZ*log(1-Z) conditional on delta==0  
      AB.R.lst[[i]][t,] = -trigamma(tmp) + ( digamma(par.post$av.post) - digamma(tmp) ) * ( digamma(par.post$bv.post) - digamma(tmp) )
      eZ0.mat[t,] = par.post$av.post/tmp
      eZStar.mat[t,] = digamma(par.post$av.post) - digamma(par.post$bv.post)
    }
    
    Del.R.lst[[i]] = eDelZ0.mat
    eDelZ0.lst.store[[i]] = c(eDelZ0.mat)
    
    beta.par = c(alpha0, beta0)
    beta.data = list(Del = c(Del.R.lst[[i]]), A = c(A.R.lst[[i]]), B = c(B.R.lst[[i]]), 
                     Xa = xa.tmp, Xb = xb.tmp)
    pseudoll.score.beta.lst.store[[i]] = .GenerateScoreFun(beta.par, beta.data)
    
    ntot = ntot + colSums((1-eDelZ0.mat)^2)
    N2 = N2 + sum(.getSubDiag( (1-eDelZ0.mat)%*%t(1-eDelZ0.mat)))
    
    # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
    eZ0 = c(eZ0.mat); mu0 = c(mu0.mat); phi0 = c(phi0.mat)
    eZStar = c(eZStar.mat); muStar = c(muStar.mat)
    eZStar.lst.store[[i]] = eZStar; muStar.lst.store[[i]] = muStar
    eZ0.lst.store[[i]] = eZ0; mu0.lst.store[[i]] = mu0; phi0.lst.store[[i]] = phi0
  }
  
  #print(paste0("Converge?", ifelse(l == max, 0, 1)))
  return(list(CONV = ifelse(l == max, 0, 1), W = W, Xa = Xa, Xb = Xb,
              gamma0 = matrix(gnew, ncol = K), alpha0 = matrix(anew, ncol = K), beta0 = matrix(bnew, ncol = K), 
              sig0 = sig.new, rho0 = rho.new, eDelZ0.lst = eDelZ0.lst.store,
              A.R.lst = A.R.lst, B.R.lst = B.R.lst,
              A2.R.lst = A2.R.lst, B2.R.lst = B2.R.lst, AB.R.lst = AB.R.lst,
              pZ0.lst = pZ0.lst.store, eZ0.lst = eZ0.lst.store, mu0.lst = mu0.lst.store, phi0.lst = phi0.lst.store,
              eZStar.lst = eZStar.lst.store, muStar.lst = muStar.lst.store,
              pseudoll.score.beta.lst = pseudoll.score.beta.lst.store))
}

.score_test_zero <- function(ID, X, X.index, Y, est.para.lst, zi.id){
  K = ncol(Y)-1
  K.zi = length(zi.id)
  stat = df = numeric(K.zi)
  for (j.zi in 1:K.zi) {
    j = zi.id[j.zi]
    Y.tmp = cbind(Y[,j], rowSums(Y[,(j+1):(K+1),drop=FALSE]))
    
    est.para = est.para.lst[[j]]
    ID.tmp = ID; W.tmp = X; W.r.tmp = est.para$W; Xa.r.tmp = est.para$Xa; Xb.r.tmp = est.para$Xb
    dw = ncol(W.tmp); dw.r = ncol(W.r.tmp); da.r = ncol(Xa.r.tmp); db.r = ncol(Xb.r.tmp)
    
    id.unique = sort(unique(ID.tmp)); N = length(id.unique); n.reps = as.vector(table(ID.tmp))
    ng = dw
    
    para.ind = sort(sapply(X.index, function(x) seq(x, ng, ng)))
    
    gamma0 = est.para$gamma0; alpha0 = est.para$alpha0; beta0 = est.para$beta0
    sig0 = est.para$sig0; rho0 = est.para$rho0
    eDelZ0.lst = est.para$eDelZ0.lst
    A.R.lst = est.para$A.R.lst; B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst; B2.R.lst = est.para$B2.R.lst; AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst; 
    eZ0.lst = est.para$eZ0.lst; mu0.lst = est.para$mu0.lst; phi0.lst = est.para$phi0.lst
    eZStar.lst = est.para$eZStar.ls; muStar.lst = est.para$muStar.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, ng, ng)
    ind.g = 1:ng
    
    for(i in 1:N){
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel,,drop=FALSE]
      w.tmp = W.tmp[ID.sel,,drop=FALSE]
      n.T = n.reps[i]
      
      pZ0 = pZ0.lst[[i]] # unlist(pZ0.lst[i,])
      eDelZ0 = eDelZ0.lst[[i]]
      cor.DelZ = diag(sig0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(sig0,each=n.T), n.T)
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      gamma0.w = matrix(0, dw, 1)
      gamma0.w[1:dw.r,] = gamma0
      tmp.g.w = t(.BlockDiag(.dpg(w.tmp, gamma0.w), 1, n.T)) %*% ginv(VdelZ)
      dvd.g.w = t(.BlockDiag(.dpg(w.tmp, gamma0.w), 1, n.T)) %*% ginv(VdelZ) %*% .BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)
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
    stat[j.zi] = c(t(U.r[para.ind,,drop=F]) %*% ginv(H.r %*% M.r %*% t(H.r)) %*% U.r[para.ind,,drop=F])
    df[j.zi] = length(para.ind)
  }
  stat.sum = sum(stat); df.sum = sum(df)
  pval = 1 - pchisq(stat.sum, df.sum)
  return(list(stat.sum = stat.sum, pval.sum = pval, stat = stat))
}

.score_test_mean <- function(ID, X, X.index, Y, est.para.lst){
  # X = Xa
  K = ncol(Y)-1
  stat = df = numeric(K)
  # ind.tmp = NULL
  for (j in 1:K) {
    # j = 1
    Y.tmp = cbind(Y[,j], rowSums(Y[,(j+1):(K+1),drop=FALSE]))
    
    est.para = est.para.lst[[j]]
    ID.tmp = ID; W.r.tmp = est.para$W; Xa.tmp = X; Xa.r.tmp = est.para$Xa; Xb.r.tmp = est.para$Xb
    dw.r = ncol(W.r.tmp); db.r = ncol(Xb.r.tmp); da = ncol(Xa.tmp); da.r = ncol(Xa.r.tmp)
    
    id.unique = sort(unique(ID.tmp)); N = length(id.unique); n.reps = as.vector(table(ID.tmp))
    na = da; nb = db.r
    
    para.ind = sort(sapply(X.index, function(x) seq(x, na, na)))
    
    gamma0 = est.para$gamma0; alpha0 = est.para$alpha0; beta0 = est.para$beta0
    sig0 = est.para$sig0; rho0 = est.para$rho0
    eDelZ0.lst = est.para$eDelZ0.lst; pZ0.lst = est.para$pZ0.lst
    A.R.lst = est.para$A.R.lst; B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst; B2.R.lst = est.para$B2.R.lst; AB.R.lst = est.para$AB.R.lst
    eZ0.lst = est.para$eZ0.lst; mu0.lst = est.para$mu0.lst; phi0.lst = est.para$phi0.lst
    eZStar.lst = est.para$eZStar.ls; muStar.lst = est.para$muStar.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, na+nb, na+nb)
    
    ind.a = 1:na; ind.b = (na+1):(na+nb)
    
    for(i in 1:N){
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel,,drop=FALSE]
      xa.r.tmp = Xa.r.tmp[ID.sel,,drop=FALSE]
      xa.tmp = Xa.tmp[ID.sel,,drop=FALSE]
      xb.r.tmp = Xb.r.tmp[ID.sel,,drop=FALSE]
      n.T = n.reps[i]
      
      pZ0 = pZ0.lst[[i]] # unlist(pZ0.lst[i,])
      
      mu0 = mu0.lst[[i]]; phi0 = phi0.lst[[i]]
      eZ0 = eZ0.lst[[i]]; eDelZ0 = eDelZ0.lst[[i]]
      eZStar = eZStar.lst[[i]]; muStar = muStar.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      alpha0.x = matrix(0, da, 1)
      alpha0.x[1:da.r,] = alpha0
      tmp.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDelZ0), length(eDelZ0))
      dvd.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDelZ0), length(eDelZ0)) %*% .BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZStar - muStar), 2, sum))
      
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.r.tmp)
      # dvd.b.x = t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.r.tmp)) %*% t(t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.r.tmp)))
      geeb.x = as.matrix(pseudoll.score.beta)
      
      dvd.ab.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDelZ0), length(eDelZ0)) %*% .BlockDiag(.dmb(xb.r.tmp, beta0, mu0, phi0), 1, n.T)
      
      ES = rbind(geea.x, geeb.x)
      U = U + ES
      
      a = (1/phi0 - 1) * mu0; a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0); b[b < 0] = 0
      A = digamma(a) - digamma(a+b)
      B = digamma(b) - digamma(a+b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1-eDelZ0)*( - A + A.post )
      two.R = (1-eDelZ0)*( - B + B.post )
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
    stat[j] = c(t(U.r[para.ind,,drop=F]) %*% ginv(H.r %*% M.r %*% t(H.r)) %*% U.r[para.ind,,drop=F])
    df[j] = length(para.ind)
    
  }
  stat.sum = sum(stat); df.sum = sum(df)
  pval = 1 - pchisq(stat.sum, df.sum)
  return(list(stat.sum = stat.sum, pval.sum = pval, stat = stat))
}

.score_test_disp <- function(ID, X, X.index, Y, est.para.lst){
  # X = Xb
  K = ncol(Y)-1
  stat = df = numeric(K)
  # ind.tmp = NULL
  for (j in 1:K) {
    Y.tmp = cbind(Y[,j], rowSums(Y[,(j+1):(K+1),drop=FALSE]))
    
    est.para = est.para.lst[[j]]
    ID.tmp = ID; W.r.tmp = est.para$W; Xb.tmp = X; Xb.r.tmp = est.para$Xb; Xa.r.tmp = est.para$Xa
    dw.r = ncol(W.r.tmp); da.r = ncol(Xa.r.tmp); db = ncol(Xb.tmp); db.r = ncol(Xb.r.tmp)
    
    id.unique = sort(unique(ID.tmp)); N = length(id.unique); n.reps = as.vector(table(ID.tmp))
    ng = no = dw.r; na = da.r; nb = db
    
    para.ind = sort(sapply(X.index, function(x) seq(x, nb, nb))) + na
    
    gamma0 = est.para$gamma0; alpha0 = est.para$alpha0; beta0 = est.para$beta0
    sig0 = est.para$sig0; rho0 = est.para$rho0
    eDelZ0.lst = est.para$eDelZ0.lst
    A.R.lst = est.para$A.R.lst; B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst; B2.R.lst = est.para$B2.R.lst; AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst
    eZ0.lst = est.para$eZ0.lst; mu0.lst = est.para$mu0.lst; phi0.lst = est.para$phi0.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, na+nb, na+nb)
    ind.a = 1:na; ind.b = (na+1):(na+nb)
    
    for(i in 1:N){
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel,,drop=FALSE]
      xa.r.tmp = Xa.r.tmp[ID.sel,,drop=FALSE]
      xb.tmp = Xb.tmp[ID.sel,,drop=FALSE]
      n.T = n.reps[i]
      
      pZ0 = pZ0.lst[[i]] # unlist(pZ0.lst[i,])

      mu0 = mu0.lst[[i]]; phi0 = phi0.lst[[i]]
      eZ0 = eZ0.lst[[i]]; eDelZ0 = eDelZ0.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      tmp.a.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDelZ0), length(eDelZ0))
      dvd.a.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDelZ0), length(eDelZ0)) %*% .BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZ0 - mu0), 2, sum))
      
      beta0.x = matrix(0, db, 1)
      beta0.x[1:db.r,] = beta0
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.tmp)
      # dvd.b.x = t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)) %*% t(t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)))
      geeb.x = as.matrix(pseudoll.score.beta)
      
      dvd.ab.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDelZ0), length(eDelZ0)) %*% .BlockDiag(.dmb(xb.tmp, beta0.x, mu0, phi0), 1, n.T)
      
      ES = rbind(geea.x, geeb.x)
      U = U + ES
      
      a = (1/phi0 - 1) * mu0; a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0); b[b < 0] = 0
      A = digamma(a) - digamma(a+b)
      B = digamma(b) - digamma(a+b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1-eDelZ0)*( - A + A.post )
      two.R = (1-eDelZ0)*( - B + B.post )
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
    stat[j] = c(t(U.r[para.ind,,drop=F]) %*% ginv(H.r %*% M.r %*% t(H.r)) %*% U.r[para.ind,,drop=F])
    df[j] = length(para.ind)
  }
  stat.sum = sum(stat); df.sum = sum(df)
  pval = 1 - pchisq(stat.sum, df.sum)
  return(list(stat.sum = stat.sum, pval.sum = pval, stat = stat))
}

# adaptive multiple cores, much faster
.score_test_perm <- function(ID, X, X.index, Y, est.para.lst, stat.sum, test, Rconf,
                                  n.perm, n.cores, R.sel, zi.id = NULL, seed){
  #ID=ID.tmp;X=W.tmp;Y=Y.tmp;stat.sum=stat.zero.sum; test="Zero";Rconf=Rconf.zero;n.perm=1000;n.cores=10
  n.perm = ceiling(n.perm/n.cores)*n.cores # change the maximal permutation times based on number of cores
  m = 0; Nexc = 0; stat.perm = numeric(n.perm)
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
  
  # The parallel package provides a way to reproducibly generate random numbers in a parallel environment via 
  # the “L’Ecuyer-CMRG” random number generator. 
  # Note that this is not the default random number generator so you 
  # will have to set it explicitly.
  while(all(Nexc < R.sel, m < n.perm)){
    tmp.subset = mclapply(1:n.cores, function(x){
      set.seed(seed+m+x)
      X.int.uniq = matrix(NA, N, length(X.index))
      for(i in 1:N){
        tmpX = X[which(ID == id.unique[i]), X.index, drop = FALSE]
        if(dim(unique(tmpX))[1]!=1){
          stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
        }
        X.int.uniq[i,] = unique(tmpX)
      }
      cluster.p = sample(1:N)
      X.p.int = NULL
      for(i in 1:N){
        n.T = n.reps[i]
        X.p.int = rbind(X.p.int, matrix(rep(X.int.uniq[cluster.p[i],], n.T), nrow = n.T, byrow = TRUE))
      }
      X[,X.index] = Rconf %*% X.p.int
      X.p = X
      if(test == "Zero")
        tmp = .score_test_zero(ID, X.p, X.index, Y, est.para.lst, zi.id)
      if(test == "Mean")
        tmp = .score_test_mean(ID, X.p, X.index, Y, est.para.lst)
      if(test == "Disp")
        tmp = .score_test_disp(ID, X.p, X.index, Y, est.para.lst)
      if(test == "Omni"){
        if(length(zi.id) > 0){
          tmp.zero = .score_test_zero(ID, X.p, X.index, Y, est.para.lst, zi.id)
        }else{
          tmp.zero = list()
          tmp.zero$stat.sum = NULL
        }
        tmp.mean = .score_test_mean(ID, X.p, X.index, Y, est.para.lst)
        tmp.disp = .score_test_disp(ID, X.p, X.index, Y, est.para.lst)
        tmp = list()
        tmp$stat.sum = sum(tmp.zero$stat.sum, tmp.mean$stat.sum, tmp.disp$stat.sum)
      }
      return(tmp$stat.sum)
    })
    Nexc = Nexc + sum(unlist(tmp.subset) >= stat.sum)
    stat.perm[m+1:n.cores] = unlist(tmp.subset)
    m = m + n.cores
  }
  
  if(m < n.perm){
    pval = (Nexc+1)/(m+1)
    cat(sprintf("# of permutations: %g\n", m))
    cat("Use ECDF approximation p-value\n")
  }else{
    if(Nexc <= 10){
      pval = tryCatch(.gpd_approx(stat.perm, 250, stat.sum), error=function(err) NA)
      cat(sprintf("# of permutations: %g\n", n.perm))
      if(is.na(pval)){
        pval = (Nexc+1)/(n.perm+1)
        cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
      }else{
        cat("Use GPD approximation p-value\n")
      }
    }else{
      pval = (Nexc+1)/(n.perm+1)
      cat(sprintf("# of permutations: %g\n", n.perm))
      cat("Use ECDF approximation p-value\n")
    }
  }
  return(pval)
}

.simData.GDM <- function(a, b, SeqDepth){
  N = nrow(a)
  K = ncol(a)
  
  # Generate Z matrix
  Z = matrix(rbeta(N * K, a, b), N, K)
  
  # Calculate P matrix
  cumprod_1_minus_Z = matrixStats::rowCumprods(cbind(1, 1 - Z))
  P = Z * cumprod_1_minus_Z[,-(K+1),drop=FALSE]
  P = cbind(P, 1 - rowSums(P))
  P[P < 0] = 0
  
  # Generate Y matrix using apply for rowwise operation
  Y = t(vapply(1:N, function(i) rmultinom(1, SeqDepth[i], P[i,]), numeric(K+1)))
  return(Y)
}


# http://motsinger-reif-lab.org/files/documents/choose-sampling-parameters-for-adaptive-permutation-1.r

# determines the number of test statistics that should 
# be sampled in adaptive permutation 
.choose_r <- function(alpha, c) {
  error <- alpha * c
  R <- 0
  foundR <- FALSE
  while(!foundR) {
    R <- R + 1
    brange <- qnbinom(c(0.1586553, 0.8413447), R, alpha)
    pvalRange <- R / (R + brange)
    diff <- max(abs(pvalRange - alpha))
    if(diff < error) {
      foundR <- TRUE
    }
  }
  return(R)
}

# use Nexc most exterme test statistic to approximate GPD 
.gpd_approx <- function(yperm, nexc, obsts){
  yperm = na.omit(yperm)
  y = sort(yperm, decreasing = TRUE)
  tmp = .gpd_params_est(nexc, y)
  a_hat = tmp[1]
  k_hat = tmp[2]
  t = tmp[3]
  ### goodness-fit test
  nexc_re = .gpd_goft(nexc, y[1:nexc]-t)
  if(nexc_re[2] == nexc){
    z0 = obsts - t
    p = .gpd_pval(a_hat, k_hat, z0, length(yperm), nexc)
    return(p)
  }else{
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
.gpd_params_est <- function(nexc, y){
  t = (y[nexc] + y[nexc+1])/2
  z = y[1:nexc] - t
  z2 = z^2
  m = mean(z)
  m2 = mean(z2)
  a_hat = m*m2/(m2-m^2)/2
  k_hat = (m^2/(m2-m^2)-1)/2
  return(c(a_hat, k_hat, t))
}

# p-value for the generalized Pareto distribution approximation
.gpd_pval <- function(a_hat, k_hat, z0, nperm, nexc){
  p = nexc/nperm*((1-k_hat*z0/a_hat)**(1/k_hat))
  return(p)
}

# iteratively reduce Nexc by 10 until a goodness-of-fit satisfy
.gpd_goft <- function(nexc, y){
  # y: the sorted test statistic - threshold t
  nexc = length(y)
  p = .gp_test(y)
  
  nexc.c = seq(0,nexc-10,10)
  z = y
  i = 0
  
  re = c()
  for(i in nexc.c) {
    z = y[1:(nexc-i)]
    p = .gp_test(z)
    re = rbind(re, c(i, p))
    if(!is.na(p) & p > 0.05) break
    i = i + 10
  }
  
  if(nrow(re) >=2) {
    nexc.c2 = seq(re[(nrow(re)-1),1]+1, re[nrow(re),1],1)
    re = c()
    for(i in nexc.c2) {
      z = y[1:(nexc-i)]
      p = .gp_test(z)
      re = rbind(re, c(i, p))
      if(!is.na(p) & p > 0.05) break
      i = i + 10
    }
  }
  
  p = re[nrow(re),2]
  len = nexc-re[nrow(re),1]
  
  return(c(p, len))
}

# Bootstrap goodness-of-fit test for the Generalized Pareto distribution
# test of fit for the GPD with unknown parameter
# refer to "goft"
.gp_test <- function(x, B = 2999){
  x = x[!is.na(x)]
  x = as.vector(x)
  n = length(x)   #  sample size without NA values
  samplerange = max(x) - min(x)
  gammap = .amle_method(x, k = ceiling(.2 * n))[1]    
  gamman = .combined_method(x)[1]
  r1 = .R1(x)     # observed value of R^-
  r2 = .R2(x)     # observed value of R^+
  # use "combined" for fitting GPD with negative shape parameter
  # use "amle" for fitting GPD with non-negative shape parameter
  p.value1 = sum(replicate(B, .R1(.rgp(n, shape = gamman))) < r1) / B  # bootstrap p-value for H_0^- 
  p.value2 = sum(replicate(B, .R2(.rgp(n, shape = gammap))) < r2) / B  # bootstrap p-value for H_0^+ 
  p.value = max(p.value1, p.value2)    # p-value of the intersection-union test
  return(p.value)
}

# Asymptotic maximum likelihood estimators 
.amle_method <- function(x, k){
  x = sort(x)
  n = length(x)
  nk = n - k
  x1 = x[(nk+1):n]
  w = log(x1)
  g = - (w[1] - sum(w) / k)    
  sigma = g * exp(w[1] + g * log(k / n))
  return(c(g, sigma))
}

# Combined estimators 
.combined_method <- function(x){
  m = mean(x)
  maxi = max(x)
  g = m / (m - maxi) 
  sigma = - g * maxi     
  return(c(g, sigma))
}

# Test statistic for H_0^-
.R1 <- function(x){
  gamma_neg = .combined_method(x)[1]
  Fn = ecdf(x)
  x1 = x[x != max(x)]
  z1 = (1 - Fn(x1))^( - gamma_neg) 
  return(abs(cor(x1, z1)))
}

# Test statistic for H_0^+
.R2  <- function(x){
  n = length(x)
  Fn = ecdf(x)
  gamma_positive = .amle_method(x, ceiling(.2 * n))[1]
  x1 = x[x != max(x)]
  y1 = (1 - Fn(x1))^( - gamma_positive) 
  x.star = log(x1)
  y.star = log( y1 -1 )
  if (gamma_positive <= 0.5)	return(cor(x1, y1))  
  if (gamma_positive > 0.5)  return((cor(x.star, y.star)))
}

# Simulation of random numbers from the gPd
.rgp  <- function (n, shape){
  if (shape != 0) 
    return((1 / shape) * (runif(n)^(-shape) - 1))
  else return(rexp(n, 1))
}

# ACAT function from ACAT package
# https://github.com/yaowuliu/ACAT
# date: 4/16/2024
.ACAT <- function(Pvals,weights=NULL,is.check=TRUE){
  Pvals<-as.matrix(Pvals)
  if (is.check){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
      stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
      stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(colSums(Pvals==0)>=1)
    is.one<-(colSums(Pvals==1)>=1)
    if (sum((is.zero+is.one)==2)>0){
      stop("Cannot have both 0 and 1 p-values in the same column!")
    }
    
    if (sum(is.zero)>0){
      warning("There are p-values that are exactly 0!")
    }
    if (sum(is.one)>0){
      warning("There are p-values that are exactly 1!")
    }
    
  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights)){
    is.weights.null<-TRUE
  }else{
    is.weights.null<-FALSE
    weights<-as.matrix(weights)
    if (sum(dim(weights)!=dim(Pvals))>0){
      stop("The dimensions of weights and Pvals must be the same!")
    }else if (is.check & (sum(weights<0)>0)){
      stop("All the weights must be nonnegative!")
    }else{
      w.sum<-colSums(weights)
      if (sum(w.sum<=0)>0){
        stop("At least one weight should be positive in each column!")
      }else{
        for (j in 1:ncol(weights)){
          weights[,j]<-weights[,j]/w.sum[j]
        }
      }
    }
    
  }
  
  #### check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small<-(Pvals<1e-15)
  if (is.weights.null){
    Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-1/Pvals[is.small]/pi
    cct.stat<-colMeans(Pvals)
  }else{
    Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
    cct.stat<-colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}

