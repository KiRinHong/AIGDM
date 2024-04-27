library(extraDistr)
AIGDM_GEE.Cluster <- function(ID, Y, Xa, Xb, X.index, W, Tax = NULL, min.depth = 0, model = "AZIGDM",
                              n.cores = 1, n.boot = NULL, fdr.alpha = 0.05, perm.type = NULL, n.perm = NULL, seed = 123){
  # ID = id; Y = Y.m; Xa = Xb = W = X.m; Tax = tax; min.depth = 0; model = "GDM"; n.cores = 10; perm.type = "BTW"
  # model could be GDM/ZIGDM/AZIGDM/AIGDM
  remove.subject = which(rowSums(Y) < min.depth)
  if(length(remove.subject) > 0){
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth))
    ID = ID[-remove.subject]
    Xa = Xa[-remove.subject, , drop=FALSE]
    Xb = Xb[-remove.subject, , drop=FALSE]
    W = W[-remove.subject, , drop=FALSE]
    Y = Y[-remove.subject, , drop=FALSE]
  }
  
  keep = which(colSums(Y)>0)
  Y = Y[, keep, drop=FALSE]
  
  Xa = cbind(1, Xa) # add the intercept term
  Xb = cbind(1, Xb) # add the intercept term
  W = cbind(1, W)
  X.index = X.index + 1
  Xa.r = Xa[, -X.index, drop=FALSE]
  Xb.r = Xb[, -X.index, drop=FALSE]
  W.r = W[, -X.index, drop=FALSE]
  
  if(is.null(Tax)){
    colOrder = order(colMeans(Y/rowSums(Y), na.rm = T), decreasing = TRUE)
    Y = Y[,colOrder]
    K = ncol(Y) - 1
    R.sel = .choose_r(fdr.alpha/K, 0.05)
    
    if(model == "GDM"){
      zi.check.pval.all.wo = oi.check.pval.all.wo = rep(1, K)
    }else if(model == "ZIGDM"){
      zi.check.pval.all.wo = rep(0, K); oi.check.pval.all.wo = rep(1, K)
    }else{
      print("Zero-inflation diagnostic test!")
      # TRUE means add zero-inflation
      zi.check.pval.all.wo = sapply(1:K, function(x){
        y.tmp = cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE]))
        y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
        return(.check_zeroinflation_fork(ID, y.tmp, n.boot = n.boot, n.cores = n.cores, R.sel))
      })
      if(model == "AIGDM"){
        R.sel = .choose_r(fdr.alpha/K/2, 0.05)
        print("One-inflation diagnostic test!")
        # TRUE means add one-inflation
        oi.check.pval.all.wo = sapply(1:K, function(x){
          y.tmp = cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE]))
          y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
          return(.check_oneinflation(ID, y.tmp, n.boot = n.boot, n.cores = n.cores, R.sel))
        })
      }else{
        oi.check.pval.all.wo = rep(1, K)
      }
    }
    
    if(model == "AZIGDM"){
      zi.check.pval.all.w = p.adjust(zi.check.pval.all.wo, method = "BH")
      oi.check.pval.all.w = rep(1, length(oi.check.pval.all.wo))
    }else if(model == "AIGDM"){
      check.pval.all = p.adjust(c(zi.check.pval.all.wo, oi.check.pval.all.wo), method = "BH")
      zi.check.pval.all.w = check.pval.all[1:length(zi.check.pval.all.wo)]
      oi.check.pval.all.w = check.pval.all[length(zi.check.pval.all.wo) + 1:length(oi.check.pval.all.wo)]
    }else if(model == "ZIGDM"){
      zi.check.pval.all.w = rep(0, length(zi.check.pval.all.wo))
      oi.check.pval.all.w = rep(1, length(oi.check.pval.all.wo))
    }else if(model == "GDM"){
      zi.check.pval.all.w = rep(1, length(zi.check.pval.all.wo))
      oi.check.pval.all.w = rep(1, length(oi.check.pval.all.wo))
    }
    
    zi.check = zi.check.pval.all.w < fdr.alpha
    oi.check = oi.check.pval.all.w < fdr.alpha
    
    keep.ind = which(rowSums(Y) != 0)
    Y.tmp = Y[keep.ind,,drop=FALSE]
    ID.tmp = ID[keep.ind]
    W.r.tmp = W.r[keep.ind,,drop=FALSE]
    Xa.r.tmp = Xa.r[keep.ind,,drop=FALSE]; Xa.tmp = Xa[keep.ind,,drop=FALSE]
    Xb.r.tmp = Xb.r[keep.ind,,drop=FALSE]; Xb.tmp = Xb[keep.ind,,drop=FALSE]
    colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
    Y.tmp = Y.tmp[,colOrder]
    
    est.para.lst = lapply(1:K, function(x){
      .est_para(ID.tmp, W.r.tmp, Xa.r.tmp, Xb.r.tmp, cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])),
                zi.check[x], oi.check[x])
    })
    
    R.sel = .choose_r(fdr.alpha, 0.05)
    asym.mean = .score_test_mean(ID.tmp, Xa.tmp, X.index, Y.tmp, est.para.lst, zi.check, oi.check)
    stat.mean.sum = asym.mean$stat.sum; stat.mean.taxonwise = asym.mean$stat
    pval.mean.a = asym.mean$pval
    pval.mean.p = .score_test_perm(ID.tmp, Xa.tmp, X.index, Y.tmp, est.para.lst, stat.mean.sum, test = "mean",
                                   n.perm, perm.type, zi.check, oi.check, n.cores, R.sel, seed)
    
    asym.disp = .score_test_disp(ID.tmp, Xb.tmp, X.index, Y.tmp, est.para.lst, zi.check, oi.check)
    stat.disp.sum = asym.disp$stat.sum; stat.disp.taxonwise = asym.disp$stat
    pval.disp.a = asym.disp$pval
    pval.disp.p = .score_test_perm(ID.tmp, Xb.tmp, X.index, Y.tmp, est.para.lst, stat.disp.sum, test = "disp",
                                   n.perm, perm.type, zi.check, oi.check, n.cores, R.sel, seed)
    
    pval = c(pval.mean.a, pval.mean.p, pval.disp.a, pval.disp.p)
    names(pval) = c("Mean-Asymptotic", "Mean-Resampling", "Disp-Asymptotic", "Disp-Resampling")
    names(zi.check) = names(oi.check) = names(zi.check.pval.all.wo) = names(oi.check.pval.all.wo) = colnames(Y.tmp)
    
    return(list(pval = pval, zi = zi.check, oi = oi.check,
                zi.pval = zi.check.pval.all.wo, oi.pval = oi.check.pval.all.wo,
                stat.mean.taxonwise = stat.mean.taxonwise,
                stat.disp.taxonwise = stat.disp.taxonwise))
  }else{
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
            subtree = c(subtree, paste(Rank.low, level.uni[m], sep = "."))
          }
        }
      }
    }
    R.sel = .choose_r(fdr.alpha/length(subtree), 0.05)
    
    zi.check.pval.lst = oi.check.pval.lst = list()
    curr.ind = 1
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
            colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
            # Y.tmp = Y.tmp[,c(colOrder[-1], colOrder[1])]
            Y.tmp = Y.tmp[,colOrder]
            K = ncol(Y.tmp) - 1

            if(model == "GDM"){
              zi.check.pval.lst[[curr.ind]] = oi.check.pval.lst[[curr.ind]] = rep(1, K)
            }else if(model == "ZIGDM"){
              zi.check.pval.lst[[curr.ind]] = rep(0, K)
              oi.check.pval.lst[[curr.ind]] = rep(1, K)
            }else{
              print(paste0("### Diagnostic Test ### Now processing Rank ", n.rank-k, ".", m))
              print("Zero-inflation diagnostic test!")
              # TRUE means add zero-inflation
              zi.check.pval = sapply(1:K, function(x){
                y.tmp = cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE]))
                id.tmp = ID[rowSums(y.tmp) != 0]
                y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
                return(.check_zeroinflation_fork(id.tmp, y.tmp, n.boot = n.boot, n.cores = n.cores, R.sel))
              })
              zi.check.pval.lst[[curr.ind]] = zi.check.pval
              if(model == "AIGDM"){
                print("One-inflation diagnostic test!")
                # TRUE means add one-inflation
                oi.check.pval = sapply(1:K, function(x){
                  y.tmp = cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE]))
                  id.tmp = ID[rowSums(y.tmp) != 0]
                  y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
                  return(.check_oneinflation(id.tmp, y.tmp, n.boot = n.boot, n.cores = n.cores, R.sel))
                })
                oi.check.pval.lst[[curr.ind]] = oi.check.pval
              }else{
                oi.check.pval.lst[[curr.ind]] = rep(1, K)
              }
            }
            curr.ind = curr.ind + 1
          }
        }
      }
    }
    
    zi.check.pval.all.wo = unlist(zi.check.pval.lst)
    oi.check.pval.all.wo = unlist(oi.check.pval.lst)
    if(model == "AZIGDM"){
      zi.check.pval.all.w = p.adjust(zi.check.pval.all.wo, method = "BH") # unlist(lapply(zi.check.pval.lst, p.adjust, method = "BH")) 
      oi.check.pval.all.w = rep(1, length(oi.check.pval.all.wo))
    }else if(model == "AIGDM"){
      check.pval.all = p.adjust(c(zi.check.pval.all.wo, oi.check.pval.all.wo), method = "BH") # need to be modified
      zi.check.pval.all.w = check.pval.all[1:length(zi.check.pval.all.wo)]
      oi.check.pval.all.w = check.pval.all[length(zi.check.pval.all.wo) + 1:length(oi.check.pval.all.wo)]
    }else if(model == "ZIGDM"){
      zi.check.pval.all.w = rep(0, length(zi.check.pval.all.wo))
      oi.check.pval.all.w = rep(1, length(oi.check.pval.all.wo))
    }else if(model == "GDM"){
      zi.check.pval.all.w = rep(1, length(zi.check.pval.all.wo))
      oi.check.pval.all.w = rep(1, length(oi.check.pval.all.wo))
    }
    
    zi.check.lst = oi.check.lst = stat.mean.taxonwise.lst = stat.disp.taxonwise.lst = list()
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
            print(paste0("### Hypothesis Test ### Now processing Rank ", n.rank-k, ".", m))
            keep.ind = which(rowSums(Y.tmp) != 0)
            
            Y.tmp = Y.tmp[keep.ind,,drop=FALSE]
            ID.tmp = ID[keep.ind]
            W.r.tmp = W.r[keep.ind,,drop=FALSE]
            Xa.r.tmp = Xa.r[keep.ind,,drop=FALSE]; Xa.tmp = Xa[keep.ind,,drop=FALSE]
            Xb.r.tmp = Xb.r[keep.ind,,drop=FALSE]; Xb.tmp = Xb[keep.ind,,drop=FALSE]
            
            colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
            # Y.tmp = Y.tmp[,c(colOrder[-1], colOrder[1])]
            Y.tmp = Y.tmp[,colOrder]
            K = ncol(Y.tmp) - 1
            
            zi.check = zi.check.pval.all.w[curr.ind+1:K] < 0.05 # fdr.alpha
            print(zi.check)
            oi.check = oi.check.pval.all.w[curr.ind+1:K] < 0.05 # fdr.alpha
            curr.ind = curr.ind+K
            
            est.para.lst = lapply(1:K, function(x){
              .est_para(ID.tmp, W.r.tmp, Xa.r.tmp, Xb.r.tmp, cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])),
                        zi.check[x], oi.check[x])
            })

            # difference between current (L-BFGS) and previous (BFGS)
            # est.para = do.call(Map, c(f = cbind, est.para.lst))
            # 
            # aigdm.reg = .AIGDM_EM(Y.tmp, W.r.tmp, Xa.r.tmp, Xb.r.tmp, 
            #                       matrix(0.001, 2, 2), matrix(-Inf, 2, 2), 
            #                       matrix(0.001, 2, 2), matrix(0.001, 2, 2))
            # aigdm.reg$gamma.est
            # 
            # est.tmp = .ZIGDM.EM.PAR2(Y.tmp, W.r.tmp, Xa.r.tmp, Xb.r.tmp,
            #                matrix(0.001, 2, 2), matrix(0.001, 2, 2), 
            #                matrix(0.001, 2, 2), 1e-4, 5000)
            # t(est.tmp$c.est)
            
            asym.mean = .score_test_mean(ID.tmp, Xa.tmp, X.index, Y.tmp, est.para.lst, zi.check, oi.check)
            stat.mean.sum = asym.mean$stat.sum; stat.mean.taxonwise = asym.mean$stat
            pval.mean.a = asym.mean$pval
            
            asym.disp = .score_test_disp(ID.tmp, Xb.tmp, X.index, Y.tmp, est.para.lst, zi.check, oi.check)
            stat.disp.sum = asym.disp$stat.sum; stat.disp.taxonwise = asym.disp$stat
            pval.disp.a = asym.disp$pval
            
            
            if(is.null(n.perm)){
              pval = cbind(pval, c(pval.mean.a, pval.disp.a))
            }else{
              # compute residual forming matrix Rconf (Smith's method)
              if(ncol(Xa.r.tmp) == 1){
                Rconf.mean = diag(nrow(Xa.r.tmp))
              }else{
                Rconf.mean = diag(nrow(Xa.r.tmp)) - Xa.r.tmp[,-1,drop=F] %*% solve(t(Xa.r.tmp[,-1,drop=F]) %*% Xa.r.tmp[,-1,drop=F]) %*% t(Xa.r.tmp[,-1,drop=F])
              }
              pval.mean.p = .score_test_perm(ID.tmp, Xa.tmp, X.index, Y.tmp, est.para.lst, stat.mean.sum, test = "mean", Rconf.mean,
                                             n.perm, perm.type, zi.check, oi.check, n.cores, R.sel, seed)
              
              # compute residual forming matrix Rconf (Smith's method)
              if(ncol(Xb.r.tmp) == 1){
                Rconf.disp = diag(nrow(Xb.r.tmp))
              }else{
                Rconf.disp = diag(nrow(Xb.r.tmp)) - Xb.r.tmp[,-1,drop=F] %*% solve(t(Xb.r.tmp[,-1,drop=F]) %*% Xb.r.tmp[,-1,drop=F]) %*% t(Xb.r.tmp[,-1,drop=F])
              }
              pval.disp.p = .score_test_perm(ID.tmp, Xb.tmp, X.index, Y.tmp, est.para.lst, stat.disp.sum, test = "disp", Rconf.disp,
                                             n.perm, perm.type, zi.check, oi.check, n.cores, R.sel, seed)
              pval = cbind(pval, c(pval.mean.a, pval.mean.p, pval.disp.a, pval.disp.p))
            }
            zi.check.lst[[ncol(pval)]] = zi.check
            oi.check.lst[[ncol(pval)]] = oi.check
            stat.mean.taxonwise.lst[[ncol(pval)]] = stat.mean.taxonwise
            stat.disp.taxonwise.lst[[ncol(pval)]] = stat.disp.taxonwise
            
          }
        }
      }# lineage loop
    }# level loop
    colnames(pval) = subtree
    names(zi.check.lst) = names(oi.check.lst) = 
      names(zi.check.pval.lst) = names(oi.check.pval.lst) = 
      names(stat.mean.taxonwise.lst) = names(stat.disp.taxonwise.lst) = subtree
    # identify significant lineages
    rslt = apply(pval, 1, .identifySigLineagesNsimesTest, fdr.alpha)
    
    if(is.null(n.perm)){
      names(rslt) = c("Mean-Asymptotic", "Disp-Asymptotic")
    }else{
      names(rslt) = c("Mean-Asymptotic", "Mean-Resampling", "Disp-Asymptotic", "Disp-Resampling")
    }
    rslt = c(rslt, lineage.zi = list(zi.check.lst), lineage.oi = list(oi.check.lst),
             lineage.zi.pval = list(zi.check.pval.lst), lineage.oi.pval = list(oi.check.pval.lst),
             stat.mean.taxonwise = list(stat.mean.taxonwise.lst),
             stat.disp.taxonwise = list(stat.disp.taxonwise.lst))
    return(rslt)
  }
  
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

# # .eZ: expectation of Z in the GD
# .eZ <- function(av, bv, Y){
#   K = length(Y)-1
#   N = sum(Y)
#   
#   av.prim = av + Y[1:K]  
#   bv.prim = bv + (N - cumsum(Y[1:K]))
#   return(list(av.post = av.prim, bv.post = bv.prim))
# }

# .ZOIeZ <- function(pZv, pOv, av, bv, Y){
#   K = length(Y)-1
#   N = sum(Y)
#   
#   pv.post = pZv.post = pOv.post = rep(0, K)
#   av.prim = av + Y[1:K]
#   bv.prim = bv
#   if(length(which(Y>0)) == 0){
#     for(j in 1:K){
#       bv.prim[j] = bv.prim[j] + (N - sum(Y[1:j]))
#     }
#   }else{
#     target.id = max(which(Y>0))
#     for(j in 1:K){
#       bv.prim[j] = bv.prim[j] + (N - sum(Y[1:j]))
#       if(j==target.id){
#         if(beta(av[j], bv[j])==0){
#           pOv.post[j] = pOv[j]/(pOv[j] + (1-pOv[j]))
#         }else{
#           tmp = pOv[j] + (1-pOv[j])*beta(av.prim[j], bv.prim[j])/beta(av[j], bv[j])
#           if(tmp==0){
#             pOv.post[j] = 1
#           }else{
#             pOv.post[j] = pOv[j]/tmp
#           }
#         }
#       }
#     }
#   } 
#   for(j in 1:K){
#     if(Y[j]==0){
#       if(beta(av[j], bv[j])==0){
#         pZv.post[j] = pZv[j]/(pZv[j] + (1-pZv[j]))
#       }else{
#         tmp = pZv[j] + (1-pZv[j])*beta(av.prim[j], bv.prim[j])/beta(av[j], bv[j])
#         if(any(tmp==0, is.nan(tmp))){
#           pZv.post[j] = 1
#         }else{
#           pZv.post[j] = pZv[j]/tmp
#         }
#       }
#     }
#   }
#   pv.post = pZv.post + pOv.post
#   return(list(pv.post = pv.post, pOv.post = pOv.post, pZv.post = pZv.post, av.post = av.prim, bv.post = bv.prim))
# }

.ZOIeZ_vec <- function(pZv, pOv, av, bv, Y){
  K = length(Y)-1
  N = sum(Y)
  
  # Vectorized operations for faster computations
  av.prim = av + Y[1:K]
  bv.prim = bv + (N - cumsum(Y[1:K]))
  pOv.post = pZv.post = rep(0, K)
  # # Index where Y > 0
  # ids_Y_positive = which(Y > 0)
  # if(length(ids_Y_positive) > 0){
  #   target.id = max(ids_Y_positive)
  #   
  #   for(j in 1:K){
  #     if(j == target.id){
  #       beta_ratio = ifelse(beta(av[j], bv[j]) == 0, 
  #                           1,
  #                           beta(av.prim[j], bv.prim[j]) / beta(av[j], bv[j]))
  #       tmp = pOv[j] + (1 - pOv[j]) * beta_ratio
  #       pOv.post[j] = ifelse(tmp == 0, 1, pOv[j] / tmp)
  #     }
  #   }
  # }
  
  # Compute pZv.post values in a vectorized way
  beta_values = ifelse(beta(av, bv) == 0, 
                       1,
                       beta(av.prim, bv.prim) / beta(av, bv))
  tmp_values = pZv + (1 - pZv) * beta_values
  pZv.post[Y[1:K] == 0] = ifelse(tmp_values[Y[1:K] == 0] == 0, 
                                 1, 
                                 pZv[Y[1:K] == 0] / tmp_values[Y[1:K] == 0])
  
  # Finalize pv.post
  pv.post = pZv.post + pOv.post
  return(list(pv.post = pv.post, pOv.post = pOv.post, pZv.post = pZv.post, av.post = av.prim, bv.post = bv.prim))
}
.ZOIeZ_mat <- function(pZv, pOv, av, bv, Y){
  n = dim(Y)[1]
  K = dim(Y)[2] - 1
  
  # Vectorized operations for faster computations
  av.prim = av + Y[,1:K, drop=FALSE]
  bv.prim = bv + (rowSums(Y) - rowCumsums(Y[,1:K,drop=FALSE]))
  pOv.post = pZv.post = matrix(0, n, K)
  
    # # Index where Y > 0
    # ids_Y_positive = which(Y > 0)
    # if(length(ids_Y_positive) > 0){
    #   target.id = max(ids_Y_positive)
    #   
    #   for(j in 1:K){
    #     if(j == target.id){
    #       beta_ratio = ifelse(beta(av[j], bv[j]) == 0, 
    #                           1,
    #                           beta(av.prim[j], bv.prim[j]) / beta(av[j], bv[j]))
    #       tmp = pOv[j] + (1 - pOv[j]) * beta_ratio
    #       pOv.post[j] = ifelse(tmp == 0, 1, pOv[j] / tmp)
    #     }
    #   }
    # }
  
  # Compute pZv.post values in a vectorized way
  beta_values = ifelse(beta(av, bv) == 0, 1, beta(av.prim, bv.prim) / beta(av, bv))
  
  tmp_values = pZv + (1 - pZv) * beta_values
  zero_indices = Y[, 1:K, drop=FALSE] == 0
  pZv.post[zero_indices] = ifelse(tmp_values[zero_indices] == 0, 1, pZv[zero_indices] / tmp_values[zero_indices])
  
  # Finalize pv.post
  pv.post = pZv.post + pOv.post
  return(list(pv.post = pv.post, pOv.post = pOv.post, pZv.post = pZv.post, av.post = av.prim, bv.post = bv.prim))
}

.LogitOptim <- function(Del, W, gamma.ini){
  # Del=DelZ.R[,1]; gamma.ini=gamma.last[,1]
  Logit.par.ini = gamma.ini
  Logit.data = list(Del=Del, W=W)
  
  return(optim(par=Logit.par.ini, fn=.LogitNegLoglik, gr=.LogitNegScore, data = Logit.data, method="L-BFGS-B")$par) 
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
  # Del=Del.R[,1];A=A.R[,1];B=B.R[,1];alpha.ini=alpha.last[,1];beta.ini=beta.last[,1]
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
  
  # zstar = data$A - data$B
  # mustar = digamma(a) - digamma(b)
  # 
  # return(-  c(colSums((1-data$Del)*(phi * mu.tmp * (1/(1+tmp)) * (zstar - mustar) * data$X)),
  #             sum((1-data$Del)*(mu.tmp * (zstar - mustar) + data$B - digamma(b) + digamma(phi))) ) )
}
# # score function, obtained by differentiating the log-likelihood function w.r.t unknown parameters
# .BetaNegScore <- function(par, data){
#   # par = Beta.par.ini; data = Beta.data
#   da = ncol(data$Xa)
#   alpha = par[1:da]
#   beta = par[-(1:da)]
#   
#   tmp.a = as.numeric(exp(data$Xa %*% alpha))
#   mu.tmp = as.numeric( tmp.a/(1+tmp.a) )
#   
#   tmp.b = as.numeric( exp( data$Xb %*% beta) )
#   sigma.tmp = as.numeric( tmp.b/(1+tmp.b) )
#   
#   a = (1/sigma.tmp - 1) * mu.tmp
#   b = (1/sigma.tmp - 1) * (1 - mu.tmp)
#   a[a < 0] = 0
#   b[b < 0] = 0
#   
#   # a.a = (1/tmp.b) * mu.tmp * (1/(1+tmp.a))
#   # a.b = -(1/tmp.b) * mu.tmp
#   # b.a = - a.a
#   # b.b = -(1/tmp.b) *  (1/(1+tmp.a))
#   # 
#   # one = digamma(a+b)-digamma(a) + data$A
#   # two = digamma(a+b)-digamma(b) + data$B
#   # 
#   # -c( colSums( (one*a.a + two*b.a) * data$Xa ), colSums( (one*a.b + two*b.b) * data$Xb ) ) 
#   zstar = data$A - data$B
#   zdagger = data$B
#   mustar = digamma(a) - digamma(b)
#   mudagger = digamma(b) - digamma(a + b)
#   
#   return(-c(colSums( 1/tmp.b * mu.tmp * (1/(1+tmp.a)) * (zstar - mustar) * data$Xa),
#             colSums( -1/tmp.b * (mu.tmp * (zstar - mustar) + (zdagger - mudagger)) * data$Xb)))
#   
#   # zstar = data$A - data$B
#   # mustar = digamma(a) - digamma(b)
#   # 
#   # return(-c(colSums(phi * mu.tmp * (1/(1+tmp)) * (zstar - mustar) * data$X),
#   #           sum(mu.tmp * (zstar - mustar) + data$B - digamma(b) + digamma(phi))))
# }


# # maximum likelihood estimators of gamma, alpha and phi
# .AIGDM_EM <- function(Y, W, Xa, Xb, gamma0, omega0, alpha0, beta0, tol = 1e-4, max.iter = 1000){
#   # Y = Y.m; W = wmat.m; X = xmat.m; 
#   # alpha0 <- matrix(0.001, nrow = ncol(xmat.m), ncol = K)
#   # CONV = 0
#   # CONV.iter = max.iter
#   n = nrow(Y)
#   K = ncol(Y) - 1
#   da = ncol(Xa); db = ncol(Xb)
# 
#   gamma.last = gamma.now = gamma0; omega.last = omega.now = omega0; alpha.last = alpha.now = alpha0; beta.last = beta.now = beta0
#   DelAB.R = matrix(NA, n, 5*K)
#   
#   for (l in 1:max.iter) {
#     # E-step
#     # print(paste("====== ", l, "th ======", sep=""))
#     for(i in 1:n){
#       tmp = exp(Xa[i,] %*% alpha.last)
#       mv = as.numeric(tmp/(1+tmp))
#       
#       tmp = exp(Xb[i,] %*% beta.last)
#       sv = as.numeric(tmp/(1+tmp))
#       
#       av = mv * (1/sv - 1)
#       bv = (1 - mv) * (1/sv - 1)
#       
#       tmp = exp(W[i,] %*% gamma.last)
#       tmp[is.na(tmp)] = 0
#       pZv = as.numeric( tmp/(1+tmp) )
#       pZv[is.infinite(tmp) & tmp>0] =1  # positive inf
#       
#       tmp = exp(W[i,] %*% omega.last)
#       tmp[is.na(tmp)] = 0
#       pOv = as.numeric( tmp/(1+tmp) )
#       pOv[is.infinite(tmp) & tmp>0] =1  # positive inf
#       par.post = .ZOIeZ(pZv, pOv, av, bv, Y[i,])
#       tmp = par.post$av.post + par.post$bv.post
#       DelAB.R[i,] = c(par.post$pv.post, par.post$pZv.post, par.post$pOv.post,
#                       digamma(par.post$av.post) - digamma(tmp),
#                       digamma(par.post$bv.post) - digamma(tmp))
#       Del.R = DelAB.R[, 1:K, drop = FALSE]
#       DelZ.R = DelAB.R[, K+1:K, drop = FALSE]
#       DelO.R = DelAB.R[, 2*K+1:K, drop = FALSE]
#       A.R = DelAB.R[, 3*K+1:K, drop = FALSE]
#       B.R = DelAB.R[, 4*K+1:K, drop = FALSE]
#     }
#     # M-step
#     for (j in 1:K) {
#       if( !is.infinite(gamma.last[1,j]) ){
#         gamma.now[,j] = .LogitOptim(DelZ.R[,j], W, gamma.last[,j])
#       }
#       if( !is.infinite(omega.last[1,j]) ){
#         omega.now[,j] = .LogitOptim(DelO.R[,j], W, omega.last[,j])
#       }
#       tmp = .AIBetaOptim(Del.R[,j], A.R[,j], B.R[,j], Xa, Xb, alpha.last[,j], beta.last[,j])
#       alpha.now[,j] = tmp[1:da]; beta.now[,j] = tmp[-(1:da)]
#     }
#     
#     diff = 0
#     if(sum(!is.infinite(gamma.now)) > 0){
#       diff = diff + sum(abs(gamma.now[!is.infinite(gamma.now)] - gamma.last[!is.infinite(gamma.now)])) 
#     }
#     if(sum(!is.infinite(omega.now)) > 0){
#       diff = diff + sum(abs(omega.now[!is.infinite(omega.now)] - omega.last[!is.infinite(omega.now)])) 
#     }
#     diff = diff + sum(abs(alpha.now - alpha.last))
#     diff = diff + sum(abs(beta.now - beta.last))
#     # print(diff)
#     if(diff < tol){
#       # CONV = 1; CONV.iter = l
#       break
#     }else{
#       gamma.last = gamma.now; omega.last = omega.now; alpha.last = alpha.now; beta.last = beta.now
#     }
#   }
#   return(list(gamma.est = gamma.now, omega.est = omega.now, alpha.est = alpha.now, beta.est = beta.now))
# }

# maximum likelihood estimators of gamma, alpha and phi
.AIGDM_EM <- function(Y, W, Xa, Xb, gamma0, omega0, alpha0, beta0, tol = 1e-4, max.iter = 1000){
  n = nrow(Y)
  K = ncol(Y) - 1
  da = ncol(Xa); db = ncol(Xb)
  
  gamma.last = gamma.now = gamma0; omega.last = omega.now = omega0; alpha.last = alpha.now = alpha0; beta.last = beta.now = beta0
  Del.R = DelZ.R = DelO.R = A.R = B.R = matrix(0, n, K)
  for (l in 1:max.iter) {
    # E-step
    # print(paste("====== ", l, "th ======", sep=""))
    expXa_alpha = exp(Xa %*% alpha.last)
    expXb_beta = exp(Xb %*% beta.last)
    expW_gamma = exp(W %*% gamma.last)
    expW_omega = exp(W %*% omega.last)
    # Compute tmpMv, tmpSv using vectorized operations
    tmpMv = expXa_alpha / (1 + expXa_alpha)
    tmpSv = expXb_beta / (1 + expXb_beta)
    # Compute tmpA and tmpB using vectorized operations
    tmpA = tmpMv * (1 / tmpSv - 1)
    tmpB = (1 - tmpMv) * (1 / tmpSv - 1)
    # Compute tmpZ and tmpO using vectorized operations
    tmpZ = expW_gamma / (1 + expW_gamma)
    tmpO = expW_omega / (1 + expW_omega)
    # Fix infinite values
    tmpZ[is.na(tmpZ)] = 0
    tmpO[is.na(tmpO)] = 0
    tmpZ[is.infinite(tmpZ) & tmpZ > 0] = 1
    tmpO[is.infinite(tmpO) & tmpO > 0] = 1
    
    par.post = .ZOIeZ_mat(tmpZ, tmpO, tmpA, tmpB, Y)
    tmp = par.post$av.post + par.post$bv.post
    
    Del.R = par.post$pv.post
    DelZ.R = par.post$pZv.post
    DelO.R = par.post$pOv.post
    A.R = digamma(par.post$av.post) - digamma(tmp)
    B.R = digamma(par.post$bv.post) - digamma(tmp)
    
    # M-step
    for (j in 1:K) {
      if( !is.infinite(gamma.last[1,j]) ){
        gamma.now[,j] = .LogitOptim(DelZ.R[,j], W, gamma.last[,j])
      }
      if( !is.infinite(omega.last[1,j]) ){
        omega.now[,j] = .LogitOptim(DelO.R[,j], W, omega.last[,j])
      }
      tmp = .AIBetaOptim(Del.R[,j], A.R[,j], B.R[,j], Xa, Xb, alpha.last[,j], beta.last[,j])
      alpha.now[,j] = tmp[1:da]; beta.now[,j] = tmp[-(1:da)]
    }
    
    diffs = abs(c(gamma.now - gamma.last, omega.now - omega.last, alpha.now - alpha.last, beta.now - beta.last))
    if (sum(diffs[!is.infinite(diffs)], na.rm=TRUE) < tol) break
    gamma.last = gamma.now; omega.last = omega.now; alpha.last = alpha.now; beta.last = beta.now
  }
  return(list(gamma.est = gamma.now, omega.est = omega.now, alpha.est = alpha.now, beta.est = beta.now))
}
# .dzigdirmn <- function(Y, ZIBlogl, Del, p, alpha, beta) {
#   # ZIBlogl = tmp.z; Del = Del.R; alpha = av; beta = bv
#   logl1 = Del * log(p) + (1-Del) * log(1-p)
#   index = which(p==0 | p==1)
#   if(length(index)>0){
#     logl1[index] = 0
#   }
#   
#   m = rowSums(Y)
#   d = ncol(Y)
#   z = t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
#   logl2 = (lgamma(m + 1) + rowSums(lgamma(Y[, -d] + alpha)) +
#             rowSums(lgamma(z[, -1] + beta)) + rowSums(lgamma(alpha + beta))) -
#     (rowSums(lgamma(Y + 1)) + rowSums(lgamma(alpha)) + rowSums(lgamma(beta)) +
#        rowSums(lgamma(alpha + beta + z[, -d])))
#   print(logl2)
#   logl = sum(logl1, na.rm = T) + sum(ZIBlogl, na.rm = T) + sum(logl2, na.rm = T)
#   return(logl)
# }

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
  # eDel = eDel0; ez = eZ0; mu = mu0; phi = phi0
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

.DVD <- function(wmat, xamat, xbmat, gamma, omega, alpha, beta, pZ, pO, mu, phi, one, two, K, n.T) {
  # wmat=w.tmp; xmat=x.tmp; gamma=gamma0.w; alpha=alpha0.x; p=p0; mu=mu0; phi=phi0
  # tmp.g = c(exp(wmat %*% gamma)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  # tmp1.g = (1 + tmp.g)^2
  # wpp = 1/(pZ*(1-pZ)) * (tmp.g/tmp1.g)^2
  wpp = -pZ*(1-pZ)
  wpp[is.infinite(wpp) | is.nan(wpp)] = 0
  wmat_expand = do.call(rbind, replicate(K, wmat, simplify = FALSE))
  KppZ = t(.BlockDiag(wmat_expand, K, n.T)) %*% diag(wpp, length(wpp)) %*% (.BlockDiag(wmat_expand, K, n.T))
  
  # tmp.o = c(exp(wmat %*% omega)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  # tmp1.o = (1 + tmp.o)^2
  # wpp = 1/(pO*(1-pO)) * (tmp.o/tmp1.o)^2
  wpp = -pO*(1-pO)
  wpp[is.infinite(wpp) | is.nan(wpp)] = 0
  KppO = t(.BlockDiag(wmat_expand, K, n.T)) %*% diag(wpp, length(wpp)) %*% (.BlockDiag(wmat_expand, K, n.T))
  
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
  a.a2 = b.tmp * a3.tmp * (1-tmp.a)  # 10/12/2017 correct
  a.ab = - b.tmp * a2.tmp
  a.b2 = b.tmp * a1.tmp
  b.a2 = - b.tmp * a3.tmp * (1-tmp.a) # 10/12/2017 correct
  b.ab = b.tmp * a2.tmp
  b.b2 = b.tmp * a0.tmp
  
  waa = - (tmp.a/tmp1.a)^2 * sigma^2 * a * (1-pZ) * (1-pO) + one * a.a2 + two * b.a2
  xamat_expand = do.call(rbind, replicate(K, xamat, simplify = FALSE))
  Kaa = t(.BlockDiag(xamat_expand, K, n.T)) %*% diag(waa, length(waa)) %*% (.BlockDiag(xamat_expand, K, n.T))
  
  wab = - tmp.a/tmp1.a * (-1/tmp.b) * sigma * (mu * a - psi2) * (1-pZ) * (1-pO) + one * a.ab + two * b.ab 
  Kab = t(xamat) %*% diag(wab, length(wab)) %*% xbmat
  
  wbb = - (1/tmp.b)^2 * b * (1-pZ) * (1-pO) + one * a.b2 + two * b.b2
  Kbb = t(xbmat) %*% diag(wbb, length(wbb)) %*% xbmat
  
  return(list(Kab = Kab, Kbb = Kbb))
}

## analytically Hessian (expected information) or covariance matrix (inverse of Hessian)
.HessFun.AIBeta.Wald <- function(wmat, xamat, xbmat, gamma, omega, alpha, beta, pZ, pO, mu, phi, one, two, K, n.T) {
  # wmat=w.tmp; xmat=x.tmp; gamma=gamma0.w; alpha=alpha0.x; p=p0; mu=mu0; phi=phi0
  # tmp.g = c(exp(wmat %*% gamma)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  # tmp1.g = (1 + tmp.g)^2
  # wpp = 1/(pZ*(1-pZ)) * (tmp.g/tmp1.g)^2
  wpp = -pZ*(1-pZ)
  wpp[is.infinite(wpp) | is.nan(wpp)] = 0
  wmat_expand = do.call(rbind, replicate(K, wmat, simplify = FALSE))
  KppZ = t(.BlockDiag(wmat_expand, K, n.T)) %*% diag(wpp, length(wpp)) %*% (.BlockDiag(wmat_expand, K, n.T))
  
  # tmp.o = c(exp(wmat %*% omega)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  # tmp1.o = (1 + tmp.o)^2
  # wpp = 1/(pO*(1-pO)) * (tmp.o/tmp1.o)^2
  wpp = -pO*(1-pO)
  wpp[is.infinite(wpp) | is.nan(wpp)] = 0
  KppO = t(.BlockDiag(wmat_expand, K, n.T)) %*% diag(wpp, length(wpp)) %*% (.BlockDiag(wmat_expand, K, n.T))
  
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
  a.a2 = b.tmp * a3.tmp * (1-tmp.a)  # 10/12/2017 correct
  a.ab = - b.tmp * a2.tmp
  a.b2 = b.tmp * a1.tmp
  b.a2 = - b.tmp * a3.tmp * (1-tmp.a) # 10/12/2017 correct
  b.ab = b.tmp * a2.tmp
  b.b2 = b.tmp * a0.tmp
  
  waa = - (tmp.a/tmp1.a)^2 * sigma^2 * a * (1-pZ) * (1-pO) + one * a.a2 + two * b.a2
  xamat_expand = do.call(rbind, replicate(K, xamat, simplify = FALSE))
  Kaa = t(.BlockDiag(xamat_expand, K, n.T)) %*% diag(waa, length(waa)) %*% (.BlockDiag(xamat_expand, K, n.T))
  
  wab = - tmp.a/tmp1.a * (-1/tmp.b) * sigma * (mu * a - psi2) * (1-pZ) * (1-pO) + one * a.ab + two * b.ab 
  Kab = t(xamat) %*% diag(wab, length(wab)) %*% xbmat
  
  wbb = - (1/tmp.b)^2 * b * (1-pZ) * (1-pO) + one * a.b2 + two * b.b2
  Kbb = t(xbmat) %*% diag(wbb, length(wbb)) %*% xbmat
  
  # print(waa)
  # print(wab)
  # print(wbb)
  
  # put together K (= expected information)
  Kab_all = rbind(cbind(Kaa, Kab), cbind(t(Kab), Kbb))
  return(.BlockDiagBind(KppZ, KppO, Kab_all))
}

.HessFun.AIBeta.Score <- function(wmat, xamat, xbmat, gamma, omega, alpha, beta, pZ, pO, mu, phi, K, n.T) {
  # wmat=w.tmp; xmat=x.tmp; gamma=gamma0.w; alpha=alpha0.x; p=p0; mu=mu0; phi=phi0
  # tmp.g = c(exp(wmat %*% gamma)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  # tmp1.g = (1 + tmp.g)^2
  # wpp = 1/(pZ*(1-pZ)) * (tmp.g/tmp1.g)^2
  wpp = -pZ*(1-pZ)
  wpp[is.infinite(wpp) | is.nan(wpp)] = 0
  wmat_expand = do.call(rbind, replicate(K, wmat, simplify = FALSE))
  KppZ = t(.BlockDiag(wmat_expand, K, n.T)) %*% diag(wpp, length(wpp)) %*% (.BlockDiag(wmat_expand, K, n.T))
  
  # tmp.o = c(exp(wmat %*% omega)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  # tmp1.o = (1 + tmp.o)^2
  # wpp = 1/(pO*(1-pO)) * (tmp.o/tmp1.o)^2
  wpp = -pO*(1-pO)
  wpp[is.infinite(wpp) | is.nan(wpp)] = 0
  KppO = t(.BlockDiag(wmat_expand, K, n.T)) %*% diag(wpp, length(wpp)) %*% (.BlockDiag(wmat_expand, K, n.T))
  
  tmp.a = c(exp(xamat %*% alpha)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
  tmp1.a = (1 + tmp.a)^2
  tmp.b = c(exp(xbmat %*% beta))
  sigma = 1/phi-1
  
  psi1 = trigamma(mu * sigma)
  psi2 = trigamma((1 - mu) * sigma)
  # auxiliary transformations
  a = psi1 + psi2
  b = psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(sigma)
  
  waa = - (tmp.a/tmp1.a)^2 * sigma^2 * a * (1-pZ) * (1-pO) 
  xamat_expand = do.call(rbind, replicate(K, xamat, simplify = FALSE))
  Kaa = t(.BlockDiag(xamat_expand, K, n.T)) %*% diag(waa, length(waa)) %*% (.BlockDiag(xamat_expand, K, n.T))
  
  wab = - tmp.a/tmp1.a * (-1/tmp.b) * sigma * (mu * a - psi2) * (1-pZ) * (1-pO) 
  Kab = t(xamat) %*% diag(wab, length(wab)) %*% xbmat
  
  wbb = - (1/tmp.b)^2 * b * (1-pZ) * (1-pO)
  Kbb = t(xbmat) %*% diag(wbb, length(wbb)) %*% xbmat
  
  # print(waa)
  # print(wab)
  # print(wbb)
  
  Kab_all = rbind(cbind(Kaa, Kab), cbind(t(Kab), Kbb))
  return(.BlockDiagBind(KppZ, KppO, Kab_all))
  
  # ## put together K (= expected information)
  # Kab_all = rbind(cbind(Kaa, Kab), cbind(t(Kab), Kbb))
  # if(!inverse) K_ab else chol2inv(chol(K_ab))
}


.ScoreNInfo.ZIBeta <- function(ymat, wmat, xamat, xbmat, 
                               wmat.r, xamat.r, xbmat.r,
                               gamma, alpha, beta){
  K = ncol(ymat) - 1; n.T = nrow(ymat)
  dw = ncol(wmat); da = ncol(xamat); db = ncol(xbmat)
  dw.r = ncol(wmat.r); da.r = ncol(xamat.r); db.r = ncol(xbmat.r)
  ### score function and information matrix 
  Kd3 = K * (dw + da + db)
  K3 = K * 3
  
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for (i in 1:n.T) {
    tmp = c(exp(wmat[i,,drop=FALSE] %*% gamma))
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = c(exp(xamat %*% alpha))
    mv = as.numeric( expa/(1+expa) )
    expb = c(exp(xbmat %*% beta))
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv; bv = (1/sv - 1) * (1-mv); abv = av + bv
    
    diga.av = digamma(av); diga.bv = digamma(bv); diga.abv = digamma(abv)
    A = diga.av - diga.abv; B = diga.bv - diga.abv
    
    triga.av = trigamma(av); triga.bv = trigamma(bv); triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv; B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .ZIeZ(pv, av, bv, Y)
    # post.expectation of delta
    Del.post = par.post$pv.post; av.post = par.post$av.post; bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa)
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      # second derivative associated with gamma
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp
    }
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    ES = ES + c(kronecker(ESi[index1], wmat[i,]),
                kronecker(ESi[index2], xamat[i,]),
                kronecker(ESi[index3], xbmat[i,]))
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    X2 = X[i,] %o% X[i,]
    I.list[[i]] = (-EI - ESS + ES.2)
    I = I + kronecker( I.list[[i]] , X2) 
  }
  
}

.shuffleInCluster <- function(var, id){
  id.unique = sort(unique(id))
  n.cluster = length(id.unique)
  var.sample = var
  for(i in 1:n.cluster){
    index = which(id==id.unique[i])
    var.sample[index,] = var.sample[sample(index),]
  }
  return(var.sample)
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

# Score test --- zero inflation diagnostic test
# .check_zeroinflation <- function(Y, n.boot = 999){
#   K = ncol(Y) - 1; X = matrix(1, nrow(Y), 1)
#   pval = stat = rep(NA, K)
#   alpha0 = matrix(0.001, nrow = 1, ncol = K); phi0 = rep(1, K)
#   gdm.reg = .GDM_EM(Y = Y, X = X, alpha0 = alpha0, phi0 = phi0)
#   alpha0 = gdm.reg$alpha.est; phi0 = gdm.reg$phi.est
#   
#   for(j in 1:K){
#     m = rowSums(Y[,j:ncol(Y)])
#     id.eff = which(m!=0)
#     y = Y[id.eff,j];  m = m[id.eff]; n = length(y)
#     x = X[id.eff,,drop=FALSE]
#     alpha = alpha0[,j,drop=FALSE]
#     phi = phi0[j]; rev.phi = 1/phi
#     w0 = which(y==0)
#     
#     tmp = exp(x %*% alpha)
#     p = as.numeric(tmp/(1+tmp))
#     a = p*phi; a[a<0] = 0
#     b = (1-p)*phi; b[b<0] = 0
#     
#     U = diag(p*(1-p), length(p)) %*% x
#     
#     tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = numeric(n)
#     for (i in 1:n) {
#       mi = m[i]
#       pi = p[i]
#       r = 0:(mi-1)
#       # probI.cdf = 1-pbetabinom.ab(r, mi, a[i], b[i]) # P(Y>r)
#       # probI.rev.cdf = pbetabinom.ab(rev(r), mi, a[i], b[i]) # P(Y<mi-r)
#       probI = dbbinom(r, mi, a[i], b[i])
#       proI.cum = cumsum(probI)
#       probI.cdf = 1-proI.cum
#       probI.rev.cdf = rev(proI.cum)
#       rphi = r*rev.phi
#       tmp1[i] = prod((1+rphi)/(1+rphi-pi), na.rm = T)
#       tmp2[i] = sum(r*pi/((1+rphi-pi)*(1+rphi)), na.rm = T)
#       tmp3[i] = sum(r^2*(probI.cdf/(pi+rphi) +
#                            probI.rev.cdf/(1+rphi-pi) +
#                            1/(1+rphi)), na.rm = T)
#       tmp4[i] = sum(probI.cdf/(pi+rphi)^2 +
#                       probI.rev.cdf/(1+rphi-pi), na.rm = T)
#       tmp5[i] = sum(r*(probI.cdf/(pi+rphi)^2 -
#                          probI.rev.cdf/(1+rphi-pi)^2), na.rm = T)
#       tmp6[i] = sum(-1/(1+rphi[-1]-pi), na.rm = T)
#     }
#     
#     Z = sum(tmp1[w0])-n
#     Igg = sum(tmp1)-n
#     Igp = sum(tmp2)
#     Ipp = sum(tmp3)
#     
#     W1 = diag(tmp4)
#     W2 = diag(tmp5)
#     W3 = diag(tmp6)
#     
#     MM  = U %*% ginv(t(U) %*% W1 %*% U) %*% t(U)
#     V = Igg - sum(t(W3) %*% MM %*% W3) - 
#       (Igp - sum(t(W3) %*% MM %*% W2))^2/(Ipp - sum(t(W2) %*% MM %*% W2))
#     
#     stat[j] = drop(Z*Z/V)
#   }
#   mod.sel = rep(0, K)
#   if(is.null(n.boot)){
#     pval = pchisq(stat, df = 1, lower.tail = FALSE)
#   }else{
#     Y.seqDepth = rowSums(Y)
#     tmp = exp(X %*% alpha0); p = tmp/(1+tmp)
#     a = p%*%diag(phi0,K); a[a<0] = 0; b = (1-p)%*%diag(phi0,K); b[b<0] = 0
#     stat.mat = lapply(1:n.boot, function(x){
#       Y.b = .simData.GDM(a, b, Y.seqDepth)
#       alpha0 = matrix(0.001, nrow = 1, ncol = K); phi0 = rep(1, K)
#       gdm.reg = .GDM_EM(Y = Y.b, X = X, alpha0 = alpha0, phi0 = phi0)
#       alpha.o = gdm.reg$alpha.est; phi.o = gdm.reg$phi.est
#       tmp = exp(X %*% alpha.o); p = tmp/(1+tmp)
#       a = p%*%diag(phi.o,K); a[a<0] = 0; b = (1-p)%*%diag(phi.o,K); b[b<0] = 0
#       stat = rep(NA, K)
#       for (j in 1:K) {
#         m = rowSums(Y.b[,j:ncol(Y.b)])
#         id.eff = which(m!=0)
#         y = Y.b[id.eff,j];  m = m[id.eff]; n = length(y)
#         x = X[id.eff,,drop=FALSE]
#         alpha = alpha.o[,j,drop=FALSE]
#         phi = phi.o[j]; rev.phi = 1/phi
#         w0 = which(y==0)
#         
#         tmp = exp(x %*% alpha)
#         p = as.numeric(tmp/(1+tmp))
#         a = p*phi; a[a<0] = 0
#         b = (1-p)*phi; b[b<0] = 0
#         
#         U = diag(p*(1-p), length(p)) %*% x
#         
#         tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = numeric(n)
#         for (i in 1:n) {
#           mi = m[i]
#           pi = p[i]
#           r = 0:(mi-1)
#           probI = dbbinom(r, mi, a[i], b[i])
#           proI.cum = cumsum(probI)
#           probI.cdf = 1-proI.cum
#           probI.rev.cdf = rev(proI.cum)
#           rphi = r*rev.phi
#           tmp1[i] = prod((1+rphi)/(1+rphi-pi), na.rm = T)
#           tmp2[i] = sum(r*pi/((1+rphi-pi)*(1+rphi)), na.rm = T)
#           tmp3[i] = sum(r^2*(probI.cdf/(pi+rphi) +
#                                probI.rev.cdf/(1+rphi-pi) +
#                                1/(1+rphi)), na.rm = T)
#           tmp4[i] = sum(probI.cdf/(pi+rphi)^2 +
#                           probI.rev.cdf/(1+rphi-pi), na.rm = T)
#           tmp5[i] = sum(r*(probI.cdf/(pi+rphi)^2 -
#                              probI.rev.cdf/(1+rphi-pi)^2), na.rm = T)
#           tmp6[i] = sum(-1/(1+rphi[-1]-pi), na.rm = T)
#         }
#         
#         Z = sum(tmp1[w0])-n
#         Igg = sum(tmp1)-n
#         Igp = sum(tmp2)
#         Ipp = sum(tmp3)
#         
#         W1 = diag(tmp4)
#         W2 = diag(tmp5)
#         W3 = diag(tmp6)
#         
#         MM  = U %*% ginv(t(U) %*% W1 %*% U) %*% t(U)
#         V = Igg - sum(t(W3) %*% MM %*% W3) - 
#           (Igp - sum(t(W3) %*% MM %*% W2))^2/(Ipp - sum(t(W2) %*% MM %*% W2))
#         
#         stat[j] = drop(Z*Z/V)
#       }
#       return(stat)
#     })
#     stat.mat = matrix(unlist(stat.mat), n.boot, K, byrow = TRUE)
#     pval = (colSums(stat.mat >= matrix(stat, n.boot, K, byrow = TRUE)) + 1) / (n.boot + 1)
#   }
#   mod.sel[which(pval < 0.05)] = 1
#   return(mod.sel)
# }

.BB.log.lik <- function(Y){ # fit beta-binomial with intercept only
  X = matrix(1, nrow(Y), 1)
  omega0 = matrix(-Inf, 1, 1)
  alpha0 = beta0 = matrix(0.001, 1, 1)
  y.ori = Y[,1]
  m.ori = rowSums(Y)
  
  # BB Model
  gamma0 = matrix(-Inf, 1, 1)
  gdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, omega0 = omega0, alpha0 = alpha0, beta0 = beta0)
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

# .ZIBB.log.lik <- function(Y){ # fit zero-inflated beta-binomial with intercept only
#   K = ncol(Y) - 1
#   omega0 = matrix(-Inf, 1, K)
#   gamma0 = alpha0 = beta0 = matrix(0.001, 1, K)
#   X = matrix(1, nrow(Y), 1)
#   zigdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, omega0 = omega0, alpha0 = alpha0, beta0 = beta0)
#   tmp = exp(X %*% zigdm.reg$alpha.est); mu = tmp/(1+tmp)
#   tmp = exp(X %*% zigdm.reg$beta.est); phi = tmp/(1+tmp)
#   a.mat = mu * (1/phi - 1); b.mat = (1 - mu) * (1/phi - 1)
#   a.mat[a.mat<0] = 0; b.mat[b.mat<0] = 0
#   
#   ans = sapply(1:K, function(x){
#     y.ori = Y[,x]
#     m.ori = rowSums(Y[,x:(K+1)])
#     a.ori = a.mat[,x]; b.ori = b.mat[,x]
#     id.sel = which(y.ori > 0)
#     y = y.ori[id.sel]; a = a.ori[id.sel]; b = b.ori[id.sel]; m = m.ori[id.sel]
#     logA = lgamma(a+m+b)+lgamma(b)
#     logB = lgamma(a+b)+lgamma(m+b)
#     return(sum( - log(1-exp(logB-logA)) - logA + lgamma(m+1) +
#                   lgamma(y+a) + lgamma(m-y+b) + lgamma(a+b) - lgamma(y+1) -
#                   lgamma(m-y+1) - lgamma(a) ))
#   })
#   
#   return(list(ans=ans, gamma0=zigdm.reg$gamma.est, alpha0=zigdm.reg$alpha.est, beta0=zigdm.reg$beta.est))
# }

.BBnZIBB.log.lik <- function(Y){ # fit both beta-binomial and zero-inflated beta-binomial with intercept only
  X = matrix(1, nrow(Y), 1)
  omega0 = matrix(-Inf, 1, 1)
  alpha0 = beta0 = matrix(0.001, 1, 1)
  y.ori = Y[,1]
  m.ori = rowSums(Y)
  
  # BB Model
  gamma0 = matrix(-Inf, 1, 1)
  gdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, omega0 = omega0, alpha0 = alpha0, beta0 = beta0)
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
  zigdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, omega0 = omega0, alpha0 = alpha0, beta0 = beta0)
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


.OIBB.log.lik <- function(Y){ # fit one-inflated beta-binomial with intercept only
  K = ncol(Y) - 1

  Y.isOI = matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  for (i in 1:nrow(Y)) {
    Y.isOI[i, max(which(Y[i,] > 0))] = 1
  }
  gamma0 = matrix(-Inf, 1, K)
  omega0 = alpha0 = beta0 = matrix(0.001, 1, K)
  X = matrix(1, nrow(Y), 1)
  oigdm.reg = .AIGDM_EM(Y = Y, W = X, Xa = X, Xb = X, gamma0 = gamma0, omega0 = omega0, alpha0 = alpha0, beta0 = beta0)
  tmp = exp(X %*% oigdm.reg$alpha.est); mu = tmp/(1+tmp)
  tmp = exp(X %*% oigdm.reg$beta.est); phi = tmp/(1+tmp)
  a.mat = mu * (1/phi - 1); b.mat = (1 - mu) * (1/phi - 1)
  a.mat[a.mat<0] = 0; b.mat[b.mat<0] = 0
  
  ans = sapply(1:K, function(x){
    y.ori = Y[,x]
    m.ori = rowSums(Y[,x:(K+1)])
    a.ori = a.mat[,x]; b.ori = b.mat[,x]
    id.sel = which(Y.isOI[,x] == 0)
    y = y.ori[id.sel]; a = a.ori[id.sel]; b = b.ori[id.sel]; m = m.ori[id.sel]
    logA = lgamma(a+m+b)+lgamma(b)
    logB = lgamma(a+b)+lgamma(m+b)
    return(sum( - log(1-exp(logB-logA)) - logA + lgamma(m+1) +
                  lgamma(y+a) + lgamma(m-y+b) + lgamma(a+b) - lgamma(y+1) -
                  lgamma(m-y+1) - lgamma(a) ))
  })
  
  return(list(ans=ans, omega0=oigdm.reg$omega.est, alpha0=oigdm.reg$alpha.est, beta0=oigdm.reg$beta.est))
}

.check_zeroinflation <- function(ID, Y, n.boot = 9999, n.cores = n.cores, R.sel = R.sel){
  n.boot = ceiling(n.boot/n.cores)*n.cores
  Y.seqDepth = rowSums(Y)
  bb.ll = .BB.log.lik(Y)
  zibb.ll = .ZIBB.log.lik(Y)
  stat = -2 * (bb.ll$ans - zibb.ll$ans)
  if(is.null(n.boot)){
    pval = pchisq(stat, df = 1, lower.tail = FALSE)
  }else{
    if(is.infinite(stat)){
      cat("Use Pseudo-ECDF approximation p-value\n")
      pval = 1/(n.boot+1)
    }else{
      cl = makeCluster(n.cores)
      clusterExport(cl, c('.simData.GDM', '.BB.log.lik', '.ZIBB.log.lik',
                          '.AIGDM_EM', '.ZOIeZ', '.AIBetaOptim', '.AIBetaNegLoglik', '.AIBetaNegScore',
                          '.LogitOptim', '.LogitNegLoglik', '.LogitNegScore',
                          'stat', 'rowProds', 'Y.seqDepth', 'ID', 'Y', 'n.boot', 'R.sel'), envir = environment())
      
      m = 0
      Nexc = 0
      stat.perm = numeric(n.boot)
      while (all(Nexc < R.sel, m < n.boot)) {
        stat.perm.subset = parSapply(cl, 1:n.cores, function(x){
          id.unique = unique(ID)
          ID.tmp = sort(sample(id.unique, replace = T))
          Y.b.tmp = NULL
          for (id in ID.tmp) {
            Y.b.tmp = rbind(Y.b.tmp, Y[which(ID == id),,drop=FALSE])
          }
          # Y.b.tmp = Y
          # id.unique = sort(unique(ID)); n.cluster = length(id.unique)
          # for(i in 1:n.cluster){
          #   index = which(id==id.unique[i])
          #   Y.b.tmp[index,] = Y[sample(index, replace = TRUE),]
          # }
          bb.ll = .BB.log.lik(Y.b.tmp)
          a = colMeans(bb.ll$a.mat); b = colMeans(bb.ll$b.mat)
          a.mat = matrix(a, nrow(Y), 1); b.mat = matrix(b, nrow(Y), 1)
          # a.mat = bb.ll$a.mat; b.mat = bb.ll$b.mat
          Y.b = .simData.GDM(a.mat, b.mat, Y.seqDepth)
          bb.ll = .BB.log.lik(Y.b)
          zibb.ll = .ZIBB.log.lik(Y.b)
          return(-2 * (bb.ll$ans - zibb.ll$ans))
        })
        Nexc = Nexc + sum(stat.perm.subset >= stat)
        stat.perm[m+1:n.cores] = stat.perm.subset
        m = m + n.cores
      }
      stopCluster(cl)
      
      if(m < n.boot){
        pval = Nexc/m
        cat(sprintf("# of bootstraps: %g\n", m))
        cat("Use ECDF approximation p-value\n")
      }else{
        if(Nexc <= 10){
          pval = tryCatch(.gpd_approx(stat.perm, 250, stat), error=function(err) NA)
          cat(sprintf("# of bootstraps: %g\n", n.boot))
          if(is.na(pval)){
            pval = (Nexc+1)/(n.boot+1)
            cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
          }else{
            cat("Use GPD approximation p-value\n")
          }
        }else{
          pval = Nexc/n.boot
          cat(sprintf("# of bootstraps: %g\n", n.boot))
          cat("Use ECDF approximation p-value\n")
        }
      }
    }
    
    # stat.mat = matrix(unlist(stat.mat), n.boot, K, byrow = TRUE)
    # pval = (colSums(stat.mat >= matrix(stat, n.boot, K, byrow = TRUE)) + 1) / (n.boot + 1)
  }
  return(pval)
}

.check_zeroinflation_modify <- function(ID, Y, n.boot = 9999, R.sel = R.sel){
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
  m = 0
  Nexc = 0
  stat.perm = numeric(n.boot)
  
  while (all(Nexc < R.sel, m < n.boot)) {
    ID.tmp = sort(sample(id.unique, replace = TRUE))
    
    # Extract the rows from Y based on ID.tmp
    Y.b.list = lapply(ID.tmp, function(id) Y[indexList[[as.character(id)]],, drop=FALSE])
    Y.b.tmp = do.call(rbind, Y.b.list)
    
    bb.ll = .BB.log.lik(Y.b.tmp)
    a.mat = matrix(colMeans(bb.ll$a.mat), nrow(Y), 1); b.mat = matrix( colMeans(bb.ll$b.mat), nrow(Y), 1)
    
    Y.b = .simData.GDM(a.mat, b.mat, Y.seqDepth)
    boot.obs = .BBnZIBB.log.lik(Y.b)
    current_stat = -2 * (boot.obs$bb.ans - boot.obs$zibb.ans)
    Nexc = Nexc + sum(current_stat >= stat)
    stat.perm[m+1] = current_stat
    m = m + 1
  }
  
  if (m < n.boot) {
    pval = Nexc/m
    cat(sprintf("# of bootstraps: %g\n", m))
    cat("Use ECDF approximation p-value\n")
  } else {
    if (Nexc <= 10) {
      pval = tryCatch(.gpd_approx(stat.perm, 250, stat), error=function(err) NA)
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
  
  return(pval)
}

.check_zeroinflation_fork <- function(ID, Y, n.boot = 9999, n.cores, R.sel = R.sel){
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
  
  return(pval)
}

# .check_zeroinflation_score <- function(Y){
#   alpha0 = matrix(0.001, nrow = 1, ncol = 1); phi0 = 1
#   m = rowSums(Y)
#   X = matrix(1, nrow(Y), 1)
#   gdm.reg = .GDM_EM(Y = Y, X = X, alpha0 = alpha0, phi0 = phi0)
#   alpha = gdm.reg$alpha.est; phi = gdm.reg$phi.est; rev.phi = 1/phi
#   id.eff = which(m!=0)
#   y = Y[id.eff,1]; m = m[id.eff]; n = length(y)
#   x = X[id.eff,,drop=FALSE]
#   
#   w0 = which(y==0)
#   
#   tmp = exp(x %*% alpha)
#   p = as.numeric(tmp/(1+tmp))
#   a = p*phi; a[a<0] = 0
#   b = (1-p)*phi; b[b<0] = 0
#   
#   U = diag(p*(1-p), length(p)) %*% x
#   
#   tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = numeric(n)
#   for (i in 1:n) {
#     mi = m[i]
#     pi = p[i]
#     r = 0:(mi-1)
#     # probI.cdf = 1-pbetabinom.ab(r, mi, a[i], b[i]) # P(Y>r)
#     # probI.rev.cdf = pbetabinom.ab(rev(r), mi, a[i], b[i]) # P(Y<mi-r)
#     probI = dbbinom(r, mi, a[i], b[i])
#     proI.cum = cumsum(probI)
#     probI.cdf = 1-proI.cum
#     probI.rev.cdf = rev(proI.cum)
#     rphi = r*rev.phi
#     tmp1[i] = prod((1+rphi)/(1+rphi-pi), na.rm = T)
#     tmp2[i] = sum(r*pi/((1+rphi-pi)*(1+rphi)), na.rm = T)
#     tmp3[i] = sum(r^2*(probI.cdf/(pi+rphi) +
#                          probI.rev.cdf/(1+rphi-pi) +
#                          1/(1+rphi)), na.rm = T)
#     tmp4[i] = sum(probI.cdf/(pi+rphi)^2 +
#                     probI.rev.cdf/(1+rphi-pi), na.rm = T)
#     tmp5[i] = sum(r*(probI.cdf/(pi+rphi)^2 -
#                        probI.rev.cdf/(1+rphi-pi)^2), na.rm = T)
#     tmp6[i] = sum(-1/(1+rphi[-1]-pi), na.rm = T)
#   }
#   
#   Z = sum(tmp1[w0])-n
#   Igg = sum(tmp1)-n
#   Igp = sum(tmp2)
#   Ipp = sum(tmp3)
#   
#   W1 = diag(tmp4)
#   W2 = diag(tmp5)
#   W3 = diag(tmp6)
#   
#   MM  = U %*% ginv(t(U) %*% W1 %*% U) %*% t(U)
#   V = Igg - sum(t(W3) %*% MM %*% W3) -
#     (Igp - sum(t(W3) %*% MM %*% W2))^2/(Ipp - sum(t(W2) %*% MM %*% W2))
#   
#   stat = drop(Z*Z/V)
#   pval = pchisq(stat, df = 1, lower.tail = FALSE)
#   return(pval)
# }

.check_oneinflation <- function(ID, Y, n.boot = 9999, n.cores = n.cores, R.sel = R.sel){
  n.boot = ceiling(n.boot/n.cores)*n.cores
  Y.seqDepth = rowSums(Y)
  bb.ll = .BB.log.lik(Y)
  oibb.ll = .OIBB.log.lik(Y)
  stat = -2 * (bb.ll$ans - oibb.ll$ans)
  if(is.null(n.boot)){
    pval = pchisq(stat, df = 1, lower.tail = FALSE)
  }else{
    if(is.infinite(stat)){
      cat("Use Pseudo-ECDF approximation p-value\n")
      pval = 1/(n.boot+1)
    }else{
      cl = makeCluster(n.cores)
      clusterExport(cl, c('.simData.GDM', '.BB.log.lik', '.OIBB.log.lik',
                          '.AIGDM_EM', '.ZOIeZ', '.AIBetaOptim', '.AIBetaNegLoglik', '.AIBetaNegScore',
                          '.LogitOptim', '.LogitNegLoglik', '.LogitNegScore',
                          'stat', 'rowProds', 'Y.seqDepth', 'ID', 'Y', 'n.boot', 'R.sel'), envir = environment())
      m = 0
      Nexc = 0
      stat.perm = numeric(n.boot)
      while (all(Nexc < R.sel, m < n.boot)) {
        stat.perm.subset = parSapply(cl, 1:n.cores, function(x){
          id.unique = unique(ID)
          ID.tmp = sort(sample(id.unique, replace = T))
          Y.b.tmp = NULL
          for (id in ID.tmp) {
            Y.b.tmp = rbind(Y.b.tmp, Y[which(ID == id),,drop=FALSE])
          }
          # Y.b.tmp = Y
          # id.unique = sort(unique(ID)); n.cluster = length(id.unique)
          # for(i in 1:n.cluster){
          #   index = which(id==id.unique[i])
          #   Y.b.tmp[index,] = Y[sample(index, replace = TRUE),]
          # }
          bb.ll = .BB.log.lik(Y.b.tmp)
          a = colMeans(bb.ll$a.mat); b = colMeans(bb.ll$b.mat)
          a.mat = matrix(a, nrow(Y), 1); b.mat = matrix(b, nrow(Y), 1)
          # a.mat = bb.ll$a.mat; b.mat = bb.ll$b.mat
          Y.b = .simData.GDM(a.mat, b.mat, Y.seqDepth)
          bb.ll = .BB.log.lik(Y.b)
          oibb.ll = .OIBB.log.lik(Y.b)
          return(-2 * (bb.ll$ans - oibb.ll$ans))
        })
        Nexc = Nexc + sum(stat.perm.subset >= stat)
        stat.perm[m+1:n.cores] = stat.perm.subset
        m = m + n.cores
      }
      stopCluster(cl)
      
      if(m < n.boot){
        pval = Nexc/m
        cat(sprintf("# of bootstraps: %g\n", m))
        cat("Use ECDF approximation p-value\n")
      }else{
        if(Nexc <= 10){
          pval = tryCatch(.gpd_approx(stat.perm, 250, stat), error=function(err) NA)
          cat(sprintf("# of bootstraps: %g\n", n.boot))
          if(is.na(pval)){
            pval = (Nexc+1)/(n.boot+1)
            cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
          }else{
            cat("Use GPD approximation p-value\n")
          }
        }else{
          pval = Nexc/n.boot
          cat(sprintf("# of bootstraps: %g\n", n.boot))
          cat("Use ECDF approximation p-value\n")
        }
      }
    }
  }
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
.est_para <- function(ID, W, Xa, Xb, Y, zi.check, oi.check){
  #ID=ID.tmp; W=W.r.tmp; Xa=Xa.r.tmp; Xb=Xb.r.tmp; Y=cbind(Y.tmp[,1], rowSums(Y.tmp[,(1+1):(K+1),drop=FALSE])); zi.check=F; oi.check=F
  da = ncol(Xa); db = ncol(Xb); dw = ncol(W); K = ncol(Y)-1
  
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
  nt = nrow(Y) # the total number of observed units
  N1 = sum(n.reps*(n.reps-1)/2)
  
  r.step = 0.05 # set up the tuning parameter for updating alpha
  max = 5e3 # maximal iterations
  
  if(zi.check){
    gamma0 = matrix(0.001, nrow = dw, ncol = K)
  }else{
    gamma0 = matrix(-Inf, nrow = dw, ncol = K)
  }
  if(oi.check){
    omega0 = matrix(0.001, nrow = dw, ncol = K)
  }else{
    omega0 = matrix(-Inf, nrow = dw, ncol = K)
  }
  
  alpha0 = matrix(0.001, nrow = da, ncol = K)
  beta0 = matrix(0.001, nrow = db, ncol = K)
  # correlation of same subject, all taxa, at different times
  rho0 = sig0 = tau0 = 0.1
  # rho0 = sig0 = tau0 = 0
  
  if(zi.check + oi.check == 2){
    mod = "ZOI"
  }else if(zi.check + oi.check == 1){
    if(zi.check)
      mod = "ZI"
    if(oi.check)
      mod = "OI"
  }else{
    mod = "NI"
  }
  
  aigdm.reg = .AIGDM_EM(Y = Y, W = W, Xa = Xa, Xb = Xb, 
                        gamma0 = gamma0, omega0 = omega0, alpha0 = alpha0, beta0 = beta0)
  gamma0 = aigdm.reg$gamma.est
  omega0 = aigdm.reg$omega.est
  alpha0 = aigdm.reg$alpha.est
  beta0 = aigdm.reg$beta.est
  
  ng = no = dw; na = da
  for (l in 1:max) {
    all.ddg = matrix(0, ng, ng) # derivative of GEE
    all.geeg = matrix(0, ng, 1) # GEE
    all.sig1 = all.sig2 = 0
    
    all.ddo = matrix(0, no, no)
    all.geeo = matrix(0, no, 1) 
    all.tau1 = all.tau2 = 0
    
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
      
      tmp = exp(w.tmp %*% omega0)
      tmp[is.na(tmp)] = 0
      pO0.mat = tmp/(1+tmp)
      pO0.mat[is.infinite(tmp) & tmp>0] = 1
      pO0 = c(pO0.mat)
      
      tmp = exp(xa.tmp %*% alpha0)
      mu0.mat = tmp/(1+tmp)
      tmp = exp(xb.tmp %*% beta0)
      phi0.mat = tmp/(1+tmp)
      
      av = mu0.mat * (1/phi0.mat - 1)
      bv = (1 - mu0.mat) * (1/phi0.mat - 1)
      
      n.T = n.reps[i]
      eDelZ0.mat = eDelO0.mat = eDel0.mat = eZ0.mat = matrix(NA, n.T, 1)
      eZStar.mat = matrix(NA, n.T, 1)
      for (t in 1:n.T) {
        par.post = .ZOIeZ_vec(pZ0.mat[t,,drop=FALSE], pO0.mat[t,,drop=FALSE], av[t,,drop=FALSE], bv[t,,drop=FALSE], y.tmp[t,])
        eDelZ0.mat[t,] = par.post$pZv.post
        eDelO0.mat[t,] = par.post$pOv.post
        eDel0.mat[t,] = par.post$pv.post
        
        tmp = par.post$av.post + par.post$bv.post
        # post.expectation of logZ conditional on delta==0
        A.R.lst[[i]][t,] = digamma(par.post$av.post) - digamma(tmp)
        # post.expectation of log(1-Z) conditional on delta==0
        B.R.lst[[i]][t,] = digamma(par.post$bv.post) - digamma(tmp)
        eZ0.mat[t,] = par.post$av.post/tmp
        # post expectation of log(Z/(1-Z)) conditional on delta==0
        eZStar.mat[t,] = digamma(par.post$av.post) - digamma(par.post$bv.post)
      }
      
      eDelZ0 = c(eDelZ0.mat); eDelO0 = c(eDelO0.mat); eDel0 = c(eDel0.mat)
      Del.R.lst[[i]] = eDel0.mat
      
      ntot = ntot + colSums((1-eDel0.mat)^2)
      N2 = N2 + sum(.getSubDiag( (1-eDel0.mat)%*%t(1-eDel0.mat)))
      
      # Vdel
      cor.DelZ = diag(sig0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(sig0,each=n.T), n.T)
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      tmp.g = t(.BlockDiag(.dpg(w.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ)
      geeg = as.matrix(apply(.RowbyRow(t(tmp.g), eDelZ0 - pZ0), 2, sum))
      temp.g1 = tmp.g %*% .BlockDiag(.dpg(w.tmp, gamma0), 1, n.T) + geeg %*% t(geeg)
      all.ddg = all.ddg + temp.g1 ## sum up to N for derivative of gee for gamma
      all.geeg = all.geeg + geeg  ## sum up to N for gee for gamma
      
      cor.DelO = diag(tau0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(tau0,each=n.T), n.T)
      VdelO = sqrt(.var_logit(pO0)) %*% cor.DelO %*% sqrt(.var_logit(pO0))
      tmp.o = t(.BlockDiag(.dpg(w.tmp, omega0), 1, n.T)) %*% ginv(VdelO)
      geeo = as.matrix(apply(.RowbyRow(t(tmp.o), eDelO0 - pO0), 2, sum))
      temp.o1 = tmp.o %*% .BlockDiag(.dpg(w.tmp, omega0), 1, n.T) + geeo %*% t(geeo)
      all.ddo = all.ddo + temp.o1 
      all.geeo = all.geeo + geeo
      
      # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
      eZ0 = c(eZ0.mat); mu0 = c(mu0.mat); phi0 = c(phi0.mat); phi0 = c(phi0.mat)
      eZStar = c(eZStar.mat); muStar = digamma(av) - digamma(bv)
      
      # tmp.a = as.numeric(exp(xa.tmp%*% alpha0))
      # mu.tmp = as.numeric( tmp.a/(1+tmp.a) )
      # tmp.b = as.numeric( exp(xb.tmp %*% beta0) )
      # phi.tmp = as.numeric( tmp.b/(1+tmp.b) )
      # a = (1/phi.tmp - 1) * mu.tmp
      # b = (1/phi.tmp - 1) * (1 - mu.tmp)
      # a[a < 0] = 0
      # b[b < 0] = 0
      # 
      # zstar = eZStar; zdagger = eZDagger
      # mustar = digamma(a) - digamma(b)
      # mudagger = digamma(b) - digamma(a + b)
      # 
      # (1/tmp.b) * mu.tmp * (1/(1+tmp.a)) * xa.tmp * (zstar - mustar) 
      # 
      # 1/tmp.b <--> 1/phi0 - 1
      # mu.tmp * (1/(1+tmp.a)) * xa.tmp <--> t(.BlockDiag(.dma(xa.tmp, alpha0), 1, n.T))
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      tmp.a = t(.BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0))
      geea = as.matrix(apply(.RowbyRow(t(tmp.a), eZStar - muStar), 2, sum))
      temp.a1 = tmp.a %*% .BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T) + geea %*% t(geea)
      # Vz = sqrt(.var_Z(mu0, phi0)) %*% cor.Z %*% sqrt(.var_Z(mu0, phi0))
      # tmp.a = t(.BlockDiag(.dma(xa.tmp, alpha0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0))
      # geea = as.matrix(apply(.RowbyRow(t(tmp.a), eZ0 - mu0), 2, sum))
      # temp.a1 = tmp.a %*% .BlockDiag(.dma(xa.tmp, alpha0), 1, n.T) + geea %*% t(geea)
      all.dda = all.dda + temp.a1 # sum up to N for derivative of gee for alpha
      all.geea = all.geea + geea # sum up to N for gee for alpha
      if(n.T != 1){
        # update sig estimate
        all.sig1 = all.sig1 + .corr_zero1(eDelZ0, pZ0, n.T)
        all.sig2 = all.sig2 + .corr_zero2(eDelZ0, pZ0, n.T)
        
        # update tau estimate
        all.tau1 = all.tau1 + .corr_zero1(eDelO0, pO0, n.T)
        all.tau2 = all.tau2 + .corr_zero2(eDelO0, pO0, n.T)
        
        # update rho estimate
        all.rho1 = all.rho1 + .corr_rho1(eDel0, eZStar, muStar, mu0, phi0, n.T)
        all.rho2 = all.rho2 + .corr_rho2(eDel0, eZStar, muStar, mu0, phi0, n.T)
        
      }
    }
    
    ########## updating scheme #############
    if(l==1){all.geeg1 <- all.geeg} else{ all.geeg1 <- ((l-1)*all.geeg1 + all.geeg)/l} 
    if(l==1){all.ddg1 <- all.ddg} else{ all.ddg1 <- ((l-1)*all.ddg1 + all.ddg)/l} 
    gnew0 = c(gamma0) + r.step * ginv(all.ddg1) %*% all.geeg
    gnew = (l * c(gamma0) + gnew0)/(l + 1)
    
    if(l==1){all.geeo1 <- all.geeo} else{ all.geeo1 <- ((l-1)*all.geeo1 + all.geeo)/l} 
    if(l==1){all.ddo1 <- all.ddo} else{ all.ddo1 <- ((l-1)*all.ddo1 + all.ddo)/l} 
    onew0 = c(omega0) + r.step * ginv(all.ddo1) %*% all.geeo
    onew = (l * c(omega0) + onew0)/(l + 1)
    
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
    
    tau.new = (all.tau1/N1)/(all.tau2/nt) #Note: p can be zero then sig.new = NaN
    nan.id = which(is.nan(tau.new) | is.infinite(tau.new) | tau.new == 1)
    if(length(nan.id) != 0)
      tau.new[nan.id] = 0
    
    rho.new = (all.rho1/N2)/(all.rho2/ntot)
    nan.id = which(is.nan(rho.new) | is.infinite(rho.new))
    if(length(nan.id) != 0)
      rho.new[nan.id] = 0
    
    # sig.new = tau.new = rho.new = 0
    
    if(mod == "ZI"){
      if(max(abs(c(gnew - gamma0, anew - alpha0, bnew - beta0, sig.new - sig0, rho.new - rho0))) <  1e-4 | l == max)
        break
    }else if(mod == "OI"){
      if(max(abs(c(onew - omega0, anew - alpha0, bnew - beta0, tau.new - tau0, rho.new - rho0))) <  1e-4 | l == max)
        break
    }else if(mod == "ZOI"){
      if(max(abs(c(gnew - gamma0, onew - omega0, anew - alpha0, bnew - beta0, sig.new - sig0, tau.new - tau0, rho.new - rho0))) <  1e-4 | l == max)
        break
    }else if(mod == "NI"){
      if(max(abs(c(anew - alpha0, bnew - beta0, rho.new - rho0))) <  1e-4 | l == max)
        break
    }
    gamma0 = matrix(gnew, ncol = K); omega0 = matrix(onew, ncol = K); sig0 = sig.new; tau0 = tau.new
    alpha0 = matrix(anew, ncol = K); beta0 = matrix(bnew, ncol = K); rho0 = rho.new
  }
  
  gamma0 = matrix(gnew, ncol = K); omega0 = matrix(onew, ncol = K); sig0 = sig.new; tau0 = tau.new
  alpha0 = matrix(anew, ncol = K); beta0 = matrix(bnew, ncol = K); rho0 = rho.new
  
  DelZ.R.lst = eDelZ0.lst.store = pZ0.lst.store = 
    DelO.R.lst = eDelO0.lst.store = pO0.lst.store = 
    Del.R.lst = eDel0.lst.store = A.R.lst = B.R.lst = 
    A2.R.lst = B2.R.lst = AB.R.lst = 
    eZ0.lst.store = mu0.lst.store = phi0.lst.store = 
    eZStar.lst.store = muStar.lst.store =
    pseudoll.score.beta.lst.store = vector("list", N)
  for (i in 1:N) {
    DelZ.R.lst[[i]] = DelO.R.lst[[i]] = Del.R.lst[[i]]  = A.R.lst[[i]] = B.R.lst[[i]] = 
      A2.R.lst[[i]] = B2.R.lst[[i]] = AB.R.lst[[i]] = matrix(NA, n.reps[i], 1)
    eDelZ0.lst.store[[i]] = pZ0.lst.store[[i]] = eDelO0.lst.store[[i]] = pO0.lst.store[[i]] = 
      eDel0.lst.store[[i]] = eZ0.lst.store[[i]] = mu0.lst.store[[i]] = phi0.lst.store[[i]] =
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
    tmp[is.na(tmp)] = 0
    pZ0.mat = tmp/(1+tmp)
    pZ0.mat[is.infinite(tmp) & tmp>0] = 1
    pZ0 = c(pZ0.mat)
    pZ0.lst.store[[i]] = pZ0
    
    tmp = exp(w.tmp %*% omega0)
    tmp[is.na(tmp)] = 0
    pO0.mat = tmp/(1+tmp)
    pO0.mat[is.infinite(tmp) & tmp>0] = 1
    pO0 = c(pO0.mat)
    pO0.lst.store[[i]] = pO0
    
    tmp = exp(xa.tmp %*% alpha0)
    mu0.mat = tmp/(1+tmp)
    tmp = exp(xb.tmp %*% beta0)
    phi0.mat = tmp/(1+tmp)
    av = mu0.mat * (1/phi0.mat - 1)
    bv = (1 - mu0.mat) * (1/phi0.mat - 1)
    muStar.mat = digamma(av) - digamma(bv)
    
    n.T = n.reps[i]
    eDelZ0.mat = eDelO0.mat = eDel0.mat = eZ0.mat = eZStar.mat = matrix(NA, n.T, 1)
    for (t in 1:n.T) {
      par.post = .ZOIeZ_vec(pZ0.mat[t,,drop=FALSE], pO0.mat[t,,drop=FALSE], av[t,,drop=FALSE], bv[t,,drop=FALSE], y.tmp[t,])
      eDelZ0.mat[t,] = par.post$pZv.post
      eDelO0.mat[t,] = par.post$pOv.post
      eDel0.mat[t,] = par.post$pv.post
      
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
    
    DelZ.R.lst[[i]] = eDelZ0.mat
    eDelZ0.lst.store[[i]] = c(eDelZ0.mat)
    
    DelO.R.lst[[i]] = eDelO0.mat
    eDelO0.lst.store[[i]] = c(eDelO0.mat)
    
    Del.R.lst[[i]] = eDel0.mat
    eDel0.lst.store[[i]] = c(eDel0.mat)
    
    beta.par = c(alpha0, beta0)
    beta.data = list(Del = c(Del.R.lst[[i]]), A = c(A.R.lst[[i]]), B = c(B.R.lst[[i]]), 
                     Xa = xa.tmp, Xb = xb.tmp)
    pseudoll.score.beta.lst.store[[i]] = .GenerateScoreFun(beta.par, beta.data)
    
    ntot = ntot + colSums((1-eDel0.mat)^2)
    N2 = N2 + sum(.getSubDiag( (1-eDel0.mat)%*%t(1-eDel0.mat)))
    
    # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
    eZ0 = c(eZ0.mat); mu0 = c(mu0.mat); phi0 = c(phi0.mat)
    eZStar = c(eZStar.mat); muStar = c(muStar.mat)
    eZStar.lst.store[[i]] = eZStar; muStar.lst.store[[i]] = muStar
    eZ0.lst.store[[i]] = eZ0; mu0.lst.store[[i]] = mu0; phi0.lst.store[[i]] = phi0
  }
  
  #print(paste0("Converge?", ifelse(l == max, 0, 1)))
  return(list(CONV = ifelse(l == max, 0, 1), W = W, Xa = Xa, Xb = Xb,
              gamma0 = matrix(gnew, ncol = K), omega0 = matrix(onew, ncol = K), 
              alpha0 = matrix(anew, ncol = K), beta0 = matrix(bnew, ncol = K), 
              sig0 = sig.new, tau0 = tau.new, rho0 = rho.new,
              eDelZ0.lst = eDelZ0.lst.store, eDelO0.lst = eDelO0.lst.store, eDel0.lst = eDel0.lst.store,
              A.R.lst = A.R.lst, B.R.lst = B.R.lst,
              A2.R.lst = A2.R.lst, B2.R.lst = B2.R.lst, AB.R.lst = AB.R.lst,
              pZ0.lst = pZ0.lst.store, pO0.lst = pO0.lst.store,
              eZ0.lst = eZ0.lst.store, mu0.lst = mu0.lst.store, phi0.lst = phi0.lst.store,
              eZStar.lst = eZStar.lst.store, muStar.lst = muStar.lst.store,
              pseudoll.score.beta.lst = pseudoll.score.beta.lst.store))
}

.score_test_mean <- function(ID, X, X.index, Y, est.para.lst, zi.check, oi.check){
  K = ncol(Y)-1
  ai.check = zi.check + oi.check
  stat = df = numeric(K)
  # ind.tmp = NULL
  for (j in 1:K) {
    # j = 1
    Y.tmp = cbind(Y[,j], rowSums(Y[,(j+1):(K+1),drop=FALSE]))
    
    est.para = est.para.lst[[j]]
    ID.tmp = ID; W.r.tmp = est.para$W; Xa.tmp = X; Xa.r.tmp = est.para$Xa; Xb.r.tmp = est.para$Xb
    dw.r = ncol(W.r.tmp); db.r = ncol(Xb.r.tmp); da = ncol(Xa.tmp); da.r = ncol(Xa.r.tmp)
    
    id.unique = sort(unique(ID.tmp)); N = length(id.unique); n.reps = as.vector(table(ID.tmp))
    ng = no = dw.r; na = da; nb = db.r
    
    if(ai.check[j]){
      if(!zi.check[j]){
        # para.ind = c(sort(sapply(W.index, function(x) seq(x, dw, dw))),
        #              sort(sapply(X.index, function(x) seq(x, dx, dx))) + dw)
        para.ind = sort(sapply(X.index, function(x) seq(x, na, na))) + no
        para.keep = c(sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
        # ind.tmp = c(ind.tmp, c(rep(1, length(W.index)), rep(2, length(X.index))))
      }else if(!oi.check[j]){
        para.ind = sort(sapply(X.index, function(x) seq(x, na, na))) + ng
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
        # ind.tmp = c(ind.tmp, c(rep(0, length(W.index)), rep(2, length(X.index))))
      }else{
        # para.ind = c(sort(sapply(W.index, function(x) seq(x, dw*2, dw))),
        #              sort(sapply(X.index, function(x) seq(x, dx, dx))) + dw*2)
        para.ind = sort(sapply(X.index, function(x) seq(x, na, na))) + ng + no
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
        # ind.tmp = c(ind.tmp, c(rep(0:1, each = length(W.index)), rep(2, length(X.index))))
      }
    }else{
      para.ind = sort(sapply(X.index, function(x) seq(x, na, na)))
      para.keep = c(sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                    sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      # ind.tmp = c(ind.tmp, rep(2, length(X.index)))
    }
    
    gamma0 = est.para$gamma0; omega0 = est.para$omega0; alpha0 = est.para$alpha0; beta0 = est.para$beta0
    sig0 = est.para$sig0; tau0 = est.para$tau0; rho0 = est.para$rho0
    eDel0.lst = est.para$eDel0.lst; eDelZ0.lst = est.para$eDelZ0.lst; eDelO0.lst = est.para$eDelO0.lst
    A.R.lst = est.para$A.R.lst; B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst; B2.R.lst = est.para$B2.R.lst; AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst; pO0.lst = est.para$pO0.lst
    eZ0.lst = est.para$eZ0.lst; mu0.lst = est.para$mu0.lst; phi0.lst = est.para$phi0.lst
    eZStar.lst = est.para$eZStar.ls; muStar.lst = est.para$muStar.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, ng+no+na+nb, ng+no+na+nb)
    ind.g = 1:ng; ind.o = (ng+1):(ng+no); ind.a = (ng+no+1):(ng+no+na); ind.b = (ng+no+na+1):(ng+no+na+nb)
    
    for(i in 1:N){
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel,,drop=FALSE]
      xa.r.tmp = Xa.r.tmp[ID.sel,,drop=FALSE]
      xa.tmp = Xa.tmp[ID.sel,,drop=FALSE]
      xb.r.tmp = Xb.r.tmp[ID.sel,,drop=FALSE]
      n.T = n.reps[i]
      
      pZ0 = pZ0.lst[[i]] # unlist(pZ0.lst[i,])
      eDelZ0 = eDelZ0.lst[[i]]
      cor.DelZ = diag(sig0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(sig0,each=n.T), n.T)
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      tmp.g.w = t(.BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ)
      dvd.g.w = t(.BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ) %*% .BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)
      geeg.w = as.matrix(apply(.RowbyRow(t(tmp.g.w), eDelZ0 - pZ0), 2, sum))
      
      pO0 = pO0.lst[[i]]
      eDelO0 = eDelO0.lst[[i]]
      cor.DelO = diag(tau0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(tau0,each=n.T), n.T)
      VdelO = sqrt(.var_logit(pO0)) %*% cor.DelO %*% sqrt(.var_logit(pO0))
      tmp.o.w = t(.BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)) %*% ginv(VdelO)
      dvd.o.w = t(.BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)) %*% ginv(VdelO) %*% .BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)
      geeo.w = as.matrix(apply(.RowbyRow(t(tmp.o.w), eDelO0 - pO0), 2, sum))
      
      mu0 = mu0.lst[[i]]; phi0 = phi0.lst[[i]]
      eZ0 = eZ0.lst[[i]]; eDel0 = eDel0.lst[[i]]
      eZStar = eZStar.lst[[i]]; muStar = muStar.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      alpha0.x = matrix(0, da, 1)
      alpha0.x[1:da.r,] = alpha0
      tmp.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0))
      dvd.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0)) %*% .BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZStar - muStar), 2, sum))
      
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.r.tmp)
      # dvd.b.x = t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.r.tmp)) %*% t(t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.r.tmp)))
      geeb.x = as.matrix(pseudoll.score.beta)
      
      dvd.ab.x = t(.BlockDiag(.dma(xa.tmp, alpha0.x, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0)) %*% .BlockDiag(.dmb(xb.r.tmp, beta0, mu0, phi0), 1, n.T)
      
      ES = rbind(geeg.w, geeo.w, geea.x, geeb.x)
      U = U + ES
      
      a = (1/phi0 - 1) * mu0; a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0); b[b < 0] = 0
      A = digamma(a) - digamma(a+b)
      B = digamma(b) - digamma(a+b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1-eDel0)*( - A + A.post )
      two.R = (1-eDel0)*( - B + B.post )
      tmp.dvd = .DVD(w.r.tmp, xa.tmp, xb.r.tmp, gamma0, omega0, alpha0.x, beta0, pZ0, pO0, mu0, phi0, one.R, two.R, 1, n.T)
      # tmp.Kab = matrix(0, nrow(tmp.dvd$Kab), ncol(tmp.dvd$Kab))
      dvd.ab.all = rbind(cbind(dvd.a.x, dvd.ab.x), cbind(t(tmp.dvd$Kab), tmp.dvd$Kbb))
      EI = .BlockDiagBind(dvd.g.w, dvd.o.w, dvd.ab.all)
      I = I + (-EI)
      
      # I = I + (-EI)
      
      # ES.2 = ES %*% t(ES)
      # ESS = ES.2
      # 
      # A2.post = c(A2.R.lst[[i]])
      # B2.post = c(B2.R.lst[[i]])
      # AB.post = c(AB.R.lst[[i]])
      # b.tmp = 1/phi0-1 # 1/expb
      # a0.tmp = 1-mu0; a1.tmp = mu0; a2.tmp = mu0 * (1-mu0)
      # # a0.tmp = 1/(1+expa); a1.tmp = expa/(1+expa); a2.tmp = expa/(1+expa)^2
      # a.a = b.tmp * a2.tmp; a.b = - b.tmp * a1.tmp
      # b.a = - b.tmp * a2.tmp; b.b = - b.tmp * a0.tmp
      # p0 = (1-pZ0) * (1-pO0)
      # 
      # ESS[ind.g, ind.g] = sum(eDelZ0 - 2 * eDelZ0 * pZ0 + pZ0 * pZ0)
      # ESS[ind.o, ind.o] = sum(eDelO0 - 2 * eDelO0 * pO0 + pO0 * pO0)
      # tmp1 = -A * a.a - B * b.a
      # diag(ESS[ind.a, ind.a]) = sum((1-eDel0)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a))
      # tmp2 = -A * a.b - B * b.b
      # diag(ESS[ind.b, ind.b]) = sum((1-eDel0)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b))
      # 
      # tmp = colSums(-t(t(pZ0)) %*% t(ES[ind.a,,drop=FALSE]))
      # ESS[ind.g, ind.a] = tmp
      # ESS[ind.a, ind.g] = tmp
      # tmp = colSums(-t(t(pO0)) %*% t(ES[ind.a,,drop=FALSE]))
      # ESS[ind.o, ind.a] = tmp
      # ESS[ind.a, ind.o] = tmp
      # 
      # tmp = colSums(-t(t(pZ0)) %*% t(ES[ind.b,,drop=FALSE]))
      # ESS[ind.g, ind.b] = tmp
      # ESS[ind.b, ind.g] = tmp
      # tmp = colSums(-t(t(pO0)) %*% t(ES[ind.b,,drop=FALSE]))
      # ESS[ind.o, ind.b] = tmp
      # ESS[ind.b, ind.o] = tmp
      # 
      # tmp = sum((1-eDel0)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post ))
      # ESS[ind.a, ind.b] = tmp
      # ESS[ind.b, ind.a] = tmp
      # 
      # I = I + (-EI-ESS+ES.2)
      
      M[ind.g, ind.g] = M[ind.g, ind.g] + geeg.w %*% t(geeg.w) * (diag(1) %x% matrix(1, ng, ng))
      M[ind.g, ind.o] = M[ind.g, ind.o] + geeg.w %*% t(geeo.w) * (diag(1) %x% matrix(1, ng, no))
      M[ind.o, ind.o] = M[ind.o, ind.o] + geeo.w %*% t(geeo.w) * (diag(1) %x% matrix(1, no, no))
      M[ind.o, ind.g] = M[ind.o, ind.g] + geeo.w %*% t(geeg.w) * (diag(1) %x% matrix(1, no, ng))
      
      M[ind.a, ind.a] = M[ind.a, ind.a] + geea.x %*% t(geea.x) * (diag(1) %x% matrix(1, na, na))
      M[ind.b, ind.b] = M[ind.b, ind.b] + geeb.x %*% t(geeb.x) * (diag(1) %x% matrix(1, nb, nb))
      
      M[ind.a, ind.b] = M[ind.a, ind.b] + geea.x %*% t(geeb.x) * (diag(1) %x% matrix(1, na, nb))
      M[ind.b, ind.a] = M[ind.b, ind.a] + geeb.x %*% t(geea.x) * (diag(1) %x% matrix(1, nb, na))
      
      M[c(ind.g, ind.o), c(ind.a, ind.b)] = M[c(ind.g, ind.o), c(ind.a, ind.b)] + rbind(geeg.w, geeo.w) %*% t(rbind(geea.x, geeb.x)) * (diag(1) %x% matrix(1, ng+no, na+nb))
      M[c(ind.a, ind.b), c(ind.g, ind.o)] = M[c(ind.a, ind.b), c(ind.g, ind.o)] + rbind(geea.x, geeb.x) %*% t(rbind(geeg.w, geeo.w)) * (diag(1) %x% matrix(1, na+nb, ng+no))
    }
    
    U.r = U[para.keep,,drop=FALSE]
    I.r = I[para.keep, para.keep]
    M.r = M[para.keep, para.keep]
    
    H.r = .colwise.cbind(-I.r[para.ind, -para.ind] %*% ginv(I.r[-para.ind, -para.ind]), diag(length(para.ind)), para.ind)
    stat[j] = c(t(U.r[para.ind,,drop=F]) %*% ginv(H.r %*% M.r %*% t(H.r)) %*% U.r[para.ind,,drop=F])
    df[j] = length(para.ind)
    
    # if(j == 1){
    #   U.tmp = U.r[para.ind,,drop=F]
    #   V.tmp = H.r %*% M.r %*% t(H.r)
    # }else{
    #   U.tmp = rbind(U.tmp, U.r[para.ind,,drop=F])
    #   V.tmp = .BlockDiagBind(V.tmp, H.r %*% M.r %*% t(H.r))
    # }
  }
  stat.sum = sum(stat); df.sum = sum(df)
  pval = 1 - pchisq(stat.sum, df.sum)
  return(list(stat.sum = stat.sum, pval.sum = pval, stat = stat))
}

.score_test_disp <- function(ID, X, X.index, Y, est.para.lst, zi.check, oi.check){
  # X = Xb
  K = ncol(Y)-1
  ai.check = zi.check + oi.check
  stat = df = numeric(K)
  # ind.tmp = NULL
  for (j in 1:K) {
    Y.tmp = cbind(Y[,j], rowSums(Y[,(j+1):(K+1),drop=FALSE]))
    
    est.para = est.para.lst[[j]]
    ID.tmp = ID; W.r.tmp = est.para$W; Xb.tmp = X; Xb.r.tmp = est.para$Xb; Xa.r.tmp = est.para$Xa
    dw.r = ncol(W.r.tmp); da.r = ncol(Xa.r.tmp); db = ncol(Xb.tmp); db.r = ncol(Xb.r.tmp)
    
    id.unique = sort(unique(ID.tmp)); N = length(id.unique); n.reps = as.vector(table(ID.tmp))
    ng = no = dw.r; na = da.r; nb = db
    
    if(ai.check[j]){
      if(!zi.check[j]){
        para.ind = sort(sapply(X.index, function(x) seq(x, nb, nb))) + no + na
        para.keep = c(sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }else if(!oi.check[j]){
        para.ind = sort(sapply(X.index, function(x) seq(x, nb, nb))) + ng + na
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }else{
        para.ind = sort(sapply(X.index, function(x) seq(x, nb, nb))) + ng + no + na
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }
    }else{
      para.ind = sort(sapply(X.index, function(x) seq(x, nb, nb))) + na
      para.keep = c(sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                    sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
    }
    
    gamma0 = est.para$gamma0; omega0 = est.para$omega0; alpha0 = est.para$alpha0; beta0 = est.para$beta0
    sig0 = est.para$sig0; tau0 = est.para$tau0; rho0 = est.para$rho0
    eDel0.lst = est.para$eDel0.lst; eDelZ0.lst = est.para$eDelZ0.lst; eDelO0.lst = est.para$eDelO0.lst
    A.R.lst = est.para$A.R.lst; B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst; B2.R.lst = est.para$B2.R.lst; AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst; pO0.lst = est.para$pO0.lst
    eZ0.lst = est.para$eZ0.lst; mu0.lst = est.para$mu0.lst; phi0.lst = est.para$phi0.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    
    # update Hessian matrix with reduced models
    U = 0
    I = 0
    M = matrix(0, ng+no+na+nb, ng+no+na+nb)
    ind.g = 1:ng; ind.o = (ng+1):(ng+no); ind.a = (ng+no+1):(ng+no+na); ind.b = (ng+no+na+1):(ng+no+na+nb)
    
    for(i in 1:N){
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel,,drop=FALSE]
      xa.r.tmp = Xa.r.tmp[ID.sel,,drop=FALSE]
      xb.tmp = Xb.tmp[ID.sel,,drop=FALSE]
      n.T = n.reps[i]
      
      pZ0 = pZ0.lst[[i]] # unlist(pZ0.lst[i,])
      eDelZ0 = eDelZ0.lst[[i]]
      cor.DelZ = diag(sig0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(sig0,each=n.T), n.T)
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      tmp.g.w = t(.BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ)
      dvd.g.w = t(.BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ) %*% .BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)
      geeg.w = as.matrix(apply(.RowbyRow(t(tmp.g.w), eDelZ0 - pZ0), 2, sum))
      
      pO0 = pO0.lst[[i]]
      eDelO0 = eDelO0.lst[[i]]
      cor.DelO = diag(tau0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(tau0,each=n.T), n.T)
      VdelO = sqrt(.var_logit(pO0)) %*% cor.DelO %*% sqrt(.var_logit(pO0))
      tmp.o.w = t(.BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)) %*% ginv(VdelO)
      dvd.o.w = t(.BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)) %*% ginv(VdelO) %*% .BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)
      geeo.w = as.matrix(apply(.RowbyRow(t(tmp.o.w), eDelO0 - pO0), 2, sum))
      
      mu0 = mu0.lst[[i]]; phi0 = phi0.lst[[i]]
      eZ0 = eZ0.lst[[i]]; eDel0 = eDel0.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      tmp.a.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0))
      dvd.a.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0)) %*% .BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZ0 - mu0), 2, sum))
      
      beta0.x = matrix(0, db, 1)
      beta0.x[1:db.r,] = beta0
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.tmp)
      # dvd.b.x = t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)) %*% t(t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)))
      geeb.x = as.matrix(pseudoll.score.beta)
      
      dvd.ab.x = t(.BlockDiag(.dma(xa.r.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0)) %*% .BlockDiag(.dmb(xb.tmp, beta0.x, mu0, phi0), 1, n.T)
      
      ES = rbind(geeg.w, geeo.w, geea.x, geeb.x)
      U = U + ES
      
      a = (1/phi0 - 1) * mu0; a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0); b[b < 0] = 0
      A = digamma(a) - digamma(a+b)
      B = digamma(b) - digamma(a+b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1-eDel0)*( - A + A.post )
      two.R = (1-eDel0)*( - B + B.post )
      tmp.dvd = .DVD(w.r.tmp, xa.r.tmp, xb.tmp, gamma0, omega0, alpha0, beta0.x, pZ0, pO0, mu0, phi0, one.R, two.R, 1, n.T)
      # tmp.Kab = matrix(0, nrow(tmp.dvd$Kab), ncol(tmp.dvd$Kab))
      dvd.ab.all = rbind(cbind(dvd.a.x, dvd.ab.x), cbind(t(tmp.dvd$Kab), tmp.dvd$Kbb))
      EI = .BlockDiagBind(dvd.g.w, dvd.o.w, dvd.ab.all)
      I = I + (-EI)
      
      # ES.2 = ES %*% t(ES)
      # ESS = ES.2
      # 
      # A2.post = c(A2.R.lst[[i]])
      # B2.post = c(B2.R.lst[[i]])
      # AB.post = c(AB.R.lst[[i]])
      # b.tmp = 1/phi0-1 # 1/expb
      # a0.tmp = 1-mu0; a1.tmp = mu0; a2.tmp = mu0 * (1-mu0)
      # # a0.tmp = 1/(1+expa); a1.tmp = expa/(1+expa); a2.tmp = expa/(1+expa)^2
      # a.a = b.tmp * a2.tmp; a.b = - b.tmp * a1.tmp
      # b.a = - b.tmp * a2.tmp; b.b = - b.tmp * a0.tmp
      # p0 = (1-pZ0) * (1-pO0)
      # 
      # ESS[ind.g, ind.g] = sum(eDelZ0 - 2 * eDelZ0 * pZ0 + pZ0 * pZ0)
      # ESS[ind.o, ind.o] = sum(eDelO0 - 2 * eDelO0 * pO0 + pO0 * pO0)
      # tmp1 = -A * a.a - B * b.a
      # diag(ESS[ind.a, ind.a]) = sum((1-eDel0)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a))
      # tmp2 = -A * a.b - B * b.b
      # diag(ESS[ind.b, ind.b]) = sum((1-eDel0)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b))
      # 
      # tmp = colSums(-t(t(pZ0)) %*% t(ES[ind.a,,drop=FALSE]))
      # ESS[ind.g, ind.a] = tmp
      # ESS[ind.a, ind.g] = tmp
      # tmp = colSums(-t(t(pO0)) %*% t(ES[ind.a,,drop=FALSE]))
      # ESS[ind.o, ind.a] = tmp
      # ESS[ind.a, ind.o] = tmp
      # 
      # tmp = colSums(-t(t(pZ0)) %*% t(ES[ind.b,,drop=FALSE]))
      # ESS[ind.g, ind.b] = tmp
      # ESS[ind.b, ind.g] = tmp
      # tmp = colSums(-t(t(pO0)) %*% t(ES[ind.b,,drop=FALSE]))
      # ESS[ind.o, ind.b] = tmp
      # ESS[ind.b, ind.o] = tmp
      # 
      # tmp = sum((1-eDel0)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post ))
      # ESS[ind.a, ind.b] = tmp
      # ESS[ind.b, ind.a] = tmp
      # 
      # I = I + (-EI-ESS+ES.2)
      
      M[ind.g, ind.g] = M[ind.g, ind.g] + geeg.w %*% t(geeg.w) * (diag(1) %x% matrix(1, ng, ng))
      M[ind.g, ind.o] = M[ind.g, ind.o] + geeg.w %*% t(geeo.w) * (diag(1) %x% matrix(1, ng, no))
      M[ind.o, ind.o] = M[ind.o, ind.o] + geeo.w %*% t(geeo.w) * (diag(1) %x% matrix(1, no, no))
      M[ind.o, ind.g] = M[ind.o, ind.g] + geeo.w %*% t(geeg.w) * (diag(1) %x% matrix(1, no, ng))
      
      M[ind.a, ind.a] = M[ind.a, ind.a] + geea.x %*% t(geea.x) * (diag(1) %x% matrix(1, na, na))
      M[ind.b, ind.b] = M[ind.b, ind.b] + geeb.x %*% t(geeb.x) * (diag(1) %x% matrix(1, nb, nb))
      
      M[ind.a, ind.b] = M[ind.a, ind.b] + geea.x %*% t(geeb.x) * (diag(1) %x% matrix(1, na, nb))
      M[ind.b, ind.a] = M[ind.b, ind.a] + geeb.x %*% t(geea.x) * (diag(1) %x% matrix(1, nb, na))
      
      M[c(ind.g, ind.o), c(ind.a, ind.b)] = M[c(ind.g, ind.o), c(ind.a, ind.b)] + rbind(geeg.w, geeo.w) %*% t(rbind(geea.x, geeb.x)) * (diag(1) %x% matrix(1, ng+no, na+nb))
      M[c(ind.a, ind.b), c(ind.g, ind.o)] = M[c(ind.a, ind.b), c(ind.g, ind.o)] + rbind(geea.x, geeb.x) %*% t(rbind(geeg.w, geeo.w)) * (diag(1) %x% matrix(1, na+nb, ng+no))
    }
    
    U.r = U[para.keep,,drop=FALSE]
    I.r = I[para.keep, para.keep]
    M.r = M[para.keep, para.keep]
    
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
                                  n.perm, perm.type, zi.check, oi.check, n.cores, R.sel, seed){
  # ID=ID.tmp; X=Xa.tmp; Y=Y.tmp; stat.sum=stat.mean.sum; test = "mean"; Rconf=Rconf.mean
  cl = makeCluster(n.cores)
  n.perm = ceiling(n.perm/n.cores)*n.cores # change the maximal permutation times based on number of cores
  clusterExport(cl, c('.shuffleInCluster', '.score_test_mean', '.score_test_disp',
                      '.var_logit', '.dpg', '.var_ZStar', '.dma', '.dmb', '.DVD',
                      '.RowbyRow', '.GenerateScoreFun', 
                      '.BlockDiag', '.BlockDiagBind', '.colwise.cbind',
                      'ID', 'X', 'X.index', 'Y', 'est.para.lst', 
                      'stat.sum', 'test', 'Rconf', 'n.perm', 'perm.type',
                      'zi.check', 'oi.check', 'R.sel', 'seed', 'ginv'), envir = environment())
  m = 0; Nexc = 0; stat.perm = numeric(n.perm)
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
  
  while(all(Nexc < R.sel, m < n.perm)){
    tmp.subset = parSapply(cl, 1:n.cores, function(x){
      set.seed(seed+m+x)
      if(perm.type == "WTH"){
        var.new = .shuffleInCluster(X[,X.index,drop=FALSE], ID)
        X[,X.index] = Rconf %*% var.new
        X.p = X
      }else if(perm.type == "BTW"){
        # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
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
      }
      if(test == "mean")
        tmp = .score_test_mean(ID, X.p, X.index, Y, est.para.lst, zi.check, oi.check)
      if(test == "disp")
        tmp = .score_test_disp(ID, X.p, X.index, Y, est.para.lst, zi.check, oi.check)
      return(tmp$stat.sum)
    })
    Nexc = Nexc + sum(tmp.subset >= stat.sum)
    stat.perm[m+1:n.cores] = tmp.subset
    m = m + n.cores
  }
  stopCluster(cl)
  
  if(m < n.perm){
    pval = Nexc/m
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
      pval = Nexc/n.perm
      cat(sprintf("# of permutations: %g\n", n.perm))
      cat("Use ECDF approximation p-value\n")
    }
  }
  return(pval = pval)
  # return(list(pval = pval, stat.perm = stat.perm))
  
}

# adaptive multiple cores, much faster
.score_test_perm_fork <- function(ID, X, X.index, Y, est.para.lst, stat.sum, test, Rconf,
                             n.perm, perm.type, zi.check, oi.check, n.cores, R.sel, seed){
  #Y = Y.tmp;
  n.perm = ceiling(n.perm/n.cores)*n.cores # change the maximal permutation times based on number of cores
  m = 0; Nexc = 0; stat.perm = numeric(n.perm)
  splitX = split(X[, X.index, drop = FALSE], ID); N = length(splitX)
  # The parallel package provides a way to reproducibly generate random numbers in a parallel environment via 
  # the “L’Ecuyer-CMRG” random number generator. 
  # Note that this is not the default random number generator so you 
  # will have to set it explicitly.
  while(all(Nexc < R.sel, m < n.perm)){
    tmp.subset = mclapply(1:n.cores, function(x){
      # set.seed(seed+m+x)
      if(perm.type == "WTH"){
        X.p.int = .shuffleInCluster(X[,X.index,drop=FALSE], ID)
      }else if(perm.type == "BTW"){
        cluster.p = sample(1:N)
        X.p.list = lapply(cluster.p, function(id) X[splitX[[as.character(id)]],X.index, drop=FALSE])
        X.p.int = do.call(rbind, X.p.list)
      }
      X[,X.index] = Rconf %*% X.p.int
      X.p = X
      if(test == "mean")
        tmp = .score_test_mean(ID, X.p, X.index, Y, est.para.lst, zi.check, oi.check)
      if(test == "disp")
        tmp = .score_test_disp(ID, X.p, X.index, Y, est.para.lst, zi.check, oi.check)
      return(tmp$stat.sum)
    })
    Nexc = Nexc + sum(unlist(tmp.subset) >= stat.sum)
    stat.perm[m+1:n.cores] = unlist(tmp.subset)
    m = m + n.cores
  }
  
  if(m < n.perm){
    pval = Nexc/m
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
      pval = Nexc/n.perm
      cat(sprintf("# of permutations: %g\n", n.perm))
      cat("Use ECDF approximation p-value\n")
    }
  }
  return(pval)
}

# .simData.GDM <- function(a, b, SeqDepth){
#   N = nrow(a); K = ncol(a)
#   Z = matrix(rbeta(N*K, a, b), N, K)
#   P = matrix(NA, N, K+1)
#   for (j in 1:K) {
#     P[,j] = Z[,j]*matrixStats::rowProds(1 - Z[,0:(j-1), drop = FALSE])
#   }
#   P[, ncol(P)] =  1 - rowSums(P[, -ncol(P), drop = FALSE])
#   P[P < 0] = 0
#   Y = c()
#   for (i in 1:N) {
#     Y = rbind(Y, t(rmultinom(1, SeqDepth[i], P[i,])))
#   }
#   return(Y)
# }

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

# .choose_b <- function(alpha, c) {
#   error <- alpha * c
#   B <- alpha*(1 - alpha) / (error^2)
#   return(B)
# }
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

.F.test <- function(x){
  
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.pchisqsum2.davies <- function(Q, lambda){
  
  delta = rep(0, length(lambda))
  acc = 1e-4
  delta <- delta[lambda > 0]
  lambda <- lambda[lambda > 0]
  
  tmp <- CompQuadForm::davies(q = Q, lambda = lambda, 
                              delta = delta, acc = acc)
  if (tmp$ifault > 0) {
    lambda <- zapsmall(lambda, digits = 2)
    delta <- delta[lambda > 0]
    lambda <- lambda[lambda > 0]
    tmp <- CompQuadForm::farebrother(q = Q, lambda = lambda, 
                                     delta = delta)
  }
  Qq <- if ("Qq" %in% names(tmp)) 
    tmp$Qq
  else tmp$res
  return(list(p = Qq, errflag = 0))
}