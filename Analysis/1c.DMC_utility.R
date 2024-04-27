DM.Cluster <- function(ID, Y, case, Tax = NULL, min.depth = 0,
                       perm.type = NULL, n.perm = NULL, fdr.alpha = 0.05, seed){
  # ID = id; Y = Y.m; case = exposure; Tax = tax; min.depth = 0; perm.type = "BTW"
  
  remove.subject = which(rowSums(Y) < min.depth)
  if(length(remove.subject) > 0){
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth))
    case = case[-remove.subject]
    Y = Y[-remove.subject, , drop=FALSE]
  }
  
  keep = which(colSums(Y)>0)
  Y = Y[, keep, drop=FALSE]
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
  
  pval = NULL
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
      Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
      remove.index = which(colSums(Y.tmp)==0)
      if(length(remove.index)==ncol(Y)){
        next
      }else{
        if(length(remove.index)>0){
          Y.tmp = Y.tmp[, -remove.index, drop=FALSE] 
        }
        if(ncol(Y.tmp)==1){
          next
        }else{
          print(paste0("Now processing Rank ", n.rank-k, ".", m))
          keep.ind = which(rowSums(Y.tmp) != 0)
          
          Y.tmp = Y.tmp[keep.ind,,drop=FALSE]
          case.tmp = case[keep.ind]
          ID.tmp = ID[keep.ind]
          
          colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)
          # Y.tmp = Y.tmp[,c(colOrder[-1], colOrder[1])]
          Y.tmp = Y.tmp[,colOrder]
          
          group.data = NULL
          for(a.case in unique(case.tmp) ){
            group.data = c(group.data, list(Y.tmp[case.tmp==a.case,,drop=FALSE]) )
          }
          
          # dm.stat.m = tryCatch( Xmcupo.sevsample(group.data)$`Xmcupo statistics`, error=function(err) NA )
          # if(length(dm.stat.m) != 1)
          #   dm.stat.m = NA
          # if(!is.na(dm.stat.m)){
          #   pval.m = .dm_perm(ID.tmp, Y.tmp, case.tmp, dm.stat.m, test = "mean", perm.type = "BTW", n.perm, R.sel, seed)
          # }else{
          #   pval.m = NA
          # }
          # 
          # dm.stat.o = .tryCatch.Xoc.sevsample(group.data)
          # if(length(dm.stat.o) != 1)
          #   dm.stat.o = NA
          # if(!is.na(dm.stat.o)){
          #   pval.o = .dm_perm(ID.tmp, Y.tmp, case.tmp, dm.stat.o, test = "disp", "BTW", n.perm, R.sel, seed)
          # }else{
          #   pval.o = NA
          # }
          # pval = cbind(pval, c(pval.m, pval.o))
          
          dm.stat.d = tryCatch( Xdc.sevsample(group.data)$`Xdc statistics`, error=function(err) NA )
          if(length(dm.stat.d) != 1)
            dm.stat.d = NA
          if(!is.na(dm.stat.d)){
            pval.d = .dm_perm(ID.tmp, Y.tmp, case.tmp, dm.stat.d, test = "omni", "BTW", n.perm, R.sel, seed)
          }else{
            pval.d = NA
          }
          
          pval = c(pval, pval.d)
        }
      }
    }# lineage loop
  }# level loop
  # colnames(pval) = subtree
  # identify significant lineages
  # rslt = apply(pval, 1, .identifySigLineagesNsimesTest, fdr.alpha)
  # names(rslt) = c("Mean-Resampling", "Disp-Resampling")
  names(pval) = subtree
  rslt = .identifySigLineagesNsimesTest(pval, fdr.alpha)
  return(rslt)
  
}

.dm_perm <- function(ID, Y, case, stat, test, perm.type, n.perm, R.sel, seed){
  # ID=ID.tmp; Y=Y.tmp; case=case.tmp; stat=dm.stat.m; test = "mean"
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID)); n = nrow(Y)
  m = 0; Nexc = 0; stat.perm = numeric(n.perm)
  while (all(Nexc < R.sel, m < n.perm)) {
    set.seed(seed+m)
    
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    case.unique = NULL
    for(i in 1:N){
      tmp = case[which(ID==id.unique[i])]
      if(length(unique(tmp))!=1)
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      case.unique = c(case.unique, unique(tmp))
    }
    case.perm = NULL
    sample.cluster.idx = sample(1:N)
    for(i in 1:N){
      n.T = n.reps[i]
      case.perm = c(case.perm, rep(case.unique[sample.cluster.idx[i]], n.T))
    }
    
    group.data.perm = NULL
    for(a.case in unique(case.perm) ){
      group.data.perm = c(group.data.perm, list(Y[case.perm==a.case,,drop=FALSE]) )
    }
    
    if(test == "mean"){
      tmp.perm = tryCatch( Xmcupo.sevsample(group.data.perm)$`Xmcupo statistics`,  error=function(err) NA)
      if(length(tmp.perm) != 1)
        tmp.perm = NA
      if(!is.na(tmp.perm) ){
        Nexc = Nexc + (tmp.perm >= stat)
      }
      stat.perm[m+1] = tmp.perm
    }
    
    if(test == "disp"){
      tmp.perm = .tryCatch.Xoc.sevsample(group.data.perm)
      if(length(tmp.perm) != 1)
        tmp.perm = NA
      if(!is.na(tmp.perm) ){
        Nexc = Nexc + (tmp.perm >= stat)
      }
      stat.perm[m+1] = tmp.perm
    }
    
    if(test == "omni"){
      tmp.perm = tryCatch( Xdc.sevsample(group.data.perm)$`Xdc statistics`,  error=function(err) NA)
      if(length(tmp.perm) != 1)
        tmp.perm = NA
      if(!is.na(tmp.perm) ){
        Nexc = Nexc + sum(tmp.perm >= stat)
      }
      stat.perm[m+1] = tmp.perm
    }
    m = m + 1
  }
  
  if(m < n.perm){
    pval = Nexc/m
    cat(sprintf("# of permutations: %g\n", m))
    cat("Use ECDF approximation p-value\n")
  }else{
    if(Nexc <= 10){
      pval = tryCatch(.gpd_approx(stat.perm, 250, stat), error=function(err) NA)
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

.tryCatch.Xoc.sevsample <- function(x){
  setTimeLimit(cpu = 120, elapsed = 120, transient = TRUE)
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  })
  
  tryCatch({
    # do some stuff
    Xoc.sevsample(x)$`Xoc statistics`
    # Sys.sleep(11)
  }, error = function(e) {
    if (grepl("reached elapsed time limit|reached CPU time limit", e$message)) {
      # we reached timeout, apply some alternative method or do something else
      NA
    } else {
      # error not related to timeout
      NA
    }
  })
}