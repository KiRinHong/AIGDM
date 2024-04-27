library(extraDistr)
# NORmal-To-Anything (NORTA) algorithm
# simulate correlated Dirichlet from Gamma 
# .corrData.Gamma <- function(n.T, r, mu){
#   # mu = mu.mat
#   N = nrow(mu)/n.T
#   K = ncol(mu)
#   mu1 = mu[1:N,]
#   mu2 = mu[(N+1):(N*n.T),]
#   gammavars = c()
#   for(j in 1:K){
#     tmp = matrix(NA, nrow = N, ncol = n.T)
#     rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
#                       Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
#     mu.j = cbind(mu1[,j], mu2[,j])
#     for(t in 1:n.T){
#       tmp[,t] = qgamma(pnorm(rawvars[,t]), mu.j[,t]) 
#     }
#     gammavars = cbind(gammavars, tmp)
#   }
#   
#   gamma.datalist = vector("list", length = n.T)
#   names(gamma.datalist) = paste0("TimePoint", 1:n.T)
#   for (t in 1:n.T) {
#     gamma.datalist[[t]] = gammavars[, seq(t, K * n.T, n.T)]
#   }
#   return(gamma.datalist)
# }
# 
# .corrData.Beta <- function(n.T, r, beta.a, beta.b){
#   N = nrow(beta.a)/n.T
#   K = ncol(beta.a)
#   beta.a1 = beta.a[1:N,]
#   beta.a2 = beta.a[(N+1):(N*n.T),]
#   beta.b1 = beta.b[1:N,]
#   beta.b2 = beta.b[(N+1):(N*n.T),]
#   
#   betavars = c()
#   for (j in 1:K) {
#     tmp = matrix(NA, nrow = N, ncol = n.T)
#     rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
#                       Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
#     beta.a.j = cbind(beta.a1[,j], beta.a2[,j])
#     beta.b.j = cbind(beta.b1[,j], beta.b2[,j])
#     for (t in 1:n.T) {
#       tmp[,t] = qbeta(pnorm(rawvars[,t]), beta.a.j[,t], beta.b.j[,t]) 
#     }
#     betavars = cbind(betavars, tmp)
#   }
#   beta.datalist = vector("list", length = n.T)
#   names(beta.datalist) = paste0("TimePoint", 1:n.T)
#   for (t in 1:n.T) {
#     beta.datalist[[t]] = betavars[, seq(t, K * n.T, n.T)]
#   }
#   return(beta.datalist)
# }

# simulate correlated ZI Log-Normal from Normal
# .corrData.ZINormal <- function(n.T, r, mu, var, pstr){
#   # mu = mu.mat; var = var.mat; cov = cov.real; pstr = pstr.mat
#   N = nrow(mu)/n.T
#   K = ncol(mu)
#   mu1 = mu[1:N,]
#   mu2 = mu[(N+1):(N*n.T),]
#   var1 = var[1:N,]
#   var2 = var[(N+1):(N*n.T),]
#   pstr.1 = pstr[1:N,]
#   pstr.2 = pstr[(N+1):(N*n.T),]
#   normvars = c()
#   
#   # tmp = exp( mvrnorm(1, mu=mu[1,], Sigma=cov)  )
#   
#   for(j in 1:K){
#     tmp = tmp_G = tmp_Delta = matrix(NA, nrow = N, ncol = n.T)
#     rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
#                       Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
#     mu.j = cbind(mu1[,j], mu2[,j])
#     sd.j = sqrt(cbind(var1[,j], var2[,j]))
#     for(t in 1:n.T){
#       tmp_G[,t] = qnorm(pnorm(rawvars[,t]), mu.j[,t], sd.j[,t]) 
#     }
#     rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
#                       Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
#     pstr.j = cbind(pstr.1[,j], pstr.2[,j])
#     for (t in 1:n.T) {
#       tmp_Delta[,t] = qbinom(pnorm(rawvars[,t]), 1, pstr.j[,t])
#       tmp[,t] = (1-tmp_Delta[,t])*exp(tmp_G[,t])
#     }
#     normvars = cbind(normvars, tmp)
#   }
#   
#   norm.datalist = vector("list", length = n.T)
#   names(norm.datalist) = paste0("TimePoint", 1:n.T)
#   
#   for (t in 1:n.T) {
#     norm.datalist[[t]] = exp(normvars[, seq(t, K * n.T, n.T)])
#   }
#   return(norm.datalist)
# }

.corrData.ZIBeta <- function(n.T, r, beta.a, beta.b, pstr){
  # beta.a = a.mat; beta.b = b.mat; pstr = pstr.mat
  N = nrow(beta.a)/n.T
  K = ncol(beta.a)
  beta.a1 = beta.a[1:N,]
  beta.a2 = beta.a[(N+1):(2*N),]
  beta.a3 = beta.a[(2*N+1):(3*N),]
  beta.a4 = beta.a[(3*N+1):(4*N),]
  beta.a5 = beta.a[(4*N+1):(5*N),]
  
  beta.b1 = beta.b[1:N,]
  beta.b2 = beta.b[(N+1):(2*N),]
  beta.b3 = beta.b[(2*N+1):(3*N),]
  beta.b4 = beta.b[(3*N+1):(4*N),]
  beta.b5 = beta.b[(4*N+1):(5*N),]
  
  pstr.1 = pstr[1:N,]
  pstr.2 = pstr[(N+1):(2*N),]
  pstr.3 = pstr[(2*N+1):(3*N),]
  pstr.4 = pstr[(3*N+1):(4*N),]
  pstr.5 = pstr[(4*N+1):(5*N),]
  betavars = c()
  
  for (j in 1:K) {
    tmp = tmp_Z = tmp_Delta = matrix(NA, nrow = N, ncol = n.T)
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    beta.a.j = cbind(beta.a1[,j], beta.a2[,j], beta.a3[,j], beta.a4[,j], beta.a5[,j])
    beta.b.j = cbind(beta.b1[,j], beta.b2[,j], beta.b3[,j], beta.b4[,j], beta.b5[,j])
    for (t in 1:n.T) {
      tmp_Z[,t] = qbeta(pnorm(rawvars[,t]), beta.a.j[,t], beta.b.j[,t]) 
    }
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    pstr.j = cbind(pstr.1[,j], pstr.2[,j], pstr.3[,j], pstr.4[,j], pstr.5[,j])
    for (t in 1:n.T) {
      tmp_Delta[,t] = qbinom(pnorm(rawvars[,t]), 1, pstr.j[,t])
      # if(sum(1 - tmp_Delta[,t]) != 0) # probably all(tmp_Delta == 1) --> all 0
      #   tmp[,t] = (1-tmp_Delta[,t])*tmp_Z[,t]
      tmp[,t] = (1-tmp_Delta[,t])*tmp_Z[,t]
    }
    betavars = cbind(betavars, tmp)
  }
  beta.datalist = vector("list", length = n.T)
  names(beta.datalist) = paste0("TimePoint", 1:n.T)
  for (t in 1:n.T) {
    beta.datalist[[t]] = betavars[, seq(t, K * n.T, n.T)]
  }
  return(beta.datalist)
}
.corrData.ZIBeta.diab <- function(n.T, r, beta.a, beta.b, pstr){
  # beta.a = a.mat; beta.b = b.mat; pstr = pstr.mat
  N = nrow(beta.a)/n.T
  K = ncol(beta.a)
  beta.a1 = beta.a[1:N,]
  beta.a2 = beta.a[(N+1):(2*N),]
  beta.a3 = beta.a[(2*N+1):(3*N),]
  beta.a4 = beta.a[(3*N+1):(4*N),]
  beta.a5 = beta.a[(4*N+1):(5*N),]
  beta.a6 = beta.a[(5*N+1):(6*N),]
  beta.a7 = beta.a[(6*N+1):(7*N),]
  beta.a8 = beta.a[(7*N+1):(8*N),]
  beta.a9 = beta.a[(8*N+1):(9*N),]
  beta.a10 = beta.a[(9*N+1):(10*N),]
  beta.a11 = beta.a[(10*N+1):(11*N),]
  beta.a12 = beta.a[(11*N+1):(12*N),]
  beta.a13 = beta.a[(12*N+1):(13*N),]
  beta.a14 = beta.a[(13*N+1):(14*N),]
  beta.a15 = beta.a[(14*N+1):(15*N),]
  beta.a16 = beta.a[(15*N+1):(16*N),]
  beta.a17 = beta.a[(16*N+1):(17*N),]
  beta.a18 = beta.a[(17*N+1):(18*N),]
  beta.a19 = beta.a[(18*N+1):(19*N),]
  beta.a20 = beta.a[(19*N+1):(20*N),]
  beta.b1 = beta.b[1:N,]
  beta.b2 = beta.b[(N+1):(2*N),]
  beta.b3 = beta.b[(2*N+1):(3*N),]
  beta.b4 = beta.b[(3*N+1):(4*N),]
  beta.b5 = beta.b[(4*N+1):(5*N),]
  beta.b6 = beta.b[(5*N+1):(6*N),]
  beta.b7 = beta.b[(6*N+1):(7*N),]
  beta.b8 = beta.b[(7*N+1):(8*N),]
  beta.b9 = beta.b[(8*N+1):(9*N),]
  beta.b10 = beta.b[(9*N+1):(10*N),]
  beta.b11 = beta.b[(10*N+1):(11*N),]
  beta.b12 = beta.b[(11*N+1):(12*N),]
  beta.b13 = beta.b[(12*N+1):(13*N),]
  beta.b14 = beta.b[(13*N+1):(14*N),]
  beta.b15 = beta.b[(14*N+1):(15*N),]
  beta.b16 = beta.b[(15*N+1):(16*N),]
  beta.b17 = beta.b[(16*N+1):(17*N),]
  beta.b18 = beta.b[(17*N+1):(18*N),]
  beta.b19 = beta.b[(18*N+1):(19*N),]
  beta.b20 = beta.b[(19*N+1):(20*N),]
  pstr.1 = pstr[1:N,]
  pstr.2 = pstr[(N+1):(2*N),]
  pstr.3 = pstr[(2*N+1):(3*N),]
  pstr.4 = pstr[(3*N+1):(4*N),]
  pstr.5 = pstr[(4*N+1):(5*N),]
  pstr.6 = pstr[(5*N+1):(6*N),]
  pstr.7 = pstr[(6*N+1):(7*N),]
  pstr.8 = pstr[(7*N+1):(8*N),]
  pstr.9 = pstr[(8*N+1):(9*N),]
  pstr.10 = pstr[(9*N+1):(10*N),]
  pstr.11 = pstr[(10*N+1):(11*N),]
  pstr.12 = pstr[(11*N+1):(12*N),]
  pstr.13 = pstr[(12*N+1):(13*N),]
  pstr.14 = pstr[(13*N+1):(14*N),]
  pstr.15 = pstr[(14*N+1):(15*N),]
  pstr.16 = pstr[(15*N+1):(16*N),]
  pstr.17 = pstr[(16*N+1):(17*N),]
  pstr.18 = pstr[(17*N+1):(18*N),]
  pstr.19 = pstr[(18*N+1):(19*N),]
  pstr.20 = pstr[(19*N+1):(20*N),]
  betavars = c()
  
  for (j in 1:K) {
    tmp = tmp_Z = tmp_Delta = matrix(NA, nrow = N, ncol = n.T)
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    beta.a.j = cbind(beta.a1[,j], beta.a2[,j], beta.a3[,j], beta.a4[,j], beta.a5[,j],
                     beta.a6[,j], beta.a7[,j], beta.a8[,j], beta.a9[,j], beta.a10[,j],
                     beta.a11[,j], beta.a12[,j], beta.a13[,j], beta.a14[,j], beta.a15[,j],
                     beta.a16[,j], beta.a17[,j], beta.a18[,j], beta.a19[,j], beta.a20[,j])
    beta.b.j = cbind(beta.b1[,j], beta.b2[,j], beta.b3[,j], beta.b4[,j], beta.b5[,j],
                     beta.b6[,j], beta.b7[,j], beta.b8[,j], beta.b9[,j], beta.b10[,j],
                     beta.b11[,j], beta.b12[,j], beta.b13[,j], beta.b14[,j], beta.b15[,j],
                     beta.b16[,j], beta.b17[,j], beta.b18[,j], beta.b19[,j], beta.b20[,j])
    for (t in 1:n.T) {
      tmp_Z[,t] = qbeta(pnorm(rawvars[,t]), beta.a.j[,t], beta.b.j[,t]) 
    }
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    pstr.j = cbind(pstr.1[,j], pstr.2[,j], pstr.3[,j], pstr.4[,j], pstr.5[,j],
                   pstr.6[,j], pstr.7[,j], pstr.8[,j], pstr.9[,j], pstr.10[,j],
                   pstr.11[,j], pstr.12[,j], pstr.13[,j], pstr.14[,j], pstr.15[,j],
                   pstr.16[,j], pstr.17[,j], pstr.18[,j], pstr.19[,j], pstr.20[,j])
    for (t in 1:n.T) {
      tmp_Delta[,t] = qbinom(pnorm(rawvars[,t]), 1, pstr.j[,t])
      # if(sum(1 - tmp_Delta[,t]) != 0) # probably all(tmp_Delta == 1) --> all 0
      #   tmp[,t] = (1-tmp_Delta[,t])*tmp_Z[,t]
      tmp[,t] = (1-tmp_Delta[,t])*tmp_Z[,t]
    }
    betavars = cbind(betavars, tmp)
  }
  beta.datalist = vector("list", length = n.T)
  names(beta.datalist) = paste0("TimePoint", 1:n.T)
  for (t in 1:n.T) {
    beta.datalist[[t]] = betavars[, seq(t, K * n.T, n.T)]
  }
  return(beta.datalist)
}

.corrData.OIBeta <- function(n.T, r, beta.a, beta.b, pstr){
  # beta.a = a.mat; beta.b = b.mat; pstr = pstr.mat
  N = nrow(beta.a)/n.T
  K = ncol(beta.a)
  beta.a1 = beta.a[1:N,]
  beta.a2 = beta.a[(N+1):(2*N),]
  beta.a3 = beta.a[(2*N+1):(3*N),]
  beta.a4 = beta.a[(3*N+1):(4*N),]
  beta.a5 = beta.a[(4*N+1):(5*N),]
  beta.b1 = beta.b[1:N,]
  beta.b2 = beta.b[(N+1):(2*N),]
  beta.b3 = beta.b[(2*N+1):(3*N),]
  beta.b4 = beta.b[(3*N+1):(4*N),]
  beta.b5 = beta.b[(4*N+1):(5*N),]
  pstr.1 = pstr[1:N,]
  pstr.2 = pstr[(N+1):(2*N),]
  pstr.3 = pstr[(2*N+1):(3*N),]
  pstr.4 = pstr[(3*N+1):(4*N),]
  pstr.5 = pstr[(4*N+1):(5*N),]
  betavars = c()
  
  for (j in 1:K) {
    tmp = tmp_Z = tmp_Delta = matrix(NA, nrow = N, ncol = n.T)
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    beta.a.j = cbind(beta.a1[,j], beta.a2[,j], beta.a3[,j], beta.a4[,j], beta.a5[,j])
    beta.b.j = cbind(beta.b1[,j], beta.b2[,j], beta.b3[,j], beta.b4[,j], beta.b5[,j])
    
    for (t in 1:n.T) {
      tmp_Z[,t] = qbeta(pnorm(rawvars[,t]), beta.a.j[,t], beta.b.j[,t]) 
    }
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    pstr.j = cbind(pstr.1[,j], pstr.2[,j], pstr.3[,j], pstr.4[,j], pstr.5[,j])
    for (t in 1:n.T) {
      tmp_Delta[,t] = qbinom(pnorm(rawvars[,t]), 1, pstr.j[,t])
      # if(sum(1 - tmp_Delta[,t]) != 0) # probably all(tmp_Delta == 1) --> all 0
      #   tmp[,t] = (1-tmp_Delta[,t])*tmp_Z[,t]
      tmp[,t] = (1-tmp_Delta[,t])*tmp_Z[,t]
    }
    tmp[tmp == 0] = 1
    betavars = cbind(betavars, tmp)
  }
  beta.datalist = vector("list", length = n.T)
  names(beta.datalist) = paste0("TimePoint", 1:n.T)
  for (t in 1:n.T) {
    beta.datalist[[t]] = betavars[, seq(t, K * n.T, n.T)]
  }
  return(beta.datalist)
}

# Z --> P --> Y
.simCount.ZIGDM <- function(Z.mat, SeqDepth){
  #Z.mat = Z.datalist$TimePoint2; SeqDepth = 1e3
  P = matrix(NA, nrow = nrow(Z.mat), ncol = ncol(Z.mat)+1)
  for (j in 1:(ncol(P)-1)) {
    P[,j] = Z.mat[,j]*rowProds(1 - Z.mat[,0:(j-1), drop = FALSE])
  }
  P[, ncol(P)] =  1 - rowSums(P[, -ncol(P), drop = FALSE])
  P[P < 0] = 0
  Y = c()
  for (i in 1:nrow(Z.mat)) {
    Y = rbind(Y, t(rmultinom(1, rpois(1, SeqDepth), P[i,]))) # high zero-inflation level --> -2.22e-16
  }
  return(Y)
}

.simCount.OIGDM <- function(Z.mat, SeqDepth){
  #Z.mat = Z.datalist$TimePoint2; SeqDepth = 1e3
  P = matrix(NA, nrow = nrow(Z.mat), ncol = ncol(Z.mat)+1)
  for (j in 1:(ncol(P)-1)) {
    P[,j] = Z.mat[,j]*rowProds(1 - Z.mat[,0:(j-1), drop = FALSE])
  }
  P[, ncol(P)] =  1 - rowSums(P[, -ncol(P), drop = FALSE])
  P[P < 0] = 0
  
  Y = c()
  for (i in 1:nrow(Z.mat)) {
    Y = rbind(Y, t(rmultinom(1, rpois(1, SeqDepth), P[i,]))) # high zero-inflation level --> -2.22e-16
  }
  return(Y)
}

# G --> P --> Y
.simCount.DM <- function(G.mat, SeqDepth){
  P = G.mat/rowSums(G.mat)
  P[P < 0] = 0
  
  Y = c()
  for (i in 1:nrow(G.mat)) {
    Y = rbind(Y, t(rmultinom(1, rpois(1, SeqDepth), P[i,]))) 
  }
  return(Y)
}

# .simCount.ZILN <- function(G.mat, SeqDepth){
#   P = G.mat/rowSums(G.mat)
#   P[P < 0] = 0
#   
#   Y = c()
#   for (i in 1:nrow(G.mat)) {
#     Y = rbind(Y, t(rmultinom(1, rpois(1, SeqDepth), P[i,]))) 
#   }
#   return(Y)
# }

# simulate correlated Multivariate normal
.cor2cov <- function(cor.mat, sd, discrepancy=1e-05){
  if (dim(cor.mat)[1] != dim(cor.mat)[2]) 
    stop("'cor.mat' should be a square matrix")
  n <- sqrt(length(cor.mat))
  if (n != length(sd)) 
    stop("The length of 'sd' should be the same as the number of rows of 'cor.mat'")
  if (length(sd[sd > 0]) != n) 
    stop("The elements in 'sd' shuold all be non-negative")
  if (isSymmetric(cor.mat)) 
    IS.symmetric <- TRUE
  else IS.symmetric <- FALSE
  p <- dim(cor.mat)[1]
  q <- p * (p - 1)/2
  if (isTRUE(all.equal(cor.mat[lower.tri(cor.mat)], 
                       rep(0, q))) || isTRUE(all.equal(cor.mat[upper.tri(cor.mat)], 
                                                       rep(0, q)))) 
    IS.triangular <- TRUE
  else IS.triangular <- FALSE
  if (!IS.symmetric & !IS.triangular) 
    stop("The object 'cor.mat' should be either a symmetric or a triangular matrix")
  cov.mat <- diag(sd) %*% cor.mat %*% diag(sd)
  colnames(cov.mat) <- rownames(cov.mat) <- colnames(cor.mat)
  return(cov.mat)
}  ## borrow code from library (MBESS)


.corrData.ZINormal <- function(zim, mm, varm, normal.corr, TT.corr){
  # zim = pstr.mat[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9,idx10,
  #                  idx11,idx12,idx13,idx14,idx15,idx16,idx17,idx18,idx19,idx20),]
  # mm = mu.mat[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9,idx10,
  #               idx11,idx12,idx13,idx14,idx15,idx16,idx17,idx18,idx19,idx20),]
  # varm = var.mat[idx1,]
  n.T = nrow(mm); K = ncol(mm) # K
  phi.mn = kronecker(TT.corr, normal.corr)
  phi.mn = .cor2cov(phi.mn, sqrt(kronecker(rep(1, n.T), varm)) )
  tmp = exp( mvrnorm(1, mu=as.numeric(t(mm)), Sigma=phi.mn )  )
  mln = matrix(tmp, nrow=n.T, byrow=TRUE)
  Delta = matrix(NA, n.T, K)
  for(j in 1:K){
    normalvars =  mvrnorm(1, mu=rep(0,n.T), Sigma=TT.corr)
    for(tt in 1:n.T){
      Delta[tt,j] = qbinom(pnorm(normalvars[tt]), 1, zim[tt,j])
    }
  }
  out = matrix(0, nrow=n.T, ncol=K)
  for(tt in 1:n.T){
    if(sum( 1-Delta[tt,])!=0 ){
      tmp = (1-Delta[tt,])*mln[tt,]
      out[tt,] = tmp/(sum(tmp)+1)
    }
  }
  return(cbind(out, 1-rowSums(out))) 
}

.simCount.ZILN <- function(mu.mat, var.mat, pstr.mat, normal.corr, n.T, r, SeqDepth){
  # var.TT = normal.var.TT; SeqDepth = 1000
  n = nrow(mu.mat)
  N = n/n.T
  
  TT.corr = matrix(r, nrow=n.T, ncol=n.T) + diag(n.T)*(1-r)
  
  idx1 = 1
  idx2 = N+1
  idx3 = 2*N+1
  idx4 = 3*N+1
  idx5 = 4*N+1
  
  Y1 = Y2 = Y3 = Y4 = Y5 = NULL
  for(i in 1:N){
    P = .corrData.ZINormal(pstr.mat[c(idx1,idx2,idx3,idx4,idx5),], 
                           mu.mat[c(idx1,idx2,idx3,idx4,idx5),], 
                           var.mat[idx1,], normal.corr, TT.corr) 
    P[P<0] = 0
    tmp = NULL
    for(t in 1:n.T){
      tmp = rbind(tmp, as.numeric( rmultinom(1, rpois(1, SeqDepth), P[t,] ) )  )
    }
    Y1 = rbind(Y1, tmp[1,]) 
    Y2 = rbind(Y2, tmp[2,])
    Y3 = rbind(Y3, tmp[3,])
    Y4 = rbind(Y4, tmp[4,])
    Y5 = rbind(Y5, tmp[5,])
    idx1 = idx1 + 1
    idx2 = idx2 + 1
    idx3 = idx3 + 1
    idx4 = idx4 + 1
    idx5 = idx5 + 1
  }
  Y = rbind(Y1, Y2, Y3, Y4, Y5)
  return(Y)  
}

.simCount.ZILN.diab <- function(mu.mat, var.mat, pstr.mat, normal.corr, n.T, r, SeqDepth){
  # normal.corr = cor.real; SeqDepth = 1000
  n = nrow(mu.mat)
  N = n/n.T
  
  TT.corr = matrix(r, nrow=n.T, ncol=n.T) + diag(n.T)*(1-r)
  
  idx1 = 1
  idx2 = N+1
  idx3 = 2*N+1
  idx4 = 3*N+1
  idx5 = 4*N+1
  idx6 = 5*N+1
  idx7 = 6*N+1
  idx8 = 7*N+1
  idx9 = 8*N+1
  idx10 = 9*N+1
  idx11 = 10*N+1
  idx12 = 11*N+1
  idx13 = 12*N+1
  idx14 = 13*N+1
  idx15 = 14*N+1
  idx16 = 15*N+1
  idx17 = 16*N+1
  idx18 = 17*N+1
  idx19 = 18*N+1
  idx20 = 19*N+1
  
  Y1 = Y2 = Y3 = Y4 = Y5 = Y6 = Y7 = Y8 = Y9 = Y10 = 
    Y11 = Y12 = Y13 = Y14 = Y15 = Y16 = Y17 = Y18 = Y19 = Y20 = NULL
  for(i in 1:N){
    P = .corrData.ZINormal(pstr.mat[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9,idx10,
                                      idx11,idx12,idx13,idx14,idx15,idx16,idx17,idx18,idx19,idx20),], 
                           mu.mat[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9,idx10,
                                    idx11,idx12,idx13,idx14,idx15,idx16,idx17,idx18,idx19,idx20),], 
                           var.mat[idx1,], normal.corr, TT.corr) 
    P[P<0] = 0
    tmp = NULL
    for(t in 1:n.T){
      tmp = rbind(tmp, as.numeric( rmultinom(1, rpois(1, SeqDepth), P[t,] ) )  )
    }
    Y1 = rbind(Y1, tmp[1,]) 
    Y2 = rbind(Y2, tmp[2,])
    Y3 = rbind(Y3, tmp[3,])
    Y4 = rbind(Y4, tmp[4,])
    Y5 = rbind(Y5, tmp[5,])
    Y6 = rbind(Y6, tmp[6,]) 
    Y7 = rbind(Y7, tmp[7,])
    Y8 = rbind(Y8, tmp[8,])
    Y9 = rbind(Y9, tmp[9,])
    Y10 = rbind(Y10, tmp[10,])
    Y11 = rbind(Y11, tmp[11,]) 
    Y12 = rbind(Y12, tmp[12,])
    Y13 = rbind(Y13, tmp[13,])
    Y14 = rbind(Y14, tmp[14,])
    Y15 = rbind(Y15, tmp[15,])
    Y16 = rbind(Y16, tmp[16,]) 
    Y17 = rbind(Y17, tmp[17,])
    Y18 = rbind(Y18, tmp[18,])
    Y19 = rbind(Y19, tmp[19,])
    Y20 = rbind(Y20, tmp[20,])
    idx1 = idx1 + 1
    idx2 = idx2 + 1
    idx3 = idx3 + 1
    idx4 = idx4 + 1
    idx5 = idx5 + 1
    idx6 = idx6 + 1
    idx7 = idx7 + 1
    idx8 = idx8 + 1
    idx9 = idx9 + 1
    idx10 = idx10 + 1
    idx11 = idx11 + 1
    idx12 = idx12 + 1
    idx13 = idx13 + 1
    idx14 = idx14 + 1
    idx15 = idx15 + 1
    idx16 = idx16 + 1
    idx17 = idx17 + 1
    idx18 = idx18 + 1
    idx19 = idx19 + 1
    idx20 = idx20 + 1
  }
  Y = rbind(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10,
            Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18, Y19, Y20)
  return(Y)  
}
AIGDM_GEE.Cluster <- function(ID, Y, Xa, Xb, X.index, W = NULL, model = "AZIGDM", zi.id = NULL, 
                              n.boot = NULL, fdr.alpha = 0.05, perm.type = NULL, n.perm = NULL, seed = 123){
  # Y = Y.m; Xa = Xb = W = xmat.m; X.index = 1; model = "GDM";
  # n.boot = NULL; fdr.alpha = 0.05; perm.type = "BTW"; n.perm = NULL; seed = 12+sim
  # model could be GDM/ZIGDM/AZIGDM/AIGDM
  
  Xa = cbind(1, Xa) # add the intercept term
  Xb = cbind(1, Xb) # add the intercept term
  W = cbind(1, W)
  X.index = X.index + 1
  Xa.r = Xa[, -X.index, drop=FALSE]
  Xb.r = Xb[, -X.index, drop=FALSE]
  W.r = W[, -X.index, drop=FALSE]
  
  K = ncol(Y) - 1
  R.sel = .choose_r(fdr.alpha, 0.05)
  
  if(model == "GDM"){
    zi.check.pval.all.wo = rep(1, K)
  }else if(model == "ZIGDM"){
    zi.check.pval.all.wo = rep(0, K)
  }else if(model == "AZIGDM"){
    print("Zero-inflation diagnostic test!")
    # TRUE means add zero-inflation
    zi.check.pval.all.wo = sapply(1:K, function(x){
      y.tmp = cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE]))
      id.tmp = ID[rowSums(y.tmp) != 0]
      y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]
      return(.check_zeroinflation(id.tmp, y.tmp, n.boot = n.boot, R.sel, type = "cluster"))
    })
  }else if(model == "AZIGDM-gold"){
    zi.check.pval.all.wo = rep(1, K)
    if(!is.null(zi.id))
      zi.check.pval.all.wo[zi.id] = 0
  }
  
  zi.check.pval.all.w = p.adjust(zi.check.pval.all.wo, method = "BH")
  zi.check = zi.check.pval.all.w < fdr.alpha
  zi.id = which(zi.check)
  
  est.para.lst = lapply(1:K, function(x){
    .est_para(ID, W.r, Xa.r, Xb.r, cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE])), zi.check[x])
  })
  
  # est.para.mean.lst = lapply(1:K, function(x){
  #   .est_para(ID, W, Xa.r, Xb, cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE])), zi.check[x])
  # })
  # est.para.disp.lst = lapply(1:K, function(x){
  #   .est_para(ID, W, Xa, Xb.r, cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE])), zi.check[x])
  # })
  
  asym.mean = .score_test_mean(ID, Xa, X.index, Y, est.para.lst)
  stat.mean.sum = asym.mean$stat.sum
  pval.mean.a = asym.mean$pval
  
  asym.disp = .score_test_disp(ID, Xb, X.index, Y, est.para.lst)
  stat.disp.sum = asym.disp$stat.sum
  pval.disp.a = asym.disp$pval
  
  if(length(zi.id) > 0){
    asym.zero = .score_test_zero(ID, W, X.index, Y, est.para.lst, zi.id)
    stat.zero.sum = asym.zero$stat.sum
    pval.zero.a = asym.zero$pval

    pval.omni.c.a = .ACAT(c(pval.zero.a, pval.mean.a, pval.disp.a))
    # pval.omni.m.a = 1 - pchisq(sum(stat.zero.sum, stat.mean.sum, stat.disp.sum), length(zi.id) + K*2)

  }else{
    pval.zero.a = NA
    pval.omni.c.a = .ACAT(c(pval.mean.a, pval.disp.a))
    # pval.omni.m.a = 1 - pchisq(sum(stat.mean.sum, stat.disp.sum), K*2)
  }
  pval.score = c(pval.zero.a, pval.mean.a, pval.disp.a, pval.omni.c.a)
  # names(pval.score) = c("Score-Zero-Asymptotic", "Score-Mean-Asymptotic", "Score-Disp-Asymptotic", "Score-Omni-Asymptotic")
  
  if(!is.null(n.perm)){
    pval.mean.p = .score_test_perm(ID, Xa, X.index, Y, est.para.lst, stat.mean.sum, test = "mean", n.perm, perm.type, R.sel, NULL, seed)
    pval.disp.p = .score_test_perm(ID, Xb, X.index, Y, est.para.lst, stat.disp.sum, test = "disp", n.perm, perm.type, R.sel, NULL, seed)
    if(length(zi.id) > 0){
      pval.zero.p = .score_test_perm(ID, W, X.index, Y, est.para.lst, stat.zero.sum, test = "zero", n.perm, perm.type, R.sel, zi.id, seed)
      pval.omni.c.p = .ACAT(c(pval.zero.p, pval.mean.p, pval.disp.p))
      # pval.omni.m.p = .score_test_perm_omni(ID, W, X.index, Y, est.para.lst, sum(stat.zero.sum, stat.mean.sum, stat.disp.sum), 
      #                                       n.perm, perm.type, R.sel, zi.id, seed)
      
    }else{
      pval.zero.p = NA
      pval.omni.c.p = .ACAT(c(pval.mean.p, pval.disp.p))
      # pval.omni.m.p = .score_test_perm_omni(ID, W, X.index, Y, est.para.lst, sum(stat.mean.sum, stat.disp.sum), 
      #                                       n.perm, perm.type, R.sel, NULL, seed)
    }
    pval.score = c(pval.zero.a, pval.zero.p, pval.mean.a, pval.mean.p, pval.disp.a, pval.disp.p, pval.omni.c.a, pval.omni.c.p)
    # names(pval.score) = c("Score-Mean-Asymptotic", "Score-Mean-Resampling", "Score-Disp-Asymptotic", "Score-Disp-Resampling")
  }
  
  # est.para.lst = lapply(1:K, function(x){
  #   .est_para(ID, W.r, Xa, Xb, cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE])),
  #             zi.check[x], oi.check[x])
  # })
  # pval.meanNdisp.sand = .wald_test_sand(ID, X.index, Y, est.para.lst, zi.check, oi.check)
  # pval.wald = c(pval.meanNdisp.sand$pval.mean, pval.meanNdisp.sand$pval.disp)
  # names(pval.wald) = c("Wald-Mean-Asymptotic", "Wald-Disp-Asymptotic")
  # pval.meanNdisp.boot = .wald_test_boot(ID, X.index, W.r, Xa, Xb, Y, est.para.lst, zi.check, oi.check)
  # pval.wald = c(pval.meanNdisp.sand$pval.mean, pval.meanNdisp.boot$pval.mean,
  #               pval.meanNdisp.sand$pval.disp, pval.meanNdisp.boot$pval.disp)
  # names(pval.wald) = c("Wald-Mean-Asymptotic", "Wald-Mean-Bootstrap", "Wald-Disp-Asymptotic", "Wald-Disp-Bootstrap")
  # pval = c(pval.wald, pval.score)
  pval = pval.score
  return(list(pval = pval, zi = zi.check, zi.pval = zi.check.pval.all.wo))
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
  Logit.par.ini = gamma.ini
  Logit.data = list(Del=Del, W=W)
  return(optim(par=Logit.par.ini, fn=.LogitNegLoglik, gr=.LogitNegScore, data = Logit.data, method = "L-BFGS-B")$par) 
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
  return(optim(par=Beta.par.ini, fn = .AIBetaNegLoglik, gr = .AIBetaNegScore, data = Beta.data, method="BFGS")$par)
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
#   phi.tmp = as.numeric( tmp.b/(1+tmp.b) )
#   
#   a = (1/phi.tmp - 1) * mu.tmp
#   b = (1/phi.tmp - 1) * (1 - mu.tmp)
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
# # .dma: dmu/dalpha, first derivative of alpha in the mean model 
# .dma <- function(xmat, alpha){ # xmat 2*2 alpha 2*5 
#   K = ncol(alpha)
#   tmp = c(exp(xmat %*% alpha)) # taxa1.t1 taxa1.t2 taxa2.t1 taxa2.t2 ... taxa5.t1 taxa5.t2
#   xmat_expand = do.call(rbind, replicate(K, xmat, simplify = FALSE))
#   tmp1 = (1 + tmp)^2
#   return(.RowbyRow(xmat_expand, tmp/tmp1))
# }
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
# # .var_Z: variance of Z in mean model (matrix Bi=Diag(mu_i*(1-mu_i)/(phi_i+1)))
# .var_Z <- function(mu, phi){
#   return(diag(mu*(1-mu)*phi, length(mu)))
# }
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

.check_zeroinflation <- function(ID, Y, n.boot = 9999, R.sel = R.sel, type = "cluster"){
  Y.seqDepth = rowSums(Y); n = row(Y)
  obs.ans = .BBnZIBB.log.lik(Y)
  stat = -2 * (obs.ans$bb.ans - obs.ans$zibb.ans)
  if (is.infinite(stat)) {
    cat("Use Pseudo-ECDF approximation p-value\n")
    return(1/(n.boot+1))
  }
  
  m = 0
  Nexc = 0
  stat.perm = numeric(n.boot)
  if(type == "cluster"){
    id.unique = unique(ID)
    # Create an index list for each unique ID
    indexList = split(1:nrow(Y), ID)
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
  }else if(type == "simple"){
    while (all(Nexc < R.sel, m < n.boot)) {
      Y.b.tmp = Y[sample(n, replace = TRUE),]
      
      bb.ll = .BB.log.lik(Y.b.tmp)
      a.mat = matrix(colMeans(bb.ll$a.mat), nrow(Y), 1); b.mat = matrix( colMeans(bb.ll$b.mat), nrow(Y), 1)
      
      Y.b = .simData.GDM(a.mat, b.mat, Y.seqDepth)
      boot.obs = .BBnZIBB.log.lik(Y.b)
      current_stat = -2 * (boot.obs$bb.ans - boot.obs$zibb.ans)
      Nexc = Nexc + sum(current_stat >= stat)
      stat.perm[m+1] = current_stat
      m = m + 1
    }
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

.check_oneinflation <- function(ID, Y, n.boot = 9999, R.sel = R.sel){
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
      # Y.b.tmp = Y
      id.unique = unique(ID)
      ID.tmp = sort(sample(id.unique, replace = T))
      Y.b.tmp = NULL
      for (id in ID.tmp) {
        Y.b.tmp = rbind(Y.b.tmp, Y[which(ID == id),,drop=FALSE])
      }
      bb.ll = .BB.log.lik(Y.b.tmp)
      a = colMeans(bb.ll$a.mat); b = colMeans(bb.ll$b.mat)
      a.mat = matrix(a, nrow(Y), 1); b.mat = matrix(b, nrow(Y), 1)
      m = 0
      Nexc = 0
      stat.perm = numeric(n.boot)
      while (all(Nexc < R.sel, m < n.boot)) {
        Y.b = .simData.GDM(a.mat, b.mat, Y.seqDepth)
        bb.ll = .BB.log.lik(Y.b)
        oibb.ll = .OIBB.log.lik(Y.b)
        stat.perm.subset = -2 * (bb.ll$ans - oibb.ll$ans)
        Nexc = Nexc + sum(stat.perm.subset >= stat)
        stat.perm[m+1] = stat.perm.subset
        m = m + 1
      }
      
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


.est_para_pl <- function(ID, W, Xa, Xb, Y, zi.check, oi.check){
  #W=W.r; Y=cbind(Y[,1], rowSums(Y[,(1+1):(5+1),drop=FALSE]));zi.check=zi.check[1];oi.check=oi.check[1]
  da = ncol(Xa); db = ncol(Xb); dw = ncol(W); K = ncol(Y)-1
  
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
  nt = nrow(Y) # the total number of observed units
  N1 = sum(n.reps*(n.reps-1)/2)
  
  r.step = 0.1 # set up the tuning parameter for updating alpha
  max = 1e3 # maximal iterations
  
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
  rho0 = sig0 = tau0 = 0.8
  
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
        par.post = .ZOIeZ(pZ0.mat[t,,drop=FALSE], pO0.mat[t,,drop=FALSE], av[t,,drop=FALSE], bv[t,,drop=FALSE], y.tmp[t,])
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
    
    sig.new = tau.new = rho.new = 0.8
    
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
      par.post = .ZOIeZ(pZ0.mat[t,,drop=FALSE], pO0.mat[t,,drop=FALSE], av[t,,drop=FALSE], bv[t,,drop=FALSE], y.tmp[t,])
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

.est_para_em <- function(ID, W, Xa, Xb, Y, zi.check, oi.check){
  # W = W.r; Xa = Xa.r; Xb = Xb.r; Y=cbind(Y[,1], rowSums(Y[,(1+1):(K+1),drop=FALSE])); zi.check=TRUE; oi.check=FALSE
  da = ncol(Xa); db = ncol(Xb); dw = ncol(W); K = ncol(Y)-1
  
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
  
  
  aigdm.reg = .AIGDM_EM(Y = Y, W = W, Xa = Xa, Xb = Xb, 
                        gamma0 = gamma0, omega0 = omega0, alpha0 = alpha0, beta0 = beta0)
  gamma0 = aigdm.reg$gamma.est
  omega0 = aigdm.reg$omega.est
  alpha0 = aigdm.reg$alpha.est
  beta0 = aigdm.reg$beta.est
  
  #print(paste0("Converge?", ifelse(l == max, 0, 1)))
  return(list(gamma0 = gamma0, omega0 = omega0, 
              alpha0 = alpha0, beta0 = beta0))
}

.wald_test_sand <- function(ID, X.index, Y, est.para.lst, zi.check, oi.check){
  K = ncol(Y)-1
  ai.check = zi.check + oi.check
  for (j in 1:K) {
    Y.tmp = cbind(Y[,j], rowSums(Y[,(j+1):(K+1),drop=FALSE]))
    
    est.para = est.para.lst[[j]]
    ID.tmp = ID; W.r.tmp = est.para$W; Xa.tmp = est.para$Xa; Xb.tmp = est.para$Xb
    dw.r = ncol(W.r.tmp); da = ncol(Xa.tmp); db = ncol(Xb.tmp)
    
    id.unique = sort(unique(ID.tmp)); N = length(id.unique); n.reps = as.vector(table(ID.tmp))
    ng = no = dw.r; na = da; nb = db
    
    if(ai.check[j]){
      if(!zi.check[j]){
        para.ind.alpha = sort(sapply(X.index, function(x) seq(x, na, na))) + no
        para.ind.beta = sort(sapply(X.index, function(x) seq(x, nb, nb))) + no + na
        para.keep = c(sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }else if(!oi.check[j]){
        para.ind.alpha = sort(sapply(X.index, function(x) seq(x, na, na))) + ng
        para.ind.beta = sort(sapply(X.index, function(x) seq(x, nb, nb))) + ng + na
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }else{
        para.ind.alpha = sort(sapply(X.index, function(x) seq(x, na, na))) + ng + no
        para.ind.beta = sort(sapply(X.index, function(x) seq(x, nb, nb))) + ng + no + na
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }
    }else{
      para.ind.alpha = sort(sapply(X.index, function(x) seq(x, na, na)))
      para.ind.beta = sort(sapply(X.index, function(x) seq(x, nb, nb))) + na
      para.keep = c(sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                    sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
    }
    
    gamma0 = est.para$gamma0; omega0 = est.para$omega0; alpha0 = est.para$alpha0; beta0 = est.para$beta0
    sig0 = est.para$sig0; tau0 = est.para$tau0; rho0 = est.para$rho0
    eDel0.lst = est.para$eDel0.lst; eDelZ0.lst = est.para$eDelZ0.lst; eDelO0.lst = est.para$eDelO0.lst
    eZStar.lst = est.para$eZStar.lst; muStar.lst = est.para$muStar.lst
    A.R.lst = est.para$A.R.lst; B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst; B2.R.lst = est.para$B2.R.lst; AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst; pO0.lst = est.para$pO0.lst
    eZ0.lst = est.para$eZ0.lst; mu0.lst = est.para$mu0.lst; phi0.lst = est.para$phi0.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    
    # update DVD matrix with full models
    U = 0
    I = 0
    M = matrix(0, ng+no+na+nb, ng+no+na+nb)
    ind.g = 1:ng; ind.o = (ng+1):(ng+no); ind.a = (ng+no+1):(ng+no+na); ind.b = (ng+no+na+1):(ng+no+na+nb)
    for(i in 1:N){
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel,,drop=FALSE]
      xa.tmp = Xa.tmp[ID.sel,,drop=FALSE]
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
      eZStar = eZStar.lst[[i]]; muStar = muStar.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      tmp.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0))
      dvd.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0)) %*% .BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZStar - muStar), 2, sum))
      
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.tmp)
      # dvd.b.x = - t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)) %*% t(t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)))
      geeb.x = as.matrix(pseudoll.score.beta)
      
      ES = rbind(geeg.w, geeo.w, geea.x, geeb.x)
      U = U + ES
      # EI = .BlockDiagBind(dvd.g.w, dvd.o.w, dvd.a.x, dvd.b.x)
      # I = I + (-EI)
      
      a = (1/phi0 - 1) * mu0; a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0); b[b < 0] = 0
      A = digamma(a) - digamma(a+b)
      B = digamma(b) - digamma(a+b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1-eDel0)*( - A + A.post )
      two.R = (1-eDel0)*( - B + B.post )
      tmp.dvd = .DVD(w.r.tmp, xa.tmp, xb.tmp, gamma0, omega0, alpha0, beta0, pZ0, pO0, mu0, phi0, one.R, two.R, 1, n.T)
      tmp.Kab = matrix(0, nrow(tmp.dvd$Kab), ncol(tmp.dvd$Kab))
      dvd.ab.all = rbind(cbind(dvd.a.x, tmp.Kab), cbind(t(tmp.dvd$Kab), tmp.dvd$Kbb))
      EI = .BlockDiagBind(dvd.g.w, dvd.o.w, dvd.ab.all)
      I = I + (-EI)
      
      M[ind.g, ind.g] = M[ind.g, ind.g] + geeg.w %*% t(geeg.w) * (diag(1) %x% matrix(1, ng, ng))
      M[ind.g, ind.o] = M[ind.g, ind.o] + geeg.w %*% t(geeo.w) * (diag(1) %x% matrix(1, ng, no))
      M[ind.o, ind.o] = M[ind.o, ind.o] + geeo.w %*% t(geeo.w) * (diag(1) %x% matrix(1, no, no))
      M[ind.o, ind.g] = M[ind.o, ind.g] + geeo.w %*% t(geeg.w) * (diag(1) %x% matrix(1, no, ng))
      
      M[ind.a, ind.a] = M[ind.a, ind.a] + geea.x %*% t(geea.x) * (diag(1) %x% matrix(1, na, na))
      M[ind.b, ind.b] = M[ind.b, ind.b] + geeb.x %*% t(geeb.x) * (diag(1) %x% matrix(1, nb, nb))
      
      M[ind.a, ind.b] = M[ind.a, ind.b] + geea.x %*% t(geeb.x) * (diag(1) %x% matrix(1, na, nb))
      M[ind.b, ind.a] = M[ind.a, ind.b] + geeb.x %*% t(geea.x) * (diag(1) %x% matrix(1, nb, na))
      
      M[c(ind.g, ind.o), c(ind.a, ind.b)] = M[c(ind.g, ind.o), c(ind.a, ind.b)] + rbind(geeg.w, geeo.w) %*% t(rbind(geea.x, geeb.x)) * (diag(1) %x% matrix(1, ng+no, na+nb))
      M[c(ind.a, ind.b), c(ind.g, ind.o)] = M[c(ind.a, ind.b), c(ind.g, ind.o)] + rbind(geea.x, geeb.x) %*% t(rbind(geeg.w, geeo.w)) * (diag(1) %x% matrix(1, na+nb, ng+no))
    }
    
    I.r = I[para.keep, para.keep]
    M.r = M[para.keep, para.keep]
    H.r = ginv(I.r)
    
    if(j == 1){
      alpha.est = alpha0[X.index,,drop=FALSE]
      beta.est = beta0[X.index,,drop=FALSE]
      sand.cov = H.r %*% M.r %*% t(H.r)
      para.alpha = para.ind.alpha
      para.beta = para.ind.beta
      para.keep.all = para.keep
    }else{
      alpha.est = rbind(alpha.est, alpha0[X.index,,drop=FALSE])
      beta.est = rbind(beta.est, beta0[X.index,,drop=FALSE])
      sand.cov = .BlockDiagBind(sand.cov, H.r %*% M.r %*% t(H.r))
      para.alpha = c(para.alpha, para.ind.alpha + length(para.keep.all))
      para.beta = c(para.beta, para.ind.beta + length(para.keep.all))
      para.keep.all = c(para.keep.all, para.keep)
    }
  }
  stat.alpha = c(t(alpha.est) %*% ginv(sand.cov[para.alpha, para.alpha]) %*% alpha.est)
  pval.alpha = 1 - pchisq(stat.alpha, nrow(alpha.est))
  # pval.alpha = 1 - pnorm(abs(alpha.est)/sqrt(diag(sand.cov)[para.alpha]))
  se.alpha.full = sqrt(diag(sand.cov[sort(c(para.alpha-1,para.alpha)), sort(c(para.alpha-1,para.alpha))]))
  stat.beta = c(t(beta.est) %*% ginv(sand.cov[para.beta, para.beta]) %*% beta.est)
  pval.beta = 1 - pchisq(stat.beta, nrow(beta.est))
  # pval.beta = 1 - pnorm(abs(beta.est)/sqrt(diag(sand.cov)[para.beta]))
  se.beta.full = sqrt(diag(sand.cov[sort(c(para.beta-1,para.beta)), sort(c(para.beta-1,para.beta))]))
  # stat.alpha = do.call(cbind, alpha.est) %*% ginv(.BlockDiagBind_list(sand.cov.alpha)) %*% t(do.call(cbind, alpha.est))
  # pval.alpha = 1 - pchisq(stat.alpha, K)
  return(list(pval.mean = pval.alpha, pval.disp = pval.beta,
              se.mean = se.alpha.full, se.disp = se.beta.full))
}

.wald_test_boot <- function(ID, X.index, W.r, Xa, Xb, Y, est.para.lst, zi.check, oi.check, n.boot = 200){
  K = ncol(Y) - 1
  da = ncol(Xa); db = ncol(Xb)
  
  est.para = do.call(Map, c(f = cbind, est.para.lst))
  est.es.est = c(est.para$alpha0, est.para$beta0)
  
  para.ind.alpha = sort(sapply(X.index, function(x) seq(x, da*K, da)))
  para.ind.beta = sort(sapply(X.index, function(x) seq(x, db*K, db))) + da*K
  
  id.unique = sort(unique(ID))
  tmp.subset.tmp1 = tmp.subset.tmp2 = vector("list", n.boot)
  for (i in 1:n.boot) {
    ID.tmp = sort(sample(id.unique, replace = T))
    print(ID.tmp)
    Y.tmp = Xa.tmp = Xb.tmp = NULL
    for (id in ID.tmp) {
      Y.tmp = rbind(Y.tmp, Y[which(ID == id),])
      Xa.tmp = rbind(Xa.tmp, Xa[which(ID == id),])
      Xb.tmp = rbind(Xb.tmp, Xb[which(ID == id),])
    }
    
    est.para.bs.lst = lapply(1:K, function(x){
      if(sum(Y.tmp[,x]==0) < 5){
        zi.check.tmp = FALSE
      }else{
        zi.check.tmp = zi.check[x]
      }
      .est_para(ID.tmp, W.r, Xa.tmp, Xb.tmp, cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])),
                zi.check.tmp, oi.check[x])
    })
    est.para.bs.tmp1 = do.call(Map, c(f = cbind, est.para.bs.lst))
    tmp.subset.tmp1[[i]] = c(est.para.bs.tmp1$alpha0, est.para.bs.tmp1$beta0)
    
    ID.tmp = sample(1:length(ID), replace = T)
    Y.tmp = Y[ID.tmp,]
    Xa.tmp = Xa[ID.tmp,]
    Xb.tmp = Xb[ID.tmp,]
    
    est.para.bs.lst = lapply(1:K, function(x){
      if(sum(Y.tmp[,x]==0) < 5){
        zi.check.tmp = FALSE
      }else{
        zi.check.tmp = zi.check[x]
      }
      .est_para(ID.tmp, W.r, Xa.tmp, Xb.tmp, cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])),
                zi.check.tmp, oi.check[x])
    })
    est.para.bs.tmp2 = do.call(Map, c(f = cbind, est.para.bs.lst))
    tmp.subset.tmp2[[i]] = c(est.para.bs.tmp2$alpha0, est.para.bs.tmp2$beta0)
  }
  
  bspara.mat = matrix(unlist(tmp.subset.tmp1), n.boot, length(est.es.est), byrow = T)
  est.es.boot.se = colSds(bspara.mat, na.rm = T)
  est.es.boot.cov = cov(bspara.mat)
  est.es.boot.cov.alpha = diag(diag(est.es.boot.cov[para.ind.alpha, para.ind.alpha, drop=FALSE]), length(para.ind.alpha))
  est.es.boot.cov.beta = diag(diag(est.es.boot.cov[para.ind.beta, para.ind.beta, drop=FALSE]), length(para.ind.beta))
  
  stat.alpha = est.es.est[para.ind.alpha] %*% ginv(est.es.boot.cov.alpha) %*% est.es.est[para.ind.alpha]
  pval.alpha.tmp1 = 1 - pchisq(stat.alpha, K)
  se.alpha.full.tmp1 = est.es.boot.se[sort(c(para.ind.alpha-1,para.ind.alpha))]
  stat.beta = est.es.est[para.ind.beta] %*% ginv(est.es.boot.cov.beta) %*% est.es.est[para.ind.beta]
  pval.beta.tmp1 = 1 - pchisq(stat.beta, K)
  se.beta.full.tmp1 = est.es.boot.se[sort(c(para.ind.beta-1,para.ind.beta))]
  
  bspara.mat = matrix(unlist(tmp.subset.tmp2), n.boot, length(est.es.est), byrow = T)
  est.es.boot.se = colSds(bspara.mat, na.rm = T)
  est.es.boot.cov = cov(bspara.mat)
  est.es.boot.cov.alpha = diag(diag(est.es.boot.cov[para.ind.alpha, para.ind.alpha, drop=FALSE]), length(para.ind.alpha))
  est.es.boot.cov.beta = diag(diag(est.es.boot.cov[para.ind.beta, para.ind.beta, drop=FALSE]), length(para.ind.beta))
  
  stat.alpha = est.es.est[para.ind.alpha] %*% ginv(est.es.boot.cov.alpha) %*% est.es.est[para.ind.alpha]
  pval.alpha.tmp2 = 1 - pchisq(stat.alpha, K)
  se.alpha.full.tmp2 = est.es.boot.se[sort(c(para.ind.alpha-1,para.ind.alpha))]
  stat.beta = est.es.est[para.ind.beta] %*% ginv(est.es.boot.cov.beta) %*% est.es.est[para.ind.beta]
  pval.beta.tmp2 = 1 - pchisq(stat.beta, K)
  se.beta.full.tmp2 = est.es.boot.se[sort(c(para.ind.beta-1,para.ind.beta))]
  return(list(pval.mean = c(pval.alpha.tmp1, pval.alpha.tmp2), pval.disp = c(pval.beta.tmp1, pval.beta.tmp2),
              se.mean = c(se.alpha.full.tmp1, se.alpha.full.tmp2), 
              se.disp = c(se.beta.full.tmp1, se.beta.full.tmp2)))
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
  # X = Xa; est.para.lst=est.para.mean.lst
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

.score_test_perm <- function(ID, X, X.index, Y, est.para.lst, stat.sum, test, 
                             n.perm, perm.type, R.sel, zi.id = NULL, seed){
  #Y = Y.tmp;
  
  m = 0; Nexc = 0; stat.perm = numeric(n.perm)
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
  
  while(all(Nexc < R.sel, m < n.perm)){
    set.seed(seed+m)
    if(perm.type == "WTH"){
      var.new = .shuffleInCluster(X[,X.index,drop=FALSE], ID)
      X[,X.index] = var.new
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
      X[,X.index] = X.p.int
      X.p = X
    }
    if(test == "zero")
      tmp = .score_test_zero(ID, X.p, X.index, Y, est.para.lst, zi.id)
    if(test == "mean")
      tmp = .score_test_mean(ID, X.p, X.index, Y, est.para.lst)
    if(test == "disp")
      tmp = .score_test_disp(ID, X.p, X.index, Y, est.para.lst)
    tmp.subset = tmp$stat.sum
    Nexc = Nexc + sum(tmp.subset >= stat.sum) # tmp.subset >= stat.sum
    stat.perm[m+1] = tmp.subset
    m = m + 1
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

.score_test_perm_omni <- function(ID, X, X.index, Y, est.para.lst, stat.sum, n.perm, perm.type, R.sel, zi.id = NULL, seed){
  #Y = Y.tmp;
  
  m = 0; Nexc = 0; stat.perm = numeric(n.perm)
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
  
  while(all(Nexc < R.sel, m < n.perm)){
    set.seed(seed+m)
    if(perm.type == "WTH"){
      var.new = .shuffleInCluster(X[,X.index,drop=FALSE], ID)
      X[,X.index] = var.new
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
      X[,X.index] = X.p.int
      X.p = X
    }
    
    tmp.mean = .score_test_mean(ID, X.p, X.index, Y, est.para.lst)
    tmp.disp = .score_test_disp(ID, X.p, X.index, Y, est.para.lst)
    if(!is.null(zi.id)){
      tmp.zero = .score_test_zero(ID, X.p, X.index, Y, est.para.lst, zi.id)
      tmp.subset = sum(tmp.zero$stat.sum, tmp.mean$stat.sum, tmp.disp$stat.sum)
    }else{
      tmp.subset = sum(tmp.mean$stat.sum, tmp.disp$stat.sum)
    }

    Nexc = Nexc + sum(tmp.subset >= stat.sum) # tmp.subset >= stat.sum
    stat.perm[m+1] = tmp.subset
    m = m + 1
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
.simData.GDM <- function(a, b, SeqDepth){
  N = nrow(a); K = ncol(a)
  Z = matrix(rbeta(N*K, a, b), N, K)
  P = matrix(NA, N, K+1)
  for (j in 1:K) {
    P[,j] = Z[,j]*matrixStats::rowProds(1 - Z[,0:(j-1), drop = FALSE])
  }
  P[, ncol(P)] =  1 - rowSums(P[, -ncol(P), drop = FALSE])
  P[P < 0] = 0
  Y = c()
  for (i in 1:N) {
    Y = rbind(Y, t(rmultinom(1, SeqDepth[i], P[i,])))
  }
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
  phi = g * exp(w[1] + g * log(k / n))
  return(c(g, phi))
}

# Combined estimators 
.combined_method <- function(x){
  m = mean(x)
  maxi = max(x)
  g = m / (m - maxi) 
  phi = - g * maxi     
  return(c(g, phi))
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

DM_resample <- function(ID, Y, case, perm.type = NULL, fdr.alpha = 0.05, n.perm = NULL){
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID)); n = nrow(Y)
  
  group.data = NULL
  for(a.case in unique(case) ){
    group.data = c(group.data, list(Y[case==a.case,,drop=FALSE]) )
  }
  # dm.stat.m = tryCatch( Xmcupo.sevsample(group.data)$`Xmcupo statistics`, error=function(err) NA )
  # dm.stat.o = tryCatch( Xoc.sevsample(group.data)$`Xoc statistics`, error=function(err) NA )
  dm.stat.d = tryCatch( Xdc.sevsample(group.data)$`Xdc statistics`, error=function(err) NA )
  pval = tryCatch( Xdc.sevsample(group.data)$`p value`, error=function(err) NA )
  R.sel = .choose_r(fdr.alpha, 0.05)
  
  # if(length(dm.stat.m) != 1)
  #   dm.stat.m = NA
  # if(!is.na(dm.stat.m)){
  #   pval.m = .dm_perm(ID, Y, case, dm.stat.m, test = "mean", "BTW", n.perm, R.sel)
  # }else{
  #   pval.m = NA
  # }
  # 
  # if(length(dm.stat.o) != 1)
  #   dm.stat.o = NA
  # if(!is.na(dm.stat.o)){
  #   pval.o = .dm_perm(ID, Y, case, dm.stat.o, test = "disp", "BTW", n.perm, R.sel)
  # }else{
  #   pval.o = NA
  # }
  
  if(length(dm.stat.d) != 1)
    dm.stat.d = NA
  if(!is.na(dm.stat.d)){
    pval.d = .dm_perm(ID, Y, case, dm.stat.d, test = "omni", "BTW", n.perm, R.sel)
  }else{
    pval.d = NA
  }
  
  # return(list(pval = c(pval.m, pval.o, pval.d)))
  return(list(pval = c(pval, pval.d)))
  
}

.dm_perm <- function(ID, Y, case, stat, test, perm.type, n.perm, R.sel){
  id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID)); n = nrow(Y)
  m = 0; Nexc = 0; stat.perm = numeric(n.perm)
  while (all(Nexc < R.sel, m < n.perm)) {
    if(perm.type == "BTW"){
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
    }else if(perm.type=="WTH"){
      sample.idx = .shuffleInCluster(1:n, ID)
      case.perm = case[sample.idx]
    }
    group.data.perm = NULL
    for(a.case in unique(case)){
      group.data.perm = c(group.data.perm, list(Y[case.perm==a.case,,drop=FALSE]) )
    }
    
    
    if(test == "mean"){
      tmp.perm = tryCatch( Xmcupo.sevsample(group.data.perm)$`Xmcupo statistics`,  error=function(err) NA)
      if(length(tmp.perm) != 1)
        tmp.perm = NA
      if(!is.na(tmp.perm) ){
        Nexc = Nexc + sum(tmp.perm >= stat)
      }
      stat.perm[m+1] = tmp.perm
    }
    
    if(test == "disp"){
      tmp.perm = .tryCatch.Xoc.sevsample(group.data.perm)
      if(length(tmp.perm) != 1)
        tmp.perm = NA
      if(!is.na(tmp.perm) ){
        Nexc = Nexc + sum(tmp.perm >= stat)
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
  setTimeLimit(cpu = 3000, elapsed = 3000, transient = TRUE)
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

.tryCatch.Xoc.sevsample.pvalue <- function(x){
  setTimeLimit(cpu = 60, elapsed = 60, transient = TRUE)
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  })
  
  tryCatch({
    # do some stuff
    Xoc.sevsample(x)$`p value`
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
# .AIGDM_MPL <- function(ID, Y){
#   # Y = Y.m
#   K = ncol(Y) - 1
#   W = Xa = Xb = matrix(1, nrow(Y.m), 1)
#   gamma0 = omega0 = matrix(-Inf, 1, K)
#   alpha0 = beta0 = matrix(0.001, 1, K)
#   tol = 1e-4; max.iter = 1000
#   
#   # CONV = 0
#   # CONV.iter = max.iter
#   n = nrow(Y)
#   K = ncol(Y) - 1
#   da = ncol(Xa); db = ncol(Xb)
#   dw = ncol(W)
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
#   
#   id.unique = sort(unique(ID)); N = length(id.unique); n.reps = as.vector(table(ID))
#   alpha.se = beta.se = numeric(K)
#   for (j in 1:K) {
#     # j = 1
#     # gamma.now = est.para.es$gamma0; omega.now = est.para.es$omega0
#     # alpha.now = est.para.es$alpha0; beta.now = est.para.es$beta0
#     gamma0 = gamma.now[,j,drop=FALSE]; omega0 = omega.now[,j,drop=FALSE]
#     alpha0 = alpha.now[,j,drop=FALSE]; beta0 = beta.now[,j,drop=FALSE]
#     U = I = M = 0
#     
#     for (i in 1:N) {
#       # i = 1
#       ID.sel = which(ID == id.unique[i])
#       w.tmp = W[ID.sel,,drop=FALSE]
#       xa.tmp = Xa[ID.sel,,drop=FALSE]
#       xb.tmp = Xb[ID.sel,,drop=FALSE]
#       n.T = n.reps[i]
#       
#       tmp.g = as.numeric(exp(w.tmp %*% gamma0))
#       pZ0 = tmp.g/(1+tmp.g)
#       score.g = colSums((DelZ.R[ID.sel,j] - pZ0) * w.tmp)
#       
#       tmp.o = as.numeric(exp(w.tmp %*% omega0))
#       pO0 = tmp.o/(1+tmp.o)
#       score.o = colSums((DelO.R[ID.sel,j] - pO0) * w.tmp)
#       
#       tmp.a = as.numeric(exp(xa.tmp %*% alpha0))
#       mu0 = as.numeric( tmp.a/(1+tmp.a) )
#       
#       tmp.b = as.numeric(exp(xb.tmp %*% beta0))
#       phi0 = as.numeric( tmp.b/(1+tmp.b) )
#       
#       a = (1/phi0 - 1) * mu0
#       b = (1/phi0 - 1) * (1 - mu0)
#       a[a < 0] = 0
#       b[b < 0] = 0
#       
#       A = digamma(a) - digamma(a+b)
#       B = digamma(b) - digamma(a+b)
#       Del.post = Del.R[ID.sel,j]
#       A.post = A.R[ID.sel,j]
#       B.post = B.R[ID.sel,j]
#       one.R = (1-Del.post)*( - A + A.post )
#       two.R = (1-Del.post)*( - B + B.post )
#       
#       zstar =  A.post - B.post
#       zdagger = B.post
#       mustar = digamma(a) - digamma(b)
#       mudagger = digamma(b) - digamma(a + b)
#       
#       score.a = colSums((1-Del.post) * (1/tmp.b) * mu0 * (1/(1+tmp.a)) * (zstar - mustar) * xa.tmp)
#       score.b = colSums((1-Del.post) * (-1/tmp.b) * (mu0 * (zstar - mustar) + (zdagger - mudagger)) * xb.tmp)
#       U = U + rbind(score.g, score.o, score.a, score.b)
#       I = I + .HessFun.AIBeta(w.tmp, xa.tmp, xb.tmp, gamma0, omega0, alpha0, beta0, pZ0, pO0, mu0, phi0, one.R, two.R, 1, n.T)
#       M = M + rbind(score.g, score.o, score.a, score.b) %*% t(rbind(score.g, score.o, score.a, score.b))
#     }
#     alpha.se[j] = sqrt(diag(ginv(I) %*% M %*% t(ginv(I))))[3]
#     beta.se[j] = sqrt(diag(ginv(I) %*% M %*% t(ginv(I))))[4]
#   }
#   return(list(# gamma.est = gamma.now, omega.est = omega.now, 
#     alpha.est = alpha.now, beta.est = beta.now,
#     alpha.se = alpha.se, beta.se = beta.se))
# }

# .AIGDM_MES <- function(ID, Y, n.boot, n.cores){
#   Xa = Xb = W.r = matrix(1, nrow(Y), 1)
#   n.boot = ceiling(n.boot/n.cores)*n.cores
#   
#   K = ncol(Y) - 1
#   zi.check = oi.check = rep(FALSE, K)
#   
#   est.para.lst = lapply(1:K, function(x){
#     .est_para(ID, W.r, Xa, Xb, cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE])),
#               zi.check[x], oi.check[x])
#   })
#   est.para = do.call(Map, c(f = cbind, est.para.lst))
#   alpha.est = est.para$alpha0; beta.est = est.para$beta0
#   
#   bspara.mat = NULL
#   id.unique = sort(unique(ID)); N = length(id.unique)
#   cl = makeCluster(n.cores)
#   clusterExport(cl, c('.est_para', '.AIGDM_EM', '.ZOIeZ', 
#                       '.AIBetaOptim', '.AIBetaNegLoglik', '.AIBetaNegScore',
#                       '.RowbyRow', '.getSubDiag', '.var_logit', '.var_ZStar', '.BlockDiag',
#                       '.dpg', '.dma', '.GenerateScoreFun',
#                       '.corr_zero1', '.corr_zero2', '.corr_rho1', '.corr_rho2', 'ginv',
#                       'ID', 'Y', 'zi.check', 'oi.check', 'n.boot'), envir = environment())
#   m = 0
#   while(m < n.boot){
#     tmp.subset = parLapply(cl, 1:n.cores, function(x){
#       ID.tmp = sort(sample(1:N, replace = T))
#       Y.tmp = NULL
#       for (id in ID.tmp) {
#         Y.tmp = rbind(Y.tmp, Y[which(ID == id),])
#       }
#       Xa = Xb = W.r = matrix(1, nrow(Y), 1)
#       est.para.bs.lst = lapply(1:K, function(x){
#         .est_para(ID.tmp, W.r, Xa, Xb, cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])),
#                   zi.check[x], oi.check[x])
#       })
#       est.para.bs = do.call(Map, c(f = cbind, est.para.bs.lst))
#       return(c(est.para.bs$alpha0, est.para.bs$beta0))
#     })
#     bspara.mat = rbind(bspara.mat, matrix(unlist(tmp.subset), n.cores, 2*K, byrow = T))
#     m = m + n.cores
#   }
#   stopCluster(cl)
#   
#   se = colSds(bspara.mat)
#   alpha.se = se[1:K]; beta.se = se[K+1:K]
#   return(list(alpha.est = alpha.est, beta.est = beta.est,
#               alpha.se = alpha.se, beta.se = beta.se))
# }

.AIGDM_BiasSE <- function(ID, Y, Xa, Xb, zi.id){
  # Y = Y.m; Xa <- Xb <- xmat.m
  
  Xa = cbind(1, Xa) # add the intercept term
  Xb = cbind(1, Xb) # add the intercept term
  
  X.index = 2
  Xa.r = Xa[, -X.index, drop=FALSE]
  Xb.r = Xb[, -X.index, drop=FALSE]
  
  W.r = matrix(1, nrow(Y), 1)
  
  K = ncol(Y) - 1
  zi.check = oi.check = rep(FALSE, K)
  if(!is.null(zi.id))
    zi.check[zi.id] = TRUE
  
  est.para.lst = lapply(1:K, function(x){
    .est_para(ID, W.r, Xa.r, Xb.r, cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE])),
              zi.check[x], oi.check[x])
  })
  asym.mean = .score_test_mean(ID, Xa, X.index, Y, est.para.lst, zi.check, oi.check)
  pval.mean.a = asym.mean$pval
  asym.disp = .score_test_disp(ID, Xb, X.index, Y, est.para.lst, zi.check, oi.check)
  pval.disp.a = asym.disp$pval
  pval.score = c(pval.mean.a, pval.disp.a)
  # names(pval.score) = c("Score-Mean-Asymptotic", "Score-Disp-Asymptotic")
  
  est.para.lst = lapply(1:K, function(x){
    .est_para(ID, W.r, Xa, Xb, cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE])),
              zi.check[x], oi.check[x])
  })
  est.para = do.call(Map, c(f = cbind, est.para.lst))
  est.es.est = c(est.para$alpha0, est.para$beta0)
  
  pval.meanNdisp.sand = .wald_test_sand(ID, X.index, Y, est.para.lst, zi.check, oi.check)
  # pval.wald = c(pval.meanNdisp.sand$pval.mean, pval.meanNdisp.sand$pval.disp)
  # names(pval.wald) = c("Wald-Mean-Asymptotic", "Wald-Disp-Asymptotic")
  pval.meanNdisp.boot = .wald_test_boot(ID, X.index, W.r, Xa, Xb, Y, est.para.lst, zi.check, oi.check, 200)
  pval.wald = c(pval.meanNdisp.sand$pval.mean, pval.meanNdisp.boot$pval.mean,
                pval.meanNdisp.sand$pval.disp, pval.meanNdisp.boot$pval.disp)
  # names(pval.wald) = c("Wald-Mean-Asymptotic", "Wald-Mean-Bootstrap", "Wald-Disp-Asymptotic", "Wald-Disp-Bootstrap")
  
  # est.es.sand.se = .generate_sand_se(ID, Y, est.para.lst, zi.check, oi.check)
  est.es.sand.se = c(pval.meanNdisp.sand$se.mean, pval.meanNdisp.sand$se.disp)
  est.es.boot.se = c(pval.meanNdisp.boot$se.mean, pval.meanNdisp.boot$se.disp)
  return(list(est.es.est = est.es.est,
              est.es.sand.se = est.es.sand.se,
              est.es.boot.se = est.es.boot.se,
              pval = c(pval.wald, pval.score)))
  
}

.generate_sand_se <- function(ID, Y, est.para.lst, zi.check, oi.check){
  K = ncol(Y)-1
  ai.check = zi.check + oi.check
  alpha.se.mat = beta.se.mat = matrix(NA, 2, K)
  for (j in 1:K) {
    Y.tmp = cbind(Y[,j], rowSums(Y[,(j+1):(K+1),drop=FALSE]))
    
    est.para = est.para.lst[[j]]
    ID.tmp = ID; W.r.tmp = est.para$W; Xa.tmp = est.para$Xa; Xb.tmp = est.para$Xb
    dw.r = ncol(W.r.tmp); da = ncol(Xa.tmp); db = ncol(Xb.tmp)
    
    id.unique = sort(unique(ID.tmp)); N = length(id.unique); n.reps = as.vector(table(ID.tmp))
    ng = no = dw.r; na = da; nb = db
    
    if(ai.check[j]){
      if(!zi.check[j]){
        para.keep = c(sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }else if(!oi.check[j]){
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }else{
        para.keep = c(sort(sapply(1:ng, function(x) seq(x, ng, ng))),
                      sort(sapply(1:no, function(x) seq(x, no, no))) + ng,
                      sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                      sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
      }
    }else{
      para.keep = c(sort(sapply(1:na, function(x) seq(x, na, na))) + ng + no,
                    sort(sapply(1:nb, function(x) seq(x, nb, nb))) + ng + no + na)
    }
    
    parakeep.len = length(para.keep)
    para.ind.alpha = (parakeep.len-3):(parakeep.len-2)
    para.ind.beta = (parakeep.len-1):parakeep.len
    # para.ind.alpha = parakeep.len-1
    # para.ind.beta = parakeep.len
    
    gamma0 = est.para$gamma0; omega0 = est.para$omega0; alpha0 = est.para$alpha0; beta0 = est.para$beta0
    sig0 = est.para$sig0; tau0 = est.para$tau0; rho0 = est.para$rho0
    eDel0.lst = est.para$eDel0.lst; eDelZ0.lst = est.para$eDelZ0.lst; eDelO0.lst = est.para$eDelO0.lst
    eZStar.lst = est.para$eZStar.lst; muStar.lst = est.para$muStar.lst
    A.R.lst = est.para$A.R.lst; B.R.lst = est.para$B.R.lst
    A2.R.lst = est.para$A2.R.lst; B2.R.lst = est.para$B2.R.lst; AB.R.lst = est.para$AB.R.lst
    pZ0.lst = est.para$pZ0.lst; pO0.lst = est.para$pO0.lst
    eZ0.lst = est.para$eZ0.lst; mu0.lst = est.para$mu0.lst; phi0.lst = est.para$phi0.lst
    pseudoll.score.beta.lst = est.para$pseudoll.score.beta.lst
    
    # update DVD matrix with full models
    U = 0
    I = 0
    M = matrix(0, ng+no+na+nb, ng+no+na+nb)
    ind.g = 1:ng; ind.o = (ng+1):(ng+no); ind.a = (ng+no+1):(ng+no+na); ind.b = (ng+no+na+1):(ng+no+na+nb)
    for(i in 1:N){
      ID.sel = which(ID.tmp == id.unique[i])
      w.r.tmp = W.r.tmp[ID.sel,,drop=FALSE]
      xa.tmp = Xa.tmp[ID.sel,,drop=FALSE]
      xb.tmp = Xb.tmp[ID.sel,,drop=FALSE]
      n.T = n.reps[i]
      
      pZ0 = pZ0.lst[[i]] # unlist(pZ0.lst[i,])
      eDelZ0 = eDelZ0.lst[[i]]
      cor.DelZ = diag(sig0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(sig0,each=n.T), n.T)
      VdelZ = sqrt(.var_logit(pZ0)) %*% cor.DelZ %*% sqrt(.var_logit(pZ0))
      tmp.g.w = t(.BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ)
      dvd.g.w = - t(.BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)) %*% ginv(VdelZ) %*% .BlockDiag(.dpg(w.r.tmp, gamma0), 1, n.T)
      geeg.w = as.matrix(apply(.RowbyRow(t(tmp.g.w), eDelZ0 - pZ0), 2, sum))
      
      pO0 = pO0.lst[[i]]
      eDelO0 = eDelO0.lst[[i]]
      cor.DelO = diag(tau0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(tau0,each=n.T), n.T)
      VdelO = sqrt(.var_logit(pO0)) %*% cor.DelO %*% sqrt(.var_logit(pO0))
      tmp.o.w = t(.BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)) %*% ginv(VdelO)
      dvd.o.w = - t(.BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)) %*% ginv(VdelO) %*% .BlockDiag(.dpg(w.r.tmp, omega0), 1, n.T)
      geeo.w = as.matrix(apply(.RowbyRow(t(tmp.o.w), eDelO0 - pO0), 2, sum))
      
      mu0 = mu0.lst[[i]]; phi0 = phi0.lst[[i]]
      eZ0 = eZ0.lst[[i]]; eDel0 = eDel0.lst[[i]]
      eZStar = eZStar.lst[[i]]; muStar = muStar.lst[[i]]
      # Vz: variance-covariance matrix of Z in the mean model
      cor.Z = diag(rho0, 1) %x% matrix(1, nrow = n.T, ncol = n.T) + diag(1-rep(rho0,each=n.T), n.T)
      Vz = sqrt(.var_ZStar(mu0, phi0)) %*% cor.Z %*% sqrt(.var_ZStar(mu0, phi0))
      tmp.a.x = t(.BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0))
      dvd.a.x = - t(.BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)) %*% ginv(Vz) %*% diag(c(1-eDel0), length(eDel0)) %*% .BlockDiag(.dma(xa.tmp, alpha0, mu0, phi0), 1, n.T)
      geea.x = as.matrix(apply(.RowbyRow(t(tmp.a.x), eZStar - muStar), 2, sum))
      
      pseudoll.score.beta = colSums(pseudoll.score.beta.lst[[i]] * xb.tmp)
      # dvd.b.x = t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)) %*% t(t(as.matrix(pseudoll.score.beta.lst[[i]] * xb.tmp)))
      geeb.x = as.matrix(pseudoll.score.beta)
      
      ES = rbind(geeg.w, geeo.w, geea.x, geeb.x)
      U = U + ES
      # EI = .BlockDiagBind(dvd.g.w, dvd.o.w, dvd.a.x, dvd.b.x)
      # I = I + (-EI)
      
      a = (1/phi0 - 1) * mu0; a[a < 0] = 0
      b = (1/phi0 - 1) * (1 - mu0); b[b < 0] = 0
      A = digamma(a) - digamma(a+b)
      B = digamma(b) - digamma(a+b)
      A.post = c(A.R.lst[[i]])
      B.post = c(B.R.lst[[i]])
      one.R = (1-eDel0)*( - A + A.post )
      two.R = (1-eDel0)*( - B + B.post )
      tmp.dvd = .DVD(w.r.tmp, xa.tmp, xb.tmp, gamma0, omega0, alpha0, beta0, pZ0, pO0, mu0, phi0, one.R, two.R, 1, n.T)
      tmp.Kab = matrix(0, nrow(tmp.dvd$Kab), ncol(tmp.dvd$Kab))
      dvd.ab.all = rbind(cbind(dvd.a.x, tmp.Kab), cbind(t(tmp.dvd$Kab), tmp.dvd$Kbb))
      EI = .BlockDiagBind(dvd.g.w, dvd.o.w, dvd.ab.all)
      I = I + (-EI)
      
      M[ind.g, ind.g] = M[ind.g, ind.g] + geeg.w %*% t(geeg.w) * (diag(1) %x% matrix(1, ng, ng))
      M[ind.g, ind.o] = M[ind.g, ind.o] + geeg.w %*% t(geeo.w) * (diag(1) %x% matrix(1, ng, no))
      M[ind.o, ind.o] = M[ind.o, ind.o] + geeo.w %*% t(geeo.w) * (diag(1) %x% matrix(1, no, no))
      M[ind.o, ind.g] = M[ind.o, ind.g] + geeo.w %*% t(geeg.w) * (diag(1) %x% matrix(1, no, ng))
      
      M[ind.a, ind.a] = M[ind.a, ind.a] + geea.x %*% t(geea.x) * (diag(1) %x% matrix(1, na, na))
      M[ind.b, ind.b] = M[ind.b, ind.b] + geeb.x %*% t(geeb.x) * (diag(1) %x% matrix(1, nb, nb))
      
      M[ind.a, ind.b] = M[ind.a, ind.b] + geea.x %*% t(geeb.x) * (diag(1) %x% matrix(1, na, nb))
      M[ind.b, ind.a] = M[ind.a, ind.b] + geeb.x %*% t(geea.x) * (diag(1) %x% matrix(1, nb, na))
      
      M[c(ind.g, ind.o), c(ind.a, ind.b)] = M[c(ind.g, ind.o), c(ind.a, ind.b)] + rbind(geeg.w, geeo.w) %*% t(rbind(geea.x, geeb.x)) * (diag(1) %x% matrix(1, ng+no, na+nb))
      M[c(ind.a, ind.b), c(ind.g, ind.o)] = M[c(ind.a, ind.b), c(ind.g, ind.o)] + rbind(geea.x, geeb.x) %*% t(rbind(geeg.w, geeo.w)) * (diag(1) %x% matrix(1, na+nb, ng+no))
    }
    
    I.r = I[para.keep, para.keep]
    M.r = M[para.keep, para.keep]
    H.r = ginv(I.r)
    tmp.sand.se = sqrt(diag(H.r %*% M.r %*% t(H.r)))
    alpha.se.mat[,j] = c(tmp.sand.se[para.ind.alpha])
    beta.se.mat[,j] = c(tmp.sand.se[para.ind.beta])
    
    # tmp.sand.se.new = sqrt(diag((ginv(I) %*% M %*% t(ginv(I)))[para.keep, para.keep]))
    # alpha.se.new.mat[,j] = c(tmp.sand.se.new[para.ind.alpha])
    # beta.se.new.mat[,j] = c(tmp.sand.se.new[para.ind.beta])
    
    # I.r11 = I[para.keep, para.keep]
    # I.r12 = I[para.keep, -para.keep]
    # I.r22 = I[-para.keep, -para.keep]
    # tmp.sand.se.new = sqrt(diag( ginv(I.r11-I.r12%*%ginv(I.r22)%*%t(I.r12)) %*% ginv((ginv(I) %*% M %*% ginv(I))[para.keep, para.keep]) %*%  ginv(I.r11-I.r12%*%ginv(I.r22)%*%t(I.r12)) ))
    # alpha.se.new.mat[,j] = c(tmp.sand.se.new[para.ind.alpha])
    # beta.se.new.mat[,j] = c(tmp.sand.se.new[para.ind.beta])
  }
  return(c(alpha.se.mat, beta.se.mat))
}

.rarefy <- function (otu.tab) {
  depth = min(rowSums(otu.tab))
  otu.tab <- as.matrix(otu.tab)
  ind <- (rowSums(otu.tab) < depth)
  sam.discard <- rownames(otu.tab)[ind]
  otu.tab <- otu.tab[!ind, ]
  rarefy <- function(x, depth) {
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y.tab))] <- y.tab
    z
  }
  otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
  rownames(otu.tab.rff) <- rownames(otu.tab)
  colnames(otu.tab.rff) <- colnames(otu.tab)
  return(otu.tab.rff)
}

# overwrite function below in miLineage
###############
# m is the total number of taxa
###############

library("geepack")

.F.test <- function(x){
  
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.simes.test <- function(x){
  
  return( min(length(x) * x/rank(x)) )
  
}

# add 06/19/2016
.diag2 <- function(x){
  
  if(length(x)>1){
    return(diag(x))
  }else{
    return(as.matrix(x))
    
  }
  
}

# assume two subjects in each strata
# last updated 05/16/2016
.sampling.strata<- function(var, strata){
  
  n.strata = length(table(strata))
  strata.uniq = unique(strata)
  var.sample = var
  for(i in 1:n.strata){
    index = which(strata==strata.uniq[i])
    if( rbinom(1, 1, 0.5)==1 ){
      tmp = var.sample[index[2]]
      var.sample[index[2]] = var.sample[index[1]]
      var.sample[index[1]] = tmp
    }
    
  }
  return(var.sample)
  
}

########################################
#                                      #
#           One Part Model             #
#                                      #
########################################

## change on 03/28/2016
.Ei.beta <- function(m, p, beta, X.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  for(j in 1:(m-1)){
    
    Ei.out[j] = exp(beta[((j-1)*p+1):(j*p)] %*% X.i)
    
  }
  
  
  Ei.out[m] = 1
  
  
  
  
  return (Ei.out)
}


## change on 03/28/2016
.fun.neg.loglik.beta <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      loglik = loglik + Y[i,] %*% log(P.i)
    }
    
  }
  
  return (-loglik)
  
}

.fun.neg.score.beta <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      Score.beta = Score.beta + kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      
    }
    
    return (-Score.beta)
  }
  
  
  
}

.fun.score.i.beta <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      
      Score.beta.i[i,] =  kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1) )
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}


.fun.hessian.beta <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; X = data$X
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  n.beta = (m-1)*p
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      Hessian.beta = Hessian.beta + kronecker( tmp.beta, ( X[i,] %o% X[i,] ) ) 
      
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }
    
    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}


.Score.test.stat <- function(Y, X, X.par.index){
  
  p = ncol(X)
  
  nY = rowSums(Y)
  
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = (m - 1)*p
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NA
    
  }else{
    
    X.reduce = X[,-X.par.index, drop=FALSE]
    p.reduce = p - length(X.par.index)    
    par.interest.index.beta =  kronecker( ((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
    
    n.par.interest.beta = length(par.interest.index.beta) 
    
    beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
    
    data.reduce.beta = list(Y=Y, X=X.reduce)
    
    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0 
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta, gr=.fun.neg.score.beta, data = data.reduce.beta, method="BFGS")$par
    
    
    data.beta = list(Y=Y, X=X)
    Score.reduce.beta = .fun.score.i.beta(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    tmp = .fun.hessian.beta(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list  
    
    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
    
    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    
    
    B2 =  matrix(0, n.beta, n.beta)
    
    # change in 04/08/2016(warning! need to change this step in the resampling function too)
    for(i in 1:n){
      B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
    }
    
    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A
    
    
  }
  
  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list )   )
  
  
}

# add 05/31/2016: allow for multiple nuisance covariate and multiple covariates of interest
# Score.test.stat.pos.4resampling is the simplier version for two-group comparison
# rename 06/16/2016 Score.test.stat.pos.4Gresampling to Score.test.stat.4Gresampling and use in both one-part model and postive part of the two-part model
.Score.test.stat.4Gresampling <- function(X.perm, X.par.index, S.beta.list, I.beta.list){
  
  n = nrow(X.perm)  
  p = ncol(X.perm)
  
  m.beta = length(S.beta.list[[1]])
  
  #n.par.interest.beta = m.beta
  n.beta = m.beta*p
  #n.beta = m.beta*2
  par.interest.index.beta =  kronecker( ((0:(m.beta-1))*p), rep(1,length(X.par.index))) + X.par.index
  #par.interest.index.beta = (1:m.beta )*2
  n.par.interest.beta = length(par.interest.index.beta) 
  
  
  Score.reduce.beta.perm = matrix(0, n, n.beta )
  Hess.reduce.beta.perm = matrix(0, n.beta, n.beta )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(X.perm[i,], ncol=1))  
    
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  X.perm[i,] %o% X.perm[i,] ) )
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  for(i in 1:n){
    B2 = B2 + Score.reduce.beta.perm.reorg[i,] %o% Score.reduce.beta.perm.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.beta.perm = A %*% ginv(B) %*% A
  
  
  return(score.stat.beta.perm)
  
  
  
}


# add 07/02/2016 for adaptive resampling
.resample.work.one <- function(X, X.par.index, score.stat.beta, S.beta.list, I.beta.list, start.nperm, end.nperm, n.one, one.acc){
  
  n = nrow(X)
  
  n.one.new = n.one
  one.acc.new = one.acc
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    X.perm = X
    X.perm[,X.par.index] = X.perm[perm.index,X.par.index]
    
    score.stat.beta.perm = try( .Score.test.stat.4Gresampling(X.perm, X.par.index, S.beta.list, I.beta.list) )
    
    if(class(score.stat.beta.perm) != "try-error"){
      
      n.one.new = n.one.new + 1
      if(score.stat.beta.perm >= score.stat.beta){
        one.acc.new = one.acc.new + 1
        
      }  
    }
  }
  
  if(one.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(one.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.one.new=n.one.new, one.acc.new=one.acc.new, flag=flag, next.end.nperm=next.end.nperm))
  
}


# Y: nxm count of microbiomes
# X: covariates
# Z: covariates 
# X.par.index: index for the parameter of interest for the X part
# Z.par.index: index for the parameter of interest for the Z part
# change 06/16/2016, change the resampling part, to accommodate multiple covariate and potential confounders

.Score.test <- function(Y, X, X.par.index, seed=11, resample=FALSE, n.replicates=NULL){
  
  
  p = ncol(X)
  
  nY = rowSums(Y)
  
  ## remove 03/28/2016
  #   nY0.index = which(nY==0)
  #   if(length(nY0.index)>0){
  #     Y = Y[-nY0.index, , drop=FALSE]
  #     X = X[-nY0.index, , drop=FALSE]
  #   }
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = (m - 1)*p
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NULL
    score.pvalue.beta = NA
    n.par.interest.beta = 0
    
  }else{
    
    tmp.one = try( .Score.test.stat(Y, X, X.par.index) )
    if(class(tmp.one) == "try-error"){
      
      score.stat.beta = NA
      score.pvalue.beta = NA
      n.par.interest.beta = NA
      
    }else{
      
      n.par.interest.beta = (m-1)*length(X.par.index)
      score.stat.beta = tmp.one$score.stat.beta
      score.pvalue.beta = 1 - pchisq(score.stat.beta,  n.par.interest.beta) 
      
    }
    
    
  }
  
  
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue.beta, df =n.par.interest.beta)
  #  print(score.stat.beta)
  
  if(resample){
    ########################################
    #                                      #
    #        Resampling Score Test         #
    # change 07/02/2016 for adaptive resampling
    #                                      #
    ########################################  
    set.seed(seed)
    if(!is.na(score.stat.beta)){
      
      n.one = 0
      one.acc = 0
      
      start.nperm = 1;
      end.nperm = min(100,n.replicates);
      flag = 1
      while(flag & end.nperm <= n.replicates){
        
        results = .resample.work.one(X, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, end.nperm, n.one, one.acc)
        n.one = results$n.one.new
        one.acc = results$one.acc.new
        flag = results$flag
        next.end.nperm = results$next.end.nperm
        
        if(flag){
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
        }
        
        if(start.nperm < n.replicates & end.nperm > n.replicates){ 
          #warning(paste( "Inaccurate pvalue with", n.replicates, "permutations"))
          results = .resample.work.one(X, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, n.replicates, n.one, one.acc)
          n.one = results$n.one.new
          one.acc = results$one.acc.new
          
        }
        
      }
      
      
      #      print(paste("Final number of resamplings: ", n.one) )
      #      if(n.one<end.nperm/2){
      #        print("Number of resamplings too small for one-part test")
      #      }
      
      tmp = (one.acc+1)/(n.one+1)
      
      #print(n.one)
      #print(one.acc)
      
    }else{
      
      tmp = NA
    }
    
    beta.results = c(beta.results, score.Rpvalue = tmp)
    
    
  }
  
  return(beta.results)
  
  
}
# Y: nxm count of microbiomes
# case (case/control status: 1 for cases, 0 for controls)
# score.stat.gamma: score statistics for gamma part
# score.stat.beta: score statistics for beta part
# S.gamma.list, S.beta.list: list with length n, the element i is the score statistics (length = m-1) of the subject i 
# I.gamma.list, I.beta.list: list with length n, the element i is the Hessian matrix ( (m-1)x(m-1) ) of the subject i 

.Tstat <- function(Y, case){
  
  m = ncol(Y)
  n = nrow(Y)
  n1 = sum(case==1)
  
  Ym = Y[,-m, drop=FALSE]
  nY = rowSums(Y)
  nY1 = sum(nY[case==1])
  nY2 = sum(nY[case==0])
  p0 = colSums(Ym)/sum(nY)
  
  A = colSums(Ym[case==1, , drop=FALSE]) - p0 * sum(nY[case==1])
  
  B = matrix(0, m-1, m-1)
  
  C = matrix(0, m-1, m-1)
  index.case = which(case==1) 
  n.case = length(index.case)
  for(i in 1:n.case){
    
    tmp = Ym[index.case[i], ] - p0 * nY[index.case[i]]
    C.tmp = tmp %o% tmp
    C = C + C.tmp
    
  }
  
  D = matrix(0, m-1, m-1)
  for(i in 1:n){
    
    tmp = Ym[i, ] - p0 * nY[i]
    D.tmp = tmp %o% tmp
    D = D + D.tmp
    
  }
  
  B = ((nY2-nY1)/sum(nY))*C + ((nY1/sum(nY))^2)*D
  
  
  score.stat = A %*% ginv(B) %*% A
  
  
  return(score.stat)
}

.Score.test.simple <- function(Y, case, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  #print(A)
  #print(B)
  m = ncol(Y)
  
  score.stat.beta = try( .Tstat(Y, case) )
  if(class(score.stat.beta) == "try-error"){
    
    score.stat.beta = NA
    score.pvalue = NA
    df = NA
    
  }else{
    
    score.pvalue = 1 - pchisq(score.stat.beta, m-1) 
    df = m-1
  }
  
  
  #print("observed score statistics:")
  #print(score.stat.beta)
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue, df = df )
  
  if(resample){
    
    if(!is.na(score.stat.beta)){
      
      n.one = 0
      one.acc = 0
      #print("simulated stat:")
      #set.seed(16)
      
      for(k in 1:n.replicates){
        
        case.perm = sample(case)
        score.stat.beta.perm = try( .Tstat(Y, case.perm) )
        
        #       if(k<10){
        #         print("case.perm"); print(case.perm)
        #         #       print("A"); print(A)
        #         #       print("Hess.reduce.beta.perm.reorg"); print(Hess.reduce.beta.perm.reorg)
        #         #       print("B1"); print(B1)
        #         #       print("B2"); print(B2)
        #         #       print("B"); print(B)
        #         print(score.stat.beta.perm)
        #         
        #       }
        if(class(score.stat.beta.perm) != "try-error"){
          
          n.one = n.one + 1
          if(score.stat.beta.perm >= score.stat.beta){
            one.acc = one.acc + 1
            
          }           
          
        }
        
        
      }
      
      if(n.one<n.replicates/2){
        print("#replicate too small for one-part test")
      }
      
      score.Rpvalue = one.acc/n.one
      
      
    }else{
      
      score.Rpvalue = NA
    }
    
    beta.results = c(beta.results, score.Rpvalue = score.Rpvalue)
    
  }
  
  
  return(beta.results)
  
}


.Score.test.simple.restrict <- function(Y, case, strata, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  
  #print(A)
  #print(B)
  m = ncol(Y)
  
  score.stat.beta = try( .Tstat(Y, case) )
  if(class(score.stat.beta) == "try-error"){
    
    score.stat.beta = NA
    score.pvalue = NA
    df = NA
    
  }else{
    
    score.pvalue = 1 - pchisq(score.stat.beta, m-1) 
    df = m-1
  }
  
  
  #print("observed score statistics:")
  #print(score.stat.beta)
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue, df = df )
  
  if(resample){
    
    if(!is.na(score.stat.beta)){
      
      n.one = 0
      one.acc = 0
      #print("simulated stat:")
      #set.seed(16)
      
      for(k in 1:n.replicates){
        
        case.perm = .sampling.strata(case, strata)
        
        score.stat.beta.perm = try( .Tstat(Y, case.perm) )
        
        
        
        #       if(k<10){
        #         print("case.perm"); print(case.perm)
        #         #       print("A"); print(A)
        #         #       print("Hess.reduce.beta.perm.reorg"); print(Hess.reduce.beta.perm.reorg)
        #         #       print("B1"); print(B1)
        #         #       print("B2"); print(B2)
        #         #       print("B"); print(B)
        #         print(score.stat.beta.perm)
        #         
        #       }
        if(class(score.stat.beta.perm) != "try-error"){
          
          n.one = n.one + 1
          if(score.stat.beta.perm >= score.stat.beta){
            one.acc = one.acc + 1
            
          }           
          
        }
        
        
      }
      
      if(n.one<n.replicates/2){
        print("#replicate too small for one-part test")
      }
      
      score.Rpvalue = one.acc/n.one
      
      
    }else{
      
      score.Rpvalue = NA
    }
    
    beta.results = c(beta.results, score.Rpvalue = score.Rpvalue)
    
  }
  
  
  return(beta.results)
  
}

########################################
#                                      #
#    Two Part Model: positive Part     #
#    ## add on 04/25/2016              #
#                                      #
########################################


.Ei.beta.pos <- function(m, p, beta, X.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  I.i = as.numeric(Y.i>0)
  
  for(j in 1:(m-1)){
    
    Ei.out[j] = I.i[j]*exp(beta[((j-1)*p+1):(j*p)] %*% X.i)
    
  }
  
  
  Ei.out[m] = I.i[m]
  
  
  return (Ei.out)
}

.fun.neg.loglik.beta.pos <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      index = which(Y[i,]>0)
      loglik = loglik + Y[i,index] %*% log(P.i[index])
      #       if(is.na(loglik)){
      #         print(i); break;
      #       }
    }
    
  }
  
  return (-loglik)
  
}

.fun.neg.score.beta.pos <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      Score.beta = Score.beta + kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      
    }
    
    return (-Score.beta)
  }
  
  
  
}

.fun.score.i.beta.pos <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      
      Score.beta.i[i,] =  kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1)) 
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}

.fun.hessian.beta.pos <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; X = data$X
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  n.beta = (m-1)*p
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      Hessian.beta = Hessian.beta + kronecker(tmp.beta, ( X[i,] %o% X[i,] ))
      
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }
    
    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}


.Score.test.stat.pos <- function(Y, X, X.par.index){
  
  p = ncol(X)
  
  nY = rowSums(Y)
  
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = (m - 1)*p
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NA
    
  }else{
    
    X.reduce = X[,-X.par.index, drop=FALSE]
    p.reduce = p - length(X.par.index)    
    par.interest.index.beta =  kronecker(((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
    
    n.par.interest.beta = length(par.interest.index.beta) 
    
    beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
    
    data.reduce.beta = list(Y=Y, X=X.reduce)
    
    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0 
    #est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos, gr=.fun.neg.score.beta.pos, data = data.reduce.beta)$par
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos, gr=.fun.neg.score.beta.pos, data = data.reduce.beta, method="BFGS")$par
    
    
    data.beta = list(Y=Y, X=X)
    Score.reduce.beta = .fun.score.i.beta.pos(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    tmp = .fun.hessian.beta.pos(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list  
    
    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
    
    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    
    
    B2 =  matrix(0, n.beta, n.beta)
    
    # change in 04/08/2016(warning! need to change this step in the resampling function too)
    for(i in 1:n){
      B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
    }
    
    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A
    
    
  }
  
  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list )   )
  
  
}

.Score.test.stat.pos.4resampling <- function(case.perm, S.beta.list, I.beta.list){
  
  n = length(case.perm)  
  
  m.beta = length(S.beta.list[[1]])
  n.par.interest.beta = m.beta
  n.beta = m.beta*2
  par.interest.index.beta = (1:m.beta )*2
  beta.acc = 0
  
  #print("simulated stat:")
  #set.seed(16)
  
  
  XZ = cbind(1,case.perm)
  
  
  Score.reduce.beta.perm = matrix(0, n, m.beta*2 )
  Hess.reduce.beta.perm = matrix(0, m.beta*2, m.beta*2 )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(XZ[i,], ncol=1))  
    
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  XZ[i,] %o% XZ[i,] ) )
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  for(i in 1:n){
    B2 = B2 + Score.reduce.beta.perm.reorg[i,] %o% Score.reduce.beta.perm.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.beta.perm = A %*% ginv(B) %*% A
  
  
  return(score.stat.beta.perm)
  
  
  
}


########################################
#                                      #
#          Two Part Model              #
#            (zero part)               #
#                                      #
########################################
.Tstat.zero <- function(Y0, case){
  
  m = ncol(Y0)
  n = nrow(Y0)
  n1 = sum(case==1)  
  n2 = n-n1
  
  A = (colSums(Y0[case==1, ,drop=FALSE])/n1  - colSums(Y0[case==0, ,drop=FALSE])/n2)
  
  BB = matrix(0, m, m)
  for(i in 1:n){
    BB = BB + Y0[i,] %o% Y0[i,]
    
  }
  
  tmp = colSums(Y0)/n
  B = (1/n1 + 1/n2)*(BB/n - tmp%o%tmp)
  
  score.stat = A %*% ginv(B) %*% A
  return(score.stat)
  
}


# add 05/31/2016: function GEE.zero to run general zero-part wald test
#Y0 = Y
#Y0[Y==0] = 1  
#Y0[Y>0] = 0
#remove.index = which(colSums(Y0)==0)
#Y0 = Y0[,-remove.index]
#Z = cbind( 1, c(rep(0, 25), rep(1, 25)), rnorm(50) )
#Z.par.index = 2
#GEE.zero(Y0, Z, Z.par.index, "independence")
.GEE.zero <- function(Y0, Z, Z.par.index, cor.stru){
  
  Z.reduce = Z[,-Z.par.index,drop=FALSE]
  n = nrow(Y0)
  m = ncol(Y0)
  p = ncol(Z)
  p.reduce = ncol(Z.reduce)
  outcome = NULL
  id = NULL
  cova = NULL
  cova.reduce = NULL
  for(i in 1:n){
    
    outcome = c(outcome, Y0[i,])
    index.start = 1
    index.end = p
    
    index.start.reduce = 1
    index.end.reduce = p.reduce
    
    for(j in 1:m){
      tmp = rep(0, m*p)
      tmp[index.start:index.end] = Z[i,]
      cova = rbind(cova, tmp )
      index.start = index.start + p
      index.end = index.end + p
      
      tmp = rep(0, m*p.reduce)
      tmp[index.start.reduce:index.end.reduce] = Z.reduce[i,]
      cova.reduce = rbind(cova.reduce, tmp )
      index.start.reduce = index.start.reduce + p.reduce
      index.end.reduce = index.end.reduce + p.reduce
      
    }
    
    
    id = c(id, rep(i, m))
  }
  
  data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  
  gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence", std.err ="jack")
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence", std.err ="jack")
  wald.test = anova(gee.full, gee.reduce)
  
  
  return(wald.test)
  
}


.Pi.alpha<-function(m, p, alpha, X.i){
  
  Pi.out = rep(NA,m)
  
  for(j in 1:m){
    
    tmp = exp(alpha[((j-1)*p+1):(j*p)] %*% X.i)
    if(is.infinite(tmp)){
      Pi.out[j] = 1
    }else{
      Pi.out[j] = tmp/(tmp + 1)
    }
    
  }
  
  
  return (Pi.out)    
}

#fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list=TRUE)
#alpha = est.reduce.alpha
#data = data.alpha 
.fun.score.i.alpha <- function(alpha, data, save.list=FALSE){
  
  Y = data$Y; Z = data$Z; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)
  
  vA.list = list()
  Vinv.list = list()
  VY.list = list()
  
  n.alpha = m*p
  
  if(length(alpha)!=n.alpha){
    
    warning("Dim of initial alpha does not match the dim of covariates")
    
  }else{
    
    Score.alpha.i = matrix(0, n, n.alpha)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha(m, p, alpha, Z[i,])
      vA.tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(vA.tmp) 
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure
      
      tmp.V.i = ginv(V.i)
      tmp.VY = tmp.V.i %*% (Y[i,] - Pi.i)
      Score.alpha.i[i,] = t.D.i %*% tmp.VY
      
      if(save.list){
        vA.list[[i]] = vA.tmp
        Vinv.list[[i]] = tmp.V.i
        VY.list[[i]] = tmp.VY
      }
    }
    
    
  }
  
  if(save.list){
    
    return ( list(Score.alpha=Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list=VY.list) )
    
  }else{
    
    return (Score.alpha.i)
  }  
  
}

#fun.hessian.alpha(est.reduce.alpha, data.alpha)
.fun.hessian.alpha <- function(alpha, data){
  
  Y = data$Y; Z = data$Z
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)
  n.alpha = m*p
  
  if(length(alpha)!=n.alpha){
    print("Waring: dim of alpha is not the same as alpha\n")
    
  }else{
    
    Hessian.alpha = matrix(0, nrow=n.alpha, ncol=n.alpha)
    nY = rowSums(Y)
    
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha(m, p, alpha, Z[i,])
      tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(tmp)
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure        
      
      
      Hessian.alpha = Hessian.alpha + t.D.i %*% ginv(V.i) %*% t(t.D.i) 
      
      
    }
    
    return (Hessian.alpha)
    
    
  }
  
  
}

# add 06/01/2016: function to run general GEE score test for the zero-part
#Y0 = Y
#Y0[Y==0] = 1  
#Y0[Y>0] = 0
#remove.index = which(colSums(Y0)==0)
#Y0 = Y0[,-remove.index]

#Z.par.index = 2
#GEE.zero(Y0, Z, Z.par.index, "independence")

.Score.test.stat.zero <- function(Y0, Z, Z.par.index, cor.stru){
  
  Z.reduce = Z[,-Z.par.index,drop=FALSE]
  n = nrow(Y0)
  m = ncol(Y0)
  p = ncol(Z)
  p.reduce = ncol(Z.reduce)
  outcome = NULL
  id = NULL
  cova.reduce = NULL
  for(i in 1:n){
    
    outcome = c(outcome, Y0[i,])
    index.start = 1
    index.end = p
    
    index.start.reduce = 1
    index.end.reduce = p.reduce
    
    for(j in 1:m){
      
      tmp = rep(0, m*p.reduce)
      tmp[index.start.reduce:index.end.reduce] = Z.reduce[i,]
      cova.reduce = rbind(cova.reduce, tmp )
      index.start.reduce = index.start.reduce + p.reduce
      index.end.reduce = index.end.reduce + p.reduce
      
    }
    
    
    id = c(id, rep(i, m))
  }
  
  #data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  #gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence")
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence")
  #wald.test = anova(gee.full, gee.reduce)
  
  
  ########### perform score test
  n.alpha = m *p
  par.interest.index.alpha =  kronecker( ((0:(m-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  n.par.interest.alpha = length(par.interest.index.alpha) 
  est.reduce.alpha = rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] = 0 
  est.reduce.alpha[-par.interest.index.alpha] = coef(gee.reduce)
  est.reduce.scale = gee.reduce
  
  data.alpha = list(Y=Y0, Z=Z)
  
  tmp = .fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list=TRUE)
  Score.reduce.alpha = tmp$Score.alpha
  # for resampling test
  vA.list = tmp$vA.list
  Vinv.list = tmp$Vinv.list
  VY.list = tmp$VY.list
  
  Hess.reduce.alpha =  .fun.hessian.alpha(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  for(i in 1:n){
    B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.alpha = A %*% ginv(B) %*% A
  score.pvalue.alpha = 1 - pchisq(score.stat.alpha, n.par.interest.alpha) 
  
  
  
  return(list(score.df.alpha=n.par.interest.alpha, score.stat.alpha = score.stat.alpha, score.pvalue.alpha=score.pvalue.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list )   )
  
}

.Score.test.stat.zero.4Gresampling <- function(Z.perm, Z.par.index, vA.list, Vinv.list, VY.list){
  
  n = nrow(Z.perm)
  p = ncol(Z.perm)
  m.alpha = length(vA.list[[1]])
  n.alpha = m.alpha*p
  
  par.interest.index.alpha =  kronecker( ((0:(m.alpha-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  n.par.interest.alpha = length(par.interest.index.alpha) 
  
  Score.reduce.alpha.perm = matrix(0, n, n.alpha )
  Hess.reduce.alpha.perm = matrix(0, n.alpha, n.alpha )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         alpha part: resampling Score test        #
    #                                                 #
    ################################################### 
    tD.tmp = kronecker(.diag2(vA.list[[i]]), as.matrix(Z.perm[i,], ncol=1))
    
    Score.reduce.alpha.perm[i,] = Score.reduce.alpha.perm[i,] + tD.tmp %*% VY.list[[i]]
    
    Hess.reduce.alpha.perm = Hess.reduce.alpha.perm + tD.tmp %*% Vinv.list[[i]] %*% t(tD.tmp)
    
    
  }
  
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha.perm[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  for(i in 1:n){
    B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.alpha = A %*% ginv(B) %*% A
  
  
  return(score.stat.alpha)
  
}

########################################
#                                      #
#          Two Part Model              #
#                                      #
#                                      #
########################################

# add 05/31/2016: function Score.test2.zerowald to run two-part test in the general setting (zero part: use wald GEE test)
# Y: nxm count of microbiomes
# X: covariates for positive part: first column is always intercept
# Z: covariates for zero part: first column is always intercept
# X.par.index: index for the parameter of interest for the X part
# Z.par.index: index for the parameter of interest for the Z part
.Score.test2.zerowald <- function(Y, X, X.par.index, Z, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  
  n = nrow(X)
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects have the same values of covariates (e.g., all case/control)
      index.cova = 1 + which(apply(X1[,-1,drop=FALSE], 2, function(x) length(table(x)) ) > 1) # index of valid covariates
      X1.par.index.ava = as.numeric( !is.na(match(X.par.index, index.cova))) # use 0/1 to indicate the index is still availiable or not
      # if no covariate left; even if have covariate left, they are not covariate of interest; no taxa
      if( length(index.cova)<1 |  sum(X1.par.index.ava)==0 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        X1 = X1[,c(1, index.cova), drop=FALSE] 
        tmp = match(X.par.index, index.cova) + 1
        X1.par.index = tmp[!is.na(tmp)]
        d1.pos = length(X1.par.index)
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, X1.par.index) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, d1.pos*(m-1) ) 
          df.pos = d1.pos * (m-1)
        }
        
      }    
      
    }
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      score.stat.zero = try( .GEE.zero(Y0, Z, Z.par.index, "independence") )
      if(class(score.stat.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        score.pvalue.zero = score.stat.zero[[3]]
        df.zero = score.stat.zero[[1]]
        score.stat.zero = score.stat.zero[[2]]
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      df.comb = df.zero + df.pos
      score.pvalue.comb = 1 - pchisq(score.stat.comb, df.comb ) 
      
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = df.pos
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = df.zero
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      #set.seed(16)
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        for(k in 1:n.replicates){
          
          perm.index = sample(1:n)
          X.perm = X
          X.perm[,X.par.index,drop=FALSE] = X.perm[perm.index,X.par.index,drop=FALSE]
          X1.perm = X.perm[index.subj.pos, , drop=FALSE]
          X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
          
          Z.perm = Z
          Z.perm[,Z.par.index,drop=FALSE] = Z.perm[perm.index,Z.par.index,drop=FALSE]
          
          score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          score.stat.zero.perm = try( .GEE.zero(Y0, Z.perm, Z.par.index, "independence")[[2]] )
          
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            } 
          }
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }  
          }
          
          if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
            
            score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
            n.comb = n.comb + 1
            if(score.stat.comb.perm >= score.stat.comb){
              comb.acc = comb.acc + 1
              
            }  
          }
          
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        if(n.comb<n.replicates/2){
          print("#replicate too small for comb test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        score.Rpvalue.zero = zero.acc/n.zero
        score.Rpvalue.comb = comb.acc/n.comb
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        for(k in 1:n.replicates){
          
          #case.perm = sample(case)
          #case1.perm = case.perm[index.subj.pos]
          perm.index = sample(1:n)
          X.perm = X
          X.perm[,X.par.index,drop=FALSE] = X.perm[perm.index,X.par.index,drop=FALSE]
          X1.perm = X.perm[index.subj.pos, , drop=FALSE]
          X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
          
          score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            }  
          }
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        for(k in 1:n.replicates){
          
          #case.perm = sample(case)
          perm.index = sample(1:n)
          Z.perm = Z
          Z.perm[,Z.par.index,drop=FALSE] = Z.perm[perm.index,Z.par.index,drop=FALSE]
          
          score.stat.zero.perm = try( .GEE.zero(Y0, Z.perm, Z.par.index, "independence")[[2]] )
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }          
          }
          
          
          
        }
        
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        
        score.Rpvalue.zero = zero.acc/n.zero
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}



# add 07/02/2016 for adaptive resampling
.resample.work.two <- function(X, X.par.index, X1.par.index, Z, Z.par.index, index.subj.pos, index.cova, score.stat.pos, score.stat.zero, pos.S.beta.list, pos.I.beta.list, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc){
  
  n = nrow(X)
  
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  n.comb.new = n.comb
  comb.acc.new = comb.acc
  
  score.stat.comb = score.stat.pos + score.stat.zero
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    X.perm = X
    X.perm[,X.par.index] = X.perm[perm.index,X.par.index]
    X1.perm = X.perm[index.subj.pos, , drop=FALSE]
    X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
    
    Z.perm = Z
    Z.perm[,Z.par.index] = Z.perm[perm.index,Z.par.index]
    
    score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, pos.S.beta.list,pos.I.beta.list) )
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling(Z.perm, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    if(class(score.stat.pos.perm) != "try-error"){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    if(class(score.stat.zero.perm) != "try-error"){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
      
      score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
      n.comb.new = n.comb.new + 1
      if(score.stat.comb.perm >= score.stat.comb){
        comb.acc.new = comb.acc.new + 1
        
      }  
    }
    
  }
  
  if(comb.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(comb.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              n.comb.new=n.comb.new, comb.acc.new=comb.acc.new,
              flag=flag, next.end.nperm=next.end.nperm))
  
}

# add 07/02/2016 for adaptive resampling
.resample.work.pos <- function(X, X.par.index, X1.par.index, index.subj.pos, index.cova, score.stat.pos, pos.S.beta.list, pos.I.beta.list, start.nperm, end.nperm, n.pos, pos.acc){
  
  n = nrow(X)
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    X.perm = X
    X.perm[,X.par.index] = X.perm[perm.index,X.par.index]
    X1.perm = X.perm[index.subj.pos, , drop=FALSE]
    X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
    
    
    score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, pos.S.beta.list,pos.I.beta.list) )
    
    
    if(class(score.stat.pos.perm) != "try-error"){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    
    
  }
  
  if(pos.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(pos.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}

# add 07/02/2016 for adaptive resampling
.resample.work.zero <- function(Z, Z.par.index, score.stat.zero, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.zero, zero.acc){
  
  n = nrow(Z)
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    
    Z.perm = Z
    Z.perm[,Z.par.index] = Z.perm[perm.index,Z.par.index]
    
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling(Z.perm, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    
    if(class(score.stat.zero.perm) != "try-error"){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    
  }
  
  if(zero.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(zero.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}

# add 06/01/2016: function Score.test2 to run two-part test in the general setting(zero part: use score GEE test)
# Y: nxm count of microbiomes
# X: covariates for positive part: first column is always intercept
# Z: covariates for zero part: first column is always intercept
# X.par.index: index for the parameter of interest for the X part
# Z.par.index: index for the parameter of interest for the Z part

.Score.test2 <- function(Y, X, X.par.index, Z, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL){
  
  
  n = nrow(X)
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects have the same values of covariates (e.g., all case/control)
      index.cova = 1 + which(apply(X1[,-1,drop=FALSE], 2, function(x) length(table(x)) ) > 1) # index of valid covariates
      X1.par.index.ava = as.numeric( !is.na(match(X.par.index, index.cova))) # use 0/1 to indicate the index is still availiable or not
      # if no covariate left; even if have covariate left, they are not covariate of interest; no taxa
      if( length(index.cova)<1 |  sum(X1.par.index.ava)==0 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        X1 = X1[,c(1, index.cova), drop=FALSE] 
        tmp = match(X.par.index, index.cova) + 1
        X1.par.index = tmp[!is.na(tmp)]
        d1.pos = length(X1.par.index)
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, X1.par.index) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, d1.pos*(m-1) ) 
          df.pos = d1.pos * (m-1)
        }
        
      }    
      
    }
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      
      tmp.zero = try( .Score.test.stat.zero(Y0, Z, Z.par.index, "independence") )
      if(class(tmp.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        df.zero = tmp.zero[[1]]
        score.stat.zero = tmp.zero[[2]]
        score.pvalue.zero = tmp.zero[[3]]
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      df.comb = df.zero + df.pos
      score.pvalue.comb = 1 - pchisq(score.stat.comb, df.comb ) 
      
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = df.pos
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = df.zero
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      set.seed(seed)
      
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        while(flag & end.nperm <= n.replicates){
          
          results = .resample.work.two(X, X.par.index, X1.par.index, Z, Z.par.index, index.subj.pos, index.cova, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc)
          
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          n.comb = results$n.comb.new
          comb.acc = results$comb.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          if(flag){
            start.nperm = end.nperm + 1;
            end.nperm = next.end.nperm;
            
          }
          
          if(start.nperm < n.replicates & end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            results = .resample.work.two(X, X.par.index, X1.par.index, Z, Z.par.index, index.subj.pos, index.cova, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            n.comb = results$n.comb.new
            comb.acc = results$comb.acc.new
            
          }
          
        }
        
        #         if(n.pos<n.replicates/2){
        #           print("#replicate too small for pos test")
        #         }
        #         if(n.zero<n.replicates/2){
        #           print("#replicate too small for zero test")
        #         }
        #         if(n.comb<n.replicates/2){
        #           print("#replicate too small for comb test")
        #         }
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        score.Rpvalue.comb = (comb.acc+1)/(n.comb+1)
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        while(flag & end.nperm <= n.replicates){
          
          results = .resample.work.pos(X, X.par.index, X1.par.index, index.subj.pos, index.cova, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, end.nperm, n.pos, pos.acc)
          
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          if(flag){
            start.nperm = end.nperm + 1;
            end.nperm = next.end.nperm;
            
          }
          
          if(start.nperm < n.replicates & end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings")) 
            results = .resample.work.pos(X, X.par.index, X1.par.index, index.subj.pos, index.cova, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, n.replicates, n.pos, pos.acc)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            
          }
          
        }
        
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        while(flag & end.nperm <= n.replicates){
          
          
          results = .resample.work.zero(Z, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.zero, zero.acc)
          
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          if(flag){
            start.nperm = end.nperm + 1;
            end.nperm = next.end.nperm;
            
          }
          
          if(start.nperm < n.replicates & end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            
            results = .resample.work.zero(Z, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.zero, zero.acc)
            
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            
          }
          
        }
        
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}

.Score.test.simple2 <- function(Y, case, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    X = cbind(1, case)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects are all case/control
      if( length(table(X1[,2]))<=1 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, 2) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, m-1) 
          df.pos = m-1
        }
        
      }
      
      
      
    }
    
    
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      score.stat.zero = try( .Tstat.zero(Y0, case) )
      if(class(score.stat.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        score.pvalue.zero = 1 - pchisq(score.stat.zero, m0) 
        df.zero = m0
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      score.pvalue.comb = 1 - pchisq(score.stat.comb, (m0 + m -1) ) 
      df.comb = m0 + m -1
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = m - 1
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = m0
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      #set.seed(16)
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = sample(case)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            } 
          }
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }  
          }
          
          if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
            
            score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
            n.comb = n.comb + 1
            if(score.stat.comb.perm >= score.stat.comb){
              comb.acc = comb.acc + 1
              
            }  
          }
          
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        if(n.comb<n.replicates/2){
          print("#replicate too small for comb test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        score.Rpvalue.zero = zero.acc/n.zero
        score.Rpvalue.comb = comb.acc/n.comb
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = sample(case)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            }  
          }
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = sample(case)
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }          
          }
          
          
          
        }
        
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        
        score.Rpvalue.zero = zero.acc/n.zero
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}

.Score.test.simple2.restrict <- function(Y, case, strata, seed=11, resample=FALSE, n.replicates=NULL){
  
  
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    X = cbind(1, case)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects are all case/control
      if( length(table(X1[,2]))<=1 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, 2) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, m-1) 
          df.pos = m-1
        }
        
      }
      
      
      
    }
    
    
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      score.stat.zero = try( .Tstat.zero(Y0, case) )
      if(class(score.stat.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        score.pvalue.zero = 1 - pchisq(score.stat.zero, m0) 
        df.zero = m0
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      score.pvalue.comb = 1 - pchisq(score.stat.comb, (m0 + m -1) ) 
      df.comb = m0 + m -1
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = m - 1
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = m0
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      set.seed(seed)
      
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = .sampling.strata(case, strata)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            } 
          }
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }  
          }
          
          if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
            
            score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
            n.comb = n.comb + 1
            if(score.stat.comb.perm >= score.stat.comb){
              comb.acc = comb.acc + 1
              
            }  
          }
          
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        if(n.comb<n.replicates/2){
          print("#replicate too small for comb test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        score.Rpvalue.zero = zero.acc/n.zero
        score.Rpvalue.comb = comb.acc/n.comb
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = .sampling.strata(case, strata)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            }  
          }
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = .sampling.strata(case, strata)
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }          
          }
          
          
          
        }
        
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        
        score.Rpvalue.zero = zero.acc/n.zero
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}

# one-part test
# no missing value is allowed in any of the inputs
#
# OTU: a matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxa. Column name is mandatory
#
# Tax: a matrix define the taxonomy ranks with each row corresponds to an OTU or a taxa and each column corresponds to a rank (start from the higher taxonomic level). Row name is mandatory and should be consistent with the column name of the OTU table,  Column name should be formated as "Rank1", "Rank2 ()"... etc 
#       If provided, tests will be performed for lineages based on the taxonomic rank. The output contains P-values for all lineages; a list of significant lineages controlling the false discovery rate (based on resampling p-value if resampling test was performed); p-values of the global tests (Fisher- and Simes-combined the p-values for testing lineages).
#       If not provided, one test will be performed with all the OTUs and one p-value will be output
#
# X: a matrix contains covariates for the positive-part test with each column pertains to one variable (pertains to the covariate of interest or the confounders)
#
# X.index: a vector indicate the columns in X for the covariate(s) of interest
#
# Z: a matrix contains covariates for the zero-part test with each column pertains to one variable (pertains to the covariate of interest or the confounders)
#
# Z.index: a vector indicate the columns in X for the covariate(s) of interest
#
# min.depth: keep samples with depths >= min.depth
#
# n.resample: perform asymptotic test is n.resample is null, other perform resampling tests using the specified number of resamplings.
#
# fdr.alpha: false discovery rate for multiple tests on the lineages.

QCAT <- function(OTU, X, X.index, Tax=NULL, min.depth=0, n.perm=NULL, fdr.alpha=0.05){
  
  n.resample = n.perm
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(missing(X.index)){
    X.index = 1:ncol(X)
  }
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ,drop=FALSE]
    OTU = OTU[-remove.subject, ,drop=FALSE]
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep, drop=FALSE]
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      pval = as.matrix( .Score.test(count, X, X.index, resample=FALSE, n.replicates=NULL)$score.pvalue )
      colnames(pval) = "Asymptotic"
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test(count, X, X.index, resample=TRUE, n.replicates=n.resample)
      pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
      names(pval) = c("Asymptotic", "Resampling")
      
    }
    
    return( list(pval=pval) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    
    tax = Tax[keep, ,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    
    n.level = n.rank-1
    
    subtree = NULL
    pval = NULL
    
    for(k in 1:n.level){
      
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              pval = cbind(pval, .Score.test(Y, X, X.index, resample=FALSE, n.replicates=NULL)$score.pvalue)
              
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test(Y, X, X.index, resample=TRUE, n.replicates=n.resample)
              pval = cbind(pval, c(tmp$score.pvalue, tmp$score.Rpvalue) )
              
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    colnames(pval) = subtree
    
    if(is.null(n.resample)){
      
      rownames(pval) = "Asymptotic"
      score.tmp = pval[1,]
      
    }else{
      
      rownames(pval) = c("Asymptotic", "Resampling")
      score.tmp = pval[2,]
    }
    
    #print(pval)
    
    # identify significant lineages
    subtree.tmp = subtree
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
    #fdr.alpha = 0.05
    
    # change 04/17/2016
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1      
    }
    
    sig.lineage = subtree.tmp[reject==1]
    
    
    # perform global test
    global.score.fisher = .F.test(score.tmp)
    global.score.min =  .simes.test(score.tmp)
    global.pval = c(global.score.fisher, global.score.min)
    names(global.pval) = c("Fisher", "Simes")
    
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
  }
  
}

#results.two = QCAT_GEE(count.rff, X, 1, X, 1, tax, n.resample=1000, fdr.alpha=0.05)
QCAT_GEE <- function(OTU, X, X.index, Z, Z.index, Tax=NULL, min.depth=0, n.perm=NULL, fdr.alpha=0.05){
  
  n.resample = n.perm
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table X is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(missing(X.index)){
    X.index = 1:ncol(X)
  }
  
  if(missing(Z)){
    Z = X
  }else{
    
    if(!is.matrix(Z)){
      warning("Covariate table Z is not a matrix")
      Z = as.matrix(Z)  
    }  
  }
  
  
  
  if(missing(Z.index)){
    Z.index = 1:ncol(Z)
  }
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ]
    OTU = OTU[-remove.subject, ,drop=FALSE]
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep,drop=FALSE]
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  Z = cbind(1, Z) # add the intercept term
  Z.index = Z.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      tmp = .Score.test2(count, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL)
      pval.comb = as.matrix( tmp$comb.results$score.pvalue )
      pval.zero = as.matrix( tmp$zero.results$score.pvalue )
      pval.pos = as.matrix( tmp$pos.results$score.pvalue )
      colnames(pval.comb) = "Asymptotic"
      colnames(pval.zero) = "Asymptotic"
      colnames(pval.pos) = "Asymptotic"
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test2(count, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample)
      pval.comb = c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue)
      pval.zero = c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue)
      pval.pos = c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue)
      names(pval.comb) = c("Asymptotic", "Resampling")
      names(pval.zero) = c("Asymptotic", "Resampling")
      names(pval.pos) = c("Asymptotic", "Resampling")
      
    }
    
    return( list(pval=pval.comb, pval.zero=pval.zero, pval.pos=pval.pos) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    tax = Tax[keep,,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    
    n.level = n.rank-1
    
    subtree = NULL
    pval.comb = NULL
    pval.zero = NULL
    pval.pos = NULL
    
    for(k in 1:n.level){
      
      #print(k)
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              tmp = .Score.test2(Y, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL)   
              
              pval.comb = cbind(pval.comb, tmp$comb.results$score.pvalue)
              pval.zero = cbind(pval.zero, tmp$zero.results$score.pvalue)
              pval.pos = cbind(pval.pos, tmp$pos.results$score.pvalue)
              
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test2(Y, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample)   
              
              
              pval.comb = cbind(pval.comb, c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue) )
              pval.zero = cbind(pval.zero, c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue) )
              pval.pos = cbind(pval.pos, c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue) )
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    colnames(pval.comb) = subtree
    colnames(pval.zero) = subtree
    colnames(pval.pos) = subtree
    
    if(is.null(n.resample)){
      
      rownames(pval.comb) = "Asymptotic"
      score.comb.tmp = pval.comb[1,]
      
      rownames(pval.zero) = "Asymptotic"
      score.zero.tmp = pval.zero[1,]
      
      rownames(pval.pos) = "Asymptotic"
      score.pos.tmp = pval.pos[1,]
      
    }else{
      
      rownames(pval.comb) = c("Asymptotic", "Resampling")
      score.comb.tmp = pval.comb[2,]
      
      rownames(pval.zero) = c("Asymptotic", "Resampling")
      score.zero.tmp = pval.zero[2,]
      
      rownames(pval.pos) = c("Asymptotic", "Resampling")
      score.pos.tmp = pval.pos[2,]
    }
    
    #print(pval.comb)
    
    # identify significant lineages
    
    global.pval = NULL
    
    sig.lineage <- vector("list",3)
    names(sig.lineage) <- c("Two-Part", "Zero-Part", "Positive-Part")
    
    for(i in 1:3){
      
      if(i==1){score.tmp = score.comb.tmp}
      if(i==2){score.tmp = score.zero.tmp}
      if(i==3){score.tmp = score.pos.tmp}
      
      subtree.tmp = subtree
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
      #fdr.alpha = 0.05
      
      # change 04/17/2016
      reject = rep(0, m.test)
      tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
      if(length(tmp)>0){
        index.reject = index.p[1:max(tmp)]
        reject[index.reject] = 1      
      }
      
      sig.lineage[[i]] = subtree.tmp[reject==1]
      
      # perform global test
      global.score.min =  .simes.test(score.tmp)
      global.pval = c(global.pval, global.score.min)
      
      
    }
    
    
    
    names(global.pval) = c("Simes_Two-Part", "Simes_Zero-Part", "Simes_Positive-Part")
    
    pval = list(pval.comb, pval.zero, pval.pos)
    names(pval) = c("Two-Part", "Zero-Part", "Positive-Part")
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
  }
  
}

library(abind)
library(CompQuadForm)
########################################
#                                      #
#   QCAT-C test for clustered data     #
#             add in v3.0              #
#                                      #
#                                      #
########################################
.diag2 <- function(x){
  
  if(length(x)>1){
    return(diag(x))
  }else{
    return(as.matrix(x))
    
  }
  
}

.F.test <- function(x){
  
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.simes.test <- function(x){
  
  return( min(length(x) * x/rank(x)) )
  
}


pchisqsum2.davies <- function(Q, lambda){
  
  delta = rep(0, length(lambda))
  acc = 1e-07
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


QCAT.Cluster <- function(ID, OTU, OTU.base=NULL, X, X.index, Tax=NULL, min.depth=0, perm.type=NULL, n.perm=NULL, fdr.alpha=0.05, test="chisq"){
  
  id = ID
  n.resample = n.perm
  if(is.null(perm.type)){
    perm.type="Simple"
  }
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(missing(X.index)){
    X.index = 1:ncol(X)
  }
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ,drop=FALSE]
    OTU = OTU[-remove.subject, ,drop=FALSE]
    if( !is.null(OTU.base) ){
      OTU.base = OTU.base[-remove.subject, ,drop=FALSE]
    }
    
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep, drop=FALSE]
  if( !is.null(OTU.base) ){
    count.base = OTU.base[,keep,drop=FALSE]
  }else{
    count.base = NULL
  }
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      pval = as.matrix( .Score.test.cluster(id, perm.type, count, Y.base=count.base, X, X.index, resample=FALSE, n.replicates=NULL, test=test)$score.pvalue )
      colnames(pval) = "Asymptotic"
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test.cluster(id, perm.type, count, Y.base=count.base, X, X.index, resample=TRUE, n.replicates=n.resample, test=test)
      pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
      names(pval) = c("Asymptotic", "Resampling")
      
    }
    
    return( list(pval=pval) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    
    tax = Tax[keep, ,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    if(!is.null(OTU.base) ){
      W.data.base = data.table(data.frame(tax, t(count.base)))
      otucols.base = names(W.data.base)[-(1:n.rank)]
    }
    
    n.level = n.rank-1
    
    subtree = NULL
    pval = NULL
    
    for(k in 1:n.level){
      
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      if(!is.null(OTU.base) ){
        tt.base = W.data.base[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols.base, by=list( get(Rank.low), get(Rank.high) )]
        W.count.base = data.matrix(tt.base[, otucols.base, with=FALSE])         
      }
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        if(!is.null(OTU.base) ){
          Y.base = t(W.count.base[which(W.tax == level.uni[j]), , drop=FALSE])
        }
        #Y = t(W.count[which(W.tax == "o__Actinomycetales"), , drop=FALSE])
        #Y.base = t(W.count.base[which(W.tax == "o__Actinomycetales"), , drop=FALSE])
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
            if(!is.null(OTU.base) ){
              Y.base = Y.base[, -remove.index, drop=FALSE]
            }
          }
          if(is.null(OTU.base) ){
            Y.base = NULL
          }
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              pval = cbind(pval, .Score.test.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, resample=FALSE, n.replicates=NULL, test=test)$score.pvalue)
              
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, resample=TRUE, n.replicates=n.resample, test=test)
              pval = cbind(pval, c(tmp$score.pvalue, tmp$score.Rpvalue) )
              
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    colnames(pval) = subtree
    
    if(is.null(n.resample)){
      
      rownames(pval) = "Asymptotic"
      score.tmp = pval[1,]
      
    }else{
      
      rownames(pval) = c("Asymptotic", "Resampling")
      score.tmp = pval[2,]
    }
    
    #print(pval)
    
    # identify significant lineages
    subtree.tmp = subtree
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
    #fdr.alpha = 0.05
    
    # change 04/17/2016
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1      
    }
    
    sig.lineage = subtree.tmp[reject==1]
    
    
    # perform global test
    global.score.fisher = .F.test(score.tmp)
    global.score.min =  .simes.test(score.tmp)
    global.pval = c(global.score.fisher, global.score.min)
    names(global.pval) = c("Fisher", "Simes")
    
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
  }
  
}


.Ei.beta.XX <- function(m, p.XX, beta, XX.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  idx1 = 0
  idx2 = 0
  
  for(j in 1:(m-1)){
    
    idx1 = idx2 + 1
    idx2 = idx2 + p.XX[j]
    
    Ei.out[j] = exp(beta[idx1:idx2] %*% XX.i[,1:p.XX[j],j])
    
    
  }
  
  
  Ei.out[m] = 1
  
  
  
  
  return (Ei.out)
}


.fun.neg.loglik.beta.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      loglik = loglik + Y[i,] %*% log(P.i)
    }
    
  }
  
  return (-loglik)
  
}


.fun.neg.score.beta.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Score.beta = Score.beta + kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta = Score.beta + tmp
    }
    
    return (-Score.beta)
  }
  
  
  
}


.fun.score.i.beta.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      
      #Score.beta.i[i,] =  kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1) )
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta.i[i,] = tmp
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}


.fun.hessian.beta.XX <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  #n.beta = (m-1)*p
  n.beta = sum(p.XX)
  
  #idx.end = (1:(m-1))*p 
  #idx.start = idx.end - (p-1)
  
  idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
  idx.end = idx.start + p.XX - 1
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      #Hessian.beta = Hessian.beta + kronecker( tmp.beta, ( X[i,] %o% X[i,] ) ) 
      Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
      for(j in 1:(m-1)){
        for(k in j:(m-1)){
          Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = tmp.beta[j,k] * (XX[i,1:p.XX[j],j] %o% XX[i,1:p.XX[k],k])
          Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
        }
      }
      Hessian.beta = Hessian.beta + Hessian.tmp
      
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }
    
    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}


# 03/28/2018 modify to allow for uneven #obs across clusters
.sampling.strata.uneven <- function(var, id){
  
  id.unique = unique(id)
  n.cluster = length(id.unique)
  var.sample = var
  for(i in 1:n.cluster){
    index = which(id==id.unique[i])
    var.sample[index] = var.sample[sample(index)]
  }
  return(var.sample)
  
}


########## functions defined for clustered data ##########
.Score.test.stat.cluster <- function(id, Y, XX, p.XX, X.par.index){
  
  #p = dim(XX)[2]
  
  nY = rowSums(Y)
  
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = sum(p.XX)
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NA
    
  }else{
    
    #XX.reduce = XX[,-X.par.index, , drop=FALSE]
    
    p.XX.reduce = p.XX - length(X.par.index) 
    XX.reduce = array(NA,dim=c(n,dim(XX)[2],m-1))
    for(j in 1:(m-1)){
      XX.reduce[,1:p.XX.reduce[j],j] = XX[,(1:p.XX[j])[-X.par.index],j]
    }    
    
    
    #par.interest.index.beta =  kronecker( ((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
    par.interest.index.beta = NULL
    tmp = 0
    for(j in 1:(m-1)){
      par.interest.index.beta = c(par.interest.index.beta, tmp + X.par.index)
      tmp = sum(p.XX[1:j])
    }
    
    n.par.interest.beta = length(par.interest.index.beta) 
    
    #beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
    beta.ini.reduce = rep(0, sum(p.XX.reduce) )    
    
    
    data.reduce.beta = list(Y=Y, XX=XX.reduce, p.XX=p.XX.reduce)
    
    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0 
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.XX, gr=.fun.neg.score.beta.XX, data = data.reduce.beta, method="BFGS")$par
    
    
    data.beta = list(Y=Y, XX=XX, p.XX=p.XX)
    Score.reduce.beta = .fun.score.i.beta.XX(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    # S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, idx.start])
    tmp = .fun.hessian.beta.XX(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list  
    
    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
    
    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    
    
    B2 =  matrix(0, n.beta, n.beta)
    
    # 03/12/2018   modify for clustered data (same #observation for each cluster for now)
    # 03/28/2018   modify to allow for uneven #obs across clusters
    id.unique = unique(id)
    n.cluster = length(id.unique)
    # TT = n/length(unique(id))  
    # n.cluster = n/TT
    # idx1 = 1
    # idx2 = TT
    for(i in 1:n.cluster){
      
      #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
      Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id==id.unique[i]),,drop=FALSE] )
      B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
      
      #idx1 = idx2 + 1
      #idx2 = idx2 + TT
      
    }
    
    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A
    
    
  }
  
  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list, U=A, V=B )   )
  
  
}


.Score.test.cluster <- function(id, perm.type, Y, Y.base=NULL, X, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, test="chisq"){
  
  n = nrow(Y)
  m = ncol(Y) 
  
  
  if(is.null(Y.base)){
    
    XX = array(NA,dim=c(n,ncol(X),m-1))
    
    
    for(j in 1:(m-1)){
      XX[,,j] = X
      
    }
    
    p.XX = rep(ncol(X), m-1)
    
  }else{
    
    XX = array(NA,dim=c(n,ncol(X)+1,m-1))
    
    p.XX = rep(0, m-1)
    
    for(j in 1:(m-1)){
      
      XX[,1:ncol(X),j] = X
      
      if(length(table(Y.base[,j]))==1){
        p.XX[j] = ncol(X)
      }else{
        p.XX[j] = ncol(X)+1
      }
      
    }
    
    XX[,(ncol(X)+1),] = Y.base[,-m,drop=FALSE]
  }
  
  #p = dim(XX)[2]
  
  nY = rowSums(Y)
  
  ## remove 03/28/2016
  #   nY0.index = which(nY==0)
  #   if(length(nY0.index)>0){
  #     Y = Y[-nY0.index, , drop=FALSE]
  #     X = X[-nY0.index, , drop=FALSE]
  #   }
  
  #n.beta = (m - 1)*p
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NULL
    score.pvalue.beta = NA
    n.par.interest.beta = 0
    
  }else{
    
    tmp.one = try( .Score.test.stat.cluster(id, Y, XX, p.XX, X.par.index) )
    #tmp.one = ( .Score.test.stat.cluster(id, Y, XX, p.XX, X.par.index) )
    if(class(tmp.one) == "try-error"){
      
      score.stat.beta = NA
      score.pvalue.beta = NA
      n.par.interest.beta = NA
      
    }else{
      
      n.par.interest.beta = (m-1)*length(X.par.index)
      
      if(test=="chisq"){
        score.stat.beta = tmp.one$score.stat.beta
        score.pvalue.beta = 1 - pchisq(score.stat.beta,  n.par.interest.beta) 
      }
      if(test=="vc"){
        inv.V = ginv(tmp.one$V)
        score.stat.beta =  tmp.one$U %*% inv.V %*% inv.V %*% tmp.one$U
        lambda = eigen(inv.V)$values
        score.pvalue.beta = pchisqsum2.davies(Q=score.stat.beta, lambda=lambda)$p 
      }
      
    }
    
    
  }
  
  
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue.beta, df =n.par.interest.beta)
  #  print(score.stat.beta)
  
  if(resample){
    ########################################
    #                                      #
    #        Resampling Score Test         #
    # change 07/02/2016 for adaptive resampling
    # change 07/17/2016 for adding vc test #
    #                                      #
    ########################################  
    set.seed(seed)
    if(!is.na(score.stat.beta)){
      
      n.one = 0
      one.acc = 0
      
      start.nperm = 1;
      end.nperm = min(100,n.replicates);
      flag = 1
      ### 07/19/2018, modify while loop
      while(flag){
        
        results = .resample.work.one.cluster(id, perm.type, XX, p.XX, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, end.nperm, n.one, one.acc, test=test)
        n.one = results$n.one.new
        one.acc = results$one.acc.new
        flag = results$flag
        next.end.nperm = results$next.end.nperm
        
        
        start.nperm = end.nperm + 1;
        end.nperm = next.end.nperm;
        
        if(start.nperm >= n.replicates){
          flag = 0
        }else if(end.nperm > n.replicates){ 
          #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
          results = .resample.work.one.cluster(id, perm.type, XX, p.XX, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, n.replicates, n.one, one.acc, test=test)
          n.one = results$n.one.new
          one.acc = results$one.acc.new
          flag = 0
        }
        
      }
      
      
      #      print(paste("Final number of resamplings: ", n.one) )
      #      if(n.one<end.nperm/2){
      #        print("Number of resamplings too small for one-part test")
      #      }
      
      tmp = (one.acc+1)/(n.one+1)
      
      #print(n.one)
      #print(one.acc)
      
    }else{
      
      tmp = NA
    }
    
    beta.results = c(beta.results, score.Rpvalue = tmp)
    
    
  }
  
  return(beta.results)
  
  
}


.resample.work.one.cluster <- function(id, perm.type, XX, p.XX, X.par.index, score.stat.beta, S.beta.list, I.beta.list, start.nperm, end.nperm, n.one, one.acc, test){
  
  n = dim(XX)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  XX.interest = XX[,X.par.index, ,drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    
    for(i in 1:n.cluster){
      tmp = XX.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmp)[1])
      if(dim(unique(tmp))[1]!=1){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        XX.interest.unique = unique(tmp)
      }else{
        XX.interest.unique = abind(XX.interest.unique, unique(tmp), along=1)
      }
      
      
    }
    
  }
  
  n.one.new = n.one
  one.acc.new = one.acc
  
  for(k in start.nperm:end.nperm){
    
    # 03/28/2018   modify to allow for uneven #obs across clusters
    # 03/29/2018   modify to allow taxon-specific covariate
    if(perm.type=="BTW"){
      
      #XX.interest.perm = NULL
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        XX.interest.perm.i = array(NA, dim=c(nobs, dim(XX.interest.unique)[2], dim(XX.interest.unique)[3]))
        for(j in 1:dim(XX.interest.unique)[3]){
          XX.interest.perm.i[,,j] =  matrix(rep(XX.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        if(i==1){
          XX.interest.perm = XX.interest.perm.i
        }else{
          XX.interest.perm = abind(XX.interest.perm, XX.interest.perm.i, along =1)
        }
        
      }
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.interest.perm
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
    }
    
    
    
    #X.perm = X
    #X.perm[,X.par.index] = X.perm[sample.idx,X.par.index]
    
    score.stat.beta.perm = try( .Score.test.stat.4Gresampling.cluster(id, XX.perm, p.XX, X.par.index, S.beta.list, I.beta.list, test=test) )
    #score.stat.beta.perm = ( .Score.test.stat.4Gresampling.cluster(id, XX.perm, p.XX, X.par.index, S.beta.list, I.beta.list) )
    
    if(class(score.stat.beta.perm)[1] != "try-error"){
      
      n.one.new = n.one.new + 1
      if(score.stat.beta.perm >= score.stat.beta){
        one.acc.new = one.acc.new + 1
        
      }  
    }
  }
  
  ### modify 07/19/2018
  if(one.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.one.new=n.one.new, one.acc.new=one.acc.new, flag=flag, next.end.nperm=next.end.nperm))
  
}


.Score.test.stat.4Gresampling.cluster <- function(id, XX.perm, p.XX, X.par.index, S.beta.list, I.beta.list, test){
  
  n = dim(XX.perm)[1]
  #p = dim(XX.perm)[2]
  
  m.beta = length(S.beta.list[[1]])
  
  #idx.end = (1:m.beta)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.XX[-m.beta]))+1
  idx.end = idx.start + p.XX - 1
  
  #n.par.interest.beta = m.beta
  #n.beta = m.beta*p
  n.beta = sum(p.XX)
  #par.interest.index.beta =  kronecker( ((0:(m.beta-1))*p), rep(1,length(X.par.index))) + X.par.index
  par.interest.index.beta = NULL
  tmp = 0
  for(j in 1:m.beta){
    par.interest.index.beta = c(par.interest.index.beta, tmp + X.par.index)
    tmp = sum(p.XX[1:j])
  }
  
  n.par.interest.beta = length(par.interest.index.beta) 
  
  
  Score.reduce.beta.perm = matrix(0, n, n.beta )
  Hess.reduce.beta.perm = matrix(0, n.beta, n.beta )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    #Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(X.perm[i,], ncol=1))  
    tmp = NULL
    for(j in 1:m.beta){
      tmp = c(tmp, S.beta.list[[i]][j] * XX.perm[i,1:p.XX[j],j] )
    }
    Score.reduce.beta.perm[i,] = tmp
    
    #Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  X.perm[i,] %o% X.perm[i,] ) )
    
    Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
    for(j in 1:m.beta){
      for(k in j:m.beta){
        Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = I.beta.list[[i]][j,k] * (XX.perm[i,1:p.XX[j],j] %o% XX.perm[i,1:p.XX[k],k])
        Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
      }
    }
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + Hessian.tmp
    
    
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id.unique = unique(id)
  n.cluster = length(id.unique)
  #TT = n/length(unique(id))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[which(id==id.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.beta.perm.reorg.cluster %o%  Score.reduce.beta.perm.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  
  B = B1 %*% B2 %*% t(B1)
  
  if(test=="chisq"){
    score.stat.beta.perm = A %*% ginv(B) %*% A
  }
  
  if(test=="vc"){
    
    inv.B = ginv(B)
    score.stat.beta.perm = A %*% inv.B %*% inv.B %*% A
  }
  
  
  return(score.stat.beta.perm)
  
  
  
}




########################################
#                                      #
#     QCAT_GEE for clustered data      #
#          two-part test               #
#                                      #
#                                      #
########################################

QCAT_GEE.Cluster <- function(ID, OTU, OTU.base=NULL, X, X.index, Z, Z.index, Tax=NULL, min.depth=0, perm.type=NULL, n.perm=NULL, fdr.alpha=0.05, test="chisq"){
  
  id = ID
  n.resample = n.perm
  if(test=="mix"){
    test="vc"
  }
  if(is.null(perm.type)){
    perm.type="Simple"
  }
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(missing(X.index)){
    X.index = 1:ncol(X)
  }
  
  if(missing(Z)){
    Z = X
  }else{
    
    if(!is.matrix(Z)){
      warning("Covariate table Z is not a matrix")
      Z = as.matrix(Z)  
    }  
  }
  
  
  
  if(missing(Z.index)){
    Z.index = 1:ncol(Z)
  }
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ]
    OTU = OTU[-remove.subject, ,drop=FALSE]
    if( !is.null(OTU.base) ){
      OTU.base = OTU.base[-remove.subject, ,drop=FALSE]
    }
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep,drop=FALSE]
  if( !is.null(OTU.base) ){
    count.base = OTU.base[,keep,drop=FALSE]
  }else{
    count.base = NULL
  }
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  Z = cbind(1, Z) # add the intercept term
  Z.index = Z.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      tmp = .Score.test2.cluster(id, perm.type, count, Y.base=count.base, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL, test=test)
      
      pval.zero = as.matrix( tmp$zero.results$score.pvalue )
      pval.pos = as.matrix( tmp$pos.results$score.pvalue )
      colnames(pval.zero) = "Asymptotic"
      colnames(pval.pos) = "Asymptotic"
      
      if(test=="chisq"){
        pval.comb = as.matrix( tmp$comb.results$score.pvalue )
        colnames(pval.comb) = "Asymptotic"
      }
      if(test=="vc"){
        #pval.comb = as.matrix( c(tmp$comb.results$scoreF.pvalue, tmp$comb.results$scoreS.pvalue), ncol=2 )
        #colnames(pval.comb) = c( "Asymptotic.Fisher", "Asymptotic.Simes" )   
        
        pval.comb = as.matrix( tmp$comb.results$score.pvalue )
        colnames(pval.comb) = "Asymptotic"  
      }
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test2.cluster(id, perm.type, count, Y.base=count.base, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, test=test)
      pval.zero = c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue)
      pval.pos = c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue)
      names(pval.zero) = c("Asymptotic", "Resampling")
      names(pval.pos) = c("Asymptotic", "Resampling")
      
      if(test=="chisq"){
        pval.comb = c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue)
        names(pval.comb) = c("Asymptotic", "Resampling")
      }
      if(test=="vc"){
        #pval.comb = c(tmp$comb.results$scoreF.pvalue, tmp$comb.results$scoreS.pvalue, tmp$comb.results$scoreF.Rpvalue, tmp$comb.results$scoreS.Rpvalue )
        #names(pval.comb) = c( "Asymptotic.Fisher", "Asymptotic.Simes", "Resampling.Fisher", "Resampling.Simes" )  
        pval.comb = c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue)
        names(pval.comb) = c("Asymptotic", "Resampling")        
        
      }
      
    }
    
    return( list(pval=pval.comb, pval.zero=pval.zero, pval.pos=pval.pos) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    tax = Tax[keep,,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    if(!is.null(OTU.base) ){
      W.data.base = data.table(data.frame(tax, t(count.base)))
      otucols.base = names(W.data.base)[-(1:n.rank)]
    }
    
    n.level = n.rank-1
    
    subtree = NULL
    pval.comb = NULL
    #pval.combF = NULL
    #pval.combS = NULL
    pval.zero = NULL
    pval.pos = NULL
    
    for(k in 1:n.level){
      
      #print(k)
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      if(!is.null(OTU.base) ){
        tt.base = W.data.base[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols.base, by=list( get(Rank.low), get(Rank.high) )]
        W.count.base = data.matrix(tt.base[, otucols.base, with=FALSE])         
      }      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        if(!is.null(OTU.base) ){
          Y.base = t(W.count.base[which(W.tax == level.uni[j]), , drop=FALSE])
        }
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
            if(!is.null(OTU.base) ){
              Y.base = Y.base[, -remove.index, drop=FALSE]
            }
          }
          if(is.null(OTU.base) ){
            Y.base = NULL
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              tmp = .Score.test2.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL, test=test)   
              
              
              pval.zero = cbind(pval.zero, tmp$zero.results$score.pvalue)
              pval.pos = cbind(pval.pos, tmp$pos.results$score.pvalue)
              if(test=="chisq"){
                pval.comb = cbind(pval.comb, tmp$comb.results$score.pvalue)
              }
              if(test=="vc"){
                pval.comb = cbind(pval.comb, tmp$comb.results$score.pvalue)
                #pval.combF = cbind(pval.combF, tmp$comb.results$scoreF.pvalue)
                #pval.combS = cbind(pval.combS, tmp$comb.results$scoreS.pvalue)
              }                
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test2.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, test=test)   
              
              pval.zero = cbind(pval.zero, c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue) )
              pval.pos = cbind(pval.pos, c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue) )
              
              if(test=="chisq"){
                pval.comb = cbind(pval.comb, c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue) )
              }
              if(test=="vc"){
                pval.comb = cbind(pval.comb, c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue) )
                #pval.combF = cbind(pval.combF, c(tmp$comb.results$scoreF.pvalue, tmp$comb.results$scoreF.Rpvalue) )
                #pval.combS = cbind(pval.combS, c(tmp$comb.results$scoreS.pvalue, tmp$comb.results$scoreS.Rpvalue) )
              }
              
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    
    colnames(pval.zero) = subtree
    colnames(pval.pos) = subtree
    if(test=="chisq"){
      colnames(pval.comb) = subtree
    }
    if(test=="vc"){
      colnames(pval.comb) = subtree
      #colnames(pval.combF) = subtree
      #colnames(pval.combS) = subtree
    }
    
    if(is.null(n.resample)){
      
      rownames(pval.zero) = "Asymptotic"
      score.zero.tmp = pval.zero[1,]
      
      rownames(pval.pos) = "Asymptotic"
      score.pos.tmp = pval.pos[1,]
      
      if(test=="chisq"){
        rownames(pval.comb) = "Asymptotic"
        score.comb.tmp = pval.comb[1,]
      }
      
      if(test=="vc"){
        rownames(pval.comb) = "Asymptotic"
        score.comb.tmp = pval.comb[1,]
        #rownames(pval.combF) = "Asymptotic"
        #score.combF.tmp = pval.combF[1,]
        #rownames(pval.combS) = "Asymptotic"
        #score.combS.tmp = pval.combS[1,]
      }
      
    }else{
      
      rownames(pval.zero) = c("Asymptotic", "Resampling")
      score.zero.tmp = pval.zero[2,]
      
      rownames(pval.pos) = c("Asymptotic", "Resampling")
      score.pos.tmp = pval.pos[2,]
      
      if(test=="chisq"){
        rownames(pval.comb) = c("Asymptotic", "Resampling")
        score.comb.tmp = pval.comb[2,]
      }
      
      if(test=="vc"){
        rownames(pval.comb) = c("Asymptotic", "Resampling")
        score.comb.tmp = pval.comb[2,]
        #rownames(pval.combF) = c("Asymptotic", "Resampling")
        #score.combF.tmp = pval.combF[2,]
        #rownames(pval.combS) = c("Asymptotic", "Resampling")
        #score.combS.tmp = pval.combS[2,]
      }  
      
    }
    
    #print(pval.comb)
    
    # identify significant lineages
    
    global.pval = NULL
    
    sig.lineage <- vector("list",3)
    names(sig.lineage) <- c("Zero-Part", "Positive-Part", "Two-Part")
    
    #if(test=="chisq"){SS = 3}
    #if(test=="vc"){SS = 4}
    SS = 3
    for(i in 1:SS){
      
      if(i==1){score.tmp = score.zero.tmp}
      if(i==2){score.tmp = score.pos.tmp}
      if(i==3){score.tmp = score.comb.tmp} 
      # if(test=="chisq" & i==3){
      #   score.tmp = score.comb.tmp
      # }
      # if(test=="vc" & i==3){
      #   score.tmp = score.combF.tmp
      # }
      # if(test=="vc" & i==4){
      #   score.tmp = score.combS.tmp
      # }      
      
      subtree.tmp = subtree
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
      #fdr.alpha = 0.05
      
      # change 04/17/2016
      reject = rep(0, m.test)
      tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
      if(length(tmp)>0){
        index.reject = index.p[1:max(tmp)]
        reject[index.reject] = 1      
      }
      
      sig.lineage[[i]] = subtree.tmp[reject==1]
      
      # perform global test
      global.score.min =  .simes.test(score.tmp)
      global.pval = c(global.pval, global.score.min)
      
      
    }
    
    
    if(test=="chisq"){
      names(global.pval) = c("Simes_Zero-Part", "Simes_Positive-Part", "Simes_Two-Part")
    }
    
    if(test=="vc"){
      #names(global.pval) = c("Simes_Zero-Part", "Simes_Positive-Part", "SimesF_Two-Part", "SimesS_Two-Part")
      names(global.pval) = c("Simes_Zero-Part", "Simes_Positive-Part", "Simes_Two-Part")
    }
    
    pval = list(pval.zero, pval.pos, pval.comb)
    names(pval) = c("Zero-Part", "Positive-Part", "Two-Part")
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
    
  }
  
}


.Ei.beta.pos.XX <- function(m, p.XX, beta, XX.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  I.i = as.numeric(Y.i>0)
  
  idx1 = 0
  idx2 = 0
  
  for(j in 1:(m-1)){
    
    idx1 = idx2 + 1
    idx2 = idx2 + p.XX[j]
    
    Ei.out[j] = I.i[j] * exp(beta[idx1:idx2] %*% XX.i[,1:p.XX[j],j])
    
  }
  
  
  Ei.out[m] = I.i[m]
  
  
  return (Ei.out)
}


.fun.neg.loglik.beta.pos.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      index = which(Y[i,]>0)
      loglik = loglik + Y[i,index] %*% log(P.i[index])
      #       if(is.na(loglik)){
      #         print(i); break;
      #       }
    }
    
  }
  
  return (-loglik)
  
}


.fun.neg.score.beta.pos.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Score.beta = Score.beta + kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta = Score.beta + tmp
      
      
    }
    
    return (-Score.beta)
  }
  
  
  
}


.fun.score.i.beta.pos.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      
      #Score.beta.i[i,] =  kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1)) 
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta.i[i,] = tmp
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}


.fun.hessian.beta.pos.XX <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  #n.beta = (m-1)*p
  n.beta = sum(p.XX)
  
  #idx.end = (1:(m-1))*p 
  #idx.start = idx.end - (p-1)
  
  idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
  idx.end = idx.start + p.XX - 1
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      #Hessian.beta = Hessian.beta + kronecker(tmp.beta, ( X[i,] %o% X[i,] ))
      Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
      for(j in 1:(m-1)){
        for(k in j:(m-1)){
          Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = tmp.beta[j,k] * (XX[i,1:p.XX[j],j] %o% XX[i,1:p.XX[k],k])
          Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
        }
      }
      Hessian.beta = Hessian.beta + Hessian.tmp
      
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }
    
    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}


.Pi.alpha.ZZ<-function(m, p.ZZ, alpha, ZZ.i){
  
  Pi.out = rep(NA,m)
  
  idx1 = 0
  idx2 = 0
  
  for(j in 1:m){
    
    idx1 = idx2 + 1
    idx2 = idx2 + p.ZZ[j]
    
    tmp = exp(alpha[idx1:idx2] %*% ZZ.i[,1:p.ZZ[j],j])
    if(is.infinite(tmp)){
      Pi.out[j] = 1
    }else{
      Pi.out[j] = tmp/(tmp + 1)
    }
    
  }
  
  
  return (Pi.out)    
}


.fun.score.i.alpha.ZZ <- function(alpha, data, save.list=FALSE){
  
  Y = data$Y; ZZ = data$ZZ; p.ZZ = data$p.ZZ
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(ZZ)[2]
  
  vA.list = list()
  Vinv.list = list()
  VY.list = list()
  
  #n.alpha = m*p
  n.alpha = sum(p.ZZ)
  
  #idx.end = (1:m)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.ZZ[-m]))+1
  idx.end = idx.start + p.ZZ - 1
  
  if(length(alpha)!=n.alpha){
    
    warning("Dim of initial alpha does not match the dim of covariates")
    
  }else{
    
    Score.alpha.i = matrix(0, n, n.alpha)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha.ZZ(m, p.ZZ, alpha, ZZ[i,,,drop=FALSE])
      vA.tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(vA.tmp) 
      #t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      t.D.i = matrix(0, nrow=n.alpha, ncol=m)
      for(j in 1:m){
        t.D.i[idx.start[j]:idx.end[j], j] = A.i[j,j] * ZZ[i,1:p.ZZ[j],j] 
      }
      
      
      V.i = A.i # independent cor structure
      
      tmp.V.i = ginv(V.i)
      tmp.VY = tmp.V.i %*% (Y[i,] - Pi.i)
      Score.alpha.i[i,] = t.D.i %*% tmp.VY
      
      if(save.list){
        vA.list[[i]] = vA.tmp
        Vinv.list[[i]] = tmp.V.i
        VY.list[[i]] = tmp.VY
      }
    }
    
    
  }
  
  if(save.list){
    
    return ( list(Score.alpha=Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list=VY.list) )
    
  }else{
    
    return (Score.alpha.i)
  }  
  
}


.fun.hessian.alpha.ZZ <- function(alpha, data){
  
  Y = data$Y; ZZ = data$ZZ; p.ZZ = data$p.ZZ
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(ZZ)[2]
  #n.alpha = m*p
  n.alpha = sum(p.ZZ)
  
  #idx.end = (1:m)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.ZZ[-m]))+1
  idx.end = idx.start + p.ZZ - 1
  
  if(length(alpha)!=n.alpha){
    print("Waring: dim of alpha is not the same as alpha\n")
    
  }else{
    
    Hessian.alpha = matrix(0, nrow=n.alpha, ncol=n.alpha)
    nY = rowSums(Y)
    
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha.ZZ(m, p.ZZ, alpha, ZZ[i,,,drop=FALSE])
      tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(tmp)
      #t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      t.D.i = matrix(0, nrow=n.alpha, ncol=m)
      for(j in 1:m){
        t.D.i[idx.start[j]:idx.end[j], j] = A.i[j,j] * ZZ[i,1:p.ZZ[j],j] 
      }      
      V.i = A.i # independent cor structure        
      
      
      Hessian.alpha = Hessian.alpha + t.D.i %*% ginv(V.i) %*% t(t.D.i) 
      
      
    }
    
    return (Hessian.alpha)
    
    
  }
  
  
}


########## functions defined for clustered data ##########

.Score.test2.cluster <- function(id, perm.type, Y, Y.base=NULL, X, X.par.index, Z, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL, test){
  
  n = nrow(Y)
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    if(test=="chisq"){
      comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    }
    if(test=="vc"){
      #comb.results = list(score.stat = NA, scoreF.pvalue = NA, scoreS.pvalue = NA, df = NA )
      comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    }    
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      if(test=="chisq"){
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      if(test=="vc"){
        #comb.results = c(comb.results, scoreF.Rpvalue = NA, scoreS.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }  
      
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE]
      if(!is.null(Y.base)){
        Y.base = Y.base[, -remove.index, drop=FALSE]
        
      }
    }
    
    m = ncol(Y)
    
    if(is.null(Y.base)){
      
      XX = array(NA,dim=c(n,ncol(X),m-1))
      
      
      for(j in 1:(m-1)){
        XX[,,j] = X
        
      }
      
      p.XX = rep(ncol(X), m-1)
      
    }else{
      
      XX = array(NA,dim=c(n,ncol(X)+1,m-1))
      
      p.XX = rep(0, m-1)
      
      for(j in 1:(m-1)){
        
        XX[,1:ncol(X),j] = X
        
        if(length(table(Y.base[,j]))==1){
          p.XX[j] = ncol(X)
        }else{
          p.XX[j] = ncol(X)+1
        }
        
      }
      
      XX[,(ncol(X)+1),] = Y.base[,-m,drop=FALSE]
    }
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      XX1 = XX[index.subj.pos, , ,drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects have the same values of covariates (e.g., all case/control)
      #index.cova = 1 + which(apply(XX1[,-1,1,drop=FALSE], 2, function(x) length(table(x)) ) > 1) # index of valid covariates
      index.cova.ls = NULL
      for(j in 1:(m-1)){
        index.cova.ls = c(index.cova.ls, list(1 + which(apply(XX1[,-1,j,drop=FALSE], 2, function(x) length(table(x)) ) > 1) ) )# index of valid covariates
      }
      X1.par.index.ava.ls = lapply(index.cova.ls, function(x) as.numeric( !is.na(match(X.par.index, x )) ) ) # use 0/1 to indicate the index is still availiable or not
      # if no covariate left; even if have covariate left, they are not covariate of interest; no taxa
      if( sum(sapply(index.cova.ls, length))<1 |   sum(sapply(X1.par.index.ava.ls, sum))<1 | m<=1 ){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        index.cova.ls = lapply(index.cova.ls, function(x) c(1, x))
        #XX1 = XX1[,c(1, index.cova), ,drop=FALSE] 
        #tmp = match(X.par.index, index.cova) + 1
        #X1.par.index = tmp[!is.na(tmp)]
        #d1.pos = length(X1.par.index)
        X1.par.index.ls = lapply(index.cova.ls, function(x) match(X.par.index, x) )
        d1.pos = sum( sapply(X1.par.index.ls, length) )
        
        tmp.pos = try( .Score.test.stat.pos.cluster(id[index.subj.pos], Y1, XX1, index.cova.ls, X1.par.index.ls) )
        #tmp.pos = ( .Score.test.stat.pos.cluster(id[index.subj.pos], Y1, XX1, index.cova.ls, X1.par.index.ls) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          df.pos = d1.pos 
          
          if(test=="chisq"){
            score.stat.pos = tmp.pos$score.stat.beta
            score.pvalue.pos = 1 - pchisq(score.stat.pos, d1.pos ) 
          }
          
          if(test=="vc"){
            inv.V = ginv(tmp.pos$V)
            score.stat.pos =  tmp.pos$U %*% inv.V %*% inv.V %*% tmp.pos$U
            lambda = eigen(inv.V)$values
            score.pvalue.pos = pchisqsum2.davies(Q=score.stat.pos, lambda=lambda)$p 
            
          }
          
        }
        
      }    
      
    }
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    if(!is.null(Y.base)){
      Y0.base = Y.base
      Y0.base[Y.base==0] = 1  
      Y0.base[Y.base>0] = 0
    }
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
        if(!is.null(Y.base)){
          Y0.base = Y0.base[, -remove.index, drop=FALSE]
          
        }
      }
      m0 = ncol(Y0)
      
      if(is.null(Y.base)){
        
        ZZ = array(NA,dim=c(n,ncol(Z),m0))
        
        
        for(j in 1:m0){
          ZZ[,,j] = Z
          
        }
        
        p.ZZ = rep(ncol(Z), m0)
        
      }else{
        
        ZZ = array(NA,dim=c(n,ncol(Z)+1,m0))
        
        p.ZZ = rep(0, m0)
        
        for(j in 1:m0){
          
          ZZ[,1:ncol(Z),j] = Z
          
          if(length(table(Y0.base[,j]))==1){
            p.ZZ[j] = ncol(Z)
          }else{
            p.ZZ[j] = ncol(Z)+1
          }
          
        }
        
        ZZ[,(ncol(Z)+1),] = Y0.base
      }
      
      
      
      tmp.zero = try( .Score.test.stat.zero.cluster(id, Y0, ZZ, p.ZZ, Z.par.index, "independence") )
      #tmp.zero = ( .Score.test.stat.zero.cluster(id, Y0, ZZ, p.ZZ, Z.par.index, "independence") )
      if(class(tmp.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        df.zero = tmp.zero[[1]]
        if(test=="chisq"){
          score.stat.zero = tmp.zero[[2]]
          score.pvalue.zero = tmp.zero[[3]]
        }
        if(test=="vc"){
          ################# 07/19/2018, for the zero part, keep the chisq test
          #inv.V = ginv(tmp.zero$V)
          #score.stat.zero =  tmp.zero$U %*% inv.V %*% inv.V %*% tmp.zero$U
          #lambda = eigen(inv.V)$values
          #score.pvalue.zero = pchisqsum2.davies(Q=score.stat.zero, lambda=lambda)$p 
          score.stat.zero = tmp.zero[[2]]
          score.pvalue.zero = tmp.zero[[3]]
        }        
        
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      df.comb = df.zero + df.pos
      
      if(test=="chisq"){
        score.stat.comb = score.stat.pos + score.stat.zero
        score.pvalue.comb = 1 - pchisq(score.stat.comb, df.comb ) 
      }
      if(test=="vc"){
        
        #score.pvalue.combF =  .F.test(c(score.pvalue.pos, score.pvalue.zero))
        #score.pvalue.combS =  .simes.test(c(score.pvalue.pos, score.pvalue.zero))
        score.pvalue.comb =  .F.test(c(score.pvalue.pos, score.pvalue.zero))
        
      }
      
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      df.comb = df.pos
      if(test=="chisq"){
        score.stat.comb = score.stat.pos
        score.pvalue.comb = score.pvalue.pos
      }
      if(test=="vc"){
        #score.pvalue.combF = score.pvalue.pos
        #score.pvalue.combS = score.pvalue.pos
        score.pvalue.comb = score.pvalue.pos
      }
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      df.comb = df.zero
      if(test=="chisq"){
        score.stat.comb = score.stat.zero
        score.pvalue.comb = score.pvalue.zero  
      }
      if(test=="vc"){
        score.pvalue.comb = score.pvalue.zero  
        #score.pvalue.combF = score.pvalue.zero
        #score.pvalue.combS = score.pvalue.zero
      }
      
    }else{
      df.comb = NA
      if(test=="chisq"){
        score.stat.comb = NA
        score.pvalue.comb =  NA
      }
      if(test=="vc"){
        score.pvalue.comb =  NA
        #score.pvalue.combF = NA
        #score.pvalue.combS = NA
      }      
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero)
    if(test=="chisq"){
      comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb )
    }
    if(test=="vc"){
      comb.results = list(score.pvalue = score.pvalue.comb )
      #comb.results = list(scoreF.pvalue = score.pvalue.combF, scoreS.pvalue = score.pvalue.combS )
    }
    
    if(resample){
      
      #print("simulated stat:")
      set.seed(seed)
      
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        
        ### 07/19/2018, modify while loop
        while(flag){
          
          results = .resample.work.two.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, ZZ, p.ZZ, Z.par.index, index.subj.pos, index.cova.ls, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc, test=test)
          
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          if(test=="chisq"){
            n.comb = results$n.comb.new
            comb.acc = results$comb.acc.new
          }
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
          if(start.nperm >= n.replicates){
            flag = 0
          }else if(end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            results = .resample.work.two.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, ZZ, p.ZZ, Z.par.index, index.subj.pos, index.cova.ls, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc, test=test)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            if(test=="chisq"){
              n.comb = results$n.comb.new
              comb.acc = results$comb.acc.new
            }
            
            flag = 0
          }
          
        }
        
        #         if(n.pos<n.replicates/2){
        #           print("#replicate too small for pos test")
        #         }
        #         if(n.zero<n.replicates/2){
        #           print("#replicate too small for zero test")
        #         }
        #         if(n.comb<n.replicates/2){
        #           print("#replicate too small for comb test")
        #         }
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        if(test=="chisq"){
          score.Rpvalue.comb = (comb.acc+1)/(n.comb+1)
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        }
        if(test=="vc"){
          score.Rpvalue.comb = .F.test(c(score.Rpvalue.pos, score.Rpvalue.zero))
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
          #score.Rpvalue.combF = .F.test(c(score.Rpvalue.pos, score.Rpvalue.zero))
          #score.Rpvalue.combS = .simes.test(c(score.Rpvalue.pos, score.Rpvalue.zero))
          #comb.results = c(comb.results, scoreF.Rpvalue = score.Rpvalue.combF, scoreS.Rpvalue = score.Rpvalue.combS)
        }
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        
        ### 07/19/2018, modify while loop
        while(flag){
          
          results = .resample.work.pos.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, index.subj.pos, index.cova.ls, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, end.nperm, n.pos, pos.acc, test)
          
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
          if(start.nperm >= n.replicates){
            flag = 0
          }else if(end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings")) 
            results = .resample.work.pos.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, index.subj.pos, index.cova.ls, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, n.replicates, n.pos, pos.acc, test)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            
            flag = 0
            
          }
          
        }
        
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        if(test=="chisq"){
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        }
        if(test=="vc"){
          #comb.results = c(comb.results, scoreF.Rpvalue = score.Rpvalue.pos, scoreS.Rpvalue = score.Rpvalue.pos)
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        }        
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        
        ### 07/19/2018, modify while loop
        while(flag){
          
          
          results = .resample.work.zero.cluster(id, perm.type, ZZ, p.ZZ, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.zero, zero.acc, test)
          
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
          if(start.nperm >= n.replicates){
            flag = 0
          }else if(end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            
            results = .resample.work.zero.cluster(id, perm.type, ZZ, p.ZZ, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.zero, zero.acc, test)
            
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            
            flag = 0
          }
          
        }
        
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        if(test=="chisq"){
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        }
        if(test=="vc"){
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
          #comb.results = c(comb.results, scoreF.Rpvalue = score.Rpvalue.zero, scoreS.Rpvalue = score.Rpvalue.zero)
        }        
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        if(test=="chisq"){
          comb.results = c(comb.results, score.Rpvalue = NA)
        }
        if(test=="vc"){
          comb.results = c(comb.results, score.Rpvalue = NA)
          #comb.results = c(comb.results, scoreF.Rpvalue = NA, scoreS.Rpvalue = NA)
        }
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}


.Score.test.stat.pos.cluster <- function(id1, Y1, XX1, index.cova.ls, X1.par.index.ls){
  
  n = nrow(Y1)
  m = ncol(Y1) 
  #p = dim(XX)[2]
  p.XX = sapply(index.cova.ls, length)
  
  XX = array(NA,dim=c(n,dim(XX1)[2],m-1))
  for(j in 1:(m-1)){
    XX[,1:p.XX[j],j] = XX1[,index.cova.ls[[j]],j]
  }
  
  nY = rowSums(Y1)
  
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if( sum( sapply(X1.par.index.ls, function(x) sum(x==1)) )>0 ){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X1.par.index.ls) || n==0){
    
    score.stat.beta = NA
    
  }else{
    
    
    #XX.reduce = XX[,-X.par.index, , drop=FALSE]
    #p.reduce = p - length(X.par.index)    
    #par.interest.index.beta =  kronecker(((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
    
    p.XX.reduce = p.XX - sapply(X1.par.index.ls, length)
    XX.reduce = array(NA,dim=c(n,dim(XX)[2],m-1))
    for(j in 1:(m-1)){
      XX.reduce[,1:p.XX.reduce[j],j] = XX[,(1:p.XX[j])[-X1.par.index.ls[[j]]] ,j]
    }    
    
    par.interest.index.beta = NULL
    tmp = 0
    for(j in 1:(m-1)){
      par.interest.index.beta = c(par.interest.index.beta, tmp + X1.par.index.ls[[j]])
      tmp = sum(p.XX[1:j])
    }
    
    
    n.par.interest.beta = length(par.interest.index.beta) 
    
    #beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
    beta.ini.reduce = rep(0, sum(p.XX.reduce) )   
    
    data.reduce.beta = list(Y=Y1, XX=XX.reduce, p.XX=p.XX.reduce)
    
    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0 
    #est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos, gr=.fun.neg.score.beta.pos, data = data.reduce.beta)$par
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos.XX, gr=.fun.neg.score.beta.pos.XX, data = data.reduce.beta, method="BFGS")$par
    
    
    data.beta = list(Y=Y1, XX=XX, p.XX=p.XX)
    Score.reduce.beta = .fun.score.i.beta.pos.XX(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    #S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, idx.start])
    tmp = .fun.hessian.beta.pos.XX(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list  
    
    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
    
    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    
    
    B2 =  matrix(0, n.beta, n.beta)
    
    # 03/12/2018   modify for clustered data
    # 03/28/2018   modify to allow for uneven #obs across clusters
    id1.unique = unique(id1)
    n.cluster = length(id1.unique)
    #TT = n/length(unique(id1))  ## same #observation for each cluster for now
    #n.cluster = n/TT
    #idx1 = 1
    #idx2 = TT
    for(i in 1:n.cluster){
      
      #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
      Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id1==id1.unique[i]),,drop=FALSE] )
      B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
      
      #idx1 = idx2 + 1
      #idx2 = idx2 + TT
      
    }
    
    
    
    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A
    
    
  }
  
  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list, U=A, V=B )   )
  
  
}


.Score.test.stat.zero.cluster <- function(id.sam, Y0, ZZ, p.ZZ, Z.par.index, cor.stru){
  
  m = ncol(Y0)
  p.ZZ.reduce = p.ZZ - length(Z.par.index)
  ZZ.reduce = array(NA, dim=c(dim(ZZ)[1], dim(ZZ)[2],m))
  for(j in 1:m){
    ZZ.reduce[,1:p.ZZ.reduce[j],j] = ZZ[,(1:p.ZZ[j])[-Z.par.index],j]
  }
  
  n = nrow(Y0)
  
  #p = dim(ZZ)[2]
  #p.reduce = dim(ZZ.reduce)[2]
  outcome = NULL
  id = NULL
  cova.reduce = NULL
  
  for(i in 1:n){
    
    outcome = c(outcome, Y0[i,])
    
    idx1 = 0
    idx2 = 0
    
    for(j in 1:m){
      
      idx1 = idx2 + 1
      idx2 = idx2 + p.ZZ.reduce[j]
      
      tmp = rep(0, sum(p.ZZ.reduce))
      tmp[idx1:idx2] = ZZ.reduce[i,1:p.ZZ.reduce[j],j]
      cova.reduce = rbind(cova.reduce, tmp )
      
    }
    
    
    id = c(id, rep(i, m))
  }
  
  #data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  #gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence")
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence")
  #wald.test = anova(gee.full, gee.reduce)
  
  
  ########### perform score test
  #n.alpha = m *p
  n.alpha = sum(p.ZZ)
  
  #par.interest.index.alpha =  kronecker( ((0:(m-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  par.interest.index.alpha = NULL
  tmp = 0
  for(j in 1:m){
    par.interest.index.alpha = c(par.interest.index.alpha, tmp + Z.par.index)
    tmp = sum(p.ZZ[1:j])
  }
  
  n.par.interest.alpha = length(par.interest.index.alpha) 
  est.reduce.alpha = rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] = 0 
  est.reduce.alpha[-par.interest.index.alpha] = coef(gee.reduce)
  est.reduce.scale = gee.reduce
  
  data.alpha = list(Y=Y0, ZZ=ZZ, p.ZZ=p.ZZ)
  
  tmp = .fun.score.i.alpha.ZZ(est.reduce.alpha, data.alpha, save.list=TRUE)
  Score.reduce.alpha = tmp$Score.alpha
  # for resampling test
  vA.list = tmp$vA.list
  Vinv.list = tmp$Vinv.list
  VY.list = tmp$VY.list
  
  Hess.reduce.alpha =  .fun.hessian.alpha.ZZ(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  # 03/12/2018   modify for clustered data
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id.unique = unique(id.sam)
  n.cluster = length(id.unique)
  #TT = n/length(unique(id.sam))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id.sam==id.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.alpha = A %*% ginv(B) %*% A
  score.pvalue.alpha = 1 - pchisq(score.stat.alpha, n.par.interest.alpha) 
  
  
  
  return(list(score.df.alpha=n.par.interest.alpha, score.stat.alpha = score.stat.alpha, score.pvalue.alpha=score.pvalue.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list, U=A, V=B )   )
  
}



# add 07/02/2016 for adaptive resampling
.resample.work.two.cluster <- function(id, perm.type, XX, X.par.index, X1.par.index.ls, ZZ, p.ZZ, Z.par.index, index.subj.pos, index.cova.ls, score.stat.pos, score.stat.zero, pos.S.beta.list, pos.I.beta.list, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc, test){
  
  n = dim(XX)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  XX.interest = XX[,X.par.index,, drop=FALSE]
  ZZ.interest = ZZ[,Z.par.index,, drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    for(i in 1:n.cluster){
      tmpX = XX.interest[which(id==id.unique[i]),,,drop=FALSE]
      tmpZ = ZZ.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmpX)[1])
      if(dim(unique(tmpX))[1]!=1 | dim(unique(tmpZ))[1]!=1){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        XX.interest.unique = unique(tmpX)
        ZZ.interest.unique = unique(tmpZ)
      }else{
        XX.interest.unique = abind(XX.interest.unique, unique(tmpX), along=1)
        ZZ.interest.unique = abind(ZZ.interest.unique, unique(tmpZ), along=1)
      }
      
      
    }
    
    
    
  }
  
  
  
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  n.comb.new = n.comb
  comb.acc.new = comb.acc
  
  score.stat.comb = score.stat.pos + score.stat.zero
  
  
  for(k in start.nperm:end.nperm){
    
    # 03/28/2018   modify to allow for uneven #obs across clusters
    if(perm.type=="BTW"){
      
      
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        XX.interest.perm.i = array(NA, dim=c(nobs, dim(XX.interest.unique)[2], dim(XX.interest.unique)[3]))
        ZZ.interest.perm.i = array(NA, dim=c(nobs, dim(ZZ.interest.unique)[2], dim(ZZ.interest.unique)[3]))
        for(j in 1:dim(XX.interest.unique)[3]){
          XX.interest.perm.i[,,j] =  matrix(rep(XX.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        for(j in 1:dim(ZZ.interest.unique)[3]){
          ZZ.interest.perm.i[,,j] =  matrix(rep(ZZ.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        if(i==1){
          XX.interest.perm = XX.interest.perm.i
          ZZ.interest.perm = ZZ.interest.perm.i
        }else{
          XX.interest.perm = abind(XX.interest.perm, XX.interest.perm.i, along =1)
          ZZ.interest.perm = abind(ZZ.interest.perm, ZZ.interest.perm.i, along =1)
        }
        
      }
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.interest.perm
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.interest.perm      
      
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
    }
    
    #X.perm = X
    #X.perm[,X.par.index] = X.perm[sample.idx,X.par.index]
    
    XX1.perm = XX.perm[index.subj.pos, , ,drop=FALSE]
    #XX1.perm = XX1.perm[,c(1, index.cova), ,drop=FALSE] 
    p.XX = sapply(index.cova.ls, length)
    XX.tmp = array(NA,dim=c(dim(XX1.perm)[1],dim(XX1.perm)[2], dim(XX1.perm)[3] ))
    for(j in 1:length(p.XX)){
      XX.tmp[,1:p.XX[j],j] = XX1.perm[,index.cova.ls[[j]],j]
    }
    XX1.perm = XX.tmp
    
    #Z.perm = Z
    #Z.perm[,Z.par.index] = Z.perm[sample.idx,Z.par.index]
    
    score.stat.pos.perm = try( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list, test) )
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list, test) )
    #score.stat.pos.perm = ( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list) )
    #score.stat.zero.perm = ( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    if(class(score.stat.pos.perm)[1] != "try-error"){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    if(class(score.stat.zero.perm)[1] != "try-error"){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    if(test=="chisq"){
      
      if(class(score.stat.pos.perm)[1] != "try-error" & class(score.stat.zero.perm)[1] != "try-error"){
        
        score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
        n.comb.new = n.comb.new + 1
        if(score.stat.comb.perm >= score.stat.comb){
          comb.acc.new = comb.acc.new + 1
          
        }  
      }
    }
    
  }
  
  # modify 07/19/2018
  if(pos.acc.new < 10 | zero.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }
  # else if(pos.acc.new<10 | zero.acc.new<10){
  #   next.end.nperm = ( end.nperm + 1) * 10 - 1;
  #   flag = 1;
  #   
  # }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              n.comb.new=n.comb.new, comb.acc.new=comb.acc.new,
              flag=flag, next.end.nperm=next.end.nperm))
  
}


# add 07/02/2016 for adaptive resampling
.resample.work.pos.cluster <- function(id, perm.type, XX, X.par.index, X1.par.index.ls, index.subj.pos, index.cova.ls, score.stat.pos, pos.S.beta.list, pos.I.beta.list, start.nperm, end.nperm, n.pos, pos.acc, test){
  
  n = dim(XX)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  XX.interest = XX[,X.par.index,, drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    for(i in 1:n.cluster){
      tmpX = XX.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmpX)[1])
      if(dim(unique(tmpX))[1]!=1 ){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        XX.interest.unique = unique(tmpX)
      }else{
        XX.interest.unique = abind(XX.interest.unique, unique(tmpX), along=1)
        
      }
      
      
    }
    
    
    
  }
  
  
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  for(k in start.nperm:end.nperm){
    
    if(perm.type=="BTW"){
      
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        XX.interest.perm.i = array(NA, dim=c(nobs, dim(XX.interest.unique)[2], dim(XX.interest.unique)[3]))
        
        for(j in 1:dim(XX.interest.unique)[3]){
          XX.interest.perm.i[,,j] =  matrix(rep(XX.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
          
        }
        if(i==1){
          XX.interest.perm = XX.interest.perm.i
          
        }else{
          XX.interest.perm = abind(XX.interest.perm, XX.interest.perm.i, along =1)
          
        }
        
      }
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.interest.perm
      
      
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
    }
    
    
    XX1.perm = XX.perm[index.subj.pos, , ,drop=FALSE]
    #XX1.perm = XX1.perm[,c(1, index.cova), ,drop=FALSE] 
    p.XX = sapply(index.cova.ls, length)
    XX.tmp = array(NA,dim=c(dim(XX1.perm)[1],dim(XX1.perm)[2],dim(XX1.perm)[3]))
    for(j in 1:length(p.XX)){
      XX.tmp[,1:p.XX[j],j] = XX1.perm[,index.cova.ls[[j]],j]
    }
    XX1.perm = XX.tmp    
    
    
    score.stat.pos.perm = try( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list, test) )
    #score.stat.pos.perm = ( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list) )
    
    
    if(class(score.stat.pos.perm)[1] != "try-error"){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    
    
  }
  
  ### modify 07/19/2018
  if(pos.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }
  # else if(pos.acc.new<10){
  #   next.end.nperm = ( end.nperm + 1) * 10 - 1;
  #   flag = 1;
  #   
  # }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}


.Score.test.stat.pos.4Gresampling.cluster <- function(id1, XX1.perm, index.cova.ls, X1.par.index.ls, S.beta.list, I.beta.list, test){
  
  n = dim(XX1.perm)[1]
  #p = dim(XX.perm)[2]
  p.XX = sapply(index.cova.ls, length)
  
  m.beta = length(S.beta.list[[1]])
  
  XX.perm = array(NA,dim=c(n,dim(XX1.perm)[2],m.beta))
  for(j in 1:m.beta){
    XX.perm[,1:p.XX[j],j] = XX1.perm[,index.cova.ls[[j]],j]
  }
  
  
  
  #idx.end = (1:m.beta)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.XX[-m.beta]))+1
  idx.end = idx.start + p.XX - 1
  
  #n.par.interest.beta = m.beta
  #n.beta = m.beta*p
  n.beta = sum(p.XX)
  #par.interest.index.beta =  kronecker( ((0:(m.beta-1))*p), rep(1,length(X.par.index))) + X.par.index
  par.interest.index.beta = NULL
  tmp = 0
  for(j in 1:m.beta){
    par.interest.index.beta = c(par.interest.index.beta, tmp + X1.par.index.ls[[j]])
    tmp = sum(p.XX[1:j])
  }
  
  n.par.interest.beta = length(par.interest.index.beta) 
  
  
  Score.reduce.beta.perm = matrix(0, n, n.beta )
  Hess.reduce.beta.perm = matrix(0, n.beta, n.beta )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    #Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(X.perm[i,], ncol=1))  
    tmp = NULL
    for(j in 1:m.beta){
      tmp = c(tmp, S.beta.list[[i]][j] * XX.perm[i,1:p.XX[j],j] )
    }
    Score.reduce.beta.perm[i,] = tmp
    
    #Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  X.perm[i,] %o% X.perm[i,] ) )
    
    Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
    for(j in 1:m.beta){
      for(k in j:m.beta){
        Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = I.beta.list[[i]][j,k] * (XX.perm[i,1:p.XX[j],j] %o% XX.perm[i,1:p.XX[k],k])
        Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
      }
    }
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + Hessian.tmp
    
    
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id1.unique = unique(id1)
  n.cluster = length(id1.unique)
  #TT = n/length(unique(id))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[which(id1==id1.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.beta.perm.reorg.cluster %o%  Score.reduce.beta.perm.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  
  B = B1 %*% B2 %*% t(B1)
  
  if(test=="chisq"){
    score.stat.beta.perm = A %*% ginv(B) %*% A
  }
  
  if(test=="vc"){
    
    inv.B = ginv(B)
    score.stat.beta.perm = A %*% inv.B %*% inv.B %*% A
  }
  
  return(score.stat.beta.perm)
  
  
  
}


# add 07/02/2016 for adaptive resampling
.resample.work.zero.cluster <- function(id, perm.type, ZZ, p.ZZ, Z.par.index, score.stat.zero, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.zero, zero.acc, test){
  
  n = dim(ZZ)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  ZZ.interest = ZZ[,Z.par.index,, drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    for(i in 1:n.cluster){
      
      tmpZ = ZZ.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmpZ)[1])
      if(dim(unique(tmpZ))[1]!=1){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        
        ZZ.interest.unique = unique(tmpZ)
      }else{
        
        ZZ.interest.unique = abind(ZZ.interest.unique, unique(tmpZ), along=1)
      }
      
      
    }
    
    
  }
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  for(k in start.nperm:end.nperm){
    
    if(perm.type=="BTW"){
      
      
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        ZZ.interest.perm.i = array(NA, dim=c(nobs, dim(ZZ.interest.unique)[2], dim(ZZ.interest.unique)[3]))
        for(j in 1:dim(ZZ.interest.unique)[3]){
          ZZ.interest.perm.i[,,j] =  matrix(rep(ZZ.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        if(i==1){
          ZZ.interest.perm = ZZ.interest.perm.i
        }else{
          ZZ.interest.perm = abind(ZZ.interest.perm, ZZ.interest.perm.i, along =1)
        }
        
      }
      
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.interest.perm   
      
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
    }
    
    
    
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list, test) )
    #score.stat.zero.perm = ( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    
    if(class(score.stat.zero.perm)[1] != "try-error"){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    
  }
  
  ### modify 07/19/2018
  if(zero.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }
  # else if(zero.acc.new<10){
  #   next.end.nperm = ( end.nperm + 1) * 10 - 1;
  #   flag = 1;
  #   
  # }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}


.Score.test.stat.zero.4Gresampling.cluster <- function(id, ZZ.perm, p.ZZ, Z.par.index, vA.list, Vinv.list, VY.list, test){
  
  n = dim(ZZ.perm)[1]
  #p = dim(ZZ.perm)[2]
  
  
  m.alpha = length(vA.list[[1]])
  #n.alpha = m.alpha*p
  n.alpha = sum(p.ZZ)
  
  #par.interest.index.alpha =  kronecker( ((0:(m.alpha-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  par.interest.index.alpha = NULL
  tmp = 0
  for(j in 1:m.alpha){
    par.interest.index.alpha = c(par.interest.index.alpha, tmp + Z.par.index)
    tmp = sum(p.ZZ[1:j])
  }
  
  
  n.par.interest.alpha = length(par.interest.index.alpha) 
  
  Score.reduce.alpha.perm = matrix(0, n, n.alpha )
  Hess.reduce.alpha.perm = matrix(0, n.alpha, n.alpha )
  
  #idx.end = (1:m.alpha)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.ZZ[-m.alpha]))+1
  idx.end = idx.start + p.ZZ - 1
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         alpha part: resampling Score test        #
    #                                                 #
    ################################################### 
    #tD.tmp = kronecker(.diag2(vA.list[[i]]), as.matrix(Z.perm[i,], ncol=1))
    tD.tmp = matrix(0, nrow=n.alpha, ncol=m.alpha)
    for(j in 1:m.alpha){
      tD.tmp[idx.start[j]:idx.end[j], j] = vA.list[[i]][j] * ZZ.perm[i,1:p.ZZ[j],j] 
    }
    
    Score.reduce.alpha.perm[i,] = Score.reduce.alpha.perm[i,] + tD.tmp %*% VY.list[[i]]
    
    Hess.reduce.alpha.perm = Hess.reduce.alpha.perm + tD.tmp %*% Vinv.list[[i]] %*% t(tD.tmp)
    
    
  }
  
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha.perm[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  # 03/12/2018   modify for clustered data
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id.unique = unique(id)
  n.cluster = length(id.unique)
  #TT = n/length(unique(id))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id==id.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  B = B1 %*% B2 %*% t(B1)
  
  if(test=="chisq"){
    score.stat.alpha.perm = A %*% ginv(B) %*% A
  }
  
  if(test=="vc"){
    
    ##################### 07/19/2018. for zero-part, still use chisq test
    #inv.B = ginv(B)
    #score.stat.alpha.perm = A %*% inv.B %*% inv.B %*% A
    
    score.stat.alpha.perm = A %*% ginv(B) %*% A
  }
  
  
  return(score.stat.alpha.perm)
  
}

