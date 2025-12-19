setwd("~/Downloads/AIGDM/LONGITUDINAL/code/Analysis/test/")
rm(list = ls())
library(phyloseq)
library(data.table)
library(matrixStats)
library(Matrix)
source("../AIGDM_GEE.R")

seed <- 12345
n.boot <- 2e3
n.perm <- 1e4
fdr.alpha <- 0.1
load("../../Data/Deriveddata/ECAM.Genus.filter10.rda")
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

ID=id; OTU=Y.m; X4mean=X4disp=X4zero=X.m;X.index=1;ZI.pos="adaptive";Tax=tax;min.depth=0;n.cores=1

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

k = 1
Rank.low = paste("Rank", n.rank-k,sep="")
Rank.high = paste("Rank", n.rank-k+1,sep="")

tmp = table(tax[,n.rank-k])
level.uni = sort( names(tmp)[which(tmp>1)] )
m.level = length(level.uni)

tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
setnames(tt, 1:2, c(Rank.low, Rank.high))
W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
W.count = data.matrix(tt[, otucols, with=FALSE])

zi.check.pval.lst = list()
curr.ind = 1

m = 1
Y.tmp = t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
colOrder = order(colMeans(Y.tmp/rowSums(Y.tmp), na.rm = T), decreasing = TRUE)

Y.tmp = Y.tmp[,colOrder]
K = ncol(Y.tmp) - 1

y.tmp = cbind(Y.tmp[,1], rowSums(Y.tmp[,(1+1):(K+1),drop=FALSE]))
id.tmp = ID[rowSums(y.tmp) != 0]
y.tmp = y.tmp[rowSums(y.tmp) != 0, , drop = FALSE]

save(y.tmp, id.tmp, file = "tmp.rdata")


Sys.setenv("PKG_CXXFLAGS"="-I/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/RcppArmadillo/include")

library(Rcpp)
sourceCpp("AIGDM_Rcpp.cpp")


fast_check_zeroinflation(y.tmp, id.tmp, nboot = 500, ncores = n.cores, Rsel = R.sel, seed = seed)

benchmark_result <- microbenchmark(
  old_check = .check_zeroinflation(ID, Y, n.boot = 500, n.cores = 1, R.sel = R.sel, seed = seed),
  new_check = fast_check_zeroinflation(Y, ID, nboot = 500, ncores = 1, Rsel = R.sel, seed = seed),
  times = 5
)
