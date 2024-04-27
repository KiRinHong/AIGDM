.fitted_model <- function(W, Xa, Xb, reg, mod, SeqDepth){
  # reg=aigdm.reg; W=W.r.tmp; Xa=X.tmp; Xb=X.tmp; mod="ZIGDM"; SeqDepth=rowSums(Y.sub)
  alpha = reg$alpha.est
  beta = reg$beta.est
  N = nrow(Xa)
  
  mu.mat = 1/(1+exp(-(Xa %*% alpha)))
  sigma.mat = 1/(1+exp(-(Xb %*% beta)))
  a.mat = mu.mat * (1/sigma.mat - 1)
  b.mat = (1 - mu.mat) * (1/sigma.mat - 1)
  
  Z.mat = NULL
  if(any(mod == "NIGDM", is.null(mod))){
    Z.mat = rbeta(N, a.mat, b.mat)
  }else if(mod == "ZIGDM"){
    gamma = reg$gamma.est
    pstr0.mat = 1/(1+exp(-(W %*% gamma)))
    tmp_Z = rbeta(N, a.mat, b.mat)
    tmp_Delta = rbinom(N, 1, pstr0.mat)
    Z.mat = (1-tmp_Delta) * tmp_Z
  }
  
  P = matrix(NA, nrow = N, ncol = 2)
  P[,1] = Z.mat
  P[,2] = 1 - P[,1]
  P[P<0] = 0
  Y = NULL
  for (i in 1:N) {
    Y = rbind(Y, t(rmultinom(1, SeqDepth[i], P[i,])))
  }
  tmp = (Y/rowSums(Y))[,1]
  return(tmp)
}

.corrData.ZIBeta <- function(n.T, r, beta.a, beta.b, pstr){
  # beta.a = a.mat; beta.b = b.mat; pstr = pstr.mat
  N = nrow(beta.a)/n.T
  K = ncol(beta.a)
  beta.a1 = beta.a[1:N,]
  beta.a2 = beta.a[(N+1):(N*n.T),]
  beta.b1 = beta.b[1:N,]
  beta.b2 = beta.b[(N+1):(N*n.T),]
  pstr.1 = pstr[1:N,]
  pstr.2 = pstr[(N+1):(N*n.T),]
  betavars = c()
  
  for (j in 1:K) {
    tmp = tmp_Z = tmp_Delta = matrix(NA, nrow = N, ncol = n.T)
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    beta.a.j = cbind(beta.a1[,j], beta.a2[,j])
    beta.b.j = cbind(beta.b1[,j], beta.b2[,j])
    for (t in 1:n.T) {
      tmp_Z[,t] = qbeta(pnorm(rawvars[,t]), beta.a.j[,t], beta.b.j[,t]) 
    }
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    pstr.j = cbind(pstr.1[,j], pstr.2[,j])
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
  beta.a2 = beta.a[(N+1):(N*n.T),]
  beta.b1 = beta.b[1:N,]
  beta.b2 = beta.b[(N+1):(N*n.T),]
  pstr.1 = pstr[1:N,]
  pstr.2 = pstr[(N+1):(N*n.T),]
  betavars = c()
  
  for (j in 1:K) {
    tmp = tmp_Z = tmp_Delta = matrix(NA, nrow = N, ncol = n.T)
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    beta.a.j = cbind(beta.a1[,j], beta.a2[,j])
    beta.b.j = cbind(beta.b1[,j], beta.b2[,j])
    for (t in 1:n.T) {
      tmp_Z[,t] = qbeta(pnorm(rawvars[,t]), beta.a.j[,t], beta.b.j[,t]) 
    }
    
    rawvars = mvrnorm(n = N, mu = rep(0, n.T), 
                      Sigma = matrix(r, nrow = n.T, ncol = n.T) + diag(n.T)*(1-r))
    pstr.j = cbind(pstr.1[,j], pstr.2[,j])
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

.fitted_model_DM <- function(mod, SeqDepth, taxa.id){
  P = mod@fitted
  P[P<0] = 0
  Y = NULL
  for (i in 1:nrow(P)) {
    Y = rbind(Y, t(rmultinom(1, SeqDepth[i], P[i,])))
  }
  tmp = (Y[,taxa.id,drop=FALSE]/rowSums(Y[,taxa.id:ncol(Y)]))[,1]
  tmp[is.nan(tmp)] = 0
  return(tmp)
}

.fitted_model_DM_sub <- function(mod, SeqDepth){
  P = mod@fitted
  Y = NULL
  for (i in 1:nrow(P)) {
    Y = rbind(Y, t(rmultinom(1, SeqDepth[i], round(P[i,],4))))
  }
  tmp = (Y/rowSums(Y))[,1]
  tmp[is.nan(tmp)] = 0
  return(tmp)
}


plot_subcomp <- function(taxa_int) {
  #taxa_int = conf_taxa[1]
  print(taxa_int)
  rank <- as.numeric(sub("Rank", "", unlist(strsplit(taxa_int, "\\."))[seq(1, 2*length(taxa_int), 2)]))
  rank_code <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")[rank]
  taxa_tax <- c("k", "p", "c", "o", "f", "g")[rank]
  taxa <- unlist(strsplit(taxa_int, "\\."))[seq(2, 2*length(taxa_int), 2)]
  
  n.rank <- 6
  k <- n.rank-rank
  Rank.low <- paste0("Rank",n.rank-k)  # change the low rank
  Rank.high <- paste0("Rank",n.rank-k+1) # change the high rank
  tt <- W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
  setnames(tt, 1:2, c(Rank.low, Rank.high))
  
  W.tax <- as.vector(unlist(tt[, Rank.low, with=FALSE]))
  W.count <- data.matrix(tt[, otucols, with=FALSE])
  Y.tmp <- t(W.count[which(W.tax == taxa), , drop=FALSE])
  id.keep <- which(rowSums(Y.tmp) != 0)
  Y.tmp <- Y.tmp[id.keep,]
  colOrder <- order(colMeans(Y.tmp/rowSums(Y.tmp)), decreasing = TRUE)
  Y.tmp <- Y.tmp[,colOrder]
  
  meta.tmp <- data.frame(expo = unlist(meta[id.keep, "delivery"]), 
                         month = factor(unlist(meta[id.keep, "life_month"]), 
                                        levels = c(0:22, 24), labels = 1:24),
                         Y.tmp/rowSums(Y.tmp))
  colnames(meta.tmp) <- c("exposure", "life_month", unlist(tt[W.tax == taxa, ..Rank.high])[colOrder]) # Rank.high
  
  # meta.count <- meta.tmp[,c(1,2)] %>%
  #   group_by(exposure, life_month) %>%
  #   dplyr::count() %>%
  #   ungroup()
  # meta.count.w <- tidyr::spread(as.data.frame(meta.count), exposure, n)
  # names(meta.count.w)[2:3] <- paste0(names(meta.count.w)[2:3], " (N) ")
  # meta.count.w <- as.matrix(meta.count.w)
  
  meta.group <- meta.tmp %>%
    group_by(exposure, life_month) %>%
    summarise_each(mean) %>%
    ungroup()
  
  meta.group.l <- tidyr::gather(meta.group, taxa_child, count, unlist(tt[W.tax == taxa, ..Rank.high])[colOrder], factor_key=FALSE)
  stackOrder <- unique(meta.group.l$taxa_child)
  colorDef <- c("#377eb8", "#e41a1c", "#4daf4a",  "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf")
  id.unc <- grep("Unclassified", stackOrder)
  newColor <- numeric(length(stackOrder))
  if (length(id.unc) > 0) {
    newColor[id.unc] <- "#999999"
    newColor[-id.unc] <- colorDef[1:(length(stackOrder)-1)]
  }else{
    newColor <- colorDef[1:length(stackOrder)]
  }
  
  if(taxa_int == "Rank5.Ruminococcaceae"){
    taxa_child_label = c("g_Faecalibacterium", "g_Oscillospira", "g_Ruminococcus", "Unclassified", "g_Gemmiger", "g_Anaerotruncus")
  }else if(taxa_int == "Rank5.Veillonellaceae"){
    taxa_child_label = c("g_Veillonella", "g_Dialister", "g_Phascolarctobacterium", "g_Megasphaera")
  }else if(taxa_int == "Rank5.Clostridiaceae"){
    taxa_child_label = c("g_Clostridium", "Unclassified", "g_SMB53")
  }else if(taxa_int == "Rank4.Pseudomonadales"){
    taxa_child_label = c("f_Moraxellaceae", "f_Pseudomonadaceae")
  }else if(taxa_int == "Rank3.Actinobacteria"){
    taxa_child_label = c("o_Bifidobacteriales", "o_Actinomycetales")
  }else if(taxa_int == "Rank2.Proteobacteria"){
    taxa_child_label = c("c_Gammaproteobacteria", "c_Betaproteobacteria", "c_Deltaproteobacteria", "c_Epsilonproteobacteria")
  }else if(taxa_int == "Rank1.Bacteria"){
    taxa_child_label = c("p_Firmicutes", "p_Bacteroidetes", "p_Proteobacteria", "p_Actinobacteria", 
                         "p_Verrucomicrobia", "p_Fusobacteria")
  }
  meta.group.l$taxa_child <- factor(meta.group.l$taxa_child, levels=stackOrder, labels = taxa_child_label)
  meta.group.l$exposure[meta.group.l$exposure == 0] <- "C-section"
  meta.group.l$exposure[meta.group.l$exposure == 1] <- "Vaginal birth"
  
  g1 <- ggplot(meta.group.l, 
               aes(fill=taxa_child, y=count, x=life_month)) +
    geom_bar(position="fill", color="black", stat="identity", width=1, linewidth=0.3) +
    scale_fill_manual(values = newColor) + 
    scale_x_discrete() +
    facet_wrap( ~exposure, scale="free_x") +
    labs(title = "", x = "Age (months)", y = "Relative abundance") + # input[5] delivery mode
    theme_bw(base_size=20) +
    theme(axis.title = element_text(face = "bold",size = rel(1.5)),
          axis.title.y = element_text(angle = 90, vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_text(size = rel(1.5)),
          axis.text.x = element_text(size = rel(1.5), angle = 90),
          # axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(),
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          strip.text.x = element_text(size = rel(1.8)),
          plot.title = element_text(face='bold', size = rel(2)),
          plot.margin = unit(c(5, 10, 5, 10), "mm"),
          legend.position = "bottom") + 
    guides(fill = guide_legend(title = paste0(taxa_tax, "_", taxa),
                               title.theme = element_text(face = "bold.italic", size = 25, angle = 0),
                               label.theme = element_text(face = "italic", size = 25, angle = 0),
                               byrow = TRUE, ncol = 4))
  
  zi.check <- rslt$AIGDMC$lineage.zi[[taxa_int]]
  K <- ncol(Y.tmp) - 1
  oi.check <- rep(FALSE, K)
  
  ID.tmp <- id[id.keep]
  W.tmp <- W[id.keep,,drop=FALSE]; W.r.tmp <- W.r[id.keep,,drop=FALSE]
  X.tmp <- X[id.keep,,drop=FALSE]; X.r.tmp <- X.r[id.keep,,drop=FALSE]
  
  # est.para.lst = lapply(1:K, function(x){
  #   .est_para(ID.tmp, W.r.tmp, X.r.tmp, X.r.tmp, cbind(Y.tmp[,x], rowSums(Y.tmp[,(x+1):(K+1),drop=FALSE])),
  #             zi.check[x], oi.check[x])
  # })
  # 
  # asym.mean = .score_test_mean(ID.tmp, X.tmp, X.index, Y.tmp, est.para.lst, zi.check, oi.check)
  # stat.mean.taxonwise = round(asym.mean$stat, 3)
  # 
  # asym.disp = .score_test_disp(ID.tmp, X.tmp, X.index, Y.tmp, est.para.lst, zi.check, oi.check)
  # stat.disp.taxonwise = round(asym.disp$stat, 3)

  meta.group.l$exposure[meta.group.l$exposure == "C-section"] <- "CS"
  meta.group.l$exposure[meta.group.l$exposure == "Vaginal birth"] <- "VB"
  
  # annotate_text_mean <- c(paste0("T[mu] ==", stat.mean.taxonwise), "")
  # annotate_text_disp <- c(paste0("T[phi] ==", stat.disp.taxonwise), "")
  # label_data <- data.frame(
  #   taxa_child = unique(meta.group.l$taxa_child),
  #   label_mean = annotate_text_mean,
  #   label_disp = annotate_text_disp
  # )
  # max_counts <- aggregate(count ~ taxa_child, data = meta.group.l, FUN = max)
  # min_counts <- aggregate(count ~ taxa_child, data = meta.group.l, FUN = min)
  # label_data <- merge(label_data, max_counts, by = "taxa_child")
  # label_data <- merge(label_data, min_counts, by = "taxa_child")
  
  if(taxa_int %in% c("Rank5.Ruminococcaceae", "Rank1.Bacteria")){
    size <- 1.5
  }else{
    size <- 1.8
  }
  g1.1 <- ggplot(meta.group.l, aes(y = count, x=factor(exposure, levels = c("CS", "VB")), 
                                   fill = taxa_child)) + 
    geom_violin(position="dodge", color = "black", trim = TRUE) +
    stat_summary(fun =function(x) c(quantile(x,0.25), quantile(x,0.75)), geom="line", color="white") +
    stat_summary(fun = "mean", geom = "point", shape = 2, size = 3, color="black") +
    scale_fill_manual(values = newColor) +
    facet_wrap(~ taxa_child, scales = "free", nrow = 1) +
    # geom_text(data = label_data, aes(x = 1.5, y = count.x - 0.1 * (count.x-count.y), label = label_mean), size = 8, color = "black", parse = T) +
    # geom_text(data = label_data, aes(x = 1.5, y = count.x - 0.2 * (count.x-count.y), label = label_disp), size = 8, color = "black", parse = T) +
    # coord_flip() +
    labs(title = "", x = "", y = "Taxon proportion") +
    theme_bw(base_size=20) +
    theme(axis.title = element_text(face = "bold",size = rel(1.5)),
          axis.title.y = element_text(angle = 90, vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_text(size = rel(1.5)),
          axis.text.x = element_text(size = rel(1.5)),
          axis.ticks.x = element_blank(),
          axis.text = element_text(),
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          strip.text = element_text(face = "italic"),
          strip.text.x = element_text(size = rel(size)),
          legend.position = "none")
  g1.1
  if(taxa_int == "Rank5.Clostridiaceae"){
    j = 1
  }else if(taxa_int == "Rank5.Ruminococcaceae"){
    j = c(1, 3)
  }else if(taxa_int == "Rank3.Actinobacteria"){
    j = 1
  }else if(taxa_int == "Rank2.Proteobacteria"){
    j = 2
  }else{
    j = NULL
  }
  
  if(!is.null(j)){
    if(length(j) == 1){
      Y.sub <- cbind(Y.tmp[,j], rowSums(Y.tmp[,(j+1):(K+1),drop=FALSE]))
      id.sel <- which(rowSums(Y.sub)!=0)
      
      dw <- ncol(W.tmp); dw.r <- ncol(W.r.tmp); dx <- ncol(X.tmp); dx.r <- ncol(X.r.tmp)
      if(zi.check[j]){
        gamma0 = matrix(0.001, nrow = dw.r, ncol = 1)
      }else{
        gamma0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      }
      omega0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      alpha0 = matrix(0.001, nrow = dx, ncol = 1)
      beta0 = matrix(0.001, nrow = dx, ncol = 1)
      aigdm.reg <- .AIGDM_EM(Y.sub, W.r.tmp, X.tmp, X.tmp, gamma0, omega0, alpha0, beta0)
      gamma0 = omega0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      gdm.reg <- .AIGDM_EM(Y.sub, W.r.tmp, X.tmp, X.tmp, gamma0, omega0, alpha0, beta0)
      dm.reg <- MGLM::MGLMreg(Y.tmp ~ 0+X.tmp, dist = "DM")
      
      obs <- (Y.sub / rowSums(Y.sub))[id.sel, 1]
      obs.qq <- obs[order(obs)]
      n <- length(obs)
      
      n.iter <- 5000
      sample.dm.qq <- sample.gdm.qq <- sample.ai.qq <- numeric(n.iter * n)
      zeromean.dm <- zeromean.gdm <- zeromean.ai <- numeric(n.iter)
      
      inc <- 1:n
      for(iter in 1:n.iter){
        tmp <- .fitted_model(W.r.tmp, X.tmp, X.tmp, aigdm.reg, "ZIGDM", rowSums(Y.sub))[id.sel]
        sample.ai.qq[inc] <- tmp[order(tmp)]
        zeromean.ai[iter] <- mean(tmp == 0)
        
        tmp <- .fitted_model(W.r.tmp, X.tmp, X.tmp, gdm.reg, NULL, rowSums(Y.sub))[id.sel]
        sample.gdm.qq[inc] <- tmp[order(tmp)]
        zeromean.gdm[iter] <- mean(tmp == 0)
        
        # tmp <- .fitted_model_DM(dm.reg, rowSums(Y.tmp[id.sel,]), j)
        tmp <- .fitted_model_DM(dm.reg, rowSums(Y.tmp), j)[id.sel]
        sample.dm.qq[inc] <- tmp[order(tmp)]
        zeromean.dm[iter] <- mean(tmp == 0)
        inc <- inc + n
      }
      
      qq_data <- data.table(obs.qq = obs.qq,
                            sample.qq = c(sample.dm.qq, sample.gdm.qq, sample.ai.qq),
                            model = rep(c("DM","GDM", "AIGDM"), each = n.iter*n))
      
      prop_data <- data.frame(
        x = factor(c("Observed", "DM", "GDM", "AIGDM"), levels = c("Observed", "DM", "GDM", "AIGDM")),
        y = c(mean(obs == 0), mean(sample.dm.qq == 0), mean(sample.gdm.qq == 0), mean(sample.ai.qq == 0)),
        low = c(NA, quantile(zeromean.dm, 0.025), quantile(zeromean.gdm, 0.025), quantile(zeromean.ai, 0.025)),
        high = c(NA, quantile(zeromean.dm, 0.975), quantile(zeromean.gdm, 0.975), quantile(zeromean.ai, 0.975))
      )
      
      g2 <- ggplot(prop_data, aes(y = y, x = x)) +
        geom_bar(stat = "identity", fill = newColor[j], color = "black", width = 0.5, position = position_dodge(0.2)) +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.3, colour = "black", alpha = 0.9) +
        labs(y = expression(bold("Taxon zero proportion"))) +
        ylim(0, max(prop_data$high) + 0.05) +
        theme_bw(base_size=20) +
        theme(axis.title = element_text(face = "bold", size = rel(1.5)),
              axis.title.y = element_text(angle = 90, vjust =2),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              axis.ticks.x = element_blank(),
              axis.text = element_text(),
              axis.line = element_line(colour="black"),
              axis.ticks = element_line(),
              strip.text.x = element_text(size = rel(1.8)),
              plot.margin = unit(c(5, 10, 5, 10), "mm"))
      
      qq_ci <- qq_data[, .(median = median(sample.qq),
                           lower = quantile(sample.qq, 0.025),
                           upper = quantile(sample.qq, 0.975)),
                       by = .(obs.qq, model)]
      qq_ci$model <- factor(qq_ci$model, levels = c("DM", "GDM", "AIGDM"))
      
      g3 <- ggplot(qq_ci, aes(x = obs.qq)) +
        facet_wrap(~ model, nrow = 1) +
        geom_abline(linetype = "dashed", color = "grey") +
        geom_point(data = qq_ci, aes(y = median), color = newColor[j], size = 2) +
        geom_ribbon(data = qq_ci,
                    aes(ymin = lower, ymax = upper, y = median, fill = "band"),
                    fill = "grey", alpha = 0.4) +
        xlab("Observed quantiles") + ylab("Theoretical quantiles") +
        scale_x_continuous(breaks = seq(0, 1, 0.2)) +
        theme_bw() +
        theme(axis.title = element_text(face = "bold", size = rel(2.6)),
              axis.title.y = element_text(angle = 90, vjust = 2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(size = rel(1.5)),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              axis.line = element_line(colour = "black"),
              strip.text.x = element_text(size = rel(2.6)),
              plot.margin = unit(c(5, 10, 5, 10), "mm"))
      
      layout <- "
AAA
BBB
CDD"
      
      pdf(paste0("../Figs/ECAM.Full.WO.", taxa_tax, "_", taxa, ".SubcompNQQ.pdf"), width = 25, height = 30)
      print(g1+g1.1+g2+g3 +plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = "bold", size = 40)))
      dev.off()
    }else{
      j1 <- j[1]
      Y.sub <- cbind(Y.tmp[,j1], rowSums(Y.tmp[,(j1+1):(K+1),drop=FALSE]))
      id.sel <- which(rowSums(Y.sub)!=0)
      
      dw <- ncol(W.tmp); dw.r <- ncol(W.r.tmp); dx <- ncol(X.tmp); dx.r <- ncol(X.r.tmp)
      if(zi.check[j1]){
        gamma0 = matrix(0.001, nrow = dw.r, ncol = 1)
      }else{
        gamma0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      }
      omega0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      alpha0 = matrix(0.001, nrow = dx, ncol = 1)
      beta0 = matrix(0.001, nrow = dx, ncol = 1)
      aigdm.reg <- .AIGDM_EM(Y.sub, W.r.tmp, X.tmp, X.tmp, gamma0, omega0, alpha0, beta0)
      gamma0 = omega0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      gdm.reg <- .AIGDM_EM(Y.sub, W.r.tmp, X.tmp, X.tmp, gamma0, omega0, alpha0, beta0)
      dm.reg <- MGLM::MGLMreg(Y.tmp ~ 0+X.tmp, dist = "DM")
      
      obs <- (Y.sub / rowSums(Y.sub))[id.sel, 1]
      obs.qq <- obs[order(obs)]
      n <- length(obs)
      
      n.iter <- 5000
      sample.dm.qq <- sample.gdm.qq <- sample.ai.qq <- numeric(n.iter * n)
      zeromean.dm <- zeromean.gdm <- zeromean.ai <- numeric(n.iter)
      
      inc <- 1:n
      for(iter in 1:n.iter){
        tmp <- .fitted_model(W.r.tmp, X.tmp, X.tmp, aigdm.reg, "ZIGDM", rowSums(Y.sub))[id.sel]
        sample.ai.qq[inc] <- tmp[order(tmp)]
        zeromean.ai[iter] <- mean(tmp == 0)
        
        tmp <- .fitted_model(W.r.tmp, X.tmp, X.tmp, gdm.reg, NULL, rowSums(Y.sub))[id.sel]
        sample.gdm.qq[inc] <- tmp[order(tmp)]
        zeromean.gdm[iter] <- mean(tmp == 0)
        
        # tmp <- .fitted_model_DM(dm.reg, rowSums(Y.tmp[id.sel,]), j)
        tmp <- .fitted_model_DM(dm.reg, rowSums(Y.tmp), j1)[id.sel]
        sample.dm.qq[inc] <- tmp[order(tmp)]
        zeromean.dm[iter] <- mean(tmp == 0)
        inc <- inc + n
      }
      
      qq_data <- data.table(obs.qq = obs.qq,
                            sample.qq = c(sample.dm.qq, sample.gdm.qq, sample.ai.qq),
                            model = rep(c("DM","GDM", "AIGDM"), each = n.iter*n))
      
      prop_data <- data.frame(
        x = factor(c("Observed", "DM", "GDM", "AIGDM"), levels = c("Observed", "DM", "GDM", "AIGDM")),
        y = c(mean(obs == 0), mean(sample.dm.qq == 0), mean(sample.gdm.qq == 0), mean(sample.ai.qq == 0)),
        low = c(NA, quantile(zeromean.dm, 0.025), quantile(zeromean.gdm, 0.025), quantile(zeromean.ai, 0.025)),
        high = c(NA, quantile(zeromean.dm, 0.975), quantile(zeromean.gdm, 0.975), quantile(zeromean.ai, 0.975))
      )
      
      g2 <- ggplot(prop_data, aes(y = y, x = x)) +
        geom_bar(stat = "identity", fill = newColor[j1], color = "black", width = 0.5, position = position_dodge(0.2)) +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.3, colour = "black", alpha = 0.9) +
        labs(y = expression(bold("Taxon zero proportion"))) +
        ylim(0, max(prop_data$high) + 0.05) +
        theme_bw(base_size=20) +
        theme(axis.title = element_text(face = "bold", size = rel(1.5)),
              axis.title.y = element_text(angle = 90, vjust =2),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              axis.ticks.x = element_blank(),
              axis.text = element_text(),
              axis.line = element_line(colour="black"),
              axis.ticks = element_line(),
              strip.text.x = element_text(size = rel(1.8)),
              plot.margin = unit(c(5, 10, 5, 10), "mm"))
      
      qq_ci <- qq_data[, .(median = median(sample.qq),
                           lower = quantile(sample.qq, 0.025),
                           upper = quantile(sample.qq, 0.975)),
                       by = .(obs.qq, model)]
      qq_ci$model <- factor(qq_ci$model, levels = c("DM", "GDM", "AIGDM"))
      
      g3 <- ggplot(qq_ci, aes(x = obs.qq)) +
        facet_wrap(~ model, nrow = 1) +
        geom_abline(linetype = "dashed", color = "grey") +
        geom_point(data = qq_ci, aes(y = median), color = newColor[j1], size = 2) +
        geom_ribbon(data = qq_ci,
                    aes(ymin = lower, ymax = upper, y = median, fill = "band"),
                    fill = "grey", alpha = 0.4) +
        xlab("Observed quantiles") + ylab("Theoretical quantiles") +
        scale_x_continuous(breaks = seq(0, 1, 0.2)) +
        theme_bw() +
        theme(axis.title = element_text(face = "bold", size = rel(2.6)),
              axis.title.y = element_text(angle = 90, vjust = 2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(size = rel(1.5)),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              axis.line = element_line(colour = "black"),
              strip.text.x = element_text(size = rel(2.6)),
              plot.margin = unit(c(5, 10, 5, 10), "mm"))
      
      j2 <- j[2]
      Y.sub <- cbind(Y.tmp[,j2], rowSums(Y.tmp[,(j2+1):(K+1),drop=FALSE]))
      id.sel <- which(rowSums(Y.sub)!=0)
      
      dw <- ncol(W.tmp); dw.r <- ncol(W.r.tmp); dx <- ncol(X.tmp); dx.r <- ncol(X.r.tmp)
      if(zi.check[j2]){
        gamma0 = matrix(0.001, nrow = dw.r, ncol = 1)
      }else{
        gamma0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      }
      omega0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      alpha0 = matrix(0.001, nrow = dx, ncol = 1)
      beta0 = matrix(0.001, nrow = dx, ncol = 1)
      aigdm.reg <- .AIGDM_EM(Y.sub, W.r.tmp, X.tmp, X.tmp, gamma0, omega0, alpha0, beta0)
      gamma0 = omega0 = matrix(-Inf, nrow = dw.r, ncol = 1)
      gdm.reg <- .AIGDM_EM(Y.sub, W.r.tmp, X.tmp, X.tmp, gamma0, omega0, alpha0, beta0)
      dm.reg <- MGLM::MGLMreg(Y.tmp ~ 0+X.tmp, dist = "DM")
      
      obs <- (Y.sub / rowSums(Y.sub))[id.sel, 1]
      obs.qq <- obs[order(obs)]
      n <- length(obs)
      
      n.iter <- 5000
      sample.dm.qq <- sample.gdm.qq <- sample.ai.qq <- numeric(n.iter * n)
      zeromean.dm <- zeromean.gdm <- zeromean.ai <- numeric(n.iter)
      
      inc <- 1:n
      for(iter in 1:n.iter){
        tmp <- .fitted_model(W.r.tmp, X.tmp, X.tmp, aigdm.reg, "ZIGDM", rowSums(Y.sub))[id.sel]
        sample.ai.qq[inc] <- tmp[order(tmp)]
        zeromean.ai[iter] <- mean(tmp == 0)
        
        tmp <- .fitted_model(W.r.tmp, X.tmp, X.tmp, gdm.reg, NULL, rowSums(Y.sub))[id.sel]
        sample.gdm.qq[inc] <- tmp[order(tmp)]
        zeromean.gdm[iter] <- mean(tmp == 0)
        
        # tmp <- .fitted_model_DM(dm.reg, rowSums(Y.tmp[id.sel,]), j)
        tmp <- .fitted_model_DM(dm.reg, rowSums(Y.tmp), j2)[id.sel]
        sample.dm.qq[inc] <- tmp[order(tmp)]
        zeromean.dm[iter] <- mean(tmp == 0)
        inc <- inc + n
      }
      
      qq_data <- data.table(obs.qq = obs.qq,
                            sample.qq = c(sample.dm.qq, sample.gdm.qq, sample.ai.qq),
                            model = rep(c("DM","GDM", "AIGDM"), each = n.iter*n))
      
      prop_data <- data.frame(
        x = factor(c("Observed", "DM", "GDM", "AIGDM"), levels = c("Observed", "DM", "GDM", "AIGDM")),
        y = c(mean(obs == 0), mean(sample.dm.qq == 0), mean(sample.gdm.qq == 0), mean(sample.ai.qq == 0)),
        low = c(NA, quantile(zeromean.dm, 0.025), quantile(zeromean.gdm, 0.025), quantile(zeromean.ai, 0.025)),
        high = c(NA, quantile(zeromean.dm, 0.975), quantile(zeromean.gdm, 0.975), quantile(zeromean.ai, 0.975))
      )
      
      g4 <- ggplot(prop_data, aes(y = y, x = x)) +
        geom_bar(stat = "identity", fill = newColor[j2], color = "black", width = 0.5, position = position_dodge(0.2)) +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.3, colour = "black", alpha = 0.9) +
        labs(y = expression(bold("Taxon zero proportion"))) +
        ylim(0, max(prop_data$high) + 0.05) +
        theme_bw(base_size=20) +
        theme(axis.title = element_text(face = "bold", size = rel(1.5)),
              axis.title.y = element_text(angle = 90, vjust =2),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              axis.ticks.x = element_blank(),
              axis.text = element_text(),
              axis.line = element_line(colour="black"),
              axis.ticks = element_line(),
              strip.text.x = element_text(size = rel(1.8)),
              plot.margin = unit(c(5, 10, 5, 10), "mm"))
      
      qq_ci <- qq_data[, .(median = median(sample.qq),
                           lower = quantile(sample.qq, 0.025),
                           upper = quantile(sample.qq, 0.975)),
                       by = .(obs.qq, model)]
      qq_ci$model <- factor(qq_ci$model, levels = c("DM", "GDM", "AIGDM"))
      
      g5 <- ggplot(qq_ci, aes(x = obs.qq)) +
        facet_wrap(~ model, nrow = 1) +
        geom_abline(linetype = "dashed", color = "grey") +
        geom_point(data = qq_ci, aes(y = median), color = newColor[j2], size = 2) +
        geom_ribbon(data = qq_ci,
                    aes(ymin = lower, ymax = upper, y = median, fill = "band"),
                    fill = "grey", alpha = 0.4) +
        xlab("Observed quantiles") + ylab("Theoretical quantiles") +
        scale_x_continuous(breaks = seq(0, 1, 0.2)) +
        theme_bw() +
        theme(axis.title = element_text(face = "bold", size = rel(2.6)),
              axis.title.y = element_text(angle = 90, vjust = 2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(size = rel(1.5)),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              axis.line = element_line(colour = "black"),
              strip.text.x = element_text(size = rel(2.6)),
              plot.margin = unit(c(5, 10, 5, 10), "mm"))
      
      layout <- "
AAA
BBB
CDD
EFF"
      
      pdf(paste0("../Figs/ECAM.Full.WO.", taxa_tax, "_", taxa, ".SubcompNQQ.pdf"), width = 25, height = 40)
      print(g1+g1.1+g2+g3+g4+g5 +plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = "bold", size = 40)))
      dev.off()
    }
    
  }else{
    layout <- "
AAA
BBB"
    pdf(paste0("../Figs/ECAM.Full.WO.", taxa_tax, "_", taxa, ".SubcompNQQ.pdf"), width = 25, height = 20)
    print(g1+g1.1+plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = "bold", size = 40)))
    dev.off()
  }
  
}



