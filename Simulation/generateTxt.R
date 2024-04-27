rm(list = ls())
setwd("~/Documents/Project/LONGITUDINAL/code/Simulation/")

tmp_power_N50 <- do.call(expand.grid, list(test=1:2, 
                                 mod=c("ZIGDM0", "ZIGDM1", "ZIGDM2",
                                       "ZILNM0", "ZILNM1", "ZILNM2"), # 
                                 sim=1:2000, r=seq(0, 0.8, 0.2),
                                 n.T=c(5,20))) 
tmp_type1_N50 <- do.call(expand.grid, list(test=3, 
                                           mod=c("ZIGDM0", "ZIGDM1", "ZIGDM2",
                                                 "ZILNM0", "ZILNM1", "ZILNM2"), 
                                           sim=1:5000, r=seq(0, 0.8, 0.2),
                                           n.T=c(5, 20)))
tmp <- rbind(tmp_power_N50, tmp_type1_N50)
tmp$outfile <- paste0("Test", tmp$test, "Mod", tmp$mod, "Sim", tmp$sim, "Corr", tmp$r, "T", tmp$n.T, ".txt")
anyDuplicated(tmp)
write.table(tmp, file = "ProposedModel/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)