rm(list = ls())
setwd("~/Documents/Project/LONGITUDINAL/code/Analysis/")
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
source("2a.utility.R")

# AIGDM proposed model ( GDM-based and LNM-based, N=50, T=5/20) ####
# test 1: alpha!=0 beta=0, power mean
# test 2: alpha=0 beta!=0, power disp
# test 3: alpha=0 beta=0, type 1 mean and disp

###### accuracy of the ZI diagnostic test (final rslts, omni test, cluster-level bootstrap LBFGSB, bootstrap reps = 500)######
tmp_gold_T5 <- read.csv("../Simulation/ProposedModel/rslts_type1_T5_omni.txt", header = FALSE,
                        col.names = c("Test", "ModelSim", "SimRep", "Correlation", "T",
                                      "AIGDM.Zero.A", "AIGDM.Zero.P", "AIGDM.Mean.A", "AIGDM.Mean.P", "AIGDM.Disp.A", "AIGDM.Disp.P", "AIGDM.Omni.A", "AIGDM.Omni.P",
                                      "ZIGDM.Zero.A", "ZIGDM.Zero.P", "ZIGDM.Mean.A", "ZIGDM.Mean.P", "ZIGDM.Disp.A", "ZIGDM.Disp.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P",
                                      "GDM.Zero.A", "GDM.Zero.P", "GDM.Mean.A", "GDM.Mean.P", "GDM.Disp.A", "GDM.Disp.P", "GDM.Omni.A", "GDM.Omni.P",
                                      "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P", 
                                      paste0("CheckZI", 1:5), paste0("ZI", 1:5)))
tmp_gold_T20 <- read.csv("../Simulation/ProposedModel/rslts_type1_T20_omni.txt", header = FALSE,
                         col.names = c("Test", "ModelSim", "SimRep", "Correlation", "T",
                                       "AIGDM.Zero.A", "AIGDM.Zero.P", "AIGDM.Mean.A", "AIGDM.Mean.P", "AIGDM.Disp.A", "AIGDM.Disp.P", "AIGDM.Omni.A", "AIGDM.Omni.P",
                                       "ZIGDM.Zero.A", "ZIGDM.Zero.P", "ZIGDM.Mean.A", "ZIGDM.Mean.P", "ZIGDM.Disp.A", "ZIGDM.Disp.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P",
                                       "GDM.Zero.A", "GDM.Zero.P", "GDM.Mean.A", "GDM.Mean.P", "GDM.Disp.A", "GDM.Disp.P", "GDM.Omni.A", "GDM.Omni.P",
                                       "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P", 
                                       paste0("CheckZI", 1:5), paste0("ZI", 1:5)))
tmp_gold <- rbind(tmp_gold_T5, tmp_gold_T20)
tmp_gold <- tmp_gold[tmp_gold$ModelSim %in% c("ZIGDM0", "ZIGDM1", "ZIGDM2"),]
tmp_wide <- tmp_gold[,c("Test","ModelSim","SimRep","Correlation","T",paste0("CheckZI",1:5),paste0("ZI",1:5))]
tmp_wide[,6:15] <- lapply(tmp_wide[,6:15], function(x) ifelse(x > 0,1,0))
names(tmp_wide)[6:10] <- paste0("ZICheck_Taxon", 1:5)
names(tmp_wide)[11:15] <- paste0("ZIReal_Taxon", 1:5)
tmp_wide <- tmp_wide[,c("Test","ModelSim","SimRep","Correlation","T",rbind(paste0("ZICheck_Taxon", 1:5),paste0("ZIReal_Taxon", 1:5)))]
tmp_long <- reshape(tmp_wide, direction='long', 
                    varying=grep("ZICheck|ZIReal", names(tmp_wide)), sep = "_",
                    v.names=c('ZICheck', 'ZIReal'),
                    timevar='TaxaID',times=paste0("Taxon",1:5),
                    idvar=c("Test","ModelSim","SimRep","Correlation","T"))
tmp_count <- tmp_long %>%
  group_by(Test, ModelSim, Correlation, T, TaxaID, ZIReal) %>%
  count() %>%
  ungroup()
tmp_rslt <- tmp_long %>%
  group_by(Test, ModelSim, Correlation, T, TaxaID, ZIReal) %>%
  summarize_all(mean, na.rm = T) %>%
  dplyr::select(-SimRep) %>%
  ungroup()
tmp_rslt$Count <- tmp_count$n
tmp_rslt_zigdm <- tmp_rslt[tmp_rslt$ModelSim %in% c("ZIGDM0", "ZIGDM1", "ZIGDM2"),]
tmp_rslt_zigdm <- tmp_rslt_zigdm[tmp_rslt_zigdm$Test == 3,]
tmp_rslt_zigdm$ModelSim <- factor(tmp_rslt_zigdm$ModelSim,
                                  labels = c("no ZI",
                                             "ZI at Taxon 1", 
                                             "ZI at Taxa 2-5"),
                                  levels = paste0("ZIGDM", 0:2))
tmp_rslt_zigdm$TaxaID <- factor(tmp_rslt_zigdm$TaxaID,
                                labels = paste0("Taxon ", 1:5),
                                levels = paste0("Taxon", 1:5))
tmp_rslt_zigdm$ZIReal <- factor(tmp_rslt_zigdm$ZIReal,
                                labels = c("classify no ZI",
                                           "classify ZI"),
                                levels = c(0, 1))
tmp_rslt_zigdm$T <- factor(tmp_rslt_zigdm$T, 
                           labels = c("T = 5", "T = 20"),
                           levels = c(5, 20))
summary(tmp_rslt_zigdm)
pdf("../Figs/ZI_diagnostic_GDM_NullModel.pdf", width = 20, height = 15)
ggplot(tmp_rslt_zigdm, aes(x=Correlation, y=ZICheck*100, group=interaction(TaxaID, ZIReal))) +
  geom_line(aes(color=TaxaID, linetype=ZIReal), linewidth = 2) +
  geom_point(aes(color = TaxaID)) +
  facet_grid(T ~ ModelSim) +
  scale_linetype_manual(values = 1:2, 
                        guide = guide_legend(order = 1, title = "Classification task",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 2, byrow= TRUE)) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#984ea3", "#ff7f00"),
                     guide = guide_legend(order = 2, title = "Taxa ID",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 2, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  labs(title = "", y="% correct classification", x="Correlation") +
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"),
        axis.title = element_text(face = "bold",size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")
dev.off()


###### Type 1 (Diff Mean or Diff Disp) Omnibus test, final rslts,  AIGDM/GDM/ZIGDM/DM/QCATC ########
tmp_T5 <- read.csv("../Simulation/ProposedModel/rslts_type1_T5_omni.txt", header = FALSE,
                   col.names = c("Test", "ModelSim", "SimRep", "Correlation", "T",
                                 "AIGDM.Zero.A", "AIGDM.Zero.P", "AIGDM.Mean.A", "AIGDM.Mean.P", "AIGDM.Disp.A", "AIGDM.Disp.P", "AIGDM.Omni.A", "AIGDM.Omni.P",
                                 "ZIGDM.Zero.A", "ZIGDM.Zero.P", "ZIGDM.Mean.A", "ZIGDM.Mean.P", "ZIGDM.Disp.A", "ZIGDM.Disp.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P",
                                 "GDM.Zero.A", "GDM.Zero.P", "GDM.Mean.A", "GDM.Mean.P", "GDM.Disp.A", "GDM.Disp.P", "GDM.Omni.A", "GDM.Omni.P",
                                 "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P", 
                                 paste0("CheckZI", 1:5), paste0("ZI", 1:5)))
head(tmp_T5)
tmp <- tmp_T5[,c("Test", "ModelSim", "SimRep", "Correlation", 
                 "AIGDM.Omni.A", "AIGDM.Omni.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P", "GDM.Omni.A", "GDM.Omni.P",
                 "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P")]

tmp[,5:16] <- tmp[,5:16] < 0.05
tmp_count <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  count() %>%
  ungroup()
tmp_rslt <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  summarize_all(mean, na.rm = T) %>%
  dplyr::select(-SimRep) %>%
  ungroup()
tmp_rslt$SimCount <- tmp_count$n

tmp_long <- gather(tmp_rslt, key = "Model", value = "Type1", AIGDM.Omni.A:QCATC2.Mean.P)
tmp_long$ModelSim <- factor(tmp_long$ModelSim,
                            labels = c("GDM no ZI", "GDM ZI at Taxon 1", "GDM ZI at Taxa 2-5",
                                       "LNM no ZI", "LNM ZI at Taxon 1", "LNM ZI at Taxa 2-5"),
                            levels = c(paste0("ZIGDM", 0:2), paste0("ZILNM", 0:2)))
tmp_long$ModelEval <- sapply(tmp_long$Model, function(x) ifelse(length(grep("AIGDM", x)) > 0, "AIGDM",
                                                                ifelse(length(grep("ZIGDM", x)) > 0, "ZIGDM", 
                                                                       ifelse(length(grep("GDM", x)) > 0, "GDM", 
                                                                              ifelse(length(grep("DM", x)) > 0, "DM", 
                                                                                     ifelse(length(grep("QCATC1", x)) > 0, "QCATC1", "QCATC2"))))))
tmp_long$TestAP <- sapply(tmp_long$Model, function(x) ifelse(length(grep("\\.A", x))>0, "Asym", "Perm"))
tmp_long$ModelEval <- factor(tmp_long$ModelEval,
                             levels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2")) # AIGDM*

range(tmp_long$Type1[tmp_long$ModelEval == "ZIGDM"]/0.05) # (0.916, 13.176)
range(tmp_long$Type1[tmp_long$ModelEval == "DM"]/0.05) # (0.096, 5.232)
tmp_long$Type1 <- ifelse(tmp_long$Type1 > 1.8*0.05, 1.8*0.05, ifelse(tmp_long$Type1 < 0.2*0.05, 0.2*0.05,  tmp_long$Type1))

p1 <- ggplot(tmp_long[tmp_long$TestAP == "Asym" ,],
             aes(x=Correlation, y=Type1/0.05, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(limits = c(0.2, 1.8), breaks = seq(0.2, 1.8, 0.2)) +
  labs(title = "Asymptotic test", y="Type 1 error / 0.05", x="Correlation") + # / 0.05
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")

p1

p2 <- ggplot(tmp_long[tmp_long$TestAP == "Perm" ,],
             aes(x=Correlation, y=Type1/0.05, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(limits = c(0.2, 1.8), breaks = seq(0.2, 1.8, 0.2)) +
  labs(title = "Permutation test", y="Type 1 error / 0.05", x="Correlation") + # / 0.05
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")
p2
pdf("../Figs/Type1_meanNdisp_T5_Omni.pdf", width = 35, height = 20)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

tmp_T20 <- read.csv("../Simulation/ProposedModel/rslts_type1_T20_omni.txt", header = FALSE,
                    col.names = c("Test", "ModelSim", "SimRep", "Correlation", "T",
                                  "AIGDM.Zero.A", "AIGDM.Zero.P", "AIGDM.Mean.A", "AIGDM.Mean.P", "AIGDM.Disp.A", "AIGDM.Disp.P", "AIGDM.Omni.A", "AIGDM.Omni.P",
                                  "ZIGDM.Zero.A", "ZIGDM.Zero.P", "ZIGDM.Mean.A", "ZIGDM.Mean.P", "ZIGDM.Disp.A", "ZIGDM.Disp.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P",
                                  "GDM.Zero.A", "GDM.Zero.P", "GDM.Mean.A", "GDM.Mean.P", "GDM.Disp.A", "GDM.Disp.P", "GDM.Omni.A", "GDM.Omni.P",
                                  "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P", 
                                  paste0("CheckZI", 1:5), paste0("ZI", 1:5)))
head(tmp_T20)
tmp <- tmp_T20[,c("Test", "ModelSim", "SimRep", "Correlation", 
                  "AIGDM.Omni.A", "AIGDM.Omni.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P", "GDM.Omni.A", "GDM.Omni.P",
                  "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P")]

tmp[,5:16] <- tmp[,5:16] < 0.05
tmp_count <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  count() %>%
  ungroup()
tmp_rslt <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  summarize_all(mean, na.rm = T) %>%
  dplyr::select(-SimRep) %>%
  ungroup()
tmp_rslt$SimCount <- tmp_count$n

tmp_long <- gather(tmp_rslt, key = "Model", value = "Type1", AIGDM.Omni.A:QCATC2.Mean.P)
tmp_long$ModelSim <- factor(tmp_long$ModelSim,
                            labels = c("GDM no ZI", "GDM ZI at Taxon 1", "GDM ZI at Taxa 2-5",
                                       "LNM no ZI", "LNM ZI at Taxon 1", "LNM ZI at Taxa 2-5"),
                            levels = c(paste0("ZIGDM", 0:2), paste0("ZILNM", 0:2)))
tmp_long$ModelEval <- sapply(tmp_long$Model, function(x) ifelse(length(grep("AIGDM", x)) > 0, "AIGDM",
                                                                ifelse(length(grep("ZIGDM", x)) > 0, "ZIGDM", 
                                                                       ifelse(length(grep("GDM", x)) > 0, "GDM", 
                                                                              ifelse(length(grep("DM", x)) > 0, "DM", 
                                                                                     ifelse(length(grep("QCATC1", x)) > 0, "QCATC1", "QCATC2"))))))
tmp_long$TestAP <- sapply(tmp_long$Model, function(x) ifelse(length(grep("\\.A", x))>0, "Asym", "Perm"))
tmp_long$ModelEval <- factor(tmp_long$ModelEval,
                             levels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2")) # AIGDM*

range(tmp_long$Type1[tmp_long$ModelEval == "ZIGDM"]/0.05) # (0.916, 18.92)
range(tmp_long$Type1[tmp_long$ModelEval == "DM"]/0.05) # (0.096, 12.1)
tmp_long$Type1 <- ifelse(tmp_long$Type1 > 1.8*0.05, 1.8*0.05, ifelse(tmp_long$Type1 < 0.2*0.05, 0.2*0.05,  tmp_long$Type1))

p1 <- ggplot(tmp_long[tmp_long$TestAP == "Asym" ,],
             aes(x=Correlation, y=Type1/0.05, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(limits = c(0.2, 1.8), breaks = seq(0.2, 1.8, 0.2)) +
  labs(title = "Asymptotic test", y="Type 1 error / 0.05", x="Correlation") + # / 0.05
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")

p1

p2 <- ggplot(tmp_long[tmp_long$TestAP == "Perm" ,],
             aes(x=Correlation, y=Type1/0.05, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(limits = c(0.2, 1.8), breaks = seq(0.2, 1.8, 0.2)) +
  labs(title = "Permutation test", y="Type 1 error / 0.05", x="Correlation") + # / 0.05
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")
p2
pdf("../Figs/Type1_meanNdisp_T20_Omni.pdf", width = 35, height = 20)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

###### Power (Diff Mean or Diff Disp) Omnibus test, final rslts,  AIGDM/GDM/ZIGDM/DM/QCATC ########
tmp_T5 <- read.csv("../Simulation/ProposedModel/rslts_power_T5_omni.txt", header = FALSE,
                   col.names = c("Test", "ModelSim", "SimRep", "Correlation", "T",
                                 "AIGDM.Zero.A", "AIGDM.Zero.P", "AIGDM.Mean.A", "AIGDM.Mean.P", "AIGDM.Disp.A", "AIGDM.Disp.P", "AIGDM.Omni.A", "AIGDM.Omni.P",
                                 "ZIGDM.Zero.A", "ZIGDM.Zero.P", "ZIGDM.Mean.A", "ZIGDM.Mean.P", "ZIGDM.Disp.A", "ZIGDM.Disp.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P",
                                 "GDM.Zero.A", "GDM.Zero.P", "GDM.Mean.A", "GDM.Mean.P", "GDM.Disp.A", "GDM.Disp.P", "GDM.Omni.A", "GDM.Omni.P",
                                 "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P", 
                                 paste0("CheckZI", 1:5), paste0("ZI", 1:5)))
head(tmp_T5)
tmp <- tmp_T5[,c("Test", "ModelSim", "SimRep", "Correlation", "AIGDM.Omni.P", "ZIGDM.Omni.P", "GDM.Omni.P", "DM.Omni.P", "QCATC1.Mean.P", "QCATC2.Mean.P")]

tmp[,5:10] <- tmp[,5:10] < 0.05
tmp_count <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  count() %>%
  ungroup()
tmp_rslt <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  summarize_all(mean, na.rm = T) %>%
  dplyr::select(-SimRep) %>%
  ungroup()
tmp_rslt$SimCount <- tmp_count$n

tmp_long <- gather(tmp_rslt, key = "Model", value = "Power", AIGDM.Omni.P:QCATC2.Mean.P)
tmp_long$ModelSim <- factor(tmp_long$ModelSim,
                            labels = c("GDM no ZI", "GDM ZI at Taxon 1", "GDM ZI at Taxa 2-5",
                                       "LNM no ZI", "LNM ZI at Taxon 1", "LNM ZI at Taxa 2-5"),
                            levels = c(paste0("ZIGDM", 0:2), paste0("ZILNM", 0:2)))
tmp_long$ModelEval <- sapply(tmp_long$Model, function(x) ifelse(length(grep("AIGDM", x)) > 0, "AIGDM",
                                                                ifelse(length(grep("ZIGDM", x)) > 0, "ZIGDM", 
                                                                       ifelse(length(grep("GDM", x)) > 0, "GDM", 
                                                                              ifelse(length(grep("DM", x)) > 0, "DM", 
                                                                                     ifelse(length(grep("QCATC1", x)) > 0, "QCATC1", "QCATC2"))))))
tmp_long$ModelEval <- factor(tmp_long$ModelEval,
                             levels = c("GDM", "ZIGDM", "AIGDM", "DM", "QCATC1", "QCATC2")) # AIGDM*
p1 <- ggplot(tmp_long[tmp_long$Test == 1 ,], # & tmp_long$TestType == "Differential mean test"
             aes(x=Correlation, y=Power, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(title = "Differential mean model", y="Power", x="Correlation") + 
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")
p1
p2 <- ggplot(tmp_long[tmp_long$Test == 2 ,], # & tmp_long$TestType == "Differential mean test"
             aes(x=Correlation, y=Power, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(title = "Differential dispersion model", y="Power", x="Correlation") + 
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")
p2

pdf("../Figs/Power_meanNdisp_T5_Omni.pdf", width = 35, height = 20)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

tmp_T20 <- read.csv("../Simulation/ProposedModel/rslts_power_T20_omni.txt", header = FALSE,
                    col.names = c("Test", "ModelSim", "SimRep", "Correlation", "T",
                                  "AIGDM.Zero.A", "AIGDM.Zero.P", "AIGDM.Mean.A", "AIGDM.Mean.P", "AIGDM.Disp.A", "AIGDM.Disp.P", "AIGDM.Omni.A", "AIGDM.Omni.P",
                                  "ZIGDM.Zero.A", "ZIGDM.Zero.P", "ZIGDM.Mean.A", "ZIGDM.Mean.P", "ZIGDM.Disp.A", "ZIGDM.Disp.P", "ZIGDM.Omni.A", "ZIGDM.Omni.P",
                                  "GDM.Zero.A", "GDM.Zero.P", "GDM.Mean.A", "GDM.Mean.P", "GDM.Disp.A", "GDM.Disp.P", "GDM.Omni.A", "GDM.Omni.P",
                                  "DM.Omni.A", "DM.Omni.P", "QCATC1.Mean.A", "QCATC1.Mean.P", "QCATC2.Mean.A", "QCATC2.Mean.P", 
                                  paste0("CheckZI", 1:5), paste0("ZI", 1:5)))
head(tmp_T20)
tmp <- tmp_T20[,c("Test", "ModelSim", "SimRep", "Correlation", "AIGDM.Omni.P", "ZIGDM.Omni.P", "GDM.Omni.P", "DM.Omni.P", "QCATC1.Mean.P", "QCATC2.Mean.P")]

tmp[,5:10] <- tmp[,5:10] < 0.05
tmp_count <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  count() %>%
  ungroup()
tmp_rslt <- tmp %>%
  group_by(Test, ModelSim, Correlation) %>%
  summarize_all(mean, na.rm = T) %>%
  dplyr::select(-SimRep) %>%
  ungroup()
tmp_rslt$SimCount <- tmp_count$n

tmp_long <- gather(tmp_rslt, key = "Model", value = "Power", AIGDM.Omni.P:QCATC2.Mean.P)
tmp_long$ModelSim <- factor(tmp_long$ModelSim,
                            labels = c("GDM no ZI", "GDM ZI at Taxon 1", "GDM ZI at Taxa 2-5",
                                       "LNM no ZI", "LNM ZI at Taxon 1", "LNM ZI at Taxa 2-5"),
                            levels = c(paste0("ZIGDM", 0:2), paste0("ZILNM", 0:2)))
tmp_long$ModelEval <- sapply(tmp_long$Model, function(x) ifelse(length(grep("AIGDM", x)) > 0, "AIGDM",
                                                                ifelse(length(grep("ZIGDM", x)) > 0, "ZIGDM", 
                                                                       ifelse(length(grep("GDM", x)) > 0, "GDM", 
                                                                              ifelse(length(grep("DM", x)) > 0, "DM", 
                                                                                     ifelse(length(grep("QCATC1", x)) > 0, "QCATC1", "QCATC2"))))))
tmp_long$ModelEval <- factor(tmp_long$ModelEval,
                             levels = c("GDM", "ZIGDM", "AIGDM", "DM", "QCATC1", "QCATC2")) # AIGDM*
p1 <- ggplot(tmp_long[tmp_long$Test == 1 ,], # & tmp_long$TestType == "Differential mean test"
             aes(x=Correlation, y=Power, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(title = "Differential mean model", y="Power", x="Correlation") + 
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")
p1
p2 <- ggplot(tmp_long[tmp_long$Test == 2 ,], # & tmp_long$TestType == "Differential mean test"
             aes(x=Correlation, y=Power, group=ModelEval)) +
  geom_line(aes(linetype=ModelEval, color=ModelEval), linewidth = 1.5) +
  geom_point(aes(shape = ModelEval, color = ModelEval),size=5) +
  facet_wrap( ~ ModelSim, ncol = 3) +
  scale_shape_manual(values = c( "AIGDM" = 1, "ZIGDM" = 1, "GDM" = 1, "DM" = 2, "QCATC1" = 0, "QCATC2" = 0), 
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method", title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_color_manual(values = c("AIGDM" = "#e41a1c", "ZIGDM" = "#984ea3", "GDM" = "#377eb8", "DM" = "#4daf4a", "QCATC1" = "#999999", "QCATC2" = "#ff7f00"),
                     breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                     labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), drop = FALSE,
                     guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                          label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_linetype_manual(values = c("AIGDM" = 1, "ZIGDM" = 3, "GDM" = 2, "DM" = 1, "QCATC1" = 1, "QCATC2" = 1),
                        breaks = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"),
                        labels = c("AIGDM", "GDM", "ZIGDM", "DM", "QCATC1", "QCATC2"), 
                        guide = guide_legend(order = 1, title = "Method",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                             label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(title = "Differential dispersion model", y="Power", x="Correlation") + 
  theme_bw(base_size = 12) +
  theme(legend.key.width = unit(2,"cm"), 
        plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5),
        axis.title = element_text(face = "bold", size = rel(2)),
        axis.title.y = element_text(angle = 90, vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)),
        strip.text = element_text(size = rel(2)),
        plot.margin = unit(c(5, 10, 5, 10), "mm"),
        legend.position="bottom")
p2

pdf("../Figs/Power_meanNdisp_T20_Omni.pdf", width = 35, height = 20)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()


