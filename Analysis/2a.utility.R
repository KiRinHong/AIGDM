
theme_Publication <- function(base_size=14, base_family="") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.5), hjust = 0.5),
            text = element_text(size = 15),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5)),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA), # key background
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1.25, "cm"),
            legend.key.height = NULL,                # key height (unit)
            legend.key.width = NULL,                 # key width (unit)
            legend.text = element_text(size=rel(1.5)),
            legend.title = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold", size = rel(1.5))
    ))
}

PlotType1 <- function(dat, TypeO.sel, TypeZI.sel, title){
  #dat = tmp_rslt; TypeO.sel = 2; TypeZI.sel = 2; title = ""
  type1.tab.w <- dat %>%
    filter(TypeO == TypeO.sel, TypeZI == TypeZI.sel) %>%
    select(N, ModelSim, Correlation, 
           DM.P, GDM.A, GDM.P, ZIGDM.A, ZIGDM.P, 
           AZIGDM_A.A, AZIGDM_A.P, AZIGDM_B.A, AZIGDM_B.P,
           QCATC1.A, QCATC1.P, QCATC2.A, QCATC2.P)
  if(TypeZI.sel == 2){
    type1.tab.w_more <- dat %>%
      filter(TypeO == TypeO.sel, TypeZI == 1, ModelSim == "ZIGDM0") %>%
      mutate(TypeZI = TypeZI.sel) %>%
      select(N, ModelSim, Correlation, 
             DM.P, GDM.A, GDM.P, ZIGDM.A, ZIGDM.P, 
             AZIGDM_A.A, AZIGDM_A.P, AZIGDM_B.A, AZIGDM_B.P,
             QCATC1.A, QCATC1.P, QCATC2.A, QCATC2.P)
    type1.tab.w <- rbind(type1.tab.w_more, type1.tab.w)
  }
  
  type1.tab.l <- gather(type1.tab.w, ModelTest, Value, DM.P:QCATC2.P, factor_key=TRUE)
  type1.tab.l$N <- ifelse(type1.tab.l$N == 50, "N = 50",
                          ifelse(type1.tab.l$N == 200, "N = 200", type1.tab.l$N))
  type1.tab.l$Model <- sub("\\.[AP]$", "", type1.tab.l$ModelTest)
  type1.tab.l$Model <- factor(type1.tab.l$Model, 
                              levels = c("AZIGDM_A", "AZIGDM_B", "ZIGDM", "GDM", "DM", "QCATC1", "QCATC2"))
  type1.tab.l$TypeT <- sub(".*\\.", "", type1.tab.l$ModelTest)
  type1.tab.l$TypeT[type1.tab.l$TypeT == "A"] <- "Asymptotic"
  type1.tab.l$TypeT[type1.tab.l$TypeT == "P"] <- "Permutation"
  type1.tab.l$TypeT <- factor(type1.tab.l$TypeT, levels = c("Asymptotic", "Permutation"))
  type1.tab.l$ModelSim <- ifelse(type1.tab.l$ModelSim == "ZIGDM0", "No ZI",
                                 ifelse(type1.tab.l$ModelSim == "ZIGDM1", "Full Low-ZI",
                                        ifelse(type1.tab.l$ModelSim == "ZIGDM2", "Full High-ZI",
                                               ifelse(type1.tab.l$ModelSim == "ZIGDM3", "Partial Low-ZI", 
                                                      ifelse(type1.tab.l$ModelSim == "ZIGDM4", "Partial High-ZI", type1.tab.l$ModelSim)))))
  type1.tab.l$Group <- paste(type1.tab.l$N, type1.tab.l$ModelSim, sep = ", ")
  type1.tab.l$Group <- factor(type1.tab.l$Group, 
                              levels = c("N = 50, No ZI", "N = 50, Partial Low-ZI", "N = 50, Partial High-ZI", 
                                         "N = 50, Full Low-ZI", "N = 50, Full High-ZI", 
                                         "N = 200, No ZI", "N = 200, Partial Low-ZI", "N = 200, Partial High-ZI", 
                                         "N = 200, Full Low-ZI", "N = 200, Full High-ZI"))
  type1.tab.l$ModelTest <- factor(type1.tab.l$ModelTest, 
                                  levels = c("AZIGDM_A.A", "AZIGDM_A.P", "AZIGDM_B.A", "AZIGDM_B.P",
                                             "ZIGDM.A", "ZIGDM.P", 
                                             "GDM.A", "GDM.P", "DM.P",
                                             "QCATC1.A", "QCATC1.P", "QCATC2.A", "QCATC2.P"))
  
  g <- ggplot(type1.tab.l, aes(x = Correlation, y = Value, group = ModelTest, 
                          color = Model, shape = Model, linetype = TypeT)) + 
    geom_hline(yintercept = 0.05, color = "#999999", size = 1) +
    geom_line(size = 0.8) + geom_point(size = 2.5) + 
    scale_size_manual(values = rep(2.5, 13), guide = guide_legend(title = "",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                                                  label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
    scale_shape_manual(values = c(1,2,rep(1, 5)), guide = guide_legend(title = "",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                                                 label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
    scale_color_manual(values = c("cyan", "#cab2d6","#377eb8", "#e41a1c", "#a65628", "#4daf4a", "#ff7f00"), 
                       guide = guide_legend(title = "",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                            label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
    scale_linetype_manual(values = c(rep(c(2,1),3), 1, rep(c(2,1),3)), guide = 
                            guide_legend(title = "", label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) + # c(1, 1, 3, 6, 1, 1, 1, 3, 3, 6, 6, 1, 1)
    facet_wrap(~ Group, ncol = 5, labeller = label_wrap_gen(multi_line = FALSE)) + #facet_grid(N~ModelSim, scales="free_y") +
    scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    labs(title = title, y="Type I error", x="Correlation") + 
    theme_bw(base_size = 14) +
    theme(legend.key.width = unit(1,"cm"), 
          axis.title = element_text(face = "bold",size = rel(1.5)),
          axis.title.y = element_text(angle = 90, vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_text(size = rel(1.5)),
          axis.text.x = element_text(size = rel(1.5)),
          strip.text = element_text(size = rel(2)),
          plot.margin = unit(c(5, 10, 5, 10), "mm"),
          legend.position="bottom")
  return(g)
}

PlotPower <- function(dat, TypeO.sel, TypeZI.sel, title){
  #dat = tmp_rslt; TypeO.sel = 1; TypeZI.sel = 1; title = "Decresing_Low2High"
  power.tab.w <- dat %>%
    filter(TypeO == TypeO.sel, TypeZI == TypeZI.sel) %>%
    select(N, ModelSim, Correlation, TypeO, TypeZI, 
           DM.P, GDM.P, ZIGDM.P, AZIGDM_A.P, AZIGDM_B.P, QCATC1.P, QCATC2.P)
  if(TypeZI.sel == 2){
    power.tab.w_more <- dat %>%
      filter(TypeO == TypeO.sel, TypeZI == 1, ModelSim == "ZIGDM0") %>%
      mutate(TypeZI = TypeZI.sel) %>%
      select(N, ModelSim, Correlation, TypeO, TypeZI, 
             DM.P, GDM.P, ZIGDM.P, AZIGDM_A.P, AZIGDM_B.P, QCATC1.P, QCATC2.P)
    power.tab.w <- as.data.frame(rbind(power.tab.w, power.tab.w_more))
  }
  power.tab.l <- gather(power.tab.w, ModelTest, Value, DM.P:QCATC2.P, factor_key=TRUE)
  power.tab.l$N <- ifelse(power.tab.l$N == 50, "N = 50",
                          ifelse(power.tab.l$N == 200, "N = 200", power.tab.l$N))
  power.tab.l$Model <- sub("\\..*", "", power.tab.l$ModelTest)
  power.tab.l$TypeT <- sub(".*\\.", "", power.tab.l$ModelTest)
  power.tab.l$ModelTest <- sub("\\.P", "", power.tab.l$ModelTest)
  power.tab.l$ModelSim <- ifelse(power.tab.l$ModelSim == "ZIGDM0", "No ZI",
                                 ifelse(power.tab.l$ModelSim == "ZIGDM1", "Full Low-ZI",
                                        ifelse(power.tab.l$ModelSim == "ZIGDM2", "Full High-ZI",
                                               ifelse(power.tab.l$ModelSim == "ZIGDM3", "Partial Low-ZI", 
                                                      ifelse(power.tab.l$ModelSim == "ZIGDM4", "Partial High-ZI", power.tab.l$ModelSim)))))
  power.tab.l$Group <- paste(power.tab.l$N, power.tab.l$ModelSim, sep = ", ")
  power.tab.l$Group <- factor(power.tab.l$Group, 
                              levels = c("N = 50, No ZI", "N = 50, Partial Low-ZI", "N = 50, Partial High-ZI", 
                                         "N = 50, Full Low-ZI", "N = 50, Full High-ZI", 
                                         "N = 200, No ZI", "N = 200, Partial Low-ZI", "N = 200, Partial High-ZI", 
                                         "N = 200, Full Low-ZI", "N = 200, Full High-ZI"))
  power.tab.l$ModelTest <- factor(power.tab.l$ModelTest, 
                                  levels = c("AZIGDM_A","AZIGDM_B","ZIGDM","GDM","DM","QCATC1","QCATC2"))
  
  g <- ggplot(power.tab.l, aes(x = Correlation, y = Value, group = ModelTest, 
                          color = ModelTest, shape = ModelTest, linetype = ModelTest)) + 
    geom_line(size = 1) + 
    geom_point(size = 2.5) + 
    scale_linetype_manual(values = rep(1, 7), guide = guide_legend(title = "",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                                                   label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) + # c(1, 1, 3, 6, 1, 1, 1, 3, 3, 6, 6, 1, 1)
    scale_size_manual(values = rep(2.5, 7), guide = guide_legend(title = "",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                                                 label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
    scale_shape_manual(values = c(1, 2, rep(1, 5)), guide = guide_legend(title = "",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                                                label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
    scale_color_manual(values = c("cyan", "#cab2d6", "#377eb8", "#e41a1c", "#a65628", "#4daf4a", "#ff7f00"), 
                       guide = guide_legend(title = "",title.theme = element_text(face = "bold", size = 25, angle = 0),
                                            label.theme = element_text(size = 25, angle = 0), nrow = 1, byrow= TRUE)) +
    facet_wrap(~ Group, ncol = 5, scales = "free_y", labeller = label_wrap_gen(multi_line = FALSE)) + #facet_grid(N~ModelSim, scales="free_y") +
    scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
    labs(title = title, y="Power", x="Correlation") + 
    theme_bw(base_size = 14) +
    theme(legend.key.width = unit(1,"cm"), 
          axis.title = element_text(face = "bold",size = rel(1.5)),
          axis.title.y = element_text(angle = 90, vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_text(size = rel(1.5)),
          axis.text.x = element_text(size = rel(1.5)),
          strip.text = element_text(size = rel(2)),
          plot.margin = unit(c(5, 10, 5, 10), "mm"),
          legend.position="bottom") 
  return(g)
}
