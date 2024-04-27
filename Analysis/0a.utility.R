get_nonrare_taxa <- function(physeq_obj, threshold = 0.1, discard_other = FALSE, other_label = "Other"){
  
  #Define a temporary physeq object
  ps_tmp <- physeq_obj
  
  #Check for 0 entries
  smpl_sms <- phyloseq::sample_sums(ps_tmp)
  if (0 %in% smpl_sms){
    stop("Error: some samples contain 0 reads. These have to be removed to avoid
         downstream problems.")
  }
  
  #Extract the otu_table as a data.frame
  otu_tbl <- phyloseq::otu_table(ps_tmp)
  if (!phyloseq::taxa_are_rows(ps_tmp)){
    otu_tbl <- t(otu_tbl)
  }
  #Get the top taxa names and discard or merge other taxa
  abun_taxa <- names(which(rowMeans(otu_tbl>0) > threshold))
  if (discard_other){
    physeq_obj <- phyloseq::prune_taxa(abun_taxa, physeq_obj)
  } else {
    to_merge <- phyloseq::taxa_names(physeq_obj)
    to_merge <- to_merge[!(to_merge %in% abun_taxa)]
    physeq_obj <- merge_taxa(physeq_obj, to_merge)
    tax_tbl <- phyloseq::tax_table(physeq_obj)
    indx <- which(row.names(tax_tbl) %in% to_merge)
    tax_tbl[indx,] <- other_label
    phyloseq::tax_table(physeq_obj) <- tax_tbl
  }
  return(physeq_obj)
}

get_summary4PTB <- function(meta, study){
  if(study == "B&J"){
    meta.tmp <- meta[meta$project %in% c("B", "J"),] # & meta$NIH.Racial.Category == "Black or African American",]
  }else{
    meta.tmp <- meta[meta$project == study,] # & (meta$NIH.Racial.Category == "Black or African American"),]
  }
 
  meta.tmp$collect_wk <- factor(meta.tmp$collect_wk, labels = as.character(1:41), levels = 1:41)
  meta.tmp$was_term[meta.tmp$was_term == "False"] <- "Term Birth"
  meta.tmp$was_term[meta.tmp$was_term == "True"] <- "Preterm Birth"
  meta.tmp.rmdup <- data.frame(as.matrix(meta.tmp[,c("participant_id", "collect_wk", "was_term")])) # remove reps for generating plot
  meta.tmp.rmdup <- meta.tmp.rmdup[!duplicated(meta.tmp.rmdup),] 
  meta.timepoint <- expand.grid(participant_id = unique(meta.tmp.rmdup$participant_id), collect_wk = 1:41, value = 0, was_term = 0)
  meta.timepoint$collect_wk <- factor(meta.timepoint$collect_wk, levels = 1:41)
  meta.timepoint$value[paste(meta.timepoint$participant_id, meta.timepoint$collect_wk, sep = ".") %in% 
                         paste(meta.tmp.rmdup$participant_id, meta.tmp.rmdup$collect_wk, sep = ".")] <- 1
  term_ids <- unique(meta.tmp.rmdup$participant_id[meta.tmp.rmdup$was_term == "Term Birth"])
  preterm_ids <- unique(meta.tmp.rmdup$participant_id[meta.tmp.rmdup$was_term == "Preterm Birth"])
  meta.timepoint$was_term[meta.timepoint$participant_id %in% term_ids] <- "Term Birth"
  meta.timepoint$was_term[meta.timepoint$participant_id %in% preterm_ids] <- "Preterm Birth"
  meta.timepoint$fill <- ifelse(meta.timepoint$value == 0, "white", 
                                ifelse(meta.timepoint$was_term == "Term Birth", "#e41a1c", "#377eb8"))
  meta.timepoint$value <- factor(meta.timepoint$value, levels = 0:1)
  
  g1 <- ggplot(meta.timepoint, aes(collect_wk, participant_id)) +
    geom_tile(aes(fill = I(fill)), color = "#2b2b2b", size=0.125) +
    scale_x_discrete() +
    scale_y_discrete() +
    theme(axis.text.x = element_text(size=14, angle = 90)) + 
    labs(title = paste0("Study ", study), x = "Collect weeks", y = "Participants in Study") + 
    theme_bw(base_size=14) +
    theme(axis.title = element_text(face = "bold",size = rel(1.5)),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = rel(1.5), angle = 90),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(),
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          strip.text.x = element_text(size = rel(1.8)),
          plot.margin = unit(c(5, 10, 5, 10), "mm"),
          legend.position = c(0.85, 0.1),
          legend.background = element_rect(fill = "white", color = "black"))
  
  for (alphaDiv in c("RootedPD", "UnrootedPD", "Shannon", "InvSimpson")) {
    if(alphaDiv == "RootedPD"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_term) %>%
        summarise(tmp.mean = mean(rooted_pd), tmp.se = sd(rooted_pd)) %>%
        ungroup()
    }else if(alphaDiv == "UnrootedPD"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_term) %>%
        summarise(tmp.mean = mean(unrooted_pd), tmp.se = sd(unrooted_pd)) %>%
        ungroup()
    }else if(alphaDiv == "Shannon"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_term) %>%
        summarise(tmp.mean = mean(shannon), tmp.se = sd(shannon)) %>%
        ungroup()
    }else if(alphaDiv == "InvSimpson"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_term) %>%
        summarise(tmp.mean = mean(inv_simpson), tmp.se = sd(inv_simpson)) %>%
        ungroup()
    }
    
    pd <- position_dodge(0.2)
    g2 <- ggplot(meta.tmp.group, aes(collect_wk, tmp.mean, 
                                     color = was_term, fill = was_term, group = was_term)) + 
      geom_errorbar(aes(ymin = tmp.mean-tmp.se, ymax = tmp.mean+tmp.se),
                    width = 0.5, linewidth = 0.4, position=pd) +
      geom_line(position=pd) +
      geom_point(position=pd, size=3, shape=21) + 
      theme(axis.text.x = element_text(size=14, angle = 90)) + 
      scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
      scale_color_manual(values = c("#377eb8", "#e41a1c"), guide = "none") +
      labs(title = paste0("Study ", study), x = "Collect weeks", y = alphaDiv)  + 
      theme_bw(base_size=14) +
      theme(axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle = 90, vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5), angle = 90),
            axis.ticks.x = element_blank(),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            strip.text.x = element_text(size = rel(1.8)),
            plot.margin = unit(c(5, 10, 5, 10), "mm"),
            legend.position = c(0.85, 0.1),
            legend.background = element_rect(fill = "white", color = "black")) + 
      guides(fill = guide_legend(title = "",
                                 title.theme = element_text(face = "bold", size = 12, angle = 0),
                                 label.theme = element_text(size = 12, angle = 0)))
    print(ggpubr::ggarrange(g1, g2, ncol=2, nrow=1))
  }
}

get_summary4ePTB <- function(meta, study){
  if(study == "B&J"){
    meta.tmp <- meta[meta$project %in% c("B", "J"),]
  }else{
    meta.tmp <- meta[meta$project == study,]
  }
  
  meta.tmp$collect_wk <- factor(meta.tmp$collect_wk, labels = as.character(1:41), levels = 1:41)
  meta.tmp$was_early_preterm[meta.tmp$was_early_preterm == "False"] <- "Not Early Preterm Birth"
  meta.tmp$was_early_preterm[meta.tmp$was_early_preterm == "True"] <- "Early Preterm Birth"
  meta.tmp.rmdup <- data.frame(as.matrix(meta.tmp[,c("participant_id", "collect_wk", "was_early_preterm")])) # remove reps for generating plot
  meta.tmp.rmdup <- meta.tmp.rmdup[!duplicated(meta.tmp.rmdup),] 
  meta.timepoint <- expand.grid(participant_id = unique(meta.tmp.rmdup$participant_id), collect_wk = 1:41, value = 0, was_early_preterm = 0)
  meta.timepoint$collect_wk <- factor(meta.timepoint$collect_wk, levels = 1:41)
  meta.timepoint$value[paste(meta.timepoint$participant_id, meta.timepoint$collect_wk, sep = ".") %in% 
                         paste(meta.tmp.rmdup$participant_id, meta.tmp.rmdup$collect_wk, sep = ".")] <- 1
  term_ids <- unique(meta.tmp.rmdup$participant_id[meta.tmp.rmdup$was_early_preterm == "Not Early Preterm Birth"])
  preterm_ids <- unique(meta.tmp.rmdup$participant_id[meta.tmp.rmdup$was_early_preterm == "Early Preterm Birth"])
  meta.timepoint$was_early_preterm[meta.timepoint$participant_id %in% term_ids] <- "Not Early Preterm Birth"
  meta.timepoint$was_early_preterm[meta.timepoint$participant_id %in% preterm_ids] <- "Early Preterm Birth"
  meta.timepoint$fill <- ifelse(meta.timepoint$value == 0, "white", 
                                ifelse(meta.timepoint$was_early_preterm == "Not Early Preterm Birth", "#e41a1c", "#377eb8"))
  meta.timepoint$value <- factor(meta.timepoint$value, levels = 0:1)
  
  g1 <- ggplot(meta.timepoint, aes(collect_wk, participant_id)) +
    geom_tile(aes(fill = I(fill)), color = "#2b2b2b", size=0.125) +
    scale_x_discrete() +
    scale_y_discrete() +
    theme(axis.text.x = element_text(size=14, angle = 90)) + 
    labs(title = paste0("Study ", study), x = "Collect weeks", y = "Participants in Study") + 
    theme_bw(base_size=14) +
    theme(axis.title = element_text(face = "bold",size = rel(1.5)),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = rel(1.5), angle = 90),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(),
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          strip.text.x = element_text(size = rel(1.8)),
          plot.margin = unit(c(5, 10, 5, 10), "mm"),
          legend.position = c(0.85, 0.1),
          legend.background = element_rect(fill = "white", color = "black"))
  
  for (alphaDiv in c("RootedPD", "UnrootedPD", "Shannon", "InvSimpson")) {
    if(alphaDiv == "RootedPD"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_early_preterm) %>%
        summarise(tmp.mean = mean(rooted_pd), tmp.se = sd(rooted_pd)) %>%
        ungroup()
    }else if(alphaDiv == "UnrootedPD"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_early_preterm) %>%
        summarise(tmp.mean = mean(unrooted_pd), tmp.se = sd(unrooted_pd)) %>%
        ungroup()
    }else if(alphaDiv == "Shannon"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_early_preterm) %>%
        summarise(tmp.mean = mean(shannon), tmp.se = sd(shannon)) %>%
        ungroup()
    }else if(alphaDiv == "InvSimpson"){
      meta.tmp.group <- meta.tmp %>%
        group_by(collect_wk, was_early_preterm) %>%
        summarise(tmp.mean = mean(inv_simpson), tmp.se = sd(inv_simpson)) %>%
        ungroup()
    }
    
    pd <- position_dodge(0.2)
    g2 <- ggplot(meta.tmp.group, aes(collect_wk, tmp.mean, 
                                     color = was_early_preterm, fill = was_early_preterm, group = was_early_preterm)) + 
      geom_errorbar(aes(ymin = tmp.mean-tmp.se, ymax = tmp.mean+tmp.se),
                    width = 0.5, linewidth = 0.4, position=pd) +
      geom_line(position=pd) +
      geom_point(position=pd, size=3, shape=21) + 
      theme(axis.text.x = element_text(size=14, angle = 90)) + 
      scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
      scale_color_manual(values = c("#377eb8", "#e41a1c"), guide = "none") +
      labs(title = paste0("Study ", study), x = "Collect weeks", y = alphaDiv)  + 
      theme_bw(base_size=14) +
      theme(axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle = 90, vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5), angle = 90),
            axis.ticks.x = element_blank(),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            strip.text.x = element_text(size = rel(1.8)),
            plot.margin = unit(c(5, 10, 5, 10), "mm"),
            legend.position = c(0.85, 0.1),
            legend.background = element_rect(fill = "white", color = "black")) + 
      guides(fill = guide_legend(title = "",
                                 title.theme = element_text(face = "bold", size = 12, angle = 0),
                                 label.theme = element_text(size = 12, angle = 0)))
    print(ggpubr::ggarrange(g1, g2, ncol=2, nrow=1))
  }
}

make_unique <- function(names) {
  # Find duplicated names
  dupes <- names[duplicated(names) | duplicated(names, fromLast = TRUE)]
  
  # For each duplicated name, add a suffix
  for (dupe in unique(dupes)) {
    indices <- which(names == dupe)
    new_names <- paste0(dupe, "_", seq_along(indices))
    names[indices] <- new_names
  }
  
  return(names)
}

# Function to check for inconsistencies in taxonomic ranks
check_inconsistencies <- function(data) {
  # Columns to check
  columns <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  
  # An empty dataframe to store rows with inconsistencies
  inconsistent_rows <- data.frame()
  
  for (i in 6:2) { # Starting from genus and moving up
    current_rank <- columns[i]
    parent_rank <- columns[i-1]
    
    # Group by current rank and summarize the unique parent ranks
    inconsistencies <- data %>%
      dplyr::group_by(.data[[current_rank]]) %>%
      dplyr::summarise(n_unique = length(unique(.data[[parent_rank]]))) %>%
      dplyr::filter(n_unique > 1) %>%
      dplyr::pull(current_rank)
    
    # If inconsistencies found, store rows in the dataframe
    if (length(inconsistencies) > 0) {
      for (inconsistency in inconsistencies) {
        inconsistent_rows <- rbind(inconsistent_rows, data[data[[current_rank]] == inconsistency, ])
      }
    }
  }
  
  return(inconsistent_rows)
}
