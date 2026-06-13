# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): Figure5E; FigureS9A; FigureS9B
# Table(s): Supplementary Table 14; Supplementary Table 15
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): none detected
# Main output(s): none detected
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
# =========================================================
# Libraries
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
})


conflicted::conflicts_prefer(ggplot2::margin)
# =========================================================
# Configuration
# =========================================================
out_dir <- "C:/Users/zqr20/Documents/tjmeta/BIG_revision"

species_single <- "Bifidobacterium pseudocatenulatum"

species_multi <- c(
  "Bifidobacterium longum",
  "Bifidobacterium pseudocatenulatum",
  "Bifidobacterium bifidum",
  "Bifidobacterium adolescentis"
)

traj_palette <- c(
  "1" = "#e95280",
  "2" = "#23b1a5",
  "3" = "#E49B0F"
)

datasets <- list(
  list(
    path = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy.rda",
    label = "medium and high quality",
    tag = "medium_high_quality"
  ),
  list(
    path = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda",
    label = "only high quality",
    tag = "only_high_quality"
  )
)

# =========================================================
# Load data
# =========================================================
load_spreadcazy_rda <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  
  if (!exists("spreadcazy", envir = e, inherits = FALSE)) {
    stop("Object 'spreadcazy' was not found in: ", path)
  }
  
  as.data.frame(e$spreadcazy)
}

standardize_cazy_data <- function(df) {
  required_cols <- c("class_1", "speciesname", "transmit")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df$class_1 <- factor(as.character(df$class_1), levels = c("1", "2", "3"))
  df$transmit <- factor(as.character(df$transmit))
  df$speciesname <- as.character(df$speciesname)
  
  df
}

# =========================================================
# Prepare plot data with fixed column names
# =========================================================
prepare_plot_data <- function(df, response, facet_col) {
  if (!response %in% names(df)) {
    stop("Response column not found: ", response)
  }
  if (!facet_col %in% names(df)) {
    stop("Facet column not found: ", facet_col)
  }
  
  out <- df
  out$value <- out[[response]]
  out$facet_group <- as.character(out[[facet_col]])
  out$class_1 <- factor(as.character(out$class_1), levels = c("1", "2", "3"))
  out <- out[!is.na(out$value) & !is.na(out$class_1) & !is.na(out$facet_group), ]
  out
}

# =========================================================
# Wilcoxon tests, manual and robust
# =========================================================
make_wilcox_table <- function(dat) {
  facet_levels <- unique(dat$facet_group)
  res_list <- list()
  
  for (fg in facet_levels) {
    sub <- dat[dat$facet_group == fg, , drop = FALSE]
    sub <- sub[!is.na(sub$value) & !is.na(sub$class_1), , drop = FALSE]
    
    obs_levels <- unique(as.character(sub$class_1))
    obs_levels <- obs_levels[!is.na(obs_levels)]
    
    if (length(obs_levels) < 2) next
    
    cmb <- t(combn(obs_levels, 2))
    if (nrow(cmb) == 0) next
    
    one_facet <- vector("list", nrow(cmb))
    
    y_max <- max(sub$value, na.rm = TRUE)
    y_min <- min(sub$value, na.rm = TRUE)
    spread <- y_max - y_min
    if (!is.finite(spread) || spread == 0) {
      spread <- max(abs(c(y_min, y_max)), 1)
    }
    
    for (i in seq_len(nrow(cmb))) {
      g1 <- cmb[i, 1]
      g2 <- cmb[i, 2]
      
      x1 <- sub$value[as.character(sub$class_1) == g1]
      x2 <- sub$value[as.character(sub$class_1) == g2]
      
      wt <- tryCatch(
        wilcox.test(x1, x2, exact = FALSE, conf.int = TRUE),
        error = function(e) NULL
      )
      
      pval <- if (is.null(wt)) NA_real_ else wt$p.value
      est  <- if (is.null(wt) || is.null(wt$estimate)) NA_real_ else unname(wt$estimate)

pval <- if (is.null(wt)) NA_real_ else wt$p.value
est  <- if (is.null(wt) || is.null(wt$estimate)) NA_real_ else unname(wt$estimate)
      
one_facet[[i]] <- data.frame(
  facet_group = fg,
  group1 = g1,
  group2 = g2,
  estimate = est,
  p = pval,
  stringsAsFactors = FALSE
)
    }
    
    one_facet <- do.call(rbind, one_facet)
    one_facet <- one_facet[!is.na(one_facet$p), , drop = FALSE]
    if (nrow(one_facet) == 0) next
    
    one_facet$p.adj <- p.adjust(one_facet$p, method = "BH")
    
    one_facet$p.adj.signif <- cut(
      one_facet$p.adj,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "ns")
    )
    
    one_facet$y.position <- y_max + 0.08 * spread * seq_len(nrow(one_facet))
    
    res_list[[length(res_list) + 1]] <- one_facet
  }
  
  if (length(res_list) == 0) return(data.frame())
  do.call(rbind, res_list)
}

# =========================================================
# Bottom labels, one label per trajectory in each facet
# =========================================================
make_bottom_labels <- function(dat) {
  facet_levels <- unique(dat$facet_group)
  out_list <- list()
  
  for (fg in facet_levels) {
    sub <- dat[dat$facet_group == fg, , drop = FALSE]
    sub <- sub[!is.na(sub$value) & !is.na(sub$class_1), , drop = FALSE]
    
    if (nrow(sub) == 0) next
    
    y_min <- min(sub$value, na.rm = TRUE)
    y_max <- max(sub$value, na.rm = TRUE)
    spread <- y_max - y_min
    if (!is.finite(spread) || spread == 0) {
      spread <- max(abs(c(y_min, y_max)), 1)
    }
    
    y_base <- y_min - 0.15 * spread
    
    levs <- unique(as.character(sub$class_1))
    levs <- levs[!is.na(levs)]
    
    tmp <- lapply(levs, function(lv) {
      x <- sub$value[as.character(sub$class_1) == lv]
      if (length(x) == 0 || all(is.na(x))) return(NULL)
      
      data.frame(
        facet_group = fg,
        class_1 = lv,
        median_val = median(x, na.rm = TRUE),
        n = sum(!is.na(x)),
        y = y_base,
        stringsAsFactors = FALSE
      )
    })
    
    tmp <- do.call(rbind, tmp)
    if (is.null(tmp) || nrow(tmp) == 0) next
    
    tmp$label <- paste0("median=", round(tmp$median_val, 2), "\n", "n=", tmp$n)
    tmp$class_1 <- factor(tmp$class_1, levels = c("1", "2", "3"))
    
    out_list[[length(out_list) + 1]] <- tmp
  }
  
  if (length(out_list) == 0) return(data.frame())
  do.call(rbind, out_list)
}

# =========================================================
# Generic plotting function
# =========================================================
plot_boxplot_panel <- function(df, response, facet_col, dataset_label, title_text, ylab_text,
                               facet_nrow = NULL, stat_filename = NULL){
  dat <- prepare_plot_data(df, response = response, facet_col = facet_col)
  stat_tbl <- make_wilcox_table(dat)
  label_tbl <- make_bottom_labels(dat)
  
  print(stat_tbl)
  if (!is.null(stat_filename) && nrow(stat_tbl) > 0) {
    write.csv(
      stat_tbl,
      file = file.path(out_dir, stat_filename),
      row.names = FALSE
    )
  }
  dat$facet_group <- factor(dat$facet_group, levels = unique(dat$facet_group))
  label_tbl$facet_group <- factor(label_tbl$facet_group, levels = levels(dat$facet_group))
  if (nrow(stat_tbl) > 0) {
    stat_tbl$facet_group <- factor(stat_tbl$facet_group, levels = levels(dat$facet_group))
  }
  
  p <- ggplot(dat, aes(x = class_1, y = value, fill = class_1)) +
    geom_boxplot(  alpha = 0.7,
                   outlier.shape = 21,
                   outlier.size = 1,
                   outlier.stroke = 0.5,
                   outlier.colour = "black",
                   outlier.fill = "white") +
    geom_text(
      data = label_tbl,
      mapping = aes(x = class_1, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.4,
      fontface = "bold",
      vjust = 1
    ) +
    facet_wrap(~facet_group, nrow = facet_nrow) +
    scale_fill_manual(values = traj_palette, drop = FALSE) +
    labs(
      x = "Trajectory",
      y = ylab_text,
      title = title_text
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 13, color = "white"),
      strip.background = element_rect(fill = "black", color = "black"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold"),
      plot.margin = margin(10, 10, 25, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.06))) +
    coord_cartesian(clip = "off")
  
  if (nrow(stat_tbl) > 0) {
    p <- p +
      stat_pvalue_manual(
        stat_tbl,
        label = "p.adj.signif",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.02,
        hide.ns = FALSE,
        inherit.aes = FALSE
      )
  }
  
  p
}

# =========================================================
# Save helper
# =========================================================
save_plot <- function(p, filename, width, height) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = p,
    width = width,
    height = height
  )
}

# =========================================================
# Main run
# =========================================================
for (ds in datasets) {
  
  print( paste0("---------------------------",ds))
  spreadcazy <- load_spreadcazy_rda(ds$path)
  spreadcazy <- standardize_cazy_data(spreadcazy)
  spreadcazy <- spreadcazy %>%
    distinct(genome,visit,.keep_all = T) 
  
  # Single species
  df_single <- spreadcazy[spreadcazy$speciesname == species_single, , drop = FALSE]
  
  
  p_single_GH <- plot_boxplot_panel(
    df = df_single,
    response = "sum_GH",
    facet_col = "transmit",
    dataset_label = ds$label,
    title_text = paste0("Stratified by Transmission Status\n", species_single, " (", ds$label, ")"),
    ylab_text = paste0("Sum of GH gene count (", ds$label, ")"),
    facet_nrow = NULL,
    stat_filename = paste0(ds$tag, "_GH_single_species_stats.csv")
  )
  
  p_single_all <- plot_boxplot_panel(
    df = df_single,
    response = "sum_all",
    facet_col = "transmit",
    dataset_label = ds$label,
    title_text = paste0("Stratified by Transmission Status\n", species_single, " (", ds$label, ")"),
    ylab_text = paste0("Richness of (", ds$label, ")"),
    facet_nrow = NULL,
    stat_filename = paste0(ds$tag, "_Richness_single_species_stats.csv")
  )
  
  # Multi-species transmitted subset
  df_multi <- spreadcazy[
    spreadcazy$speciesname %in% species_multi &
      spreadcazy$transmit == "transmitted between dyads",
    ,
    drop = FALSE
  ]
  

  
  p_multi_GH <- plot_boxplot_panel(
    df = df_multi,
    response = "sum_GH",
    facet_col = "speciesname",
    dataset_label = ds$label,
    title_text = paste0("Transmitted genomes (", ds$label, ")"),
    ylab_text = "Gene count of GH",
    facet_nrow = 1,
    stat_filename = paste0(ds$tag, "_GH_multi_species_stats.csv")
  )
  
  p_multi_all <- plot_boxplot_panel(
    df = df_multi,
    response = "sum_all",
    facet_col = "speciesname",
    dataset_label = ds$label,
    title_text = paste0("Transmitted genomes (", ds$label, ")"),
    ylab_text = "Richessness",
    facet_nrow = 1,
    stat_filename = paste0(ds$tag, "_Richessness_multi_species_stats.csv")
  )


  

  save_plot(
    p_single_GH,
    paste0(ds$tag, "_GH_", gsub(" ", "_", species_single), "_transmit_status.pdf"),
    width = 9, height = 8
  )
  
  save_plot(
    p_single_all,
    paste0(ds$tag, "_all_", gsub(" ", "_", species_single), "_transmit_status.pdf"),
    width = 9, height = 8
  )
  

  
  save_plot(
    p_multi_GH,
    paste0(ds$tag, "_GH_Bifidobacterium_transmitted.pdf"),
    width = 14, height =8
  )
  
  save_plot(
    p_multi_all,
    paste0(ds$tag, "_richness_", gsub(" ", "_", species_single), "_transmit_status.pdf"),
    width = 14, height = 8
  )
}