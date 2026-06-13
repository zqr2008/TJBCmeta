# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 03_differential_abundance_Differential_abundance_trajectory
# Figure(s): Figure3B
# Table(s): No direct submitted table matched from filename
# Purpose: Processes MaAsLin2 differential-abundance results and heatmaps.
# Main input(s): none detected
# Main output(s): none detected
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END


rm(list = ls())

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(conflicted)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::first)


# ============================================================
# 1. Set working directory
# ============================================================

setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")


# ============================================================
# 2. Read MaAsLin2 significant results
# ============================================================

significant_results1 <- read.delim(
  "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/offspring_outputSpecies_comparsion1/significant_results.tsv"
)

significant_results2 <- read.delim(
  "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/offspring_outputSpecies_comparsion2/significant_results.tsv"
)

significant_results3 <- read.delim(
  "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/offspring_outputSpecies_comparsion3/significant_results.tsv"
)


# ============================================================
# 3. Select trajectory-associated species
# ============================================================

significant_results_select <- bind_rows(
  significant_results1,
  significant_results2,
  significant_results3
) %>%
  filter(str_detect(metadata, "comparsion")) %>%
  filter(qval < 0.1)


write.table(
  significant_results_select,
  file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/offspring_speices_significant_results_select.csv",
  fileEncoding = "GBK",
  sep = ",",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)


# ============================================================
# 4. Extract phenotype-related annotation
#    Feeding type and delivery mode
# ============================================================

significant_results_phenotype <- significant_results2 %>%
  filter(
    metadata == "feeding_type_28days" |
      metadata == "delivery_mode"
  ) %>%
  pivot_wider(
    id_cols = feature,
    names_from = value,
    values_from = coef,
    values_fn = first
  ) %>%
  mutate(
    `Enriched in mixed feeding` = case_when(
      `mixed_feeding` > 0 ~ "Yes",
      TRUE ~ NA_character_
    ),
    
    `Enriched in exclusive formula fed` = case_when(
      `exclusive_formula fed` > 0 ~ "Yes",
      TRUE ~ NA_character_
    ),
    
    `Enriched in exclusive breastfeeding` = case_when(
      `mixed_feeding` < 0 | `exclusive_formula fed` < 0 ~ "Yes",
      TRUE ~ NA_character_
    ),
    
    `Delivery mode` = case_when(
      `Vaginal delivery` > 0 ~ "Enriched in Vaginal delivery",
      `Vaginal delivery` < 0 ~ "Enriched in C-section",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    feature,
    `Enriched in mixed feeding`,
    `Enriched in exclusive formula fed`,
    `Enriched in exclusive breastfeeding`,
    `Delivery mode`
  )


# ============================================================
# 5. Build heatmap coefficient matrix
#    Rows are ordered by breastfeeding/formula-feeding enrichment
# ============================================================

significant_results_traj_all <- significant_results_select %>%
  dplyr::select(value, feature, coef) %>%
  pivot_wider(
    id_cols = feature,
    names_from = value,
    values_from = coef
  ) %>%
  left_join(
    significant_results_phenotype,
    by = "feature"
  ) %>%
  mutate(
    row_group = case_when(
      `Enriched in exclusive breastfeeding` == "Yes" ~ 1,
      `Enriched in exclusive formula fed` == "Yes" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(row_group, feature) %>%
  dplyr::select(-row_group)


# Store ordered feature names before converting to rownames
ordered_features <- significant_results_traj_all$feature


# Convert to rownames
significant_results_traj_all <- significant_results_traj_all %>%
  column_to_rownames("feature")


# ============================================================
# 6. Separate row annotations from heatmap matrix
# ============================================================

significant_results_traj_cate <- significant_results_traj_all %>%
  dplyr::select(
    `Enriched in mixed feeding`,
    `Enriched in exclusive formula fed`,
    `Enriched in exclusive breastfeeding`,
    `Delivery mode`
  )


significant_results_traj <- significant_results_traj_all %>%
  dplyr::select(
    -`Enriched in mixed feeding`,
    -`Enriched in exclusive formula fed`,
    -`Enriched in exclusive breastfeeding`,
    -`Delivery mode`
  )


# Replace underscores in species names for display
rownames(significant_results_traj) <- str_replace_all(
  rownames(significant_results_traj),
  "_",
  " "
)

rownames(significant_results_traj_cate) <- rownames(significant_results_traj)


# ============================================================
# 7. Build q-value star matrix
# ============================================================

significant_results_traj_dot <- significant_results_select %>%
  dplyr::select(value, feature, qval) %>%
  pivot_wider(
    id_cols = feature,
    names_from = value,
    values_from = qval
  ) %>%
  column_to_rownames("feature")


rownames(significant_results_traj_dot) <- str_replace_all(
  rownames(significant_results_traj_dot),
  "_",
  " "
)


# Match the heatmap row order exactly
significant_results_traj_dot <- significant_results_traj_dot[
  rownames(significant_results_traj),
]


# Convert q-values to significance stars
significant_results_traj_dot <- as.data.frame(significant_results_traj_dot)

significant_results_traj_dot[] <- lapply(
  significant_results_traj_dot,
  function(x) {
    case_when(
      is.na(x) ~ "",
      x < 0.001 ~ "***",
      x < 0.01 ~ "**",
      x < 0.1 ~ "*",
      TRUE ~ ""
    )
  }
)


# ============================================================
# 8. Rename trajectory columns
# ============================================================

colnames(significant_results_traj) <- recode(
  colnames(significant_results_traj),
  traj1 = "Trajectory 1 vs. non-Trajectory 1",
  traj2 = "Trajectory 2 vs. non-Trajectory 2",
  traj3 = "Trajectory 3 vs. non-Trajectory 3"
)

colnames(significant_results_traj_dot) <- recode(
  colnames(significant_results_traj_dot),
  traj1 = "Trajectory 1 vs. non-Trajectory 1",
  traj2 = "Trajectory 2 vs. non-Trajectory 2",
  traj3 = "Trajectory 3 vs. non-Trajectory 3"
)


# Make sure column order is the same between coefficient matrix and star matrix
significant_results_traj_dot <- significant_results_traj_dot[
  ,
  colnames(significant_results_traj),
  drop = FALSE
]


# ============================================================
# 9. Heatmap colors and row annotation colors
# ============================================================

coul <- colorRampPalette(
  brewer.pal(12, "RdBu")
)(50)


ann_colors <- list(
  `Delivery mode` = c(
    `Enriched in Vaginal delivery` = "#16317d",
    `Enriched in C-section` = "#B0B0B0"
  ),
  
  `Enriched in mixed feeding` = c(
    `Yes` = "#60c6db"
  ),
  
  `Enriched in exclusive formula fed` = c(
    `Yes` = "#911f27"
  ),
  
  `Enriched in exclusive breastfeeding` = c(
    `Yes` = "#006400"
  )
)


# Italicized species names
newnames <- lapply(
  rownames(significant_results_traj),
  function(x) bquote(italic(.(x)))
)


# ============================================================
# 10. Draw heatmap
# ============================================================

pdf(
  file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure3B.pdf",
  width = 8,
  height = 12
)

pheatmap(
  significant_results_traj,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  border_color = "white",
  legend = TRUE,
  color = rev(coul),
  annotation_colors = ann_colors,
  annotation_row = significant_results_traj_cate,
  breaks = seq(-1.0, 1.0, length.out = 50),
  display_numbers = significant_results_traj_dot,
  display_numbers_color = "white",
  fontsize_number = 15,
  cellwidth = 20,
  cellheight = 13,
  labels_row = as.expression(newnames)
)

dev.off()
