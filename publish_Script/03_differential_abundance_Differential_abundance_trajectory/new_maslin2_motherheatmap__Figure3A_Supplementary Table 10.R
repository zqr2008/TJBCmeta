# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 03_differential_abundance_Differential_abundance_trajectory
# Figure(s): Figure3A
# Table(s): Supplementary Table 10
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


setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")


# ============================================================
# 1. Read MaAsLin2 significant results
# ============================================================

significant_results1 <- read.delim(
  "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/mother_output_Species_comparsion1/significant_results.tsv"
)

significant_results2 <- read.delim(
  "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/mother_output_Species_comparsion2/significant_results.tsv"
)

significant_results3 <- read.delim(
  "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/mother_output_Species_comparsion3/significant_results.tsv"
)


# ============================================================
# 2. Select trajectory-associated maternal species
# ============================================================

significant_results_select <- bind_rows(
  significant_results1,
  significant_results2,
  significant_results3
) %>%
  distinct(feature,metadata,value,.keep_all = T) %>%
  filter(str_detect(metadata, "comparsion|visit")) %>%

  filter(qval < 0.1)


write.table(
  significant_results_select,
  file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 10 Mother_base_model_species.csv",
  fileEncoding = "GBK",
  sep = ",",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)


# ============================================================
# 3. BMI annotation
# ============================================================
significant_results_select <- significant_results_select %>%
  filter(str_detect(metadata, "comparsion")) %>%
  
  filter(qval < 0.1)
significant_results_BMI <- significant_results2 %>%
  filter(metadata == "BMI_mo") %>%
  filter(qval < 0.1) %>%
  pivot_wider(
    id_cols = feature,
    names_from = value,
    values_from = coef
  ) %>%
  mutate(
    `BMI status` = case_when(
      BMI_mo < 0 ~ "Negatively associated with mother's BMI",
      BMI_mo > 0 ~ "Positively associated with mother's BMI",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(feature, `BMI status`)


# ============================================================
# 4. Trimester annotation
# ============================================================

significant_results_visit <- significant_results2 %>%
  filter(metadata == "visit") %>%
  filter(qval < 0.1) %>%
  pivot_wider(
    id_cols = feature,
    names_from = value,
    values_from = coef
  ) %>%
  mutate(
    `trimester status` = case_when(
      MG2 > 0 & MG3 > 0 ~ "Positively associated with mother's trimesters",
      MG2 < 0 & MG3 < 0 ~ "Negatively associated with mother's trimesters",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(feature, `trimester status`)


# ============================================================
# 5. UA annotation
# ============================================================

significant_results_UA <- significant_results2 %>%
  filter(metadata == "UA.x") %>%
  filter(qval < 0.1) %>%
  pivot_wider(
    id_cols = feature,
    names_from = value,
    values_from = coef
  ) %>%
  mutate(
    `UA status` = case_when(
      UA.x > 0 ~ "Positively associated with mother's UA",
      UA.x < 0 ~ "Negatively associated with mother's UA",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(feature, `UA status`)


# ============================================================
# 6. TG annotation
# ============================================================

significant_results_TG <- significant_results2 %>%
  filter(metadata == "TG") %>%
  filter(qval < 0.1) %>%
  pivot_wider(
    id_cols = feature,
    names_from = value,
    values_from = coef
  ) %>%
  mutate(
    `TG status` = case_when(
      TG > 0 ~ "Positively associated with mother's TG",
      TG < 0 ~ "Negatively associated with mother's TG",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(feature, `TG status`)


# ============================================================
# 7. Build heatmap matrix and row annotations
#    Rows are ordered according to BMI status
# ============================================================

significant_results_traj_all <- significant_results_select %>%
  mutate(trajectory = paste0(value)) %>%
  dplyr::select(trajectory, feature, coef) %>%
  pivot_wider(
    id_cols = feature,
    names_from = trajectory,
    values_from = coef
  ) %>%
  left_join(significant_results_visit, by = "feature") %>%
  left_join(significant_results_BMI, by = "feature") %>%
  left_join(significant_results_UA, by = "feature") %>%
  left_join(significant_results_TG, by = "feature") %>%
  
  # Explicit BMI-based row order
  mutate(
    BMI_order = case_when(
      `BMI status` == "Positively associated with mother's BMI" ~ 1,
      `BMI status` == "Negatively associated with mother's BMI" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(`TG status`,BMI_order,`UA status`,`trimester status`,feature) %>%
  dplyr::select(-BMI_order)


# Store row order before converting to rownames
ordered_features <- significant_results_traj_all$feature


# Row annotation table
significant_results_traj_cate <- significant_results_traj_all %>%
  dplyr::select(feature, ends_with("status")) %>%
  column_to_rownames("feature")


# Heatmap coefficient matrix
significant_results_traj <- significant_results_traj_all %>%
  dplyr::select(feature, traj1, traj2, traj3) %>%
  column_to_rownames("feature")


# Replace underscores for display
rownames(significant_results_traj) <- str_replace_all(
  rownames(significant_results_traj),
  "_",
  " "
)

rownames(significant_results_traj_cate) <- str_replace_all(
  rownames(significant_results_traj_cate),
  "_",
  " "
)


# ============================================================
# 8. Build q-value star matrix
# ============================================================

significant_results_traj_dot <- significant_results_select %>%
  mutate(trajectory = paste0(value)) %>%
  dplyr::select(trajectory, feature, qval) %>%
  pivot_wider(
    id_cols = feature,
    names_from = trajectory,
    values_from = qval
  ) %>%
  column_to_rownames("feature")


rownames(significant_results_traj_dot) <- str_replace_all(
  rownames(significant_results_traj_dot),
  "_",
  " "
)


# Match the exact heatmap row order
significant_results_traj_dot <- significant_results_traj_dot[
  rownames(significant_results_traj),
  ,
  drop = FALSE
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
# 9. Rename columns
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


# Ensure star matrix columns match heatmap matrix columns
significant_results_traj_dot <- significant_results_traj_dot[
  ,
  colnames(significant_results_traj),
  drop = FALSE
]


# ============================================================
# 10. Row labels and colors
# ============================================================

newnames <- lapply(
  rownames(significant_results_traj),
  function(x) bquote(italic(.(x)))
)


coul <- colorRampPalette(
  brewer.pal(12, "RdBu")
)(50)


# Optional but recommended: annotation colors
ann_colors <- list(
  `BMI status` = c(
    `Positively associated with mother's BMI` = "#B22222",
    `Negatively associated with mother's BMI` = "#2166AC"
  ),
  `trimester status` = c(
    `Positively associated with mother's trimesters` = "#E69F00",
    `Negatively associated with mother's trimesters` = "#56B4E9"
  ),
  `UA status` = c(
    `Positively associated with mother's UA` = "#CC79A7",
    `Negatively associated with mother's UA` = "#009E73"
  ),
  `TG status` = c(
    `Positively associated with mother's TG` = "#D55E00",
    `Negatively associated with mother's TG` = "#0072B2"
  )
)


# ============================================================
# 11. Draw heatmap
# ============================================================

pdf(
  file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure3A motherde_species.pdf",
  width = 10,
  height = 16
)

pheatmap(
  significant_results_traj,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  border_color = "white",
  legend = TRUE,
  color = rev(coul),
  breaks = seq(-1.0, 1.0, length.out = 50),
  annotation_row = significant_results_traj_cate,
  annotation_colors = ann_colors,
  display_numbers = significant_results_traj_dot,
  fontsize_number = 20,
  cellwidth = 20,
  cellheight = 20,
  labels_row = as.expression(newnames)
)

dev.off()