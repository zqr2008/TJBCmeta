# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 01_metadata_qc_baseline_Metadata_QC_baseline
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Generates quality-control summaries and QC figure/table outputs.
# Main input(s): newdata_phylo_NMrevision_0305.rds
# Main output(s): sample_count_by_visit_trajectory.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
library(tidyverse)
library(dplyr)
library(ggplot2)

data <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")
sam_data <- as.data.frame(data@sam_data)
class(sam_data) <- "data.frame"

visit_level=c("MG1", "MG2","MG3", "BM1.5", "BM6", "BM12", "BM24", "BM36")
bar_df <- sam_data %>%
  dplyr::select(class_1, visit) %>%
  group_by(class_1, visit) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(visit) %>%
  mutate(Total_visit = sum(Count),
         Percent = Count / Total_visit * 100) %>%
  mutate(visit = factor(visit, levels = visit_level))


forstatsitic <- bar_df %>%
  filter(visit == "BM36" | visit == "BM24")
  

# Make a small helper function
test_one <- function(df, k) {
  sub <- df[df$class_1 == k, ]
  
  x <- sub$Count
  n <- sub$Total_visit
  
  prop.test(x = x, n = n)
}


# Run for each trajectory
test1 <- test_one(forstatsitic, 1)
test2 <- test_one(forstatsitic, 2)
test3 <- test_one(forstatsitic, 3)

test1
test2
test3


fig <- bar_df %>%
  ggplot(aes(x = visit, y = Count, fill = as.factor(class_1))) +
  geom_bar(stat = "identity",
           position = "stack",
           color = "grey",
           width = 0.7,
           size = 0.25) +
  
  # --- percentage labels ---
  geom_text(aes(label = paste0(sprintf("%.1f", Percent), "%")),
            size = 3,
            color = "#4d4545",
            position = position_stack(vjust = 0.5)) +
  #geom_text(aes(label = paste0(Count)),
  #          size = 3,
  #          color = "#4d4545",
  #          position = position_stack(vjust = 0.2))+
  
  scale_fill_manual(name = "Trajectory",
                    values = c("#e95280", "#23b1a5", "#ffdd7e"),
                    labels = c("Trajectory 1", "Trajectory 2", "Trajectory 3")) +
  xlab("Visit") +
  ylab("Sample Size") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text.x = element_text(size = 10)
  )

pdf("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/sample_count_by_visit_trajectory.pdf", 6, 4)
print(fig)
dev.off()
