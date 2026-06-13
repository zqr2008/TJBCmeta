# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 03_differential_abundance_Differential_abundance_trajectory
# Figure(s): no direct figure
# Table(s): Supplementary Table 8
# Purpose: Processes MaAsLin2 differential-abundance results and heatmaps.
# Main input(s): significant_results.tsv
# Main output(s): Supplementary Table 8 Mother_base_model_genus.tsv; motherde_genus.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")
significant_results1 <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/mother_output_Genus_comparsion1/significant_results.tsv")
significant_results2 <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/mother_output_Genus_comparsion2/significant_results.tsv")
significant_results3 <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/mother_output_Genus_comparsion3/significant_results.tsv")

significant_results_select <- rbind(significant_results1,significant_results2)
significant_results_select <- rbind(significant_results_select,significant_results3) %>%
  filter(str_detect(metadata,"comparsion")) %>%
  filter(qval<0.1)


write.table(significant_results_select,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 8 Mother_base_model_genus.tsv",
            fileEncoding = "GBK",sep = ",",quote = F,col.names = T,row.names = F)


significant_results_BMI <- significant_results2 %>%
  filter(metadata == "BMI_mo") %>%
  filter(qval < 0.1) %>%
  mutate(trajectory = paste0(metadata,value)) %>%
  pivot_wider(id_cols = feature,
              names_from = value,
              values_from = coef ) %>%
  mutate(`BMI status` = case_when(BMI_mo < 0 ~ "Negatively associated with mother's BMI",
                                  BMI_mo > 0 ~ "Positively associated with mother's BMI"))



significant_results_visit  <- significant_results2 %>%
  filter(metadata == "visit") %>%
  filter(qval < 0.1) %>%
  mutate(trajectory = paste0(metadata,value)) %>%
  pivot_wider(id_cols = feature,
              names_from = value,
              values_from = coef ) %>%
  mutate(`trimester status` = case_when(MG2 > 0 & MG3 > 0 ~ "Positively associated with mother's trimesters",
                                        MG2 < 0 & MG3 < 0 ~ "Negatively associated with mother's trimesters"))




significant_results_UA  <- significant_results2 %>%
  filter(metadata == "UA.x") %>%
  filter(qval < 0.1) %>%
  mutate(trajectory = paste0(metadata,value)) %>%
  pivot_wider(id_cols = feature,
              names_from = value,
              values_from = coef ) %>%
  mutate(`UA status` = case_when(UA.x > 0 ~ "Positively associated with mother's UA",
                                 UA.x < 0 ~ "Negatively associated with mother's UA"))





significant_results_TG  <- significant_results2 %>%
  filter(metadata == "TG") %>%
  filter(qval < 0.1) %>%
  mutate(trajectory = paste0(metadata,value)) %>%
  pivot_wider(id_cols = feature,
              names_from = value,
              values_from = coef ) %>%
  mutate(`TG status` = case_when(TG > 0 ~ "Positively associated with mother's TG",
                                 TG < 0 ~ "Negatively associated with mother's TG"))



significant_results_traj <- significant_results_select %>%
  mutate(trajectory = paste0(value)) %>%
  dplyr::select(trajectory,feature,coef) %>%
  pivot_wider(id_cols = feature,
              names_from = trajectory,
              values_from = coef ) %>%
  left_join(significant_results_visit,by = "feature") %>%
  left_join(significant_results_BMI, by = "feature") %>%
  left_join(significant_results_UA, by = "feature") %>%
  left_join(significant_results_TG, by = "feature")


#significant_results_traj <- significant_results_traj %>%
#  arrange(`BMI status`) %>%


significant_results_traj_cate <- significant_results_traj %>%
  #dplyr::select(`Trimester associated changes`,`BMI status`)
  column_to_rownames("feature") %>%
  dplyr::select(ends_with('status'))

significant_results_traj <- significant_results_traj %>%
  dplyr::select(feature,traj1,traj2,traj3) %>%
  column_to_rownames("feature")

rownames(significant_results_traj) <- str_replace(rownames(significant_results_traj),"_"," ")
rownames(significant_results_traj_cate) <-  str_replace(rownames(significant_results_traj_cate),"_"," ")



significant_results_traj_dot <-  significant_results_select %>%
  mutate(trajectory = paste0(value)) %>%
  dplyr::select(trajectory,feature,qval) %>%
  pivot_wider(id_cols = feature,
              names_from = trajectory,
              values_from = qval ) %>%
  column_to_rownames("feature") 


rownames(significant_results_traj_dot) <- str_replace(rownames(significant_results_traj_dot),"_"," ")

significant_results_traj_dot <- significant_results_traj_dot[rownames(significant_results_traj), ]
significant_results_traj_dot[significant_results_traj_dot < 0.001] <- "***" 
significant_results_traj_dot[significant_results_traj_dot < 0.01 & significant_results_traj_dot >= 0.001] <- "**" 
significant_results_traj_dot[significant_results_traj_dot < 0.1 & significant_results_traj_dot >= 0.01] <- "*" 
significant_results_traj_dot[is.na(significant_results_traj_dot)] <- ""



newnames <- lapply(
  rownames(significant_results_traj),
  function(x) bquote(italic(.(x))))


colnames(significant_results_traj)[colnames(significant_results_traj) == "traj1"] <- "Trajectory 1 vs. non-Trajectory 1"
colnames(significant_results_traj)[colnames(significant_results_traj) == "traj2"] <- "Trajectory 2 vs. non-Trajectory 2"
colnames(significant_results_traj)[colnames(significant_results_traj) == "traj3"] <- "Trajectory 3 vs. non-Trajectory 3"

colnames(significant_results_traj_dot)[colnames(significant_results_traj_dot) == "traj1"] <- "Trajectory 1 vs. non-Trajectory 1"
colnames(significant_results_traj_dot)[colnames(significant_results_traj_dot) == "traj2"] <- "Trajectory 2 vs. non-Trajectory 2"
colnames(significant_results_traj_dot)[colnames(significant_results_traj_dot) == "traj3"] <- "Trajectory 3 vs. non-Trajectory 3"





coul <- colorRampPalette(brewer.pal(12, "RdBu"))(50)

pdf(file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/motherde_genus.pdf",width = 8, height = 16)
pheatmap(significant_results_traj,
         cluster_cols = F,
         cluster_rows = F,
         border_color = "white",
         legend = TRUE,
         color = rev(coul),
         breaks = seq(-1.0, 1.0, length.out = 50), 
         annotation_row = significant_results_traj_cate,
         display_numbers = significant_results_traj_dot,
         fontsize_number=20, 
         cellwidth = 20,
         cellheight = 20,
         labels_row = as.expression(newnames))
dev.off()
