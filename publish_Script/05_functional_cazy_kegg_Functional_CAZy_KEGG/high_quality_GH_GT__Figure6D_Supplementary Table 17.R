# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): Figure6D
# Table(s): Supplementary Table 17
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): spreadcazy_high_qual.rda; newdata_phylo_NMrevision_0305.rds
# Main output(s): Supplementary Table 17 highquality_GHstat_summary.csv; Figure6D cazy_heatmap_t.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)
library(rstatix)


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")

complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")



metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)


spreadcazy <- spreadcazy %>%
  distinct(genome,visit,.keep_all = T)

spreadcazy$class_1 <- as.factor(spreadcazy$class_1)


genenames <- c(colnames(spreadcazy)[str_detect(colnames(spreadcazy),"GH")])

genenames <- genenames[str_detect(genenames,"sum")==FALSE]

results_list <- list()

for (gg in genenames) {
  
  test_data <- spreadcazy %>%
    filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
    filter(transmit == "transmitted between dyads") %>%
    select(class_1, all_of(gg)) %>%
    filter(!is.na(.data[[gg]])) 
  
  # ---- safety checks ----
  if (!is.numeric(test_data[[gg]])) next
  if (n_distinct(test_data$class_1) < 2) next
  if (any(table(test_data$class_1) == 0)) next
  if (nrow(test_data) == 0) next
  
  statcompare <- tryCatch(
    wilcox_test(test_data,
                reformulate("class_1", gg),
                detailed = TRUE,
                p.adjust.method = "BH"),
    error = function(e) NULL
  )
  
  if (is.null(statcompare)) next
  
  # ---- compute medians ----
  med_table <- test_data %>%
    group_by(class_1) %>%
    summarise(med = mean(.data[[gg]], na.rm = TRUE),
              .groups = "drop")
  
  # ---- join medians to pairwise results ----
  statcompare <- statcompare %>%
    left_join(med_table, by = c("group1" = "class_1")) %>%
    rename(med1 = med) %>%
    left_join(med_table, by = c("group2" = "class_1")) %>%
    rename(med2 = med) %>%
    mutate(
      log2FoldChange = ifelse(med1 > 0 & med2 > 0,
                              log2(med2 / med1),
                              NA_real_),
      gene = gg
    )
  
  results_list[[gg]] <- statcompare
}

final_results <- bind_rows(results_list) 


write.table(final_results, file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 17 highquality_GHstat_summary.csv",
            fileEncoding = "GBK",col.names = T,row.names = F,
            sep = ",",quote = F)


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")


gh_vars <- c("GH2_18","GH1","GH172",
             "GH43_22","GH172","GH13_3","GH13_44")





heatmap_data <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  filter(transmit == "transmitted between dyads") %>%
  distinct(genome,visit,.keep_all = T) %>%
  dplyr::select(genome,pair_visit,transmit,visit,class_1, unique(xx$gene)) %>%
  mutate(across(unique(xx$gene), ~ as.numeric(. > 0))) %>%  # ensure presence/absence
  pivot_longer(
               cols = all_of(unique(xx$gene)),
               names_to = "GH",
               values_to = "presence") 




heatmap_n <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  filter(transmit == "transmitted between dyads") %>%
  distinct(genome,visit,.keep_all = T)  %>%
  group_by(class_1) %>%
  summarise(n=n()) 


heatmap_summary <- heatmap_data %>%
  group_by(class_1,GH) %>%
  mutate(n=n()) %>%
  summarise( 
    n = sum(presence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(heatmap_n,by="class_1")%>%
  mutate(prevalences= n.x/n.y)

heatmap_summary %>%
  mutate(Trajectory = paste0("Trajectory ",class_1)) %>%
ggplot(aes(x = GH, y = Trajectory, fill = prevalences)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", prevalences)), size = 4) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#C44E52"
  )+
  labs(
    x = "CAZy",
    y = "Trajectory",
    fill = "Prevalence"
  ) +
  theme_minimal(base_size = 16)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure6D cazy_heatmap_t.pdf",width = 6, height = 4)
