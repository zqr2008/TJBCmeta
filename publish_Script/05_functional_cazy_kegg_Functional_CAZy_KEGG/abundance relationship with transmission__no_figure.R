# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Analyzes maternal-infant strain transmission patterns and related outputs.
# Main input(s): spreadcazy.rda; Bifidobacterium.rda; ps_filtered.rds
# Main output(s): none detected
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggpubr)

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy.rda")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/Bifidobacterium.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/ps_filtered.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)

Bifidobacterium_pseudocatenulatum_list <- Bifidobacterium %>%
  filter(species == "Bifidobacterium pseudocatenulatum")


Bifidobacterium_bifidum_list <- Bifidobacterium %>%
  filter(species == "Bifidobacterium bifidum")


Bifidobacterium_adolescentis_list <- Bifidobacterium %>%
  filter(species == "Bifidobacterium adolescentis")


otu <- as.data.frame(complete_phylo@otu_table)
class(otu) <- "data.frame"
otu  <-  otu %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,Bifidobacterium_pseudocatenulatum,Bifidobacterium_bifidum) %>%
  mutate(transmit_pseudocatenulatum = case_when(sampleid %in% Bifidobacterium_pseudocatenulatum_list$Samp1 ~ "Transmitted",
                                                sampleid %in% Bifidobacterium_pseudocatenulatum_list$Samp2 ~ "Transmitted",
                                                TRUE ~ "No")) %>%
  mutate(transmit_bifidum = case_when(sampleid %in% Bifidobacterium_bifidum_list$Samp1 ~ "Transmitted",
                                                sampleid %in% Bifidobacterium_bifidum_list$Samp2 ~ "Transmitted",
                                                TRUE ~ "No")) %>%
  mutate(transmit_adolescentis = case_when(sampleid %in% Bifidobacterium_adolescentis_list$Samp1 ~ "Transmitted",
                                      sampleid %in% Bifidobacterium_adolescentis_list$Samp2 ~ "Transmitted",
                                      TRUE ~ "No")) %>%
  left_join(metadata,by="sampleid")



# Run Wilcoxon test
stat_test <- otu %>%

  group_by(visit) %>%
  wilcox_test(Bifidobacterium_pseudocatenulatum ~ transmit_pseudocatenulatum, detailed = TRUE, p.adjust.method = "BH") %>%
  add_xy_position(x= "visit")# adjust y position for p-values





# Add significance category
volcano_data <- stat_test %>%
  mutate(
    neglog10p = -log10(p),
    signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Volcano plot
ggplot(volcano_data, aes(x = -estimate, y = neglog10p, color = signif)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = visit), size = 4, max.overlaps = 100) +
  scale_color_manual(
    values = c("***" = "#d73027", "**" = "#fc8d59", "*" = "#fee090", "ns" = "grey70"),
    name = "Significance"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40") +
  theme_classic(base_size = 14) +
  labs(
    title = expression(italic("Bifidobacterium pseudocatenulatum")~" Transmission Effect by Visit"),
    x = "Effect size (Wilcoxon estimate)",
    y = expression(-log[10](p~value))
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )






# Run Wilcoxon test
stat_test <- otu %>%
  
  group_by(visit) %>%
  wilcox_test(Bifidobacterium_bifidum ~ transmit_bifidum, detailed = TRUE, p.adjust.method = "BH") %>%
  add_xy_position(x= "visit")# adjust y position for p-values





# Add significance category
volcano_data <- stat_test %>%
  mutate(
    neglog10p = -log10(p),
    signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Volcano plot
ggplot(volcano_data, aes(x = -estimate, y = neglog10p, color = signif)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = visit), size = 4, max.overlaps = 100) +
  scale_color_manual(
    values = c("***" = "#d73027", "**" = "#fc8d59", "*" = "#fee090", "ns" = "grey70"),
    name = "Significance"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40") +
  theme_classic(base_size = 14) +
  labs(
    title = expression(italic("Bifidobacterium bifidum")~" Transmission Effect by Visit"),
    x = "Effect size (Wilcoxon estimate)",
    y = expression(-log[10](p~value))
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )