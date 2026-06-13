# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): Figure5A; Figure5B
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): spreadcazy.rda; newdata_phylo_NMrevision_0305.rds; spreadcazy_high_qual.rda
# Main output(s): medium_GT_Bifidobacterium_pseudocatenulatum_transmit_status.pdf; high_GT_Bifidobacterium_pseudocatenulatum_transmit_status.pdf; high_GH_Bifidobacterium_pseudocatenulatum_transmit_status.pdf; high_GT_Bifidobacterium_transmitted.pdf; plus 4 more
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
library(pheatmap)
library(patchwork)
###########################################################
#clustering heatmap for cazy 
###########################################################
#load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit,zwhz_1:zwhz_7)


otu <- as.data.frame(complete_phylo@otu_table)
class(otu) <- "data.frame"
otu  <-  otu %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,Bifidobacterium_pseudocatenulatum,Bifidobacterium_longum,Bifidobacterium_bifidum,Bifidobacterium_breve)



###############################################################################
#starting!!!
#comparsion of Bifidobacterium pseudocatenulatum transmit vs. non-transmit
###############################################################################

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")


Bifidobacterium_pseudocatenulatum_spreadcazy <- spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum")

Bifidobacterium_pseudocatenulatum_highGH <- Bifidobacterium_pseudocatenulatum_spreadcazy %>%
  filter(sum_GH>=26)

Bifidobacterium_pseudocatenulatum_highGHotu <- otu %>%
  filter(sampleid %in% Bifidobacterium_pseudocatenulatum_spreadcazy$sampleid) %>%
  mutate(Category = case_when(sampleid %in% Bifidobacterium_pseudocatenulatum_highGH$sampleid~"high GH (>=26)",
                           TRUE ~"low GH (<26)"))

Bifidobacterium_pseudocatenulatum_highGHotu <- Bifidobacterium_pseudocatenulatum_highGHotu %>%
  group_by(Category) %>%
  mutate(Category_Label = paste0(Category, "\n(n=", n(), ")"))

Bifidobacterium_pseudocatenulatum_highGT <- Bifidobacterium_pseudocatenulatum_spreadcazy %>%
  filter(sum_GT>0)

Bifidobacterium_pseudocatenulatum_highGTotu <- otu %>%
  filter(sampleid %in% Bifidobacterium_pseudocatenulatum_spreadcazy$sampleid) %>%
  mutate(Category = case_when(sampleid %in% Bifidobacterium_pseudocatenulatum_highGT$sampleid~"high GT (>0)",
                              TRUE ~"low GT (<=0)"))

Bifidobacterium_pseudocatenulatum_highGTotu <- Bifidobacterium_pseudocatenulatum_highGTotu %>%
  group_by(Category) %>%
  mutate(Category_Label = paste0(Category, "\n(n=", n(), ")"))

# Updated pic1 (removed title, changed x to Category_Label)
p1 <- ggplot(Bifidobacterium_pseudocatenulatum_highGHotu, 
             aes(x = Category_Label, y = log10(Bifidobacterium_pseudocatenulatum + 1e-6), fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_classic(base_size = 18) +
  scale_fill_manual(values = c("#6173f4", "#9d55ad")) +
  labs(x = "GH Count Category", y = expression(log[10] ~ "Relative Abundance")) +
  theme(legend.position = "none") # Hide individual legend

# Updated pic2 (removed title, changed x to Category_Label)
p2 <- ggplot(Bifidobacterium_pseudocatenulatum_highGTotu, 
             aes(x = Category_Label, y = log10(Bifidobacterium_pseudocatenulatum + 1e-6), fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_classic(base_size = 18) +
  scale_fill_manual(values = c("#6173f4", "#9d55ad")) +
  labs(x = "GT Count Category", y = NULL) + # Remove Y label for second plot to save space
  theme(legend.position = "none")



# Combine with Patchwork
p1 + p2 + 
  plot_annotation(
    title = expression(italic("Bifidobacterium pseudocatenulatum")),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 20, vjust = 1))
  )





#############################################################################


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy.rda")
#load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")

# Step 1: Wilcoxon test
stat_test <- spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit) %>%
  wilcox_test(sum_all~ class_1, detailed = TRUE, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# Step 2: Compute max y per facet
y_pos <- spreadcazy %>%
  group_by(transmit) %>%
  summarise(y.position = max(sum_GT, na.rm = TRUE) * 1.05)

# Step 3: Merge y.position into stat_test
stat_test <- stat_test %>%
  left_join(y_pos, by = "transmit")

# Step 4: Create a column `class_1` for ggpubr to find
# We'll set it to group1 for each comparison
stat_test <- stat_test %>%
  mutate(class_1 = group1)


# Add small vertical dodge to y.position to avoid overlap
stat_test <- stat_test %>%
  group_by(transmit) %>%
  arrange(group1, group2) %>%  # optional: ensure consistent order
  mutate(y.position = y.position + (row_number() - 1) * 0.07 * max(y.position)) %>%
  ungroup()

spreadcazy$class_1 <- as.factor(spreadcazy$class_1)

median_df <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit, class_1) %>%
  summarise(median_val = median(sum_GT, na.rm = TRUE), .groups = "drop")


n_df <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit, class_1) %>%
  summarise(n = n(), .groups = "drop")

median_n_df <- median_df %>%
  left_join(n_df, by =c("transmit","class_1"))

spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
  ggplot(aes(x = class_1, y = sum_GT, fill = class_1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(aes(color = class_1), width = 0.2,height = 0, alpha = 0.7, size = 1.5) +
  geom_text(
    data = median_n_df,
    aes(x = class_1, y =median_val, label =paste0("median=",round(median_val, 2),"\n","n=",n)),
    vjust = 3,
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  )+
  facet_wrap(~transmit) +
  scale_fill_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  scale_color_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  labs(
    x = "Trajectory",
    y = "Sum of GT gene count (high-quality genomes)",
    title = "Stratified by Transmission Status\nBifidobacterium pseudocatenulatum") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  theme (plot.title = element_text (size = 16, face = "bold" ))+
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.02,
    hide.ns = FALSE,
    y.position = "y.position"
  )+
  scale_y_continuous(
    limits = c(-3, NA),
    expand = expansion(mult = c(0, 0.05))
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/medium_GT_Bifidobacterium_pseudocatenulatum_transmit_status.pdf",
       width = 6, height = 6)





####################################################################
rm(list = ls())

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")



# Step 1: Wilcoxon test
stat_test <- spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(transmit) %>%
  wilcox_test(sum_GH ~ class_1, detailed = TRUE, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# Step 2: Compute max y per facet
y_pos <- spreadcazy %>%
  group_by(transmit) %>%
  summarise(y.position = max(sum_GT, na.rm = TRUE) * 1.05)

# Step 3: Merge y.position into stat_test
stat_test <- stat_test %>%
  left_join(y_pos, by = "transmit")

# Step 4: Create a column `class_1` for ggpubr to find
# We'll set it to group1 for each comparison
stat_test <- stat_test %>%
  mutate(class_1 = group1)


# Add small vertical dodge to y.position to avoid overlap
stat_test <- stat_test %>%
  group_by(transmit) %>%
  arrange(group1, group2) %>%  # optional: ensure consistent order
  mutate(y.position = y.position + (row_number() - 1) * 0.07 * max(y.position)) %>%
  ungroup()



median_df <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit, class_1) %>%
  summarise(median_val = median(sum_GT, na.rm = TRUE), .groups = "drop")


n_df <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit, class_1) %>%
  summarise(n = n(), .groups = "drop")

median_n_df <- median_df %>%
  left_join(n_df, by =c("transmit","class_1"))

spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
  ggplot(aes(x = as.factor(class_1), y = sum_GT, fill = as.factor(class_1))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~transmit) +
  geom_text(
    data = median_n_df,
    aes(x = class_1, y = median_val, label =paste0("median=",round(median_val, 2),"\n","n=",n)),
    vjust = 2,
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  scale_color_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  labs(
    x = "Trajectory",
    y = "Sum of GT gene count (high-quality genomes)",
    title = "Stratified by Transmission Status\nBifidobacterium pseudocatenulatum"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  theme (plot.title = element_text (size = 16, face = "bold" ))+
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.02,
    hide.ns = FALSE,
    y.position = "y.position"
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/high_GT_Bifidobacterium_pseudocatenulatum_transmit_status.pdf",
       width = 6, height = 6)


# Step 1: Wilcoxon test
stat_test <- spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit) %>%
  wilcox_test(sum_GH ~ class_1, detailed = TRUE, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# Step 2: Compute max y per facet
y_pos <- spreadcazy %>%
  group_by(transmit) %>%
  summarise(y.position = max(sum_GH, na.rm = TRUE) * 1.05)

# Step 3: Merge y.position into stat_test
stat_test <- stat_test %>%
  left_join(y_pos, by = "transmit")

# Step 4: Create a column `class_1` for ggpubr to find
# We'll set it to group1 for each comparison
stat_test <- stat_test %>%
  mutate(class_1 = group1)


# Add small vertical dodge to y.position to avoid overlap
stat_test <- stat_test %>%
  group_by(transmit) %>%
  arrange(group1, group2) %>%  # optional: ensure consistent order
  mutate(y.position = y.position + (row_number() - 1) * 0.07 * max(y.position)) %>%
  ungroup()



median_df <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit, class_1) %>%
  summarise(median_val = median(sum_GH, na.rm = TRUE), .groups = "drop")


n_df <- spreadcazy %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  group_by(transmit, class_1) %>%
  summarise(n = n(), .groups = "drop")

median_n_df <- median_df %>%
  left_join(n_df, by =c("transmit","class_1"))

spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
  ggplot(aes(x = as.factor(class_1), y = sum_GH, fill = as.factor(class_1))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~transmit) +
  #geom_jitter(aes(color = class_1), width = 0.2, alpha = 0.7, size = 1.5) +
  geom_text(
    data = median_n_df,
    aes(x = class_1, y = median_val, label =paste0("median=",round(median_val, 2),"\n","n=",n)),
    vjust = 3,
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  scale_color_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  labs(
    x = "Trajectory",
    y = "Sum of GH gene count (high-quality genomes)",
    title = "Stratified by Transmission Status\nBifidobacterium pseudocatenulatum"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  theme (plot.title = element_text (size = 16, face = "bold" ))+
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.02,
    hide.ns = FALSE,
    y.position = "y.position"
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/high_GH_Bifidobacterium_pseudocatenulatum_transmit_status.pdf",
       width = 6, height = 6)


###################################################################
#Bifidobacterium 4 specis transmittted
#
###################################################################
# Ensure class_1 is factor

# Wilcoxon test per visit
sum_GT_stat <- spreadcazy %>%  
  filter(speciesname =="Bifidobacterium longum" | speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>% 
  group_by(speciesname) %>%
  filter(transmit == "transmitted between dyads") %>%
  #group_by(visit) %>%
  wilcox_test(sum_GT ~ class_1, detailed = TRUE, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup() %>%
  mutate(
    y.position = max(spreadcazy$sum_GT, na.rm = TRUE) * 1.05,
    group1 = factor(group1, levels = levels(spreadcazy$class_1)),
    group2 = factor(group2, levels = levels(spreadcazy$class_1))
  ) %>%
  add_xy_position(x = "class_1")  



median_df <- spreadcazy %>%
  filter(speciesname =="Bifidobacterium longum" |speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(speciesname, class_1) %>%
  summarise(median_val = median(sum_GT, na.rm = TRUE), .groups = "drop")


n_df <- spreadcazy %>%
  filter(speciesname =="Bifidobacterium longum" |speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(speciesname, class_1) %>%
  summarise(n = n(), .groups = "drop")

median_n_df <- median_df %>%
  left_join(n_df, by =c("speciesname","class_1"))



spreadcazy %>%
  filter(speciesname =="Bifidobacterium longum" | speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>%
  filter(transmit == "transmitted between dyads") %>%
  ggplot(aes(x = class_1, y = sum_GT, fill = class_1)) +
  facet_wrap(~speciesname,nrow = 1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_text(
    data = median_n_df,
    aes(x = class_1, y = median_val, label =paste0("median=",round(median_val, 2),"\n","n=",n)),
    vjust = 3,
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  )+  
  facet_wrap(~speciesname,nrow = 1)+
  scale_fill_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  scale_color_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  stat_pvalue_manual(
    sum_GT_stat,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    hide.ns = FALSE,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Trajectory",
    y = "Gene count of GT",
    title = "Transmitted genomes (high quality)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )+
  scale_y_continuous(
    limits = c(-2, NA),
    expand = expansion(mult = c(0, 0.05))
  )



ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/high_GT_Bifidobacterium_transmitted.pdf",
       width = 12, height = 6)



# Wilcoxon test per visit
sum_GH_stat <- spreadcazy %>%  
  filter(speciesname =="Bifidobacterium longum" | speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>%
  group_by(speciesname) %>%
  filter(transmit == "transmitted between dyads") %>%
  #group_by(visit) %>%
  wilcox_test(sum_GH ~ class_1, detailed = TRUE, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup() %>%
  mutate(
    y.position = max(spreadcazy$sum_GH, na.rm = TRUE) * 1.05,
    group1 = factor(group1, levels = levels(spreadcazy$class_1)),
    group2 = factor(group2, levels = levels(spreadcazy$class_1))
  ) %>%
  add_xy_position(x = "class_1")  




median_df <- spreadcazy %>%
  filter(speciesname =="Bifidobacterium longum" | speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(speciesname, class_1) %>%
  summarise(median_val = median(sum_GH, na.rm = TRUE), .groups = "drop")


n_df <- spreadcazy %>%
  filter(speciesname =="Bifidobacterium longum" |speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(speciesname, class_1) %>%
  summarise(n = n(), .groups = "drop")

median_n_df <- median_df %>%
  left_join(n_df, by =c("speciesname","class_1"))


spreadcazy %>%
  filter(speciesname =="Bifidobacterium longum" | speciesname=="Bifidobacterium pseudocatenulatum" | speciesname=="Bifidobacterium bifidum" | speciesname =="Bifidobacterium adolescentis")%>%
  filter(transmit == "transmitted between dyads") %>%
  ggplot(aes(x = class_1, y = sum_GH, fill = class_1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~speciesname,nrow = 1)+
  geom_text(
    data = median_n_df,
    aes(x = class_1, y = median_val, label =paste0("median=",round(median_val, 2),"\n","n=",n)),
    vjust = 3,
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  scale_color_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  stat_pvalue_manual(
    sum_GH_stat,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    hide.ns = FALSE,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Trajectory",
    y = "Gene count of GH",
    title = "Transmitted genomes (high quality)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )+
  scale_y_continuous(
    limits = c(-2, NA),
    expand = expansion(mult = c(0, 0.05))
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/high_GH_Bifidobacterium_transmitted.pdf",
       width =12, height = 6)







###############################################
#time course 
###############################################

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")


spreadcazy <- spreadcazy %>%
  distinct(genome,visit,.keep_all = T)

lastpic <- spreadcazy%>%   
  filter(speciesname=="Bifidobacterium pseudocatenulatum" )%>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(visit) %>%
  wilcox_test(sum_GH~ class_1, detailed = TRUE, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup() %>%
  mutate(
    y.position = max(spreadcazy$sum_GH, na.rm = TRUE) * 1.05
  ) %>%
  add_xy_position(x = "class_1")  




median_df <- spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum" )%>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(visit, class_1) %>%
  summarise(median_val = median(sum_GH, na.rm = TRUE), .groups = "drop")

n_df <- spreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum")%>%
  filter(transmit == "transmitted between dyads") %>%
  group_by(visit, class_1) %>%
  summarise(n = n(), .groups = "drop")

median_n_df <- median_df %>%
  left_join(n_df, by =c("visit","class_1"))


sepspreadcazy<- spreadcazy

sepspreadcazy$class_1 <- factor(sepspreadcazy$class_1)

sepspreadcazy %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum")%>%
  filter(transmit == "transmitted between dyads") %>%
  ggplot(aes(x = class_1, y = sum_GH, fill = class_1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_text(
    data = median_n_df,
    aes(x = class_1, y = median_val, label =paste0("median=",round(median_val, 2),"\n","n=",n)),
    vjust = 3,
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  )+
  facet_wrap(~visit)+
  scale_fill_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  scale_color_manual(values = c("1" = "#e95280", "2" = "#23b1a5", "3" = "#E49B0F")) +
  stat_pvalue_manual(
    lastpic,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    hide.ns = FALSE,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Trajectory",
    y = "Gene count of GH in Bifidobacterium pseudocatenulatum(only high quality)",
    title = "Transmitted genomes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )



ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/timecourseGH.pdf",
       width = 8, height =14)



###########################################################################
#for rare Bifidobacterium
###########################################################################
rm(list = ls())
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")


spreadcazy <- spreadcazy %>%
  distinct(genome,visit,.keep_all = T)
  
  
  
countrare <- spreadcazy %>%
  #filter(speciesname != "Bifidobacterium longum" & speciesname!="Bifidobacterium pseudocatenulatum" & speciesname!="Bifidobacterium bifidum" & speciesname !="Bifidobacterium adolescentis")%>%
  group_by(transmit,speciesname)  %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(speciesname) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  mutate(speciesname = factor(speciesname, levels = unique(speciesname)))


countrare2 <- spreadcazy %>%
  #filter(speciesname != "Bifidobacterium longum" & speciesname!="Bifidobacterium pseudocatenulatum" & speciesname!="Bifidobacterium bifidum" & speciesname !="Bifidobacterium adolescentis")%>%
  filter(transmit=="transmitted between dyads") %>%
  group_by(visit,speciesname)  %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(speciesname) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  mutate(speciesname = factor(speciesname, levels = unique(speciesname)))

ggplot(countrare, aes(x = speciesname, y = n, fill = transmit)) +
  geom_col(
    position = "stack",
    color = "grey30",
    linewidth = 0.1,
    width = 0.85
  ) +
  
  geom_text(
    aes(label = ifelse(n >= 3, n, "")),
    position = position_stack(vjust = 0.5),
    size = 4.5,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "transmitted between dyads" = "#D73027",
      "No" = "white"
    )
  ) +
  
  coord_flip() +
  
  labs(
    x = NULL,
    y = "Count",
    fill = NULL
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 18, face = "italic"),
    axis.text.x = element_text(size = 18, face = "bold")
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure5A transtedBifidobacterium.pdf",width =14, height = 6)

# Plot
ggplot(countrare2, aes(x = speciesname, y = n, fill = visit)) +
  geom_col(
    position = "stack",
    color = "grey30",
    linewidth = 0.1,
    width = 0.85
  ) +
  
  # labels for each segment
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  scale_fill_manual(values = c("#e1f6f4","#7db9b3","#166678",
                               "#8B8000","#ffd2a5","#ffa8b8",
                               "#d988bc","#66429b"))  +
  
  coord_flip() +
  labs(
    x = NULL,
    y = "Count",
    fill = NULL
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 18,face = "italic"),
    axis.text.x = element_text(size = 18,face = "bold")
  )


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure5B transtedBifidobacteriumbyvisit.pdf",width =13, height = 6)


