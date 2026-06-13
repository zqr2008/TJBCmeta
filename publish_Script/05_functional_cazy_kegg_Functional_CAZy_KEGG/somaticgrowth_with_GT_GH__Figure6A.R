# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): Figure6A
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): spreadcazy_high_qual.rda; newdata_phylo_NMrevision_0305.rds
# Main output(s): Figure6A.pdf; GH_count_by_transmission_status_candidate.pdf
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
library(broom)

conflicted::conflicts_prefer(stats::cor)
conflicted::conflicts_prefer(dplyr::filter)
###########################################################
#clustering heatmap for cazy 
###########################################################
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")
#load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,zwhz_1:zwhz_7) %>%
  mutate(derta = zwhz_2-zwhz_1 )


otu <- as.data.frame(complete_phylo@otu_table)
class(otu) <- "data.frame"
otu  <-  otu %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,Bifidobacterium_pseudocatenulatum,Bifidobacterium_bifidum)




spreadcazyforspecies <- spreadcazy %>%
  distinct(visit,genome,.keep_all = T) %>%
  inner_join(otu,by="sampleid") %>%
  filter(speciesname == "Bifidobacterium pseudocatenulatum") %>%
  inner_join(metadata, by ="sampleid")  %>%
  filter(transmit=="transmitted between dyads") 
  


valid_pairs <- tibble::tribble(
  ~visit,   ~zwhz,
  "BM1.5",  "zwhz_1",
  "BM1.5",  "zwhz_2",
  "BM1.5",  "zwhz_3",
  
  "BM6",    "zwhz_2",
  "BM6",    "zwhz_3",
  "BM6",    "zwhz_4",
  
  "BM12",   "zwhz_3",
  "BM12",   "zwhz_4",
  "BM12",   "zwhz_5",
  
  "BM24",   "zwhz_4",
  "BM24",   "zwhz_5",
  "BM24",   "zwhz_6",
  
  "BM36",   "zwhz_5",
  "BM36",   "zwhz_6",
  "BM36",   "zwhz_7"
)

df <- spreadcazyforspecies 


df_plot <- df %>%
  filter(visit %in% c("BM6"))

# calculate corrected p-values
cor_df <- df_plot %>%
  group_by(visit) %>%
  summarise(
    r = cor(derta, sum_GH, method = "pearson", use = "complete.obs"),
    p = cor.test(derta, sum_GH,
                 method = "pearson")$p.value
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    label = paste0(
      "rho = ", round(r, 3),
      "\nP = ", signif(p, 3)
    )
  )

ggplot(df_plot,
       aes(x = derta,
           y = sum_GH)) +
  
  geom_point(alpha = 0.9, size = 2.5) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    inherit.aes = FALSE,
    aes(x = derta, y = sum_GH),
    color = "black"
  ) +
  geom_text(
    data = cor_df,
    aes(label = label),
    x = Inf,
    y = Inf,
    hjust = 1,
    vjust = 4,
    inherit.aes = FALSE,
    size = 12
  ) +
  xlab("WFL/H z-scores change \n(from 0 month to 6 months)")+ 
  ylab("GH gene counts")+
  theme_classic(base_size = 18)


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure6A.pdf",
       width =6, height = 6)



df <- spreadcazyforspecies 


df_plot <- df %>%
  filter(visit %in% c("BM1.5","BM6"))


n_df <- df_plot %>%
  group_by(GH42) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0("n = ", n))

y_min <- min(df_plot$derta, na.rm = TRUE)
y_max <- max(df_plot$derta, na.rm = TRUE)
y_range <- y_max - y_min

p_y <- y_max + 0.12 * y_range
n_y <- y_min - 0.10 * y_range

my_comparison <- list(levels(df_plot$GH42))

ggplot(df_plot, aes(x = as.factor(GH42), y = zwhz_2)) +
  
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "grey90",
    color = "black"
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 2,
    alpha = 0.55
  ) +
  
  stat_compare_means(
    comparisons = list(c("2","3")),
    method = "wilcox.test",
    label = "p.format",
    label.y = p_y,
    tip.length = 0.03,
    bracket.size = 0.7,
    size = 5
  ) +
  
  geom_text(
    data = n_df,
    aes(
      x = GH42,
      y = n_y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  scale_y_continuous(
    limits = c(n_y - 0.05 * y_range, p_y + 0.12 * y_range)
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "GH42",
    y = "WFL/H z-score change\n(0 to 6 months)"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.margin = margin(20, 20, 20, 20)
  )




ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/GH_count_by_transmission_status_candidate.pdf",
       width =5, height = 6)
