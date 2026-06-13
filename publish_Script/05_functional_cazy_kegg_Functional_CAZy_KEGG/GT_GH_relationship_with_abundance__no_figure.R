# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): spreadcazy_high_qual.rda; Bifidobacterium.rda; newdata_phylo_NMrevision_0305.rds
# Main output(s): correlationheatmap.pdf
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
library(rstatix)
library(lme4)
library(lmerTest)

#load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy.rda")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/Bifidobacterium.rda")

complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,familyid, delivery_mode,feeding_type_28days)



spreadcazy <- spreadcazy %>%
  filter(transmit=="transmitted between dyads")
  #filter(transmit!="transmitted between dyads")

otu <- as.data.frame(complete_phylo@otu_table)
class(otu) <- "data.frame"
otu  <-  otu %>%
  rownames_to_column("sampleid") %>%
  #dplyr::select(sampleid,Bifidobacterium_pseudocatenulatum,Bifidobacterium_longum,Bifidobacterium_bifidum,Bifidobacterium_adolescentis) %>%
  left_join(metadata,by="sampleid") %>%
  inner_join(spreadcazy,by ="sampleid") %>%
  filter(is.na(visit)==FALSE) 



Bifidobacterium_pseudocatenulatum_relationship <- otu %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") 



Bifidobacterium_bifidum_relationship <- otu %>%
  filter(speciesname=="Bifidobacterium bifidum")


Bifidobacterium_adolescentis_relationship <- otu %>%
  filter(speciesname=="Bifidobacterium adolescentis")


Bifidobacterium_longum_relationship <- otu %>%
  filter(speciesname=="Bifidobacterium longum")

# Run correlations
cor_results <- tibble::tibble(
  Comparison = c("GT vs B. pseudocatenulatum", 
                 "GH vs B. pseudocatenulatum", 
                 "GT vs B. bifidum", 
                 "GH vs B. bifidum",
                 "GT vs B. adolescentis", 
                 "GH vs B. adolescentis",
                 "GT vs B. longum", 
                 "GH vs B. longum"),
  r = c(
    cor(Bifidobacterium_pseudocatenulatum_relationship$sum_GT, Bifidobacterium_pseudocatenulatum_relationship$Bifidobacterium_pseudocatenulatum, use="complete.obs", method = "spearman"),
    cor(Bifidobacterium_pseudocatenulatum_relationship$sum_GH, Bifidobacterium_pseudocatenulatum_relationship$Bifidobacterium_pseudocatenulatum, use="complete.obs", method = "spearman"),
    cor(Bifidobacterium_bifidum_relationship$sum_GT, Bifidobacterium_bifidum_relationship$Bifidobacterium_bifidum, use="complete.obs"),
    cor(Bifidobacterium_bifidum_relationship$sum_GH, Bifidobacterium_bifidum_relationship$Bifidobacterium_bifidum, use="complete.obs"),
    cor(Bifidobacterium_adolescentis_relationship$sum_GT, Bifidobacterium_adolescentis_relationship$Bifidobacterium_adolescentis, use="complete.obs"),
    cor(Bifidobacterium_adolescentis_relationship$sum_GH, Bifidobacterium_adolescentis_relationship$Bifidobacterium_adolescentis, use="complete.obs"),
    cor(Bifidobacterium_longum_relationship$sum_GT, Bifidobacterium_longum_relationship$Bifidobacterium_longum, use="complete.obs"),
    cor(Bifidobacterium_longum_relationship$sum_GH, Bifidobacterium_longum_relationship$Bifidobacterium_longum, use="complete.obs")
  ),
  p = c(
    cor.test(Bifidobacterium_pseudocatenulatum_relationship$sum_GT, Bifidobacterium_pseudocatenulatum_relationship$Bifidobacterium_pseudocatenulatum, method = "spearman")$p.value,
    cor.test(Bifidobacterium_pseudocatenulatum_relationship$sum_GH, Bifidobacterium_pseudocatenulatum_relationship$Bifidobacterium_pseudocatenulatum, method = "spearman")$p.value,
    cor.test(Bifidobacterium_bifidum_relationship$sum_GT, Bifidobacterium_bifidum_relationship$Bifidobacterium_bifidum)$p.value,
    cor.test(Bifidobacterium_bifidum_relationship$sum_GH, Bifidobacterium_bifidum_relationship$Bifidobacterium_bifidum)$p.value,
    cor.test(Bifidobacterium_adolescentis_relationship$sum_GT, Bifidobacterium_adolescentis_relationship$Bifidobacterium_adolescentis)$p.value,
    cor.test(Bifidobacterium_adolescentis_relationship$sum_GH, Bifidobacterium_adolescentis_relationship$Bifidobacterium_adolescentis)$p.value,
    cor.test(Bifidobacterium_longum_relationship$sum_GT, Bifidobacterium_longum_relationship$Bifidobacterium_longum)$p.value,
    cor.test(Bifidobacterium_longum_relationship$sum_GH, Bifidobacterium_longum_relationship$Bifidobacterium_longum)$p.value
  )
)


cor_results2 <- cor_results %>%
  mutate(
    CAZy = str_extract(Comparison, "GT|GH"),
    Species = str_extract(Comparison, "B\\. pseudocatenulatum|B\\. bifidum|B\\. adolescentis|B\\. longum"),
    signif_label = paste0("(",sprintf("p = %.3e", p),")"),
    Species = fct_relevel(Species, "B. pseudocatenulatum", "B. bifidum", "B. adolescentis","B. longum")
  )

# 2x2 heatmap
# Create named vector for italicized labels
species_labels <- c(
  "B. pseudocatenulatum" = expression(italic("B. pseudocatenulatum")),
  "B. bifidum" = expression(italic("B. bifidum")),
  "B. adolescentis" = expression(italic("B. adolescentis")),
  "B. longum" = expression(italic("B. longum"))
)

# 2x2 heatmap with italic species names
ggplot(cor_results2, aes(x = CAZy, y = Species, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(r, 3), signif_label)), size = 3) +
  scale_fill_gradient2(low = "#6D9EC1", mid = "white", high = "#E46726", midpoint = 0) +
  scale_y_discrete(labels = species_labels) +
  theme_minimal(base_size = 13) +
  labs(
    x = "CAZy gene count type (high quality)",
    y = "Species",
    fill = "Spearman r"
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank()
  )


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/correlationheatmap.pdf",
       width =7, height = 7)
