# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 07_supporting_tables_data_Supporting_tables_data
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Calculates prevalence tables from microbiome profile inputs.
# Main input(s): newdata_phylo_NMrevision_0305.rds
# Main output(s): allpre.tsv
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(Maaslin2)
library(microViz)
library(tidyverse)
library(microbiome)


data_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")
  
otu <- as.data.frame(data_phylo@otu_table)
class(otu) <- "data.frame"

phenotype <- as.data.frame(data_phylo@sam_data)
class(phenotype) <- "data.frame"


otu_pa <- otu > 0
prev_overall <- colMeans(otu_pa, na.rm = TRUE)

prev_overall_df <- data.frame(
  visit="overall",
  species = names(prev_overall),
  prevalence = prev_overall
)


prev_by_visit <- otu %>%
  mutate(visit = phenotype$visit) %>%
  group_by(visit) %>%
  summarise(across(everything(), ~ mean(. > 0))) %>%
  tidyr::pivot_longer(
    cols = -visit,
    names_to = "species",
    values_to = "prevalence"
  )


prevalencemerge<- rbind(prev_overall_df,prev_by_visit)



write.table(prevalencemerge,file = "allpre.tsv",
            fileEncoding = "GBK",sep = "\t",row.names = F,col.names = F,
            quote = F)