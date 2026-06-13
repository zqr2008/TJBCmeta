# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): Bifidobacterium.rda; Bifidobacterium_same.rda; newdata_phylo_NMrevision_0305.rds; 1.inStrain_C_genome.profile; plus 1 more
# Main output(s): spreadcazy_high_qual.rda
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/Bifidobacterium.rda")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/Bifidobacterium_same.rda")

complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)


inStrain_C_genome <- read_delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/1.inStrain_C_genome.profile", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE) %>%
  mutate(
    genus = str_extract(classification, "g__[^;]+") %>% str_remove("^g__"),
    species = str_extract(classification, "s__[^;]+") %>% str_remove("^s__")
  )  %>%
  filter(genus=="Bifidobacterium") %>%
  dplyr::select(genome,Samp1,Samp2) %>%
  filter(Samp1 %in% metadata$sampleid & Samp2 %in% metadata$sampleid)






cazy_annotation_genome_bifidoB <- read_delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/cazy_kegg_inputs/8.cazy_annotation_genome.Bifido.reformat.tsv", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE) %>%
  filter(`#ofTools`==3) %>%
  #filter(Quality_class == "HIGH" | Quality_class == "MEDIUM") %>%
  filter(Quality_class == "HIGH" ) %>%
  dplyr::select(genome,DIAMOND,classification) %>%
  mutate(countthis = 1) %>%
  mutate(speciesname = str_split_fixed(classification,"s__",n=2)[,2])


Bifidobacterium <- Bifidobacterium %>%
  dplyr::select(genome,pair_visit) 

Bifidobacterium_same <- Bifidobacterium_same %>%
  dplyr::select(genome,pair_visit) 

spreadcazy <- cazy_annotation_genome_bifidoB %>%
  pivot_wider(id_cols = c(genome,speciesname),
              names_from = DIAMOND,
              values_from =  c(countthis),
              values_fn = function(x)sum(x,na.rm = T),
              values_fill = 0) %>%
  mutate(sum_all = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  mutate(transmit = case_when(genome %in% Bifidobacterium$genome ~"transmitted between dyads",
                              TRUE ~ "No")) %>%
  mutate(persist = case_when(genome %in% Bifidobacterium_same$genome ~"persist in offspring",
                              TRUE ~ "No")) %>%
  mutate(sampleid = sub("\\..*", "", genome)) %>%
  mutate(sampleid= str_remove(sampleid,"^D")) %>%
  mutate(sampleid= str_remove(sampleid,"-1$")) %>% 
  mutate(sampleid= str_remove(sampleid,"-1-1$")) %>% 
  left_join(metadata,by ="sampleid") %>%
  filter(is.na(class_1)==FALSE) %>%
  left_join(Bifidobacterium,by="genome") %>%
  separate(pair_visit, into = c("BM", "MG"), sep = "_", remove = FALSE) %>%
  pivot_longer(cols = c("BM","MG")) %>%
  mutate(value = coalesce(value,visit))  %>%
  distinct(genome,pair_visit,value,.keep_all = T) %>%
  dplyr::select(-visit) %>%
  dplyr::rename(visit = value) 


spreadcazy$visit <- factor(spreadcazy$visit,
                        levels = c("BM1.5","BM6","BM12","BM24","BM36",
                                   "MG1","MG2","MG3"))


cazy_cols <- colnames(spreadcazy)[3:394]
# CAZyme classes to summarize
cazy_classes <- c("GH", "CBM", "GT", "CE", "AA", "PL")

# calculate row-wise sum for columns starting with each CAZyme class
for (cls in cazy_classes) {
  
  cols_to_sum <- cazy_cols[str_detect(cazy_cols, paste0("^", cls))]
  
  if (length(cols_to_sum) > 0) {
    spreadcazy[[paste0("sum_", cls)]] <- rowSums(
      spreadcazy[, cols_to_sum, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    spreadcazy[[paste0("sum_", cls)]] <- 0
  }
}

save(spreadcazy,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")

#save(spreadcazy,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy.rda")

