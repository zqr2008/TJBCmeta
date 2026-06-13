# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): Figure5C
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): spreadcazy_high_qual.rda; newdata_phylo_NMrevision_0305.rds; res_kegg_enrichment.rda
# Main output(s): mother.rda; overallko.rda; ko_by_trajectory.rda; ko_by_transmit.rda; plus 1 more
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
library(data.table)
library(tidyverse)
library(sjmisc)


setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/")

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")
bpfa <- spreadcazy %>%
   filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
   distinct(genome,transmit,visit) %>%
   mutate(genome = str_split_fixed(genome,"\\.fa",n=2)[,1] ) %>%
   dplyr::select(genome,transmit,visit) %>%
   mutate(filename = paste0(genome,".kofam.txt"))
  
# step 1: read file
df <- fread("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/cazy_kegg_inputs/bpmerged_with_filename.txt",
            sep = "\t", header = FALSE, fill = TRUE)

# step 2: split ONLY first 6 fields, keep rest as annotation
df2 <- df[, {
  parts <- strsplit(V1, "\\s+")[[1]]
  
  list(
    flag    = parts[1],
    contig  = parts[2],
    KO      = parts[3],
    score1  = as.numeric(parts[4]),
    score2  = as.numeric(parts[5]),
    evalue  = as.numeric(parts[6]),
    annotation = paste(parts[-(1:6)], collapse = " ")
  )
}, by = 1:nrow(df)]

# step 3: add filename
df_final <- cbind(df2[, -1], filename = df$V2)


complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")


metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,trajectory)

overallko <- df_final %>% 
  inner_join(bpfa,by ="filename") %>%
  mutate(sampleid = sub("\\..*", "", filename)) %>%
  mutate(sampleid= str_remove(sampleid,"^D")) %>%
  mutate(sampleid= str_remove(sampleid,"-1$")) %>% 
  mutate(sampleid= str_remove(sampleid,"-1-1$")) %>%
  inner_join(metadata,by ="sampleid") %>%
  distinct(filename, KO,.keep_all = T) %>%
  mutate(countnum =1) %>%
  pivot_wider(id_cols = c(filename,transmit,visit,trajectory),
              names_from = KO,
              values_from =  countnum,
              values_fill = 0)


mothercomparemeta <- overallko %>%
  filter(str_detect(visit,"MG")) %>%
  dplyr::select(1:4)
mothercompare <- overallko %>%
  filter(str_detect(visit,"MG")) %>%
  dplyr::select(-c(2:4)) %>%
  column_to_rownames("filename") %>%
  rotate_df()

save(mothercomparemeta,mothercompare,file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/differential_abundance_objects/mother.rda")

save(overallko,file = "overallko.rda")

ko_cols <- setdiff(colnames(overallko), c("filename", "transmit", "trajectory"))

ko_counts_df <- colSums(overallko[,ko_cols]) %>%
  as.data.frame() %>%
  rownames_to_column("KO") %>%
  rename(count = ".")

ko_by_transmit <- overallko %>%
  group_by(transmit) %>%
  summarise(across(all_of(ko_cols), sum), .groups = "drop")

ko_by_trajectory <- overallko %>%
  group_by(trajectory) %>%
  summarise(across(all_of(ko_cols), sum), .groups = "drop")


save(ko_by_trajectory,file ="ko_by_trajectory.rda")
save(ko_by_transmit,file ="ko_by_transmit.rda")


sepko <-  overallko  %>%
  mutate(sampleid = sub("\\..*", "", filename)) %>%
  mutate(sampleid= str_remove(sampleid,"^D")) %>%
  mutate(sampleid= str_remove(sampleid,"-1$")) %>% 
  mutate(sampleid= str_remove(sampleid,"-1-1$")) %>%
  left_join(metadata,by ="sampleid")


rm(list = ls())
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/res_kegg_enrichment.rda")

res_all$GeneRatio_num <- sapply(strsplit(res_all$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})

plot_df <- res_all %>%
  arrange(qvalue) %>%
  slice(1:20) %>%
  mutate(
    Description = str_replace_all(Description, "_", " "),
    Description = str_trunc(Description, 50)
  )

plot_df <- plot_df %>%
  mutate(
    Description = factor(Description, levels = rev(Description))
  )

library(ggplot2)

p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description)) +
  
  geom_point(aes(size = Count, color = -log10(qvalue)), alpha = 0.8) +
  
  scale_color_gradient(
    low = "#23b1a5",   # consistent teal (from your palette)
    high = "#e95280"   # consistent pink/red
  ) +
  
  scale_size(range = c(2.5, 6)) +
  
  labs(
    x = "Gene ratio",
    y = NULL,
    color = expression(-log[10](qvalue)),
    size = "Gene count"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

p

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure5C BPoverallkegg.pdf",width =8, height = 6)
