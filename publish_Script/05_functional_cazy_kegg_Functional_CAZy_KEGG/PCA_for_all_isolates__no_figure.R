# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): spreadcazy_high_qual.rda; experiement.cazy_summary.tsv
# Main output(s): none detected
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())


library(factoextra)
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")


experiement_cazy_summary <- read_delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/experiement.cazy_summary.tsv", 
                                       delim = "\t", escape_double = FALSE, 
                                       col_names = FALSE, trim_ws = TRUE, skip = 1)

colnames(experiement_cazy_summary)[1:7] <- c("strain","Gene.ID","EC.","HMMER",
                                             "dbCAN_sub","DIAMOND","X.ofTools")






mergestrain <- experiement_cazy_summary %>%
  filter(X.ofTools==3) %>%
  dplyr::select(Gene.ID,DIAMOND,strain) %>%
  mutate(countthis = 1) 

experiementspreadcazy <- mergestrain %>%

  pivot_wider(id_cols = c(strain),
              names_from = DIAMOND,
              values_from =  c(countthis),
              values_fn = function(x)sum(x,na.rm = T),
              values_fill = 0) %>%
  dplyr::select(strain,starts_with("GH")) %>%
  mutate(strain = paste0(strain,"_experiemental"))
  




setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/cazy_kegg_inputs/allfna.dbcan_out")
base_dir <- "C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/cazy_kegg_inputs/allfna.dbcan_out"

mergestrain <- list.files(base_dir, full.names = TRUE) %>%
  map_dfr(function(dir_path) {
    
    file_path <- file.path(dir_path, "overview.txt")
    
    if (!file.exists(file_path)) return(NULL)
    
    # extract strain name (remove ".dbcan_out")
    strain <- basename(dir_path) %>%
      str_remove("\\.dbcan_out$")
    
    read.delim(file_path) %>%
      mutate(strain = strain)
  })

mergestrain <- mergestrain %>%
  filter(strain != "ATCC_15697") %>%
  filter(X.ofTools==3) %>%
  dplyr::select(Gene.ID,DIAMOND,strain) %>%
  mutate(countthis = 1) %>%
  mutate(strain = paste0(strain,"_external"))

spreadcazyexternal <- mergestrain %>%
  pivot_wider(id_cols = c(strain),
              names_from = DIAMOND,
              values_from =  c(countthis),
              values_fn = function(x)sum(x),
              values_fill = 0)



combined <- bind_rows(
  spreadcazyexternal,
  experiementspreadcazy
)

combined[is.na(combined)]<- 0

mat <- combined %>%
  column_to_rownames("strain") 
mat2 <- mat[, apply(mat, 2, sd) > 0]

pca <- prcomp(mat2, scale. = TRUE)



pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)



colors15 <- c(
  "#D94E7A",  # rose
  "#2AA89C",  # teal
  "#D99A1A",  # amber
  "#4C72B0",  # blue
  "#55A868",  # green
  "#C44E52",  # red
  "#8172B2",  # purple
  "#937860",  # brown
  "#64B5CD",  # cyan
  "#CCB974",  # olive
  "#8C8C8C",  # gray
  "#E17C05",  # orange
  "#76B7B2",  # turquoise
  "#B07AA1",  # mauve
  "#F28E2B"   # light orange
)

ggplot(pca_df, aes(x=PC1, y=PC2,color=Sample,fill=Sample)) +
  geom_point(size = 5) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0(
      "PC1 (",
      round(100 * summary(pca)$importance[2, 1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(100 * summary(pca)$importance[2, 2], 1),
      "%)"
    ),
    title = "Principal Component Analysis"
  )+
  scale_color_manual(values = colors15)

