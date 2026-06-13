# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): FigureS9C
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): experiement.cazy_summary.tsv
# Main output(s): FigureS9C.pdf; experimental_evidence_OD.pdf; experimental_evidence_PCA.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
# ---- Source Rmd chunk: setup (permonva_revise.Rmd lines 16-25) ----
rm(list = ls())

knitr::opts_chunk$set(echo = TRUE)
setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")

set.seed(1000)

# ---- Source Rmd chunk: loading pkgs (permonva_revise.Rmd lines 28-38) ----
library(tidyverse)
library(readxl)
library(phyloseq)
library(sjmisc)
library(lubridate)
library(labelled)
library(conflicted)
library(factoextra)

conflicts_prefer(base::setdiff)
conflicted::conflicts_prefer(dplyr::filter)


# ---- Source Rmd chunk: external BPstrain (permonva_revise.Rmd lines 1179-1272) ----
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
  filter(X.ofTools==3) %>%
  dplyr::select(Gene.ID,DIAMOND,strain) %>%
  mutate(countthis = 1) 

spreadcazy <- mergestrain %>%
  pivot_wider(id_cols = c(strain),
              names_from = DIAMOND,
              values_from =  c(countthis),
              values_fn = function(x)sum(x,na.rm = T),
              values_fill = 0)


cazy_cols <- setdiff(colnames(spreadcazy), "Gene.ID")

# get main CAZyme class for each column
cazy_class <- sapply(cazy_cols, function(x){
  str_extract(x, "^(GH|CBM|GT|CE|AA|PL|GH)")
})



# keep only columns with a detected class
valid_cols <- cazy_cols[!is.na(cazy_class)]
valid_class <- cazy_class[!is.na(cazy_class)]

# sum per row for each class
for(cls in unique(valid_class)){
  cols_to_sum <- valid_cols[valid_class == cls]
  spreadcazy[[paste0("sum_", cls)]] <- rowSums(spreadcazy[, cols_to_sum], na.rm = TRUE)
}

summaryspreadcazy <- spreadcazy %>%
  group_by(strain) %>%
  summarise(GH_count=sum(sum_GH))

custom_colors <- c(
  "#D94E7A",  # muted rose (trajectory_1-like)
  "#2AA89C",  # teal (trajectory_2-like)
  "#D99A1A",  # amber (trajectory_3-like)
  "#4C72B0",  # soft blue
  "#55A868",  # muted green
  "#C44E52",  # brick red
  "#8172B2"   # desaturated purple
)

summaryspreadcazy %>%
  filter(strain != "ATCC_15697") %>%
  ggplot(aes(x = strain, y = GH_count, fill = strain)) +
  geom_col(width = 0.7, alpha = 0.85) +
  
  geom_text(aes(label = GH_count),
            vjust = -0.4,
            size = 5) +
  
  scale_fill_manual(values = custom_colors) +
  
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  
  labs(
    title = "GH gene counts across isolated \nB. pseudocatenulatum strains",
    x = "Strain",
    y = "GH gene count"
  ) +
  
  expand_limits(y = max(summaryspreadcazy$GH_count) + 5)


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS9C.pdf",width = 8, height = 8)




strain_GH <- spreadcazy %>%
  filter(strain != "ATCC_15697") %>%
  group_by(strain) %>%
  summarise(
    across(
      where(is.numeric),
      ~ as.integer(sum(.x, na.rm = TRUE) > 0)
    ),
    .groups = "drop"
  )



mat <- strain_GH %>%
  column_to_rownames("strain") %>%
  dplyr::select(-contains("sum"))

mat2 <- mat[, apply(mat, 2, sd) > 0]

pca <- prcomp(mat2, scale. = TRUE)


fviz_pca_biplot(
  pca,
  label = "all",
  repel = TRUE,
  select.var = list(contrib = 18)
)+theme_classic(base_size = 16)+
  ggtitle("Principal Component Analysis of CAZy for \nexperimental B. pseudocatenulatum")




################################################################################
#experimental
################################################################################
rm(list = ls())
library(readr)

experiement_cazy_summary <- read_delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/experiement.cazy_summary.tsv", 
                                       delim = "\t", escape_double = FALSE, 
                                       col_names = FALSE, trim_ws = TRUE, skip = 1)

colnames(experiement_cazy_summary)[1:7] <- c("strain","Gene.ID","EC.","HMMER",
                                        "dbCAN_sub","DIAMOND","X.ofTools")
                                        





mergestrain <- experiement_cazy_summary %>%
  filter(X.ofTools==3) %>%
  dplyr::select(Gene.ID,DIAMOND,strain) %>%
  mutate(countthis = 1) 

spreadcazy <- mergestrain %>%
  pivot_wider(id_cols = c(strain),
              names_from = DIAMOND,
              values_from =  c(countthis),
              values_fn = function(x)sum(x,na.rm = T),
              values_fill = 0) %>%
  filter(strain != "ZHB-PYG-1") %>%
  filter(strain != "TQ-M-4")


cazy_cols <- setdiff(colnames(spreadcazy), "Gene.ID")

# get main CAZyme class for each column
cazy_class <- sapply(cazy_cols, function(x){
  str_extract(x, "^(GH|CBM|GT|CE|AA|PL|GH)")
})



# keep only columns with a detected class
valid_cols <- cazy_cols[!is.na(cazy_class)]
valid_class <- cazy_class[!is.na(cazy_class)]

# sum per row for each class
for(cls in unique(valid_class)){
  cols_to_sum <- valid_cols[valid_class == cls]
  spreadcazy[[paste0("sum_", cls)]] <- rowSums(spreadcazy[, cols_to_sum], na.rm = TRUE)
}

summaryspreadcazy <- spreadcazy %>%
  group_by(strain) %>%
  summarise(GH_count=sum(sum_GH))

strain_colors <- c(
  "GBWPFCBAD1"  = "#D94E7A",  # muted rose
  "GBXQCBAD1"   = "#2AA89C",  # teal
  "Y0808-CBA-1" = "#4C72B0",  # soft blue
  "YCRPC1-N-11" = "#55A868",  # muted green
  "YCRPC9-N-11" = "#C44E52",  # brick red
  "YF2-M136-4"  = "#8172B2"

)

summaryspreadcazy %>%
  ggplot(aes(
    x = reorder(strain, GH_count),
    y = GH_count,
    fill = strain
  )) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = GH_count),
            vjust = -0.4,
            size = 5) +
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Comparison of GH gene counts across  \nexperimental B. pseudocatenulatum strains",
    x = "Strain",
    y = "GH Gene Count"
  ) +
  expand_limits(y = max(summaryspreadcazy$GH_count) + 5)

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/experimental_evidence_OD.pdf",width = 7, height = 8)




strain_GH <- spreadcazy %>%
  group_by(strain) %>%
  summarise(
    across(
      where(is.numeric),
      ~ as.integer(sum(.x, na.rm = TRUE) > 0)
    ),
    .groups = "drop"
  )



mat <- strain_GH %>%
  column_to_rownames("strain") %>%
  dplyr::select(-contains("sum"))

mat2 <- mat[, apply(mat, 2, sd) > 0]

pca <- prcomp(mat2, scale. = TRUE)


fviz_pca_biplot(
  pca,
  label = "all",
  repel = TRUE,
  select.var = list(contrib = 20)
)+theme_classic(base_size = 16)+
  ggtitle("Principal Component Analysis of CAZy for \nexperimental B. pseudocatenulatum")



ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/experimental_evidence_PCA.pdf",width = 6, height = 6)



table(strain_GH$strain,strain_GH$GH13_14)
table(strain_GH$strain,strain_GH$GH43_22)
table(strain_GH$strain,strain_GH$GH23)
