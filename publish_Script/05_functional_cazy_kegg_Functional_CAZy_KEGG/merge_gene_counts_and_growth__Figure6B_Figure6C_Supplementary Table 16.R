# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): Figure6B; Figure6C
# Table(s): Supplementary Table 16
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): validation.cazy.39_8.tsv; allgrow.rda
# Main output(s): Figure6B.pdf; Supplementary Table 16 isolates grow.csv; Figure6C.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(ggrepel)

validation.cazy.39_8 <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/cazy_kegg_inputs/validation.cazy.39_8.tsv", header=FALSE)

colnames(validation.cazy.39_8) <- c("strain","Gene.ID","EC.","HMMER","dbCAN_sub","DIAMOND","X.ofTools")  
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/strain_growth_objects/allgrow.rda")
allstrain_plot <- allstrain_plot %>%
  mutate(strainname = toupper(strainname))


mergestrain <- validation.cazy.39_8 %>%
  filter(X.ofTools==3) %>%
  dplyr::select(Gene.ID,DIAMOND,strain) %>%
  mutate(countthis = 1) 





spreadcazy <- mergestrain %>%
  pivot_wider(id_cols = c(strain),
              names_from = DIAMOND,
              values_from =  c(countthis),
              values_fn = function(x)sum(x,na.rm = T),
              values_fill = 0)   %>%
  filter(!strain %in% c(
                "ZHB-PYG-1.scaffolds.fasta",
                "strain19.fna",
                "ATCC_15697.fna"
              ))    

#ZHB-PYG-1 is removed because the SD too large(measurement error)
#strain19 is removed because contamination

strain_rename_key <- tibble::tribble(
  ~strain,                         ~strain_internal,
  "GBWPFCBAD1.scaffolds.fasta",    "BP_001",
  "GBXQCBAD1.scaffolds.fasta",     "BP_002",
  "Strain5",                       "BP_003",
  "Strain6",                       "BP_004",
  "Strain7",                       "BP_005",
  "TQ-M-4.scaffolds.fasta",        "BP_006",
  "Y0808-CBA-1.scaffolds.fasta",   "BP_007",
  "YCRPC1-N-11.scaffolds.fasta",   "BP_008",
  "YCRPC9-N-11.scaffolds.fasta",   "BP_009",
  "YF2-M136-4.scaffolds.fasta",    "BP_010",
  "strain1",                       "BP_011",
  "strain10.fna",                  "BP_012",
  "strain11.fna",                  "BP_013",
  "strain12.fna",                  "BP_014",
  "strain13.fna",                  "BP_015",
  "strain14.fna",                  "BP_016",
  "strain15.fna",                  "BP_017",
  "strain16.fna",                  "BP_018",
  "strain17.fna",                  "BP_019",
  "strain18.fna",                  "BP_020",
  "strain2.fna",                   "BP_021",
  "strain20.fna",                  "BP_022",
  "strain21.fna",                  "BP_023",
  "strain22.fna",                  "BP_024",
  "strain23.fna",                  "BP_025",
  "strain24.fna",                  "BP_026",
  "strain25.fna",                  "BP_027",
  "strain26.fna",                  "BP_028",
  "strain27.fna",                  "BP_029",
  "strain28.fna",                  "BP_030",
  "strain29.fna",                  "BP_031",
  "strain3.fna",                   "BP_032",
  "strain30.fna",                  "BP_033",
  "strain31.fna",                  "BP_034",
  "strain32.fna",                  "BP_035",
  "strain4.fna",                   "BP_036",
  "strain9.fna",                   "BP_037"
)


spreadcazy <- spreadcazy %>%
  left_join(strain_rename_key, by = "strain") %>%
  mutate(strainname_plot = strain_internal) %>%
  mutate(strain = toupper(str_split_fixed(strain, "\\.", n = 2)[, 1]))







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






#ZHB-PYG-1 is removed because the SD too large(measurement error)
#strain19 is removed because contamination
summaryspreadcazy <- spreadcazy %>%
  group_by(strainname_plot,strain) %>%
  summarise(GH_count = sum(sum_GH), .groups = "drop") %>%
  mutate(strainname= strain)%>%
  inner_join(allstrain_plot, by = "strainname") %>%
  mutate(
    genecoutcate = case_when(
      GH_count >= 50 ~ "High gene counts (>=50)",
      GH_count < 50  ~ "Low gene counts (<50)"
    ),
    genecoutcate = factor(genecoutcate, levels = c("Low gene counts (<50)", "High gene counts (>=50)"))
  )




n_df <- summaryspreadcazy %>%
  count(condition_formal, genecoutcate)



summaryspreadcazy <- summaryspreadcazy %>%
  dplyr::select(condition_formal,strainname_plot,GH_count,genecoutcate,mean_OD600_blank_subtracted,sd_OD600_blank_subtracted)


summaryspreadcazy %>%
  ggplot(aes(
    x = genecoutcate,
    y = mean_OD600_blank_subtracted,
    fill = genecoutcate
  )) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.75,
    outlier.shape = NA,
    color = "grey25",
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(fill = genecoutcate),
    shape = 21,
    width = 0.12,
    size = 2.3,
    alpha = 0.85,
    color = "grey20",
    stroke = 0.35
  ) +
  geom_text(
    data = n_df,
    aes(
      x = genecoutcate,
      y = -Inf,
      label = paste0("n = ", n)
    ),
    inherit.aes = FALSE,
    vjust = -0.1,
    size = 4
  ) +
  stat_compare_means(
    comparisons = list(c("Low gene counts (<50)", "High gene counts (>=50)")),
    method = "t.test",
    label = "p.format",
    size = 3.5,
    bracket.size = 0.35,
    tip.length = 0.01
  ) +
  facet_wrap(~ condition_formal, scales = "free_y", ncol = 4) +
  scale_fill_manual(
    values = c(
      "Low gene counts (<50)" = "#8DA0CB",
      "High gene counts (>=50)" = "#FC8D62"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Low gene counts (<50)" = "Low gene \ncounts(<50)",
      "High gene counts (>=50)" = "High gene \ncounts(>=50)"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.18))) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    
    strip.text = element_text(
      size = 18,
      face = "bold",
      margin = margin(t = 2, r = 2, b = 2, l = 2)
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = "grey50",
      linewidth = 0.3
    ),
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5
    ),
    
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  labs(
    title = "Growth comparison of B. pseudocatenulatum isolates \nstratified by GH gene counts",
    x = "GH gene counts category",
    y = "Mean blank-subtracted OD600 across 3 technical replicates"
  )


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure6B.pdf",width = 11, height = 9)

write.table(summaryspreadcazy, file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 16 isolates grow.csv",
            fileEncoding = "GBK",col.names = T,row.names = F,
            sep = ",",quote = F)



mat <- spreadcazy %>%
  filter(strain %in% allstrain_plot$strainname) %>%
  column_to_rownames("strainname_plot") %>%
  dplyr::select(starts_with("GH"))


mat2 <- mat[, apply(mat, 2, sd) > 0]

pca <- prcomp(mat2, scale. = TRUE)

pca_meta <- summaryspreadcazy %>%
  distinct(strainname_plot, genecoutcate) %>%
  filter(strainname_plot %in% rownames(mat2)) %>%
  arrange(match(strainname_plot, rownames(mat2)))


ind_lab <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("strainname_plot") %>%
  rename(Dim.1 = PC1, Dim.2 = PC2)

fviz_pca_biplot(
  pca,
  label = "var",              # only label arrows/GH variables
  repel = TRUE,
  labelsize = 4,              # bigger arrow names
  pointsize = 2,
  select.var = list(contrib = 18),
  habillage = pca_meta$genecoutcate,
  addEllipses = TRUE,
  ellipse.level = 0.95,
  arrowsize = 0.5
) +
  geom_text_repel(
    data = ind_lab,
    aes(x = Dim.1, y = Dim.2, label = strainname_plot),
    inherit.aes = FALSE,
    size = 1.8,                 # small dot/strain names
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.15,
    segment.size = 0.2
  ) +
  scale_color_manual(
    values = c(
      "Low gene counts (<50)" = "#8DA0CB",
      "High gene counts (>=50)" = "#FC8D62"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Low gene counts (<50)" = "#8DA0CB",
      "High gene counts (>=50)" = "#FC8D62"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  ) +
  xlab("PC 1 (15.0%)")+
  ylab("PC 2 (2.4%)")+
  labs(
    color = "GH gene count category",
    fill = "GH gene count category"
  ) +
  ggtitle("PCA of B. pseudocatenulatum isolates GH genes\nstratified by GH gene counts")

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure6C.pdf",width = 7, height = 6)
