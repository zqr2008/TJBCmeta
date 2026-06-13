# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): Figure5D
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): newdata_phylo_NMrevision_0305.rds; spreadcazy_high_qual.rda
# Main output(s): Figure5D PCA of CAZy profiles.pdf; PCA of CAZy profiles trajectory.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(tidyverse)
library(vegan)


complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/cazy_kegg_objects/spreadcazy_high_qual.rda")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid")



spreadcazy <- spreadcazy %>%
  distinct(genome,visit,.keep_all = T) %>%
  filter(speciesname=="Bifidobacterium pseudocatenulatum") %>%
  inner_join(metadata,by ="sampleid") %>%
  mutate(type = case_when(str_detect(visit.x,"MG") & transmit =="No"~ "mother specific",
                          str_detect(visit.x,"BM") & transmit =="No"~ "offspring specific",
                          transmit =="transmitted between dyads" ~ "transmitted"))
  
 

# 1. 提取子集
cazy_subset <- spreadcazy[, 3:339]
richness <- as.data.frame(specnumber(cazy_subset))
colnames(richness)[1] <- "index_richness"


cazy_subset[is.na(cazy_subset)] <- 0

# 2. 检查每一列的方差，找出方差不为 0 的列
# apply(..., 2, sd) 计算每一列的标准差
nonzero_cols <- apply(cazy_subset, 2, sd) > 0

# 3. 仅保留有变异的列
cazy_filtered <- cazy_subset[, nonzero_cols]


# 5. 重新运行 PCA
pca_res <- prcomp(cazy_filtered, center = TRUE, scale. = TRUE)

# 6. 查看结果
summary(pca_res)

# Create a data frame for plotting
pca_data <- data.frame(
  SampleID = spreadcazy[, "sampleid"],
  Group = spreadcazy[, "GH2_18"], 
  Group = spreadcazy[, "CE20"], 
  Group = spreadcazy[, "GH127"],
  Group = spreadcazy[, "GH43_22"],
  Group = spreadcazy[, "sum_GH"],
  Group = spreadcazy[, "trajectory"],
  Group = spreadcazy[, "transmit"],
  Group = spreadcazy[, "visit.x"],
  Group = spreadcazy[, "type"],
  index_richness = richness[, "index_richness"],
  pca_res$x[, 1:2]            # Extract PC1 and PC2
)

# Calculate variance for axis labels
var_exp <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)



type_cols <- c(
  "mother specific" = "#66429b",   # muted blue
  "offspring specific"  = "#166678",   # muted green
  "transmitted"  = "#d73027"    # muted red
)


# make convex hull for each type
hull_data <- pca_data %>%
  group_by(type) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

type_fill_colors <- c(
  "offspring specific"  = "#9ECAE1",
  "transmitted" = "#F4A6A6",
  "mother specific"= "#A1D99B"
)

pca_data %>%
  ggplot(aes(x = PC1, y = PC2)) +
  
  geom_polygon(
    data = hull_data,
    aes(x = PC1, y = PC2, fill = type, group = type),
    alpha = 0.18,
    color = NA
  ) +
  
  geom_point(
    aes(color = index_richness, shape = type),
    size = 2,
    alpha = 0.45
  ) +
  
  scale_fill_manual(values = type_fill_colors) +
  scale_color_viridis_c(option = "magma", end = 0.95) +
  
  theme_classic(base_size = 18) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "PCA of CAZy profiles in \nBifidobacterium pseudocatenulatum genomes",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    fill = "Type",
    shape = "Type",
    color = "Richness"
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure5D PCA of CAZy profiles.pdf",width =7, height = 6)





pca_data %>%
  mutate(GH2_18 = as.factor(GH2_18)) %>%
  ggplot(aes(x = PC1, y = PC2,
             color = trajectory,
             shape = GH2_18)) +
  
  geom_point(size = 4, alpha = 0.3) +
  scale_color_manual(values =  c(
    "trajectory_1" = "#e95280",
    "trajectory_2" = "#23b1a5",
    "trajectory_3" = "#E49B0F"
  ))+
  theme_classic(base_size = 18) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  
  labs(
    title = "PCA of CAZy profiles in \nBifidobacterium pseudocatenulatum genomes",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  )

ggsave("PCA of CAZy profiles trajectory.pdf",width =6, height = 4)

