# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 03_differential_abundance_Differential_abundance_trajectory
# Figure(s): FigureS5
# Table(s): No direct submitted table matched from filename
# Purpose: Generates or supports submitted figure(s): FigureS5.
# Main input(s): newdata_phylo_NMrevision_0305.rds
# Main output(s): exposure_specific_strains.csv; FigureS5 specific_strain_exposure.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(tidyverse)
library(broom)
library(phyloseq)
library(rstatix)
# -----------------------------
# 1. YOUR ORIGINAL INPUT (unchanged)
# -----------------------------
complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

otu <- as.data.frame(complete_phylo@otu_table)
class(otu) <- "data.frame"

# -----------------------------
# 2. Add sampleID (same as you did)
# -----------------------------
metadata$sampleID <- rownames(metadata)
otu$sampleID <- rownames(otu)

# -----------------------------
# 3. Define variables of interest
# -----------------------------
meta_vars <- c(
  "feeding_type_28days",
  "delivery_mode",
  "has_siblings_early_preg",
  "is_threegenerations_early_preg",
  "secondhand_smoke_exposure_42day",
  "bedroom_floor_tile_stone_early_preg",
  "GDM_group2023",
  "coffee_freq_prepreg"
)

# keep only needed metadata
metadata_sub <- metadata %>%
  select(sampleID,visit, all_of(meta_vars))


metadata_sub[meta_vars] <- lapply(metadata_sub[meta_vars], as.factor)

# -----------------------------
# 4. Define OTU columns CORRECTLY
# -----------------------------
otu_cols <- setdiff(colnames(otu), "sampleID")

# ensure numeric before log transform
otu[otu_cols] <- lapply(otu[otu_cols], as.numeric)

# log transform ONLY OTU columns
otu[otu_cols] <- log1p(otu[otu_cols])

# optional: remove very sparse taxa
otu_cols <- otu_cols[colSums(otu[otu_cols] > 0, na.rm = TRUE) > 10]


allotuwithinfo <-  metadata_sub %>%
  inner_join(otu,by = "sampleID")

# -----------------------------
# 5. Variable type detection
# -----------------------------
get_var_type <- function(x) {
  if (is.numeric(x)) {
    return("continuous")
  } else {
    x <- as.factor(x)
    if (nlevels(x) == 2) return("binary")
    else return("categorical")
  }
}

# -----------------------------
# 6. Association test
# -----------------------------
run_test <- function(x, y) {
  
  df <- data.frame(x = x, y = y) %>%
    drop_na()
  
  if (nrow(df) < 20) return(NULL)
  
  type <- get_var_type(df$y)
  
  tryCatch({
    
    if (type == "continuous") {
      
      res <- cor.test(df$x, df$y, method = "spearman")
      
      tibble(
        test = "spearman",
        estimate = unname(res$estimate),
        p = res$p.value
      )
      
    } else if (type == "binary") {
      
      df$y <- as.factor(df$y)
      lev <- levels(df$y)
      
      res <- df %>%
        wilcox_test(x ~ y,detailed = TRUE,
                             p.adjust.method = "BH")
      print(res)

    
      
    } else {
      
      df$y <- as.factor(df$y)
  
      
      res <- df %>%
        pairwise_wilcox_test(
          x ~ y,
          p.adjust.method = "BH",  # or "bonferroni"
          detailed = TRUE
        )
      
      print(res)
}
    
  }, error = function(e) NULL)
}

# -----------------------------
# 7. Run WITHOUT merging (key fix)
# -----------------------------

result_two <- data.frame()
result_pairwise <- data.frame()
for (i in otu_cols){
  for (j in meta_vars){
    for (tt in c("BM","MG")){
    
  suballotuwithinfo <- allotuwithinfo %>%
    filter(str_detect(visit,tt))
  
  if(sum(suballotuwithinfo[[i]]==0)/(sum(suballotuwithinfo[[i]]==0)+sum(suballotuwithinfo[[i]]!=0)) > 0.9){
   next  
}
  
  x <- suballotuwithinfo[[i]]
  y <- suballotuwithinfo[[j]]
  res <- run_test(x, y) %>%
    mutate(species = i) %>%
    mutate(pheno = j) %>%
    mutate(type =tt)
   
  if(dim(res)[2]==15){
    result_two <- rbind(result_two,res)
  }
  if(dim(res)[2]==17){
    result_pairwise <- rbind(result_pairwise,res)
  }
} 
}
}


exposure_specific_strains <- result_pairwise %>%
  dplyr::select(-c(p.adj,p.adj.signif)) 
 
exposure_specific_strains <- rbind(exposure_specific_strains,result_two) %>%
  filter(pheno != "GDM_group2023") %>%
  mutate(multipletest.q = p.adjust(p,method = "BH")) 


write.csv(exposure_specific_strains, "C:/Users/zqr20/Documents/tjmeta/BIG_revision/06_working_tables/misc_tables/exposure_specific_strains.csv", row.names = FALSE)
##################################################################################
#plot
##################################################################################

top_df <- exposure_specific_strains %>%
  mutate(group1= case_when(group1==1~"Yes",
                           group1==0~"No",
                           TRUE~ group1),
         group2= case_when(group2==1~"Yes",
                           group2==0~"No",
                           TRUE~ group2)) %>%
  filter(group1 != "2-3 times a week" & group1 != "once a week" & group1 != "sometimes")  %>%
  filter(group2 != "2-3 times a week" & group2 != "once a week" & group2 != "sometimes")  %>%
  mutate(species= str_replace(species,"_"," ")) %>%
  filter(!is.na(multipletest.q)) %>%
  group_by(pheno, type) %>%
  arrange(multipletest.q,abs(estimate), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()  

top_df <- top_df %>%
  mutate(
    enriched_group = ifelse(estimate > 0, group1, group2),
    direction_label = paste0(enriched_group, " enriched")
  )

top_df <- top_df %>%
  group_by(pheno, type) %>%
  mutate(species = fct_reorder(species, estimate)) %>%
  ungroup()

top_df <- top_df %>%
  mutate(
    q_label = paste0("q=", formatC(multipletest.q, format = "e", digits = 2))
  )

pheno_labels <- c(
  "is_threegenerations_early_preg" = "Grandparent in household",
  "bedroom_floor_tile_stone_early_preg" = "Tile/stone bedroom floor",
  "has_siblings_early_preg" = "Presence of siblings",
  "secondhand_smoke_exposure_42day" = "Second-hand smoke exposure",
  "delivery_mode" = "Delivery mode",
  "feeding_type_28days" = "Feeding type",
  "coffee_freq_prepreg" = "Coffee drinking frequency"
)

legend_order <- c(
  "Vaginal delivery", "Caesarean section",
  "exclusive breastfeeding", "exclusive formula fed", "mixed feeding",
  "above 3 times a week","2-3 times a week","once a week","never",
  "Yes", "No"
)


top_df$enriched_group <- factor(top_df$enriched_group, levels = legend_order)
my_colors <- c(
  "Vaginal delivery" = "#1b9e77",
  "Caesarean section" = "#d95f02",
  "exclusive breastfeeding" = "#7570b3",
  "exclusive formula fed" = "#e7298a",
  "mixed feeding" = "#66a61e",
  
  "above 3 times a week" = "#a6761d",
  "never" = "#666666",
  
  "Yes" = "#4daf4a",
  "No" = "#e41a1c"
)

p <- ggplot(top_df, aes(x = species, y =  estimate, fill = enriched_group)) +
  geom_col(width = 1, color = "black", linewidth = 0)+
  geom_point(
    aes(y = estimate, color = enriched_group, show.legend = FALSE),
    size = 2
  )+ 
  coord_flip() +
  
  facet_grid(
    type ~ pheno,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(pheno = pheno_labels)
  ) +
  
  scale_fill_manual(values = my_colors, drop = FALSE) +
  scale_color_manual(values = my_colors, drop = FALSE,guide = "none") +
  labs(
    x = NULL,
    y = "Effect size",
    fill = "Enriched in"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "italic"),  # italic species
    axis.text.x = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.spacing = unit(0.8, "lines")
  )

p + geom_text(
  aes(label =q_label),
  hjust = ifelse(top_df$estimate > 0, -0.04, 1),
  size =2.5
) +
  expand_limits(    y = c(min(top_df$estimate) * 2.5,
                          max(top_df$estimate) * 2.6))

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS5 specific_strain_exposure.pdf",
       width = 21, height = 12)

