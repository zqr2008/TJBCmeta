# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 02_diversity_permanova_Diversity_PERMANOVA
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Builds PCoA ordination analysis outputs.
# Main input(s): ps_filtered.rds; dist_aitchison_pc_rda.rda
# Main output(s): none detected
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(microViz)
library(tidyverse)
library(microbiome)
library(cluster)

setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")


data_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/ps_filtered.rds")


supple <- phyloseq_validate(data_phylo, remove_undetected = TRUE)
supple <- tax_fix(supple)



offspring_1.5 <- supple %>%
  ps_filter(
    str_detect(visit,"BM6") | str_detect(visit,"MG3") 
  ) %>%
  ps_filter(
    is.na(class_1) == FALSE
  ) 


input_data_offspring <- as.data.frame(t(offspring_1.5@otu_table))
Metadata_offspring <- as.data.frame((offspring_1.5@sam_data))
class(Metadata_offspring) <- "data.frame"
class(input_data_offspring) <- "data.frame"


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/dist_aitchison_pc_rda.rda")


bmpcoa <- dist_aitchison_pc_rda[rownames(Metadata_offspring),rownames(Metadata_offspring)]

bmpcoa <- cmdscale(bmpcoa, k=2, eig=TRUE) 

pcoa_scores <- bmpcoa$points  # sample x axis

# Make sure sample order matches
abund <- input_data_offspring[,rownames(pcoa_scores)]

# 2️⃣ Correlate species (rows) with PCoA axes
corr_axis1 <- apply(abund, 1, function(x) cor(x, pcoa_scores[, 1], method = "spearman"))
corr_axis2 <- apply(abund, 1, function(x) cor(x, pcoa_scores[, 2], method = "spearman"))

# 3️⃣ Combine results
corr_df <- data.frame(
  species = rownames(abund),
  axis1 = corr_axis1,
  axis2 = corr_axis2
)

# 4️⃣ Identify top drivers
top_axis1 <- corr_df[order(abs(corr_df$axis1), decreasing = TRUE), ][1:15, ]
top_axis2 <- corr_df[order(abs(corr_df$axis2), decreasing = TRUE), ][1:15, ]




df_pcoa <- data.frame(
  SampleID = rownames(bmpcoa$points),
  Axis1 = bmpcoa$points[,1],
  Axis2 = bmpcoa$points[,2]
)

joinpcoa <- Metadata_offspring %>%
  rownames_to_column("SampleID") %>%
  inner_join(df_pcoa,by = "SampleID") %>%
  dplyr::select(SampleID,familyid,visit,class_1,Axis1,Axis2) %>%
  group_by(familyid) %>%
  mutate(n=n())%>%
  filter(n==2)%>%
  dplyr::select(-n) %>%
  pivot_wider(id_cols =familyid,
              names_from = visit,
              values_from = c(Axis1,Axis2,class_1) ) %>%
  mutate(Axis1diff = Axis1_MG3 - Axis1_BM6,
         Axis2diff = Axis2_MG3 - Axis2_BM6)

class(joinpcoa) <- "data.frame"


stat.test_1 <- joinpcoa %>%
  t_test(reformulate("class_1_MG3", "Axis2diff" ),detailed = TRUE,p.adjust.method ="BH") %>%
  add_significance("p.adj")   %>%
  filter(p.adj.signif != "ns") %>%
  add_xy_position(x = "class_1_MG3",step.increase = 1,fun = "median_iqr",
                  dodge = 0.5) %>%
  mutate(color = case_when(group1 == "1" & group2 == "2" ~ "1",
                           group1 == "2" & group2 == "3" ~ "2",
                           group1 == "1" & group2 == "3" ~ "3")) 

joinpcoa %>% 
  ggplot(aes( x = class_1_MG3 ,y = Axis2diff, fill =as.factor(class_1_MG3))) +
  geom_boxplot() + 
  theme_classic()+
  xlab(paste0("Trajectories")) +
  ylab("Axis2 differences with their mother ") +
  stat_ellipse(linetype = 2,
               lwd = 1.0)+
  scale_fill_manual(values = c("#e95280","#23b1a5","#E49B0F","grey")) +
  scale_color_manual(values = c("#e95280","#23b1a5","#E49B0F","grey")) +
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        plot.title = element_text(size=22, face="bold"),
        plot.subtitle=element_text(size=18,  color="black"))+
  stat_pvalue_manual(
    stat.test_1, label = "p.adj.signif", tip.length = 0,
    color = "color",fill = "color",hide.ns = F, size= 6
  )
