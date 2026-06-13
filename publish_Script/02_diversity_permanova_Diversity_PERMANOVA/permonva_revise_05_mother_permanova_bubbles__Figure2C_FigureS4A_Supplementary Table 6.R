# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 02_diversity_permanova_Diversity_PERMANOVA
# Figure(s): Figure2C; FigureS4A
# Table(s): Supplementary Table 6
# Purpose: Runs or organizes PERMANOVA association summaries.
# Main input(s): none detected
# Main output(s): Figure2C.pdf; FigureS4A_technical_and_biological.pdf
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
conflicts_prefer(base::setdiff)

# ---- Source Rmd chunk: mother bubble (permonva_revise.Rmd lines 689-1001) ----
attentionlist_MG<- c(
  "Age_mo",
  "BMI_mo",
  "coffee_freq_prepreg",
  "edu_mo",
  "GDM_group2023",
  "HDLC",
  "hypothyroidism_history_mo",
  "TG",
  "UA.y",
  "weight_before_delivery",
  "weight_mo_first_check",
  "weight_mo_pre_preg",
  "trajectory",
  "realsampleday",
  "antibiotic_midpregnancy",
  "antibiotic_latepregnancy",
    "is_threegenerations_early_preg",
  "bedroom_floor_tile_stone_early_preg",
  "has_siblings_early_preg",  
  "secondhand_smoke_exposure_42day",
  "coffee_freq_prepreg"
)

setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/permanova_results/final_permonva")

MG_hold <- list()

patterns <- c("^MG1.*\\.rda$","^MG2.*\\.rda$","^MG3.*\\.rda$")

for (kk in patterns) {
  
  myFiles <- list.files(pattern = kk)
  
  temp_list <- list()
  
  for (i in myFiles) {
    
    e <- new.env()
    load(i, envir = e)
    
    obj_name <- ls(e)[1]   # object stored in rda
    temp <- as.data.frame(e[[obj_name]])
    
    temp$file <- i
    temp$visittype <- kk
    
    temp_list[[i]] <- temp
  }
  
  MG_hold[[kk]] <- bind_rows(temp_list)
}


holdresult <- data.frame()

MG_all_adjust_technical <- bind_rows(MG_hold) %>%
  filter(R2>0.00001) %>%
  filter(mainfocus ==term) %>%
  filter(term %in% attentionlist_MG ) %>% 
  group_by(visit)  %>%
  filter(covariate=="technical(depth and batch)" )%>%

  mutate(q = p.adjust(`Pr(>F)`,method="BH",n = n())) %>%
  mutate(Significance = case_when(
    q < 0.01  ~ "Multiple-testing–corrected q-values <0.01",
    q < 0.05 & q >= 0.01 ~ "Multiple-testing–corrected q-values <0.05",
    q < 0.1 & q > 0.05 ~ "Multiple-testing–corrected q-values <0.1",
    TRUE ~ "Not significant") )%>%
  mutate( Significance = factor(Significance, levels = c(
      "Multiple-testing–corrected q-values <0.01",
      "Multiple-testing–corrected q-values <0.05",
      "Multiple-testing–corrected q-values <0.1",
      "Not significant")
  ))





MG_all_adjust_technical$visit <- factor(MG_all_adjust_technical$visit,levels = c("MG1","MG2","MG3"))



vars_to_keep <- c(
  "trajectory",
  "realsampleday",
  "weight_mo_pre_preg",
  "weight_mo_first_check",
  "weight_before_delivery",
  "BMI_mo",
  "GDM_group2023",
  "hypothyroidism_history_mo",
  "Age_mo",
  "edu_mo",
  "TG",
  "HDLC",
  "UA.y",
  "coffee_freq_prepreg",

  "realsampleday",
  "antibiotic_midpregnancy",
  "antibiotic_latepregnancy",
  "is_threegenerations_early_preg",
  "bedroom_floor_tile_stone_early_preg" ,
  "has_siblings_early_preg",
  "secondhand_smoke_exposure_42day"

)

# Corresponding readable names
readable_names <- c(
  "Trajectory",
  "gestational age",
  "Prepregnancy weight",
  "Weight at first prenatal visit",
  "Weight before delivery",
  "Prepregnancy BMI",
  "GDM",
  "History of maternal hypothyroidism",
  "Age of mother",
  "Education level of mother",
  "TG",
  "HDLC",
  "UA",
  "Coffee drinking frequency before pregnancy",
 
  "gestational age",
  "Mid-pregnancy antibiotic exposure",
  "Late-pregnancy antibiotic exposure",
  "Grandparent in household",
  "Tile/stone bedroom floor",
  "Presence of siblings",
  "Second-hand smoke exposure"
)

# Create a named vector for renaming
name_map <- setNames(readable_names, vars_to_keep)

# Apply to your dataset
MG_filter <- MG_all_adjust_technical %>%
  dplyr::mutate(newname = name_map[term])

MG_filter <- MG_filter %>%
  dplyr::select(newname,Df,SumOfSqs,R2,`F`,`Pr(>F)`,visit,covariate,q,Significance)

holdresult <- rbind(holdresult,MG_filter)




target <- c(
  "Trajectory",
  # Maternal anthropometrics
  "Prepregnancy weight",
  "Weight at first prenatal visit",
  "Weight before delivery",
  "Prepregnancy BMI",
  "Age of mother",
  "Education level of mother",
  
  # Maternal health / conditions
  "GDM",
  "History of maternal hypothyroidism",
  
  # Laboratory measures
  "TG",
  "HDLC",
  "UA",
  
  # Lifestyle / exposures
  "Coffee drinking frequency before pregnancy",
  "Grandparent in household",
  "Tile/stone bedroom floor",
  "Presence of siblings",
  "Second-hand smoke exposure",
  
  # Child trajectory

  "gestational age",
  "Mid-pregnancy antibiotic exposure",
  "Late-pregnancy antibiotic exposure"
)

MG_filter %>%
  ggplot(aes(x = factor(newname,level = target), y = visit,
             color = R2*100, shape=Significance)) +
  geom_point(size = 5) +
  scale_color_gradientn(
    colours = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45"),
    values  = scales::rescale(c(0, 1, 5, 10)),  # expand low end
    name = "Variation explained \n (partial R-squared)"
  )+
  ylab("Visit") +
  xlab("Phenotype") +
  ggtitle(label = "")+
  scale_shape_manual(values = c(
    "Multiple-testing–corrected q-values <0.01" = 16,
    "Multiple-testing–corrected q-values <0.05" = 17,
    "Multiple-testing–corrected q-values <0.1"  = 15,
    "Not significant" = 4
  ))+
  theme (plot.title = element_text (size = 14, face = "bold" ))+
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=18,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.text.x = element_text(size=18,face="bold",angle = 60, hjust = 1))+
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+ 
  guides(fill = guide_legend(override.aes = list(size=8)))  



ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure2C.pdf",width = 12, height = 8)






MG_technical_and_biologic <- bind_rows(MG_hold) %>%
  filter(R2>0.00001) %>%
  filter(mainfocus ==term) %>%
  filter(term %in% attentionlist_MG ) %>% 
  group_by(visit)  %>%
 
  filter(covariate=="technical(depth and batch) and biologic(BMI_mo,TG,UA)") %>%
  mutate(q = p.adjust(`Pr(>F)`,method="BH",n = n())) %>%
  mutate(Significance = case_when(
    q < 0.01  ~ "Multiple-testing–corrected q-values <0.01",
    q < 0.05 & q >= 0.01 ~ "Multiple-testing–corrected q-values <0.05",
    q < 0.1 & q > 0.05 ~ "Multiple-testing–corrected q-values <0.1",
    TRUE ~ "Not significant") )%>%
  mutate( Significance = factor(Significance, levels = c(
      "Multiple-testing–corrected q-values <0.01",
      "Multiple-testing–corrected q-values <0.05",
      "Multiple-testing–corrected q-values <0.1",
      "Not significant")
  ))





MG_technical_and_biologic$visit <- factor(MG_technical_and_biologic$visit,levels = c("MG1","MG2","MG3"))


# Create a named vector for renaming
name_map <- setNames(readable_names, vars_to_keep)

# Apply to your dataset
MG_filter <- MG_technical_and_biologic %>%
  dplyr::mutate(newname = name_map[term])

MG_filter <- MG_filter %>%
  dplyr::select(newname,Df,SumOfSqs,R2,`F`,`Pr(>F)`,visit,covariate,q,Significance)

holdresult <- rbind(holdresult,MG_filter)

MG_filter %>%
  ggplot(aes(x = factor(newname,level = target), y = visit,
             color = R2*100, shape=Significance)) +
  geom_point(size = 5) +
  scale_color_gradientn(
    colours = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45"),
    values  = scales::rescale(c(0, 1, 5, 10)),  # expand low end
    name = "Variation explained \n (partial R-squared)"
  )+
  ylab("Visit") +
  xlab("Phenotype") +
  ggtitle(label = "")+
  theme (plot.title = element_text (size = 14, face = "bold" ))+
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=18,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.text.x = element_text(size=18,face="bold",angle = 60, hjust = 1))+
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+ 
  scale_shape_manual(values = c(
    "Multiple-testing–corrected q-values <0.01" = 16,
    "Multiple-testing–corrected q-values <0.05" = 17,
    "Multiple-testing–corrected q-values <0.1"  = 15,
    "Not significant" = 4
  ))+
  guides(fill = guide_legend(override.aes = list(size=8)))  

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS4A_technical_and_biological.pdf",width = 12, height = 8)


write.table(holdresult,
            file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 6 MG_all.tsv",fileEncoding = "GBK",row.names = F,
            col.names = T,quote = F,sep = "\t"
            )
