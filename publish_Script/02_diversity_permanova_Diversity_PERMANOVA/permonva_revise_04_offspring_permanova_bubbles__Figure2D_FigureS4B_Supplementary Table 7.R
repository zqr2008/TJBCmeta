# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 02_diversity_permanova_Diversity_PERMANOVA
# Figure(s): Figure2D; FigureS4B
# Table(s): Supplementary Table 7
# Purpose: Runs or organizes PERMANOVA association summaries.
# Main input(s): none detected
# Main output(s): Figure2D.pdf; FigureS4B_technical_and_biological.pdf
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

# ---- Source Rmd chunk: offspring bubble (permonva_revise.Rmd lines 390-688) ----
conflicted::conflicts_prefer(dplyr::filter)

attentionlist_BM <- c(
  "Age_fa",
  "Age_mo",
  "AST",
  "BMI_mo",
  "trajectory",
  "BUN",
  "delivery_mode",
  "edu_fa",
  "edu_mo",
  "feeding_type_28days",
  "GDM_group2023",
  "gender_kid",
  "antibiotic_usage_180day",
  "HDLC",
  "preterm",
  "TG",
  "UA.y",
  "is_threegenerations_early_preg",
  "bedroom_floor_tile_stone_early_preg",
  "has_siblings_early_preg",  
  "secondhand_smoke_exposure_42day",
  "coffee_freq_prepreg"
)




setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/permanova_results/final_permonva")

BM_hold <- list()

patterns <- c("BM1\\.5.*\\.rda",
              "^BM6.*\\.rda$",
              "^BM12.*\\.rda$",
              "^BM24.*\\.rda$",
              "^BM36.*\\.rda$")

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
  
  BM_hold[[kk]] <- bind_rows(temp_list)
}


holdresult <- data.frame()

BM_all_adjust_technical <- bind_rows(BM_hold) %>%
  filter(R2>0.00001) %>%
  filter(mainfocus ==term) %>%
  filter(term %in% attentionlist_BM ) %>% 
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


BM_all_adjust_technical$visit <- factor(
  BM_all_adjust_technical$visit,
  levels = c("BM1.5", "BM6", "BM12", "BM24", "BM36")
)

BM_filter <- BM_all_adjust_technical %>% 
  mutate(
    newname = case_when(
      term == "trajectory" ~ "Trajectory",
      term == "delivery_mode" ~ "Delivery mode",
      term == "feeding_type_28days" ~ "Feeding type",
      term == "BMI_mo" ~ "Prepregnancy BMI",
      term == "preterm" ~ "Preterm",
      term == "gender_kid" ~ "Gender",
      term == "edu_fa" ~ "Education level of father",
      term == "edu_mo" ~ "Education level of mother",
      term == "family_income" ~ "Family income",
      term == "gravida" ~ "Gravida",
      term == "Age_mo" ~ "Age of mother",
      term == "Age_fa" ~ "Age of father",
      term == "TG" ~ "TG",
      term == "UA.y" ~ "UA", 
      term == "HDLC" ~ "HDLC",
      term == "AST" ~ "AST",
      term == "BUN" ~ "BUN",
      term == "GDM_group2023" ~ "GDM",
      term == "antibiotic_usage_180day" ~ "Antibiotic use within 180 days",
      term == "secondhand_smoke_exposure_42day" ~ "Second-hand smoke exposure",
      term == "is_threegenerations_early_preg" ~ "Grandparent in household",
      term == "bedroom_floor_tile_stone_early_preg" ~ "Tile/stone bedroom floor",
      term == "has_siblings_early_preg" ~ "Presence of siblings",
      term == "coffee_freq_prepreg" ~ "Coffee drinking frequency",
      TRUE ~ term
    )
  )

target <- c(
  "Trajectory",
  "Delivery mode", 
  "Feeding type",
  "Prepregnancy BMI",
  "Preterm",
  "Gender",
  "Education level of father",
  "Education level of mother",
  "Age of mother",
  "Age of father",
  "HDLC",
  "UA",
  "TG",
  "AST",
  "BUN",
  "GDM",
  "Antibiotic use within 180 days",
  "Second-hand smoke exposure",
  "Grandparent in household",
  "Tile/stone bedroom floor",
  "Presence of siblings",
  "Coffee drinking frequency"
)

BM_filter <- BM_filter %>%
  mutate(
    newname = factor(newname, levels = target)
  ) %>%
  dplyr::select(
    newname, Df, SumOfSqs, R2, `F`, `Pr(>F)`,
    visit, covariate, q, Significance
  )

holdresult <- rbind(holdresult, BM_filter)

BM_pic <- BM_filter %>%
  ggplot(
    aes(
      x = newname,
      y = visit,
      color = R2 * 100,
      shape = Significance
    )
  ) +
  geom_point(size = 5) +
  scale_shape_manual(
    values = c(
      "Multiple-testing–corrected q-values <0.01" = 16,
      "Multiple-testing–corrected q-values <0.05" = 17,
      "Multiple-testing–corrected q-values <0.1"  = 15,
      "Not significant" = 4
    )
  ) +
  scale_color_gradientn(
    colours = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45"),
    values  = scales::rescale(c(0, 1, 5, 10)),
    name = "Variation explained\n(partial R-squared)"
  ) +
  labs(
    x = "Phenotype",
    y = "Visit",
    title = ""
  ) +
  guides(
    color = guide_colorbar(
      barheight = unit(4, "cm")
    ),
    shape = guide_legend(
      override.aes = list(size = 5)
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    legend.text = element_text(
      size = 14,
      face = "bold"
    ),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    )
  )

BM_pic

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure2D.pdf",width = 12, height = 6)




#################################################################
# additional adjust for covariates
#################################################################





BM_all_technical_and_biologic <- bind_rows(BM_hold) %>%
  filter(R2>0.00001) %>%
  filter(mainfocus ==term) %>%
  filter(term %in% attentionlist_BM ) %>% 
  group_by(visit)  %>%
  filter(covariate=="technical(depth and batch) and biologic(delivery_mode,feeding_type_28days,antibiotic_usage_180day") %>%
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




BM_all_technical_and_biologic$visit <- factor(BM_all_technical_and_biologic$visit,levels = c("BM1.5","BM6","BM12","BM24","BM36"))

BM_filter <- BM_all_technical_and_biologic %>% 
  mutate(newname = case_when( term == "trajectory" ~ "Trajectory",
                              term ==  "delivery_mode" ~ "Delivery mode",
                              term == "feeding_type_28days" ~ "Feeding type",
                              term == "BMI_mo" ~ "Prepregnancy BMI",
                              term == "preterm" ~ "Preterm",
                              term == "gender_kid" ~ "Gender",
                              term == "edu_fa" ~ "Education level of father",
                              term == "edu_mo" ~ "Education level of mother",
                              term == "family_income" ~ "Family income",
                              term == "gravida" ~ "Gravida",
                              term == "Age_mo" ~ "Age of mother",
                              term == "Age_fa" ~ "Age of father",
                              term == "TG" ~ "TG",
                              term == "UA.y" ~ "UA", 
                              term == "HDLC" ~ "HDLC",
                              term == "AST" ~ "AST",
                              term == "BUN" ~ "BUN",
                              term == "GDM_group2023" ~ "GDM",
                             
                              term =="antibiotic_usage_180day"~"Antibiotic use within 180 days",
                              term =="secondhand_smoke_exposure_42day"~"Second-hand smoke exposure",
                              term =="is_threegenerations_early_preg" ~"Grandparent in household",
                              term =="bedroom_floor_tile_stone_early_preg" ~"Tile/stone bedroom floor",
                              term =="has_siblings_early_preg" ~"Presence of siblings",
                              term =="coffee_freq_prepreg" ~"Coffee drinking frequency",
                              
                              TRUE ~ term
  ))


target <- c("Trajectory","Delivery mode", "Feeding type","Prepregnancy BMI","Preterm","Gender",
            "Education level of father","Education level of mother","Age of mother","Age of father",
            "HDLC","UA","TG","AST","BUN","GDM","Antibiotic use within 180 days","Second-hand smoke exposure",
            "Grandparent in household","Tile/stone bedroom floor","Presence of siblings","Coffee drinking frequency")



BM_filter <- BM_filter %>%  dplyr::select(newname,Df,SumOfSqs,R2,`F`,`Pr(>F)`,visit,covariate,q,Significance)

holdresult <- rbind(holdresult,BM_filter)



BM_filter %>%
  
  ggplot(aes(x = factor(newname,level = target), 
             y = visit,color =R2*100,
             shape = Significance)) +
  geom_point(size=5) +
  scale_shape_manual(values = c(
    "Multiple-testing–corrected q-values <0.01" = 16,
    "Multiple-testing–corrected q-values <0.05" = 17,
    "Multiple-testing–corrected q-values <0.1"  = 15,
    "Not significant" = 4
  ))+
  ylab("Visit") +
  xlab("Phenotype") +
  scale_color_gradientn(
    colours = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45"),
    values  = scales::rescale(c(0, 1, 5, 10)),  # expand low end
    name = "Variation explained \n (partial R-squared)"
  )+
  ggtitle(label = "")+
  theme (plot.title = element_text (size = 14, face = "bold" ))+
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=18,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        
        plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.text.x = element_text(size=18,face="bold",angle = 45, hjust = 1))+
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


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS4B_technical_and_biological.pdf",width = 12, height = 6)


holdresult <- rbind(holdresult,BM_filter) %>%
  dplyr::select(newname,Df,SumOfSqs,R2,`F`,`Pr(>F)`,visit,covariate,q,Significance)

write.table(holdresult,
            file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 7 BM_all.tsv",fileEncoding = "GBK",row.names = F,
            col.names = T,quote = F,sep = "\t"
            )
