# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 02_diversity_permanova_Diversity_PERMANOVA
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Runs or organizes PERMANOVA association summaries.
# Main input(s): none detected
# Main output(s): none detected
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

# ---- Source Rmd chunk: other exposure MG (permonva_revise.Rmd lines 1003-1120) ----
attentionlist_exposure<- c(
  "has_siblings_early_preg",
  "is_threegenerations_early_preg",
  "current_pet_ownership_early_preg",
  "mother_pregnancy_secondhand_smoke_early_preg",
  "secondhand_smoke_exposure_42day",
  "bedroom_floor_tile_stone_early_preg",
  "gourd_vegetables_prepregnancy_consumed",
  "leafy_vegetables_prepregnancy_consumed",
  "freshwater_fish_midpregnancy_consumed",
  "sugar_free_non_carbonated_prepregnancy_consumed",
  "probiotic_drinks_prepregnancy_consumed",
  "shellfish_seafood_prepregnancy_consumed",
  "candy_chocolate_prepregnancy_consumed",
  "rice_prepregnancy_consumed"
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

MG_all_adjust_technical <- bind_rows(MG_hold) %>%
  filter(R2>0.00001) %>%
  filter(mainfocus ==term) %>%
  filter(term %in% attentionlist_exposure ) 

setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/permanova_results/final_permonva")


attentionlist_exposure<- c(
  "has_siblings_early_preg",
  "is_threegenerations_early_preg",
  "current_pet_ownership_early_preg",
  "mother_pregnancy_secondhand_smoke_early_preg",
  "secondhand_smoke_exposure_42day",
  "bedroom_floor_tile_stone_early_preg",
  "gourd_vegetables_prepregnancy_consumed",
  "leafy_vegetables_prepregnancy_consumed",
  "freshwater_fish_midpregnancy_consumed",
  "sugar_free_non_carbonated_prepregnancy_consumed",
  "probiotic_drinks_prepregnancy_consumed",
  "shellfish_seafood_prepregnancy_consumed",
  "candy_chocolate_prepregnancy_consumed",
  "rice_prepregnancy_consumed"
)

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
  filter(term %in% attentionlist_exposure ) 


explore <- rbind(BM_all_adjust_technical,MG_all_adjust_technical) %>%
  dplyr::select(1:9)
