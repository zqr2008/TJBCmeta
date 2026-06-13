# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 03_differential_abundance_Differential_abundance_trajectory
# Figure(s): no direct figure
# Table(s): Supplementary Table 9
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

# ---- Source Rmd chunk: organize maaslin2 table be one table (permonva_revise.Rmd lines 1126-1173) ----
rm(list = ls())
setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/")
conflicted::conflicts_prefer(dplyr::filter)
base_dir <- "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final"
out_dir <- "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results"

final_table <- list.files(base_dir, full.names = TRUE) %>%
  map_dfr(function(dir_path) {
    
    file_path <- file.path(dir_path, "significant_results.tsv")
    if (!file.exists(file_path)) return(NULL)
    
    df <- read.delim(file_path)
    folder <- basename(dir_path)
    
    # ---- filter FIRST (your logic) ----
    df <- df %>%
      filter(str_detect(metadata, "comparsion")) 
    
    # ---- annotate ----
    df %>%
      mutate(
        type = ifelse(str_detect(folder, "^mother"), "mother", "offspring"),
        taxonomy = ifelse(str_detect(folder, "Genus"), "Genus", "Species"),
        covariate_setting = case_when(
          str_detect(folder, "addagemo") ~ "base model+maternal age",
          str_detect(folder, "addantibiotic_latepregnancy") ~ "base model+late-pregnancy antibiotic exposure",
          str_detect(folder, "addantibiotic_midpregnancy") ~ "base model+mid-pregnancy antibiotic exposure",
          str_detect(folder, "addgestionalage") ~ "base model+gestational age",
          str_detect(folder, "addBMI_mo") ~ "base model+pre-pregnancy BMI",
          str_detect(folder, "addgender_kid") ~ "base model+infant gender",
          str_detect(folder, "addpreterm") ~ "base model+preterm",
          str_detect(folder, "addantibiotic_usage_180day") ~ "base model+antibiotic exposure during the first 180 days of life",
          TRUE ~ "only base model"
        ),
        comparsion = str_extract(folder, "comparsion\\d+"),
        source_folder = folder
      )
    
  })


write.table(final_table,
            file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 9 allde.csv",fileEncoding = "GBK",row.names = F,
            col.names = T,quote = F,sep = ","
            )
