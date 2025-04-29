rm(list = ls())
library(Maaslin2)
library(microViz)
library(tidyverse)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta")
load("t_level.rda")

data_phylo <- readRDS("complete_phylo.rds")


supple <- phyloseq_validate(data_phylo, remove_undetected = TRUE)
supple <- tax_fix(supple)

mother_ps <- supple %>%
  ps_filter(
    type == "mother"
  ) %>%
  ps_filter(
    is.na(trajcluster) == FALSE
  )

offspring_ps <- supple %>%
  ps_filter(
    type == "offspring"
  ) %>%
  ps_filter(
    is.na(trajcluster) == FALSE
  )


input_data_mother <- as.data.frame(t(mother_ps@otu_table))
Metadata_mother <- as.data.frame((mother_ps@sam_data))
class(Metadata_mother) <- "data.frame"


input_data_offspring <- as.data.frame(t(offspring_ps@otu_table))
Metadata_offspring <- as.data.frame((offspring_ps@sam_data))
class(Metadata_offspring) <- "data.frame"


Metadata_mother <- Metadata_mother %>%
  mutate(trajcluster = paste0("traj",trajcluster))


Metadata_offspring <- Metadata_offspring %>%
  mutate(trajcluster = paste0("traj",trajcluster))


fit_data_mother = Maaslin2(input_data  = input_data_mother, 
                    input_metadata = Metadata_mother, 
                    normalization  = "NONE",
                    output         = "statistics_results_final/mother_output", 
                    fixed_effects  = c("visit","trajcluster","BMI_mo","TG","UA.x"),
                    random_effects = c('familyid'),
                    reference      = c("visit,MG1",
                                       "trajcluster,traj2"),
                    min_abundance = 0.05,
                    min_prevalence = 0.05,
                    correction = "BH",
                    plot_scatter = FALSE)



fit_data_offspring = Maaslin2(input_data = input_data_offspring, 
                    input_metadata = Metadata_offspring, 
                    normalization  = "NONE",
                    output         = "statistics_results_final/offspring_output", 
                    fixed_effects  = c("visit",
                                       "delivery_mode_kid", 
                                       "feeding_type_28days","BMI_mo","trajcluster"),
                    random_effects = c('sjid_kid'),
                    reference      = c("visit,BM1.5",
                                       "feeding_type_28days",
                                       "trajcluster,traj2"),
                    min_abundance = 0.05,
                    min_prevalence = 0.1,
                    correction = "BH",
                    plot_scatter = FALSE)
