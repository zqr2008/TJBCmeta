# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 03_differential_abundance_Differential_abundance_trajectory
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Processes MaAsLin2 differential-abundance results and heatmaps.
# Main input(s): newdata_phylo_NMrevision_0305.rds
# Main output(s): none detected
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(Maaslin2)
library(microViz)
library(tidyverse)
library(microbiome)

setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")



maternalco <- c("QC__RawBasesCountGp","realsampleday",
                "GDM_group2023","Age_mo","antibiotic_midpregnancy",
                "antibiotic_latepregnancy")

offspringco<-c("QC__RawBasesCountGp","BMI_mo",
               "GDM_group2023","gender_kid","preterm")




for (level in c("Genus","Species")){
  
data_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")
  
  
if (level == "Genus"){
data_phylo <- aggregate_taxa(data_phylo, level)
}


supple <- phyloseq_validate(data_phylo, remove_undetected = TRUE)
supple <- tax_fix(supple)

mother_ps <- supple %>%
  ps_filter(
     str_detect(visit,"MG")
  ) %>%
  ps_filter(
    is.na(class_1) == FALSE
  ) 

offspring_ps <- supple %>%
  ps_filter(
    str_detect(visit,"BM")
  ) %>%
  ps_filter(
    is.na(class_1) == FALSE
  ) 



#saveRDS(mother_ps,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/mother_ps.rds")
#saveRDS(offspring_ps,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/offspring_ps.rds")

input_data_mother <- as.data.frame(t(mother_ps@otu_table))
Metadata_mother <- as.data.frame((mother_ps@sam_data))
class(Metadata_mother) <- "data.frame"


input_data_offspring <- as.data.frame(t(offspring_ps@otu_table))
Metadata_offspring <- as.data.frame((offspring_ps@sam_data))
class(Metadata_offspring) <- "data.frame"
Metadata_offspring$antibiotic_usage_180day <- as.factor(Metadata_offspring$antibiotic_usage_180day )

Metadata_mother <- Metadata_mother %>%
  mutate(trajcluster = paste0("traj",class_1)) %>%
  mutate(comparsion1 = case_when(trajcluster=="traj1"~"traj1",
                                 TRUE ~"ref"),
         comparsion2 =case_when(trajcluster=="traj2"~"traj2",
                                TRUE ~"ref"),
         comparsion3 = case_when(trajcluster=="traj3"~"traj3",
                                 TRUE ~"ref"))


Metadata_offspring <- Metadata_offspring %>%
  mutate(trajcluster = paste0("traj",class_1)) %>%
  mutate(feeding_type_28days=str_replace(feeding_type_28days," ","_")) %>%
  mutate(comparsion1 = case_when(trajcluster=="traj1"~"traj1",
                                 TRUE ~"ref"),
         comparsion2 =case_when(trajcluster=="traj2"~"traj2",
                                TRUE ~"ref"),
         comparsion3 = case_when(trajcluster=="traj3"~"traj3",
                                 TRUE ~"ref"))


for (comparelevel in c("comparsion1","comparsion2","comparsion3")){

fit_data_mother = Maaslin2(input_data  = input_data_mother, 
                           input_metadata = Metadata_mother, 
                           output         = paste0("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/mother_output_addantibiotic_latepregnancy",level,"_",comparelevel), 
                           fixed_effects  = c("visit",comparelevel,"BMI_mo","TG","UA.x","QC__RawBasesCountGp","antibiotic_latepregnancy"),
                           random_effects = c('familyid',"Sequencing"),
                           reference      = c("visit,MG1",
                                              paste0(comparelevel,"ref")),
                           min_abundance = 0.05,
                           normalization = "TSS",
                           min_prevalence = 0.05,
                           correction = "BH",
                           plot_scatter = FALSE)



fit_data_offspring = Maaslin2(input_data = input_data_offspring, 
                              input_metadata = Metadata_offspring, 
                              output         =paste0("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/statistics_results_final/offspring_output_addpreterm",level,"_",comparelevel), 
                              fixed_effects  = c(comparelevel,
                                                 "delivery_mode","feeding_type_28days","visit","QC__RawBasesCountGp","preterm"),
                              random_effects = c('familyid',"Sequencing"),
                              reference      = c("visit,BM1.5","preterm,Term delivery",
                                                 "feeding_type_28days,exclusive_breastfeeding",
                                                 "delivery_mode,Vaginal delivery",
                                                 paste0(comparelevel,"ref")),
                              min_abundance = 0.01,
                              normalization = "TSS",
                              correction = "BH",
                              plot_scatter = FALSE)

}
}