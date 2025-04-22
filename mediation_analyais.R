rm(list = ls())
library(tidyverse)
library(flextable)
library(knitr)
library(mediation)
library(broom)
library(flextable)
library(factoextra)
library(readxl)
library(ggpubr)
library(MetBrewer)


`%ni%` <- Negate(`%in%`)
set.seed(1993)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/")
complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/complete_phylo.rds")
otudata <- as.data.frame(complete_phylo@otu_table)
class(otudata) <- "data.frame"

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid")


load("DEsig.transmission_event.rda") 


organizedmetadata <- metadata %>%
  dplyr::select(familyid, sampleid,visit,type,zwhz_5,zwhz_4,
                zwhz_3,zwhz_2,zwhz_1,feeding_type_28days,delivery_mode_kid,
                BMI_mo,gender_kid) 


organizedmetadata$visit[organizedmetadata$visit=="BM1.5"] <- "BM1_5"
DEsig.transmission_event$visit.y[DEsig.transmission_event$visit.y=="BM1.5"] <- "BM1_5"

organizedmetadata$feeding_type_28days[is.na(organizedmetadata$feeding_type_28days)] <- "unknown"

organizedmetadata$delivery_mode_kid[is.na(organizedmetadata$delivery_mode_kid)] <- "unknown"

significant_results_select_offspring <- read.delim("significant_results_select_offspring")


saveresult <- data.frame()

visitname <- c("BM1_5","BM6","BM12","BM24","BM36")

  
for (vv in c(1:length(visitname))){
    
    selected_organizedmetadata  <- organizedmetadata %>%
      dplyr::filter(visit == visitname[vv] ) 
    
    for (name in significant_results_select_offspring$feature){
      
      selected_DEsig.transmission_event <- DEsig.transmission_event %>%
        dplyr::filter(speciesname == name)  %>%
        dplyr::filter(visit.y == visitname[vv] )
      
      otudata_single <- otudata %>%
        dplyr::select(any_of(significant_results_select_offspring$feature)) %>%
        dplyr::select(name) %>%
        rownames_to_column("sampleid")  %>%
        mutate(transmissiondetect = case_when(sampleid %in% selected_DEsig.transmission_event$sampleid ~ 1,
                                              TRUE ~ 0))
    
      
      organizedmetadata_single <- otudata_single %>%
        right_join(selected_organizedmetadata,by = "sampleid") 
      
      

      organizedmetadata_single <- organizedmetadata_single[complete.cases(organizedmetadata_single), ]
      
      
      if(dim(organizedmetadata_single)[1] < 50) {
        next
      }
      
      
      med.fit <- glm(reformulate(paste0("delivery_mode_kid +  feeding_type_28days"),"transmissiondetect"),
                     data = organizedmetadata_single,family = binomial(link = "logit"))
      
    
      
      out.fit <- glm(reformulate(paste0("transmissiondetect + delivery_mode_kid + feeding_type_28days + 
                gender_kid"),name),
                     data = organizedmetadata_single)
      
      
      try(med.out <- mediate(med.fit, out.fit, treat = "delivery_mode_kid", mediator = "transmissiondetect",
                             robustSE = TRUE, sims = 500)
      )
      
      if(med.out[["d0.p"]] < 0.05){
        
        print(paste0("result is significant for ",name, "at" ,visitname[vv])) 
        print(summary(med.out))
        
        
        singleresult <- tidy(med.out) %>%
          mutate(visit = visitname[vv],
                 species = name,
                 samplesize =length(med.fit[["fitted.values"]]),
                 totaleffectp = med.out[["tau.p"]],
                 totaleffectcoef =med.out[["tau.coef"]],
                 prop =med.out[["n.avg"]])
        
        saveresult <- rbind(saveresult,singleresult)
        
        
        
      }
      else{
        print(paste0("result is not significant for ",name, " at " ,visitname[vv])) 
        print(name)
        
      }
    }
  }
