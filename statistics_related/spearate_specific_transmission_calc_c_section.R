rm(list = ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggh4x)
library(patchwork)
library(sjmisc)
library(rstatix)
library(Maaslin2)
library(microViz)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/")

load("joinanalysis.rda")
load("DEsig.transmission_event.rda")
load("t_level.rda")

`%nin%` = Negate(`%in%`)

#############################################################################################################################################
#Some species were deliberately excluded because the analysis presented in Figure 3a indicated that they originated primarily from the mother.
#############################################################################################################################################
DEsig.transmission_event <- DEsig.transmission_event %>%
  filter(familyid.x %in% joinanalysis$familyid) %>%
  dplyr::filter(speciesname != "Blautia_hansenii" &
                speciesname != "Enterococcus_faecium" &
                speciesname != "Enterococcus_casseliflavus" &
                speciesname != "Klebsiella_oxytoca"  ) %>%
  #dplyr::filter(delivery_mode_kid.x == "Vaginal delivery")  %>%
  dplyr::filter(delivery_mode_kid.x == "Caesarean section") %>%
  mutate(speciesname = case_when(speciesname == "Bacteroides_fragilis" ~V1,
                                 TRUE ~ speciesname)) 


joinanalysis <- joinanalysis %>%
  left_join(t_level, by = "sampleid") %>%
  relocate(c(t__SGB1853,t__SGB1855),.before = familyid) %>%
  select(-Bacteroides_fragilis) %>%
  #dplyr::filter(delivery_mode_kid == "Vaginal delivery") %>%
  dplyr::filter(delivery_mode_kid == "Caesarean section") %>%
  dplyr::select(-c(Blautia_hansenii,Enterococcus_faecium,Enterococcus_casseliflavus,Klebsiella_oxytoca)) 



countn <- length(levels(factor(DEsig.transmission_event$speciesname)))




transmissibillity_overall <- data.frame()
x = 1


trajectorylevel <- c("1","2","3")

for (speciesnumber in c(2:(countn+2))){
  for (tt in trajectorylevel){
    
    
    #detected shared number of pairs for a particular species in MetaPhlAn4  
    #specieswithmother <- joinanalysis %>%
    #  filter(trajcluster == tt) %>%
    #  dplyr::filter( get(colnames(joinanalysis)[speciesnumber])!=0 ) %>%
    #  distinct(familyid,type,.keep_all = T) %>%
    #  group_by(familyid) %>%
    #  mutate(n = n()) %>%    
    #  filter(n >= 2) %>%
    #  ungroup()%>%
    #  distinct(familyid,.keep_all = T)
    
    
    specieswithmother <- joinanalysis %>%
      filter(trajcluster == tt) %>%
      dplyr::filter( get(colnames(joinanalysis)[speciesnumber])!=0 ) %>%
      filter(type == "mother") %>%
      #filter(type == "offspring") %>%
      distinct(familyid,.keep_all = T)
    
    
    #number of species transmission events for a species
    DEsig.transmission_event_species <-  DEsig.transmission_event %>%
      filter(speciesname == colnames(joinanalysis)[speciesnumber]) %>%
      distinct(familyid.x,.keep_all = T) %>%
      filter(trajcluster.x == tt) 
      
      
    
    
    transmissibillity_overall[x,1] <- colnames(joinanalysis)[speciesnumber]
    transmissibillity_overall[x,2] <- "Overall"
    transmissibillity_overall[x,3] <- tt
    transmissibillity_overall[x,4] <- dim(DEsig.transmission_event_species)[1]
    transmissibillity_overall[x,5] <- dim(DEsig.transmission_event_species)[1]/dim(specieswithmother)[1]
    transmissibillity_overall[x,6] <- dim(specieswithmother)[1] -  dim(DEsig.transmission_event_species)[1]
    
    x = x + 1
  }
}


transmissibillity_overall <- transmissibillity_overall %>%
  dplyr::rename(Visit = V2,
                Trajectory = V3,
                `Transmitted number` = V4,
                `Transmissibility` = V5,
                `Not transmitted number` = V6) 
 

transmissibillity <- data.frame()
x = 1

offspringvisit <- c("BM1.5","BM6","BM12","BM24","BM36")
mm <- c("MG1","MG2","MG3")


trajectorylevel <- c("1","2","3")

for (speciesnumber in c(2:(countn+2))){
  for (oo in offspringvisit){
    for (tt in trajectorylevel){
      
      #specieswithmother <- joinanalysis %>%
      #  filter(visit %in% mm | visit == oo) %>%
      #  filter(trajcluster == tt) %>%
      #  dplyr::filter( get(colnames(joinanalysis)[speciesnumber])!=0 ) %>%
      #  group_by(familyid) %>%
      #  mutate(n = n()) %>%    
      #  filter(n >= 2) %>%
      #  ungroup()%>%
      #  distinct(familyid,.keep_all = T)
      
  
      specieswithmother <- joinanalysis %>%
        filter(visit %in% mm | visit == oo) %>%
        filter(trajcluster == tt) %>%
        filter(type == "mother") %>%
        #filter(type == "offspring") %>%
        dplyr::filter( get(colnames(joinanalysis)[speciesnumber])!=0 ) %>%
        distinct(familyid,.keep_all = T)
      
      
      
      DEsig.transmission_event_species <-  DEsig.transmission_event %>%
        filter(visit.x %in% mm & visit.y == oo) %>% 
        filter(speciesname ==  colnames(joinanalysis)[speciesnumber]) %>%
        distinct(familyid.x,.keep_all = T) %>%
        filter(trajcluster.x == tt) 
        
      
      transmissibillity[x,1] <- colnames(joinanalysis)[speciesnumber]
      transmissibillity[x,2] <- oo
      transmissibillity[x,3] <- tt
      transmissibillity[x,4] <- dim(DEsig.transmission_event_species)[1]
      transmissibillity[x,5] <- dim(DEsig.transmission_event_species)[1]/dim(specieswithmother)[1]
      transmissibillity[x,6] <- dim(specieswithmother)[1] -  dim(DEsig.transmission_event_species)[1]
      
      x = x + 1
    }
  } 
}



transmissibillity <- transmissibillity %>%
  dplyr::rename(Visit = V2,
                Trajectory = V3,
                `Transmitted number` = V4,
                `Transmissibility` = V5,
                `Not transmitted number` = V6)


transmissibillity_merge <- rbind(transmissibillity_overall,transmissibillity) %>%
  group_by(V1) %>%
  arrange(Trajectory)

#############################################################################################################################################
#write to supplementary
#############################################################################################################################################
write.table(transmissibillity_merge,file = "transmissibillity_merge_c_section.tsv",fileEncoding = "GBK",
            row.names = F, col.names = T, sep = "\t", quote = F)


#############################################################################################################################################
#pairwise_fisher_test for transmissibillity 
#############################################################################################################################################
storedf <- data.frame()
y = 1

for (species in levels(factor(transmissibillity_merge$V1))){
  for (visit in levels(factor(transmissibillity_merge$Visit))){

    
    subtable <- transmissibillity_merge %>%
      filter(V1 == species & Visit == visit)
    
   
    xtable  <- cbind(as.matrix(subtable$`Transmitted number`),as.matrix(subtable$ `Not transmitted number`)) %>%
    pairwise_fisher_test() %>%
      mutate(group = paste0(group1,group2)) %>%
      filter(group != "row1row3") %>%
      mutate(group = case_when(group == "row1row2" ~"Trajectory 1 vs. Trajectory 2",
                               group == "row2row3" ~"Trajectory 3 vs. Trajectory 2")) %>%
      mutate(contrast =  NA) %>%
      mutate(newp = case_when(p < 0.05 ~ "*")) %>%
      mutate(speciesname = species) %>%
      mutate(visitvisit = visit) %>%
      mutate(newgroup = paste0(visitvisit,group)) 
      
    xtable[1,"contrast"] <- subtable[1,"Transmissibility"] - subtable[2,"Transmissibility"] 
    xtable[2,"contrast"] <- subtable[3,"Transmissibility"] - subtable[2,"Transmissibility"]   
    
    
    storedf <- rbind(storedf,xtable)
    
  
}
}



storedf <- storedf %>%
  mutate(sign = case_when(contrast > 0 ~ "Increase",
                          contrast < 0 ~ "Decrease")) %>%
  filter(is.na(sign) == FALSE)

storedf$newgroup <- factor(storedf$newgroup,
                           levels = c("OverallTrajectory 1 vs. Trajectory 2",
                                      "OverallTrajectory 3 vs. Trajectory 2",
                                      "BM1.5Trajectory 1 vs. Trajectory 2",
                                      "BM1.5Trajectory 3 vs. Trajectory 2",
                                      "BM6Trajectory 1 vs. Trajectory 2",
                                      "BM6Trajectory 3 vs. Trajectory 2",
                                      "BM12Trajectory 1 vs. Trajectory 2",
                                      "BM12Trajectory 3 vs. Trajectory 2",
                                      "BM24Trajectory 1 vs. Trajectory 2",
                                      "BM24Trajectory 3 vs. Trajectory 2",
                                      "BM36Trajectory 1 vs. Trajectory 2",
                                      "BM36Trajectory 3 vs. Trajectory 2"))


storedf$speciesname <-factor(storedf$speciesname)

#save(storedf,file = "storedf_vaginal.rda")
#write.table(storedf,file = "storedf_vaginal.tsv",fileEncoding = "GBK",
#            row.names = F, col.names = T, sep = "\t", quote = F)




save(storedf,file = "storedf_cection.rda")
write.table(storedf,file = "storedf_cection.tsv",fileEncoding = "GBK",
            row.names = F, col.names = T, sep = "\t", quote = F)
