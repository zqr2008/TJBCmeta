# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 04_transmission_and_strain_Transmission_strain_persistence
# Figure(s): FigureS7A
# Table(s): No direct submitted table matched from filename
# Purpose: Analyzes maternal-infant strain transmission patterns and related outputs.
# Main input(s): ps_filtered.rds; 3.mother2kid_share_genome.filter_0.99999_0.99.profile; instrain_transmission_overall_deliverymode.rda; instrain_transmission_deliverymode.rda
# Main output(s): instrain_transmission_deliverymode.rda; instrain_transmission_overall_deliverymode.rda; FigureS7A Instrain_transmission_dividedbydeliverymode.pdf; Instrain_transmission_byvisitanddeliverymode.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(readr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(ggh4x)
setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/")
complete_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/ps_filtered.rds")

#complete_phylo <- aggregate_taxa(complete_phylo, "Genus")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"


metadata_instrain <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit,delivery_mode) %>%
  mutate(Samp1 = sampleid)




inStrain_C_genome <- read_delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/3.mother2kid_share_genome.filter_0.99999_0.99.profile", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE) %>%
  mutate(
    genus = str_extract(classification, "g__[^;]+") %>% str_remove("^g__"),
    species = str_extract(classification, "s__[^;]+") %>% str_remove("^s__")
  ) %>%
  inner_join(metadata_instrain,by ="Samp1") %>%
  mutate(
    pair_visit = paste(pmin(Visit1, Visit2), pmax(Visit1, Visit2), sep = "_")
  ) 




offspringvisit <- c("BM1.5","BM6","BM12","BM24","BM36")

instrain_transmission <- data.frame()
x = 1


for (speciesnumber in levels(factor(inStrain_C_genome$species))){
  for (mm in c(1:(length(offspringvisit)))){
    for (tt in levels(factor(inStrain_C_genome$delivery_mode))){
      
      
      specieswithnext <- inStrain_C_genome %>%
        filter(str_detect(pair_visit,offspringvisit[mm]) & str_detect(pair_visit,"MG3") ) %>%
        filter(delivery_mode == tt) %>%
        filter(species == speciesnumber)
      
      
      specieswithmother <- metadata %>%
        mutate(type = case_when(str_detect(visit,"MG")~"mother",
                                str_detect(visit,"BM")~"offspring"
        )) %>%
        filter( str_detect(visit,"MG3") | str_detect(visit,offspringvisit[mm])) %>%
        filter(delivery_mode == tt) %>%
        group_by(familyid) %>%
        summarise(
          has_mother = any(type == "mother"),
          has_offspring = any(type == "offspring"),
          .groups = "drop"
        )   %>%
        filter(has_mother & has_offspring) %>%     # require both mother and offspring non-zero
        summarise(family_count = n())
      
      
      
      
      instrain_transmission[x,1] <- speciesnumber
      instrain_transmission[x,2] <- offspringvisit[mm]
      instrain_transmission[x,3] <- "MG3"
      instrain_transmission[x,4] <- tt
      instrain_transmission[x,5] <- dim(specieswithnext)[1]
      instrain_transmission[x,6] <- dim(specieswithnext)[1]/as.numeric(specieswithmother[1,1])
      instrain_transmission[x,7] <- as.numeric(specieswithmother[1,1])- dim(specieswithnext)[1] 
      
      print(paste0(offspringvisit[mm],"MG3",speciesnumber))
      x = x + 1
    }
  } 
}


save(instrain_transmission,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/instrain_transmission_deliverymode.rda")




instrain_transmission_overall <- data.frame()
x = 1


for (speciesnumber in levels(factor(inStrain_C_genome$species))){
  for (tt in levels(factor(inStrain_C_genome$delivery_mode))){
    
    
    specieswithnext <- inStrain_C_genome %>%
      filter(str_detect(pair_visit,"MG3") & str_detect(pair_visit,"BM") ) %>%
      filter(delivery_mode == tt) %>%
      filter(species == speciesnumber)
    
    
    specieswithmother <- metadata %>%
      filter(delivery_mode ==tt) %>%
      mutate(type = case_when(str_detect(visit,"MG")~"mother",
                              str_detect(visit,"BM")~"offspring"
      )) %>%
      filter( str_detect(visit,"MG3") | str_detect(visit,"BM")) %>%

      group_by(familyid) %>%
      summarise(
        has_mother = any(type == "mother"),
        has_offspring = any(type == "offspring"),
        .groups = "drop"
      )   %>%
      filter(has_mother & has_offspring) %>%     # require both mother and offspring non-zero
      summarise(family_count = n())
    
    
    
    
    instrain_transmission_overall[x,1] <- speciesnumber
    instrain_transmission_overall[x,2] <- "overall"
    instrain_transmission_overall[x,3] <- "MG3"
    instrain_transmission_overall[x,4] <- tt
    instrain_transmission_overall[x,5] <- dim(specieswithnext)[1]
    instrain_transmission_overall[x,6] <- dim(specieswithnext)[1]/as.numeric(specieswithmother[1,1])
    instrain_transmission_overall[x,7] <- as.numeric(specieswithmother[1,1])- dim(specieswithnext)[1] 
    
    
    print(paste0("overall","MG3",speciesnumber))
    x = x + 1
  }
  
}



save(instrain_transmission_overall,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/instrain_transmission_overall_deliverymode.rda")




rm(list = ls())
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/instrain_transmission_overall_deliverymode.rda")


transmissibillity_merge <-instrain_transmission_overall %>%
  dplyr::rename(Visit = V2,
                materanlvisit = V3,
                Trajectory = V4,
                `Transmitted number` = V5,
                `Transmissibility` = V6,
                `Not transmitted number` = V7) %>%
  group_by(Visit,V1) %>%
  mutate(Transmissibility = as.numeric(Transmissibility)) %>%
  filter(sum(Transmissibility, na.rm = TRUE) > 0) %>%
  filter(sum(`Transmitted number`, na.rm = TRUE) > 50) %>%
  filter(!is.na(Transmissibility),
         Transmissibility >= 0, Transmissibility <= 1) %>%
  ungroup()

#write.table(instrain_transmission_overall, file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/instrain_transmission_overall.csv",
#            fileEncoding = "GBK",row.names = F,col.names = T,
#            quote = F,sep = ",")




storedf <- data.frame()
y = 1

for (species in levels(factor(transmissibillity_merge$V1))){
  for (visit in levels(factor(transmissibillity_merge$Visit))){
    
    
    subtable <- transmissibillity_merge %>%
      filter(V1 == species & Visit == visit)
    if(dim(subtable)[1]==0){
      next
    }
    
    

    
    
    mat <- matrix(
      c(
        subtable$`Transmitted number`,
        subtable$`Not transmitted number`
      ),
      nrow = 2,
      byrow = FALSE
    )
    
    rownames(mat) <- subtable$Trajectory
    colnames(mat) <- c("Transmitted", "Not_transmitted")
    
    mat
    
    res <- fisher.test(mat)
    
    fisher_table <- tibble(
      Comparison = paste(subtable$Trajectory, collapse = " vs "),
      odds_ratio = unname(res$estimate),
      ci_lower   = res$conf.int[1],
      ci_upper   = res$conf.int[2],
      p_value    = res$p.value
    ) %>%
      mutate(speciname = species)
    
    
    
    storedf <- rbind(storedf,fisher_table)
    
  }
}


signforstoredf <- storedf %>% 
  filter(p_value < 0.05)







# 1️⃣ Prepare data
df <- transmissibillity_merge %>%
  filter(V1 %in% signforstoredf$speciname) %>%
  rename(
    Species = V1,
    transmitted = `Transmitted number`,
    not_transmitted = `Not transmitted number`
  )

df_long <- df %>%
  pivot_longer(
    cols = c(transmitted, not_transmitted),
    names_to = "status",
    values_to = "count"
  ) %>%
  mutate(
    fill_group = ifelse(
      status == "transmitted",
      paste0(Trajectory),   # T1, T2, T3
      "NT"
    )
  )


df_long$fill_group <- factor(df_long$fill_group,
                             levels = c( "NT",
                                         "Vaginal delivery",
                                        "Caesarean section" ))
# 2️⃣ Plot
ggplot(df_long, aes(x = Species, y = count, fill = fill_group)) +
  
  geom_bar(stat = "identity", position = "stack", color = "black") +
  
  geom_text(
    data = df,
    aes(
      x = Species,
      y = transmitted + not_transmitted,
      label = paste0(round(transmitted/(transmitted + not_transmitted),3))
    ),
    hjust =1,
    size = 5,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~Trajectory, scales = "free_x") +
  
  
  scale_fill_manual(
    values = c(
      "Caesarean section" = "#D7191C",
      "Vaginal delivery" ="#2C7BB6",
      "NT"="white"

    ),
    labels = c(
      "T1" = "Transmitted (Traj 1)",
      "T2" = "Transmitted (Traj 2)",
      "T3" = "Transmitted (Traj 3)",
      "NT" = "Not transmitted"
    )
  ) +
  
  coord_flip() +
  
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic")
  ) +
  
  labs(
    x = "Species",
    y = "Number of dyads",
    fill = "",
    title = "Strain transmissibility by species and delivery mode"
  )



ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS7A Instrain_transmission_dividedbydeliverymode.pdf",width = 14, height =7)

  
  
  
  rm(list = ls())
  load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/instrain_transmission_deliverymode.rda")
  
  transmissibillity_merge <-instrain_transmission %>%
    dplyr::rename(Visit = V2,
                  materanlvisit = V3,
                  Trajectory = V4,
                  `Transmitted number` = V5,
                  `Transmissibility` = V6,
                  `Not transmitted number` = V7) %>%
    group_by(Visit,V1) %>%
    mutate(Transmissibility = as.numeric(Transmissibility)) %>%
    filter(sum(Transmissibility, na.rm = TRUE) > 0) %>%
    filter(sum(`Transmitted number`, na.rm = TRUE) > 5) %>%
    filter(!is.na(Transmissibility),
           Transmissibility >= 0, Transmissibility <= 1) %>%
    ungroup()
  
  
  
  # 1️⃣ Prepare data
  df <- transmissibillity_merge %>%
    rename(
      Species = V1,
      transmitted = `Transmitted number`,
      not_transmitted = `Not transmitted number`
    )
  
  df_long <- df %>%
    pivot_longer(
      cols = c(transmitted, not_transmitted),
      names_to = "status",
      values_to = "count"
    ) %>%
    mutate(
      fill_group = ifelse(
        status == "transmitted",
        paste0(Trajectory),   # T1, T2, T3
        "NT"
      )
    )
  
  
  df$Visit <- factor(df$Visit,
                     levels = c("BM1.5","BM6","BM12","BM24","BM36"))
  
  df_long$fill_group <- factor(df_long$fill_group,
                               levels = c( "NT",
                                           "Vaginal delivery",
                                           "Caesarean section" ))
  
  df_long$Visit <- factor(df_long$Visit,
                               levels = c( "BM1.5","BM6","BM12",
                                           "BM24","BM36" ))
  
  A<- df_long %>%
    filter(fill_group=="Vaginal delivery") %>%
    
    ggplot(aes(x = Species, y = count, fill = fill_group)) +
    
    geom_bar(stat = "identity", position = "stack", color = "black") +
    
    geom_text(
      data = df %>% filter(Trajectory=="Vaginal delivery"),
      aes(
        x = Species,
        y = transmitted + not_transmitted,
        label = paste0(round(transmitted/(transmitted + not_transmitted),3))
      ),
      hjust =1,
      size = 5,
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~Visit, scales = "free_x", nrow = 1) +
    
    scale_fill_manual(
      values = c(
        "Caesarean section" = "#D7191C",
        "Vaginal delivery" = "#2C7BB6",
        "NT" = "white"
      )
    ) +
    
    coord_flip() +
    
    theme_bw(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(face = "italic")
    ) +
    
    labs(
      x = "Species",
      y = "Number of dyads",
      fill = "",
      title = "Strain transmissibility by species and delivery mode"
    )
  
  
  
  
  B <- df_long %>%
    filter(fill_group=="Caesarean section") %>%
    
    ggplot(aes(x = Species, y = count, fill = fill_group)) +
    
    geom_bar(stat = "identity", position = "stack", color = "black") +
    
    geom_text(
      data = df %>% filter(Trajectory=="Caesarean section"),
      aes(
        x = Species,
        y = transmitted + not_transmitted,
        label = paste0(round(transmitted/(transmitted + not_transmitted),3))
      ),
      hjust =1,
      size = 5,
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~Visit, scales = "free_x", nrow = 1) +
    
    scale_fill_manual(
      values = c(
        "Caesarean section" = "#D7191C",
        "Vaginal delivery" = "#2C7BB6",
        "NT" = "white"
      )
    ) +
    
    coord_flip() +
    
    theme_bw(base_size = 19) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(face = "italic")
    ) +
    
    labs(
      x = "Species",
      y = "Number of dyads",
      fill = "",
      title = "Strain transmissibility by species and delivery mode"
    )
  
  
  library(patchwork)
  
  A/B
  ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/05_working_figures/exploratory_figures/Instrain_transmission_byvisitanddeliverymode.pdf",width = 17, height = 17)
  
