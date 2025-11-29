  rm(list = ls())
  library(readr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(rstatix)
  library(ggh4x)
  setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/")
  complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
  
  #complete_phylo <- aggregate_taxa(complete_phylo, "Genus")
  
  metadata <- as.data.frame(complete_phylo@sam_data)
  class(metadata) <- "data.frame"
  
  
  metadata_instrain <- metadata %>%
    rownames_to_column("sampleid") %>%
    dplyr::select(sampleid,class_1,visit,delivery_mode) %>%
    mutate(Samp1 = sampleid)
  
  
  
  
  inStrain_C_genome <- read_delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/3.mother2kid_share_genome.filter_0.99999_0.99.profile", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE) %>%
    mutate(
      genus = str_extract(classification, "g__[^;]+") %>% str_remove("^g__"),
      species = str_extract(classification, "s__[^;]+") %>% str_remove("^s__")
    ) %>%
    left_join(metadata_instrain,by ="Samp1") %>%
    mutate(
      pair_visit = paste(pmin(Visit1, Visit2), pmax(Visit1, Visit2), sep = "_")
    ) 
  

  Bifidobacterium  <- inStrain_C_genome %>%
    filter(genus=="Bifidobacterium")
  
  Phocaeicola_vulgatus <- inStrain_C_genome %>%
    filter(species=="Phocaeicola vulgatus")
  
  
  save(Bifidobacterium,file = "C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium.rda")
  
  save(Phocaeicola_vulgatus,file = "C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Phocaeicola_vulgatus.rda")
  
  
  offspringvisit <- c("BM1.5","BM6","BM12","BM24","BM36")
  
  instrain_transmission <- data.frame()
  x = 1
  
  
  for (speciesnumber in levels(factor(inStrain_C_genome$species))){
    for (mm in c(1:(length(offspringvisit)))){
      for (tt in levels(factor(inStrain_C_genome$class_1))){
        
        
        specieswithnext <- inStrain_C_genome %>%
          filter(str_detect(pair_visit,offspringvisit[mm]) & str_detect(pair_visit,"MG3") ) %>%
          filter(class_1 == tt) %>%
          filter(species == speciesnumber)
        
        
        specieswithmother <- metadata %>%
          mutate(type = case_when(str_detect(visit,"MG")~"mother",
                                  str_detect(visit,"BM")~"offspring"
                                  )) %>%
          filter( str_detect(visit,"MG3") | str_detect(visit,offspringvisit[mm])) %>%
          filter(class_1 == tt) %>%
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
  
  
  save(instrain_transmission,file = "C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/instrain_transmission.rda")


load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/instrain_transmission.rda")
  

instrain_transmission_overall <- data.frame()
x = 1


for (speciesnumber in levels(factor(inStrain_C_genome$species))){
    for (tt in levels(factor(inStrain_C_genome$class_1))){
      
      
      specieswithnext <- inStrain_C_genome %>%
        filter(str_detect(pair_visit,"MG3") & str_detect(pair_visit,"BM") ) %>%
        filter(class_1 == tt) %>%
        filter(species == speciesnumber)
      
      
      specieswithmother <- metadata %>%
        mutate(type = case_when(str_detect(visit,"MG")~"mother",
                                str_detect(visit,"BM")~"offspring"
        )) %>%
        filter( str_detect(visit,"MG3") | str_detect(visit,"BM")) %>%
        filter(class_1 == tt) %>%
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



save(instrain_transmission_overall,file = "C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/instrain_transmission_overall.rda")



##################################################################################
#handle overall!
################################################################################
rm(list = ls())
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/instrain_transmission_overall.rda")

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

write.table(instrain_transmission_overall, file = "C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/instrain_transmission_overall.csv",
            fileEncoding = "GBK",row.names = F,col.names = T,
            quote = F,sep = ",")



storedf <- data.frame()
y = 1

for (species in levels(factor(transmissibillity_merge$V1))){
  for (visit in levels(factor(transmissibillity_merge$Visit))){
    
    
    subtable <- transmissibillity_merge %>%
      filter(V1 == species & Visit == visit)
    if(dim(subtable)[1]==0){
      next
    }
    
    
    subtable2 <- subtable %>%
      mutate(Trajectory = as.character(Trajectory)) %>%
      
      # Create long format
      pivot_longer(cols = c(`Transmitted number`, `Not transmitted number`),
                   names_to = "Type", values_to = "Count") %>%
      
      # For each trajectory, create comparison vs the others
      group_by(Type) %>%
      summarise(
        Comparison = c(
          "Trajectory 1 vs non-1",
          "Trajectory 2 vs non-2",
          "Trajectory 3 vs non-3"
        ),
        Focal = c(
          Count[Trajectory == "1"],
          Count[Trajectory == "2"],
          Count[Trajectory == "3"]
        ),
        Others = c(
          sum(Count[Trajectory != "1"]),
          sum(Count[Trajectory != "2"]),
          sum(Count[Trajectory != "3"])
        ),
        .groups = "drop"
      )
    
    
    # 1. reshape into wide format
    fish_input <- subtable2 %>%
      pivot_wider(
        names_from = Type,
        values_from = c(Focal, Others)
      )
    # Columns created:
    # Focal_Transmitted number
    # Others_Transmitted number
    # Focal_Not transmitted number
    # Others_Not transmitted number
    
    
    # 2. Apply Fisher's test with CI
    stat_table <- fish_input %>%
      rowwise() %>%
      mutate(
        fisher_res = list(
          fisher.test(
            matrix(
              c(
                `Focal_Transmitted number`,
                `Others_Transmitted number`,
                `Focal_Not transmitted number`,
                `Others_Not transmitted number`
              ),
              nrow = 2,
              byrow = TRUE
            ),
            alternative = "two.sided",
            conf.int = TRUE
          )
        ),
        p_value = fisher_res$p.value,
        odds_ratio = unname(fisher_res$estimate),
        ci_lower = fisher_res$conf.int[1],
        ci_upper = fisher_res$conf.int[2]
      ) %>%
      ungroup() %>%
      select(
        Comparison,
        odds_ratio,
        ci_lower,
        ci_upper,
        p_value,
        `Focal_Transmitted number`,
        `Others_Transmitted number`,
        `Focal_Not transmitted number`,
        `Others_Not transmitted number`
      ) %>%
      mutate(speciesname = species) %>%
      mutate(visitvisit = visit)  %>%
      mutate(p.adj= p.adjust(p_value,method="BH"))
    
  
    storedf <- rbind(storedf,stat_table)
    
  }
}


# storedf should contain:
# speciesname, group1, group2, estimate, conf.low, conf.high, p

# Prepare table
storedf2 <- storedf %>%
  mutate(
    contrast = Comparison,
    speciesname = factor(speciesname)
  ) %>%
  group_by(contrast) %>%
  arrange(odds_ratio, .by_group = TRUE) %>%
  mutate(speciesname = factor(speciesname, levels = speciesname)) %>%
  ungroup()

# Moderate expansion to avoid over-compression & overshoot
xmin_global <- min(storedf2$ci_lower) - 0.1
xmax_global <- max(storedf2$ci_upper) + 0.6


storedf2 %>%
ggplot( aes(x = odds_ratio, y = speciesname)) +
  
  # --- 95% CI error bars ---
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0.15,
    linewidth = 0.4,
    color = "black"
  ) +
  
  # --- Point Estimate ---
  geom_point(size = 1.3, color = "black") +
  
  # --- p-values ---
  geom_text(
    aes(label = paste0("p.adj=", signif(p.adj, 3)),
        x = pmax(ci_upper, odds_ratio) + 0.05),
    hjust = 0, size = 3
  ) +
  
  # --- Reference line at OR = 1 ---
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  # --- Dashed x-axis style ---
  theme(axis.line.x = element_line(color = "grey60", linetype = "dashed", linewidth = 0.4)) +
  
  # --- Facet by contrast, KEEP titles ---
  facet_wrap(~ contrast, ncol = 3, scales = "free_y") +
  
  # Labels
  labs(
    x = "Odds Ratio (95% CI)",
    y = "Transmission rate comparisons"
  ) +
  
  # --- Cell Host & Microbe-like theme ---
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text = element_text(size = 12, face = "bold"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.text.x = element_text(size = 10),
    
    legend.position = "none",
    panel.border = element_rect(color = "black", linewidth = 0.3)
  ) +
  
  # --- Clip points but avoid squeezing error bars ---
  coord_cartesian(
    #xlim = c(xmin_global, xmax_global),
    xlim = c(0, 8),
    clip = "on"
  ) +
  
  # --- Keep y-axis only for first facet column ---
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_discrete(),                           # 1st facet → keep y-axis
      scale_y_discrete(labels = NULL, breaks = NULL), # 2nd facet → no y-axis
      scale_y_discrete(labels = NULL, breaks = NULL)  # 3rd facet → no y-axis
    )
  ) + 
  theme(axis.line.x = element_line(linetype = "dashed"))



ggsave("transmission_spearate.pdf",width = 16, height = 6)

write.table(storedf2, file = "C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/transmissibillity_statistic.csv",
            fileEncoding = "GBK",row.names = F,col.names = T,
            quote = F,sep = ",")

