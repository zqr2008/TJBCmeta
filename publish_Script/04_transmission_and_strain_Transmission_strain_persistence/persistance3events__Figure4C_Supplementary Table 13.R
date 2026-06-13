# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 04_transmission_and_strain_Transmission_strain_persistence
# Figure(s): Figure4C
# Table(s): Supplementary Table 13
# Purpose: Analyzes maternal-infant strain transmission patterns and related outputs.
# Main input(s): newdata_phylo_NMrevision_0305.rds; 1.inStrain_C_genome.profile
# Main output(s): Supplementary Table 13 persistancedifferencestat.csv; Figure4C persistance3event.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
library(readr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(data.table)


complete_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit) %>%
  mutate(Samp1 = sampleid)



Persistance <- read_delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/1.inStrain_C_genome.profile", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE) %>%
  mutate(
    genus = str_extract(classification, "g__[^;]+") %>% str_remove("^g__"),
    species = str_extract(classification, "s__[^;]+") %>% str_remove("^s__")
  ) %>%
  left_join(metadata,by ="Samp1") %>%
  filter(Samp1 %in% metadata$sampleid & Samp2 %in% metadata$sampleid ) %>%
  mutate(
    pair_visit = paste(pmin(Visit1, Visit2), pmax(Visit1, Visit2), sep = "_")
  ) %>%
  filter(popANI>0.99999 & percent_compared >0.5)



newcal <- Persistance %>%
  group_by(Assembly,classification,FID) %>%
  filter(is.na(class_1)==FALSE) %>%
  mutate(n=n()) %>%
  filter(n>=3)



instrain_release_rate <- data.frame()
x = 1


for (speciesnumber in levels(factor(newcal$species))){
    for (tt in levels(factor(Persistance$class_1))){

      specieswithnext <- newcal %>%
        filter(class_1 == tt) %>%
        filter(species == speciesnumber) %>%
        distinct(FID)
      


      basnumber <- Persistance %>%
       filter(class_1 == tt ) %>%
        filter(species == speciesnumber)  %>%
        distinct(FID)
      
    
      instrain_release_rate[x,1] <- speciesnumber
      instrain_release_rate[x,2] <- tt
      instrain_release_rate[x,3] <- dim(specieswithnext)[1]
      instrain_release_rate[x,4] <- dim(specieswithnext)[1]/dim(basnumber)[1]
      instrain_release_rate[x,5] <- dim(basnumber)[1] - dim(specieswithnext)[1]
      
      print(paste0(speciesnumber))
      x = x + 1
    }
  } 




instrain_release_rate_df <- instrain_release_rate %>%
  dplyr::rename( Species = V1,
                Trajectory = V2,
                `Persistance >= 3 number` = V3,
                `Persistance` = V4,
                `Persistance < 3 number` = V5)  %>%
  mutate(Persistance = as.numeric(Persistance))



storedf <- data.frame()
y = 1

for (species in levels(factor(instrain_release_rate_df$Species))){

    #species <- "Bifidobacterium pseudocatenulatum"
    
    subtable <- instrain_release_rate_df %>%
      filter(Species == species )
    if(dim(subtable)[1]==0){
      next
    }
    
    
    subtable2 <- subtable %>%
      mutate(Trajectory = as.character(Trajectory)) %>%
      
      # Create long format
      pivot_longer(cols = c(`Persistance >= 3 number`, `Persistance < 3 number`),
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

    # 2. Apply Fisher's test with CI
    stat_table <- fish_input %>%
      rowwise() %>%
      mutate(
        fisher_res = list(
          fisher.test(
            matrix(
              c(
                `Focal_Persistance >= 3 number`,
                `Focal_Persistance < 3 number`,
                `Others_Persistance >= 3 number`,
                `Others_Persistance < 3 number`
      
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
        `Focal_Persistance >= 3 number`,
        `Focal_Persistance < 3 number`,
        `Others_Persistance >= 3 number`,
        `Others_Persistance < 3 number`
      ) %>%
      mutate(speciesname = species) %>%
      mutate(p.adj= p.adjust(p_value,method="BH"))
    
    
    storedf <- rbind(storedf,stat_table)
    
}



write.table(storedf, file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 13 persistancedifferencestat.csv",
            fileEncoding = "GBK",row.names = F,col.names = T,
            quote = F,sep = ",")













signforstoredf <- storedf %>% 
  filter(p.adj < 0.05)






# 1️⃣ Prepare data
df <- instrain_release_rate_df %>%
  filter(Species %in% unique(signforstoredf$speciesname))

df_long <- df %>%
  pivot_longer(
    cols = c(`Persistance < 3 number`, `Persistance >= 3 number`),
    names_to = "status",
    values_to = "count"
  ) %>%
  mutate(
    fill_group = ifelse(
      status == "Persistance >= 3 number",
      paste0("T", Trajectory),   # T1, T2, T3
      "Persistance < 3"
    )
  )



df_long$plabel <- " "

for (i in 1:nrow(signforstoredf)) {
  
  target_species <- as.character(signforstoredf[i, "speciesname"])
  
  target_traj <- str_extract(
    as.character(signforstoredf[i, "Comparison"]),
    "\\d+"
  )
  
  df_long <- df_long %>%
    mutate(
      plabel = ifelse(
        Species == target_species &
          Trajectory == target_traj &
          fill_group !="Persistance < 3"
        ,
        "*",
        plabel
      )
    )
}

# 2️⃣ Plot
ggplot(df_long, aes(x = Species, y = count, fill = fill_group)) +
  
  geom_bar(stat = "identity", position = "stack", color = "black") +
  
  geom_text(
    data = df,
    aes(
      x = Species,
      y = `Persistance >= 3 number` + `Persistance < 3 number`,
      label = paste0(round(`Persistance >= 3 number`/(`Persistance >= 3 number` + `Persistance < 3 number`),3))
    ),
    hjust =1,
    size = 7,
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(
      x = Species,
      y = count,
      label = plabel
    ),
    hjust =0.5,
    size = 14)+
  facet_wrap(~Trajectory, scales = "free_x") +
  
  
  scale_fill_manual(
    values = c(
      "T1" = "#e95280",
      "T2" = "#23b1a5",
      "T3" = "#E49B0F",
      "Persistance < 3" = "white"
    ),
    labels = c(
      "T1" = "Persistance > 3 (Traj 1)",
      "T2" = "Persistance > 3 (Traj 2)",
      "T3" = "Persistance > 3 (Traj 3)",
      "NT" = "Not transmitted"
    )
  ) +
  
  coord_flip() +
  
  theme_bw(base_size = 22) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic")
  ) +
  
  labs(
    x = "Species",
    y = "Number of offspring with Persistance detected",
    fill = "",
    title = "Strain Persistance proportion by trajectory"
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure4C persistance3event.pdf",width = 17, height = 5)




