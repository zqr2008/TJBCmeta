# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 04_transmission_and_strain_Transmission_strain_persistence
# Figure(s): FigureS6
# Table(s): No direct submitted table matched from filename
# Purpose: Analyzes maternal-infant strain transmission patterns and related outputs.
# Main input(s): ps_filtered.rds; data_index.rda; allvisittransmission.rda; shareratedf.rda; plus 1 more
# Main output(s): shareratedf.rda; acquisitionrate.rda; FigureS6 oursharerate.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(viridis)
library(rstatix)
library(ggpubr)
library(data.table)
library(patchwork)
###############################################
#otu table
###############################################
  
#setwd("C:/Users/zqr20/Documents/tjmeta/")

subcomplete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/ps_filtered.rds")
otudata <- as.data.frame(subcomplete_phylo@otu_table)
class(otudata) <- "data.frame"


metadata <- as.data.frame(subcomplete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid") 


joinfamily <- otudata %>%
  rownames_to_column("sampleid") %>%
  left_join(metadata,by = "sampleid") %>%
  dplyr::select(sampleid,class_1,familyid,visit,any_of(colnames(otudata)))


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/data_index.rda")

data_index <- data_index %>%
  filter(sampleid %in% rownames(subcomplete_phylo@sam_data))

richness <- data_index %>%
  group_by(visit) %>%
  t_test(data_richness~class_1)


shannon <- data_index %>%
  group_by(visit) %>%
  t_test(data_shannon~class_1)

evenness <- data_index %>%
  group_by(visit) %>%
  t_test(data_evenness~class_1)

indexde <- rbind(richness,shannon)
indexde <- rbind(indexde,evenness)



##################################
#handle transmission table
##################################

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/allvisittransmission.rda")

#######################################################
#speciesname %in% colnames(otudata) filter species KQ analysis(no filtering before straplan) 
#######################################################
TJ_1k <- allvisittransmission %>%
  mutate(species = str_remove(Species,"s__")) %>%
  dplyr::filter(species %in% colnames(otudata)) %>%
  group_by(SAMP1,SAMP2) %>%
  mutate(n = n()) %>%
  filter(n>5) %>%
  mutate(key = paste0(SAMP1,SAMP2)) %>%
  distinct(key,Strain,.keep_all = T)

 
TJ_1k_index <-  TJ_1k %>%
  distinct(SAMP1,SAMP2,.keep_all = T)


############################################
#different rate calculated
############################################

shareratedf <- data.frame()
TJ_1k <- as.data.table(TJ_1k)

print(Sys.time())
for (familybind in  (1:nrow(TJ_1k_index))){

  
  samp1 = as.character(TJ_1k_index[familybind,"SAMP1"])
  samp2 = as.character(TJ_1k_index[familybind,"SAMP2"])                    
  
  makespeciesunique <- TJ_1k %>%
    filter(SAMP1 == samp1 & SAMP2== samp2) 
  
  
  overallspecies <- joinfamily %>%
    filter(sampleid == samp1 | sampleid == samp2) %>%
    select_if(is.numeric)   
  
  overallspecies <- overallspecies[,!sapply(overallspecies, function(x) any(x <0.01 ))]
  
  strain_sharing_rate <- dim(makespeciesunique)[1]/ (dim(overallspecies)[2]-1) 
  shareratedf[familybind,1] <- strain_sharing_rate
  shareratedf[familybind,2] <- TJ_1k_index[familybind,"SAMP1"]
  shareratedf[familybind,3] <- TJ_1k_index[familybind,"SAMP2"]
  print(paste0(Sys.time(),"row number:",familybind))

}


save(shareratedf,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/shareratedf.rda")


acquisitionrate <- data.frame()

joinfamilyoffspring <- joinfamily %>%
  filter(str_detect(visit,"BM")) 


for (familybind in  (1:nrow(TJ_1k_index))){
  
  
  samp1 = as.character(TJ_1k_index[familybind,"SAMP1"])
  samp2 = as.character(TJ_1k_index[familybind,"SAMP2"])         
  
  makespeciesunique <- TJ_1k %>%
    filter(SAMP1 == samp1 & SAMP2== samp2) 
  
  offspringspecies <- joinfamilyoffspring %>%
    filter(sampleid == samp1 | sampleid == samp2) %>%
    select_if(is.numeric)   
  
  offspringspecies <- offspringspecies[,!sapply(offspringspecies, function(x) any(x <0.01))]
  
  acquisitionrate_rate <- dim(makespeciesunique)[1]/  (dim(offspringspecies)[2] -1)
  
  acquisitionrate[familybind,1] <- acquisitionrate_rate
  acquisitionrate[familybind,2] <- TJ_1k_index[familybind,"SAMP1"]
  acquisitionrate[familybind,3] <- TJ_1k_index[familybind,"SAMP2"]
  print(paste0(Sys.time(),"row number:",familybind))
}

save(acquisitionrate,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/acquisitionrate.rda")






rm(list = ls())
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/shareratedf.rda")

metaforplot  <- metadata %>%
  dplyr::select(sampleid:visit,class_1,delivery_mode) %>%
  dplyr::rename(SAMP1 = sampleid)







metaforplot2 <- metadata %>%
  dplyr::select(sampleid:visit,class_1) %>%
  dplyr::rename(SAMP2 = sampleid)

shareratedfplot <- shareratedf %>%
  left_join(metaforplot, by = "SAMP1") %>%
  left_join(metaforplot2, by = "SAMP2") %>%
  mutate(
    pair_visit = paste(pmin(visit.x, visit.y), pmax(visit.x, visit.y), sep = "_")
  )%>%
  mutate(
    pair_visit = factor(
      pair_visit,
      levels = c(
        paste0("BM1.5_", c("MG1", "MG2", "MG3")),
        paste0("BM6_",   c("MG1", "MG2", "MG3")),
        paste0("BM12_",  c("MG1", "MG2", "MG3")),
        paste0("BM24_",  c("MG1", "MG2", "MG3")),
        paste0("BM36_",  c("MG1", "MG2", "MG3"))
      )
    )
  )




visit_order <- c("BM1.5","BM6","BM12","BM24","BM36")

# base filtered data
df_main <- shareratedfplot %>%
  filter(str_detect(visit.x, "MG3")) %>%
  mutate(visit.y = factor(visit.y, levels = visit_order))

# global y range (for stable annotation placement)
y_min <- min(df_main$V1, na.rm = TRUE)
y_range <- diff(range(df_main$V1, na.rm = TRUE))
y_pos <- y_min - 0.08 * y_range


visit_order <- c("BM1.5","BM6","BM12","BM24","BM36")

# base filtered data
df_main <- shareratedfplot %>%
  filter(str_detect(visit.x, "MG3")) %>%
  mutate(visit.y = factor(visit.y, levels = visit_order))

# global y range (for stable annotation placement)
y_min <- min(df_main$V1, na.rm = TRUE)
y_range <- diff(range(df_main$V1, na.rm = TRUE))
y_pos <- y_min - 0.08 * y_range

df_combined <- df_main %>%
  mutate(
    visit_group = case_when(
      visit.y %in% c("BM1.5","BM6","BM12") ~ "BM1.5–BM12",
      TRUE ~ as.character(visit.y)
    ),
    visit_group = factor(visit_group, levels = c("BM1.5–BM12","BM24","BM36"))
  )

df_stat2 <- df_combined %>%
  group_by(visit_group) %>%
  summarise(
    median = median(V1, na.rm = TRUE),
    Q1 = quantile(V1, 0.25, na.rm = TRUE),
    Q3 = quantile(V1, 0.75, na.rm = TRUE),
    label = sprintf("M=%.2f\n%.2f–%.2f", median, Q1, Q3),
    .groups = "drop"
  ) %>%
  mutate(y_pos = y_pos)

p2 <- ggplot(df_combined, aes(x = visit_group, y = V1, fill = visit_group)) +
  geom_boxplot() +
  geom_text(
    data = df_stat2,
    aes(x = visit_group, y = y_pos, label = label),
    size = 4,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Visits of offspring", y = "Sharing rate", fill = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


dodge_width <- 0.85

df_stat3 <- df_main %>%
  group_by(visit.y, delivery_mode) %>%
  summarise(
    median = median(V1, na.rm = TRUE),
    Q1 = quantile(V1, 0.25, na.rm = TRUE),
    Q3 = quantile(V1, 0.75, na.rm = TRUE),
    label = sprintf("M=%.2f\n%.2f–%.2f", median, Q1, Q3),
    .groups = "drop"
  ) %>%
  mutate(y_pos = y_pos)

p3 <- ggplot(df_main, aes(x = visit.y, y = V1, fill = delivery_mode)) + 
  geom_boxplot(position = position_dodge(width = dodge_width)) +
  geom_text(
    data = df_stat3,
    aes(x = visit.y, y = y_pos, label = label, group = delivery_mode),
    position = position_dodge(width = dodge_width),
    size = 2,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Visits of offspring", y = "Sharing rate", fill = "Delivery mode") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



(p1 / p2) | p3
ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS6 oursharerate.pdf",width = 14, height = 10)
  



comparestat <- shareratedfplot %>%
  filter(is.na(familyid.x)==FALSE) %>%
  #dplyr::filter(str_detect(pair_visit,"MG3")) %>%
  group_by(pair_visit) %>%
  t_test(V1~class_1.x)



stat.test <- shareratedfplot %>%
  filter(is.na(familyid.x)==FALSE) %>%
  group_by(pair_visit) %>%
  t_test(V1~class_1.x,detailed = TRUE) %>%
  add_significance("p.adj") 

stat.test <- stat.test %>%
  add_xy_position(x = "pair_visit",step.increase = 0.05,fun = "median_iqr",
                  dodge = 0.5) %>%
  mutate(color = case_when(group1 == "1" & group2 == "2" ~ "1",
                           group1 == "2" & group2 == "3" ~ "2",
                           group1 == "1" & group2 == "3" ~ "3")) 

shareratedfplot %>%
  filter(is.na(familyid.x)==FALSE) %>%
  #left_join(trajinfor,by = "V3") %>%
  #drop_na(trajcluster) %>%
  ggplot(aes(fill=as.factor(class_1.x), y=V1, x=pair_visit)) + 
  geom_boxplot()+
  xlab("visits of offspring")+
  ylab("sharing rate")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_viridis(discrete = T) +
  labs(fill="visits of mother") +
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=18,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.text.x = element_text(size=20,face="bold",angle = 315, hjust = -0.01))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0,
    color = "color",fill = "color",hide.ns = F, size= 6,y.position = -0.03
  )






load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/transmission_objects/acquisitionrate.rda")



acquisitionrateplot <- acquisitionrate %>%
  left_join(metaforplot, by = "SAMP1") %>%
  left_join(metaforplot2, by = "SAMP2") %>%
  mutate(
    pair_visit = paste(pmin(visit.x, visit.y), pmax(visit.x, visit.y), sep = "_")
  )%>%
  mutate(
    pair_visit = factor(
      pair_visit,
      levels = c(
        paste0("BM1.5_", c("MG1", "MG2", "MG3")),
        paste0("BM6_",   c("MG1", "MG2", "MG3")),
        paste0("BM12_",  c("MG1", "MG2", "MG3")),
        paste0("BM24_",  c("MG1", "MG2", "MG3")),
        paste0("BM36_",  c("MG1", "MG2", "MG3"))
      )
    )
  )


comparestat <- acquisitionrateplot %>%
  filter(is.na(familyid.x)==FALSE) %>%
  group_by(pair_visit) %>%
  t_test(V1~class_1.x)



stat.test <- acquisitionrateplot %>%
  filter(is.na(familyid.x)==FALSE) %>%
  group_by(pair_visit) %>%
  t_test(V1~class_1.x,detailed = TRUE) %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  add_xy_position(x = "pair_visit",step.increase = 0.05,fun = "median_iqr",
                  dodge = 0.5) %>%
  mutate(color = case_when(group1 == "1" & group2 == "2" ~ "1 v.s. 2",
                           group1 == "2" & group2 == "3" ~ "2 v.s. 3",
                           group1 == "1" & group2 == "3" ~ "1 v.s. 3")) 

acquisitionrateplot %>%
  filter(is.na(familyid.x)==FALSE) %>%
  ggplot(aes(fill=as.factor(class_1.x), y=V1, x=pair_visit)) + 
  geom_boxplot()+
  xlab("visits of offspring")+
  ylab("acquisition rate")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_viridis(discrete = T) +
  theme_classic()+
  labs(fill="visits of mother") +
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=18,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.text.x = element_text(size=20,face="bold",angle = 315, hjust = -0.01))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0,
    color = "color",fill = "color",hide.ns = F, size= 6,y.position = -0.03
  )







