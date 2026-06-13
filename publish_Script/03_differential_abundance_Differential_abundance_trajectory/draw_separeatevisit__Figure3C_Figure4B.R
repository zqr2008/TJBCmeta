# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 03_differential_abundance_Differential_abundance_trajectory
# Figure(s): Figure3C; Figure4B
# Table(s): No direct submitted table matched from filename
# Purpose: Generates or supports submitted figure(s): Figure3C; Figure4B.
# Main input(s): newdata_phylo_NMrevision_0305.rds; ps_filtered.rds
# Main output(s): none detected
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(Maaslin2)
library(microViz)
library(tidyverse)
library(microbiome)
library(sjmisc)
library(tibble)
library(MetBrewer)
library(labelled)
library(sjmisc)
library(rstatix)
library(ggrepel)
library(ggpubr)
library(haven)

##########################################
#handle data
##########################################


setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")


data_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

#data_phylo <- aggregate_taxa(data_phylo, "Genus")


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




complete_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/ps_filtered.rds")
m6 <-read_sav("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/metadata/metadata/喂养数据250617/m6.sav")
m6 <- remove_labels(m6)


feeding6month <- m6  %>%
  dplyr::select(sjid_kid= IDchild,
                feeding_m6) 

feeding6month$sjid_kid <- as.character(feeding6month$sjid_kid)

otudf <- as.data.frame(complete_phylo@otu_table)

trajectory <- as.data.frame(complete_phylo@sam_data)

class(trajectory) <- "data.frame"
class(otudf) <- "data.frame"

otudf <- otudf %>%
  rownames_to_column("sampleid")



trajectory <- trajectory %>%
  rownames_to_column("sampleid") %>%
  left_join(otudf,by = "sampleid")   %>%
  mutate(trajcluster= as.factor(class_1))


trajectory <- trajectory %>%
  mutate(across(
    .cols = intersect(colnames(otudf), names(.)[sapply(., is.numeric)]),
    .fns = ~ log10(. + 1e-6)
  ))

trajectorymaternal <- trajectory %>%
  filter(str_detect(visit,"MG"))

trajectory <- trajectory  %>%
  left_join(feeding6month,by ="sjid_kid")

trajectory$visit <- factor(trajectory$visit,levels = c("BM1.5",
                                                       "BM6",
                                                       "BM12",
                                                       "BM24",
                                                       "BM36",
                                                       "MG1",
                                                       "MG2",
                                                       "MG3"))






##########################################
#plot de for offspring
##########################################


setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/differentiatingspecies_for_offspring/")

trajectory<- trajectory %>%
  mutate(comparisontype1 = case_when(trajcluster==1~1,
                                   TRUE~0)) %>%
  mutate(comparisontype2 = case_when(trajcluster==2~1,
                                     TRUE~0)) %>%
  mutate(comparisontype3 = case_when(trajcluster==3~1,
                                     TRUE~0))

# Create the vector
demother <- c("GGB9712_SGB15244", 
               "Lacrimispora_celerecrescens", 
               "Bacteroides_eggerthii", 
               "Clostridium_sp_AT4", 
               "Eubacterium_ramulus", 
               "GGB51884_SGB49168", 
               "Firmicutes_bacterium_AF16_15", 
               "GGB3746_SGB5089", 
               "GGB9593_SGB15015", 
               "Clostridiaceae_bacterium_Marseille_Q3526", 
               "Butyricimonas_faecalis", 
               "GGB9747_SGB15356", 
               "Odoribacter_splanchnicus", 
               "Clostridiales_bacterium_KLE1615", 
               "GGB6518_SGB9205", 
               "Flavonifractor_plautii", 
               "Coprobacter_fastidiosus", 
               "GGB33469_SGB15236", 
               "Parabacteroides_johnsonii", 
               "Blautia_glucerasea", 
               "Bacteroides_thetaiotaomicron", 
               "GGB13493_SGB15238", 
               "Clostridium_sp_AF27_2AA", 
               "GGB9708_SGB15233", 
               "GGB9634_SGB15093", 
               "Coprobacter_secundus", 
               "GGB3175_SGB4191", 
               "Clostridium_sp_AF36_4", 
               "Bacteroides_salyersiae", 
               "Clostridium_sp_AM49_4BH", 
               "GGB9345_SGB14311", 
               "Coprococcus_comes", 
               "Clostridiales_bacterium_KLE1615", 
               "Bifidobacterium_pseudocatenulatum", 
               "GGB9712_SGB15244", 
               "Ruminococcus_torques", 
               "Odoribacter_splanchnicus", 
               "Anaerobutyricum_soehngenii", 
               "Firmicutes_bacterium_AF16_15", 
               "GGB3746_SGB5089", 
               "Bacteroides_faecis", 
               "Ruminococcus_SGB4421", 
               "Bacteroides_uniformis", 
               "Bifidobacterium_pseudocatenulatum", 
               "GGB51884_SGB49168", 
               "Eggerthella_lenta", 
               "Clostridium_sp_AT4", 
               "Gemmiger_SGB15295", 
               "Phocaeicola_coprocola", 
               "Clostridiaceae_unclassified_SGB4771", 
               "Clostridiaceae_bacterium_Marseille_Q3526", 
               "GGB3619_SGB4894", 
               "Eubacterium_SGB6276", 
               "Blautia_stercoris", 
               "GGB13493_SGB15238", 
               "GGB9595_SGB15019", 
               "Prevotella_copri_clade_A")



demother <-c("Bifidobacterium_pseudocatenulatum")

holddf<-data.frame()

for (j in demother){
  
  stat.test <- trajectory %>%
    group_by(visit) %>%
    filter(str_detect(visit,"BM")) %>%
    wilcox_test(reformulate("comparisontype1", j ),p.adjust.method = "BH") %>%
    mutate(comparsion = "Trajectory 1 vs. non-Trajectory 1")
  
  stat.test2 <- trajectory %>%
    group_by(visit) %>%
    filter(str_detect(visit,"BM")) %>%
    wilcox_test(reformulate("comparisontype2", j ),p.adjust.method = "BH") %>%
    mutate(comparsion = "Trajectory 2 vs. non-Trajectory 2")
  
  
  stat.test3 <- trajectory %>%
    group_by(visit) %>%
    filter(str_detect(visit,"BM")) %>%
    wilcox_test(reformulate("comparisontype3", j),p.adjust.method = "BH") %>%
    mutate(comparsion = "Trajectory 3 vs. non-Trajectory 3")
  
  stat.test <-rbind(stat.test,stat.test2)
  stat.test <-rbind(stat.test,stat.test3) 
  stat.test$p_scientific <- format(stat.test$p, scientific = TRUE, digits = 3)
  
  holddf<-rbind(holddf,stat.test)
  stat.test <- stat.test %>%
    add_significance("p") %>%
    add_xy_position(x = "visit",step.increase = 0.05,fun = "median_iqr",
                    dodge = 0.5) %>%
    mutate(color = case_when(comparsion == "Trajectory 1 vs. non-Trajectory 1" ~ "1",
                             comparsion == "Trajectory 2 vs. non-Trajectory 2" ~ "2",
                             comparsion == "Trajectory 3 vs. non-Trajectory 3" ~ "3"))
  
  newname <- str_replace(j,"_"," ")
  
  pic <- trajectory %>% 
    filter(str_detect(visit,"BM")) %>%
    ggplot(aes(x = as.numeric(visit), y = get(j),
               color = trajcluster))+ 
    geom_smooth(method = "loess", se = TRUE)  +
    theme_classic() +
    xlab("Visit") +
    ylab(expression(log[10] ~ "Relative Abundance")) +
    scale_fill_manual(values = c("#e95280","#23b1a5","#e49b0f")) +
    scale_color_manual(values = c("#e95280","#23b1a5","#e49b0f"),
                       labels = c("1" = "Trajectory 1 vs. non-Trajectory 1", 
                                  "2" = "Trajectory 2 vs. non-Trajectory 2", 
                                  "3" = "Trajectory 3 vs. non-Trajectory 3")) +
    ggtitle(label = paste0(newname,"\n","Overall"))+
    theme (plot.title = element_text (size = 14, face = "bold" ))+
    guides(color=guide_legend(title="Trajectory"))+
    scale_x_discrete(limits=c("BM1.5","BM6","BM12","BM24","BM36"))+
    theme(axis.text=element_text(size=20,face="bold", color = "black"), 
          axis.title=element_text(size=24,face="bold"), 
          legend.text = element_text(size=14,face="bold"), 
          legend.title = element_text(size=18,face="bold"), 
          plot.title = element_text(size=26, face="italic", hjust = 0.5), 
          axis.text.x = element_text(size=20,face="bold",angle = 315, hjust = -0.01))  +
    stat_pvalue_manual(
      stat.test,  
      y.position = -4.5, 
      size = 11, 
      label = "p.signif",
      tip.length = 0, 
      color = "color", 
      hide.ns = TRUE,
      vjust = -0.5, 
      fontface = "bold",  # Make the p-value bold
      hjust = 0.5,  # Adjust horizontally for better spacing if needed
      position = position_dodge(width = 2)
    )
  ggsave(paste0(j,".pdf"),width =25,height = 20,units = "cm")
}





##########################################
#plot in vaginal delivery 
##########################################


for (j in c("Bifidobacterium_pseudocatenulatum")){
  

    
    stat.test <- trajectory %>%
      filter(delivery_mode=="Vaginal delivery")  %>%
      group_by(visit) %>%
      filter(str_detect(visit,"BM")) %>%
      wilcox_test(reformulate("comparisontype1", j ),p.adjust.method = "BH") %>%
      mutate(comparsion = "Trajectory 1 vs. non-Trajectory 1")
    
    stat.test2 <- trajectory %>%
      filter(delivery_mode=="Vaginal delivery")  %>%
      group_by(visit) %>%
      filter(str_detect(visit,"BM")) %>%
      wilcox_test(reformulate("comparisontype2", j ),p.adjust.method = "BH") %>%
      mutate(comparsion = "Trajectory 2 vs. non-Trajectory 2")
    
    
    stat.test3 <- trajectory %>%
      filter(delivery_mode=="Vaginal delivery")  %>%
      group_by(visit) %>%
      filter(str_detect(visit,"BM")) %>%
      wilcox_test(reformulate("comparisontype3", j),p.adjust.method = "BH") %>%
      mutate(comparsion = "Trajectory 3 vs. non-Trajectory 3")
    
    stat.test <-rbind(stat.test,stat.test2)
    stat.test <-rbind(stat.test,stat.test3) 
    stat.test$p_scientific <- format(stat.test$p, scientific = TRUE, digits = 3)
    
    holddf<-rbind(holddf,stat.test)
    
    stat.test <- stat.test %>%
      add_significance("p") %>%
      add_xy_position(x = "visit",step.increase = 0.05,fun = "median_iqr",
                      dodge = 0.5) %>%
      mutate(color = case_when(comparsion == "Trajectory 1 vs. non-Trajectory 1" ~ "1",
                               comparsion == "Trajectory 2 vs. non-Trajectory 2" ~ "2",
                               comparsion == "Trajectory 3 vs. non-Trajectory 3" ~ "3"))
    
    newname <- str_replace(j,"_"," ")
    
    pic <- trajectory %>% 
      filter(delivery_mode=="Vaginal delivery")  %>%
      filter(str_detect(visit,"BM")) %>%
      ggplot(aes(x = as.numeric(visit), y = get(j),
                 color = trajcluster))+ 
      geom_smooth(method = "loess", se = TRUE) +
      theme_classic() +
      xlab("Visit") +
      ylab(expression(log[10] ~ "Relative Abundance")) +
      scale_fill_manual(values = c("#e95280","#23b1a5","#e49b0f")) +
      scale_color_manual(values = c("#e95280","#23b1a5","#e49b0f"),
                         labels = c("1" = "Trajectory 1 vs. non-Trajectory 1", 
                                    "2" = "Trajectory 2 vs. non-Trajectory 2", 
                                    "3" = "Trajectory 3 vs. non-Trajectory 3")) +
      ggtitle(label = paste0(newname,"\n","Vaginal delivery"))+
      theme (plot.title = element_text (size = 14, face = "bold" ))+
      guides(color=guide_legend(title="Trajectory"))+
      scale_x_discrete(limits=c("BM1.5","BM6","BM12","BM24","BM36"))+
      theme(axis.text=element_text(size=20,face="bold", color = "black"), 
            axis.title=element_text(size=24,face="bold"), 
            legend.text = element_text(size=14,face="bold"), 
            legend.title = element_text(size=18,face="bold"), 
            plot.title = element_text(size=26, face="italic", hjust = 0.5), 
            axis.text.x = element_text(size=20,face="bold",angle = 315, hjust = -0.01))  +
      stat_pvalue_manual(
        stat.test,  
        y.position = -4.0, 
        size = 11, 
        label = "p.signif",
        tip.length = 0, 
        color = "color", 
        hide.ns = TRUE,
        vjust = -0.5, 
        fontface = "bold",  # Make the p-value bold
        hjust = 0.5,  # Adjust horizontally for better spacing if needed
        position = position_dodge(width = 2)
      )
  ggsave(paste0(j,"Vaginal delivery.pdf"),width =25,height = 20,units = "cm")
}






##########################################
#plot in exclusively breastfeeding
##########################################

for (j in c("Bifidobacterium_pseudocatenulatum")){
  
  
  
  stat.test <- trajectory %>%
    filter(feeding_m6==1)  %>%
    group_by(visit) %>%
    filter(str_detect(visit,"BM")) %>%
    wilcox_test(reformulate("comparisontype1", j ),p.adjust.method = "BH") %>%
    mutate(comparsion = "Trajectory 1 vs. non-Trajectory 1")
  
  stat.test2 <- trajectory %>%
    filter(feeding_m6==1)  %>%
    group_by(visit) %>%
    filter(str_detect(visit,"BM")) %>%
    wilcox_test(reformulate("comparisontype2", j ),p.adjust.method = "BH") %>%
    mutate(comparsion = "Trajectory 2 vs. non-Trajectory 2")
  
  
  stat.test3 <- trajectory %>%
    filter(feeding_m6==1)  %>%
    group_by(visit) %>%
    filter(str_detect(visit,"BM")) %>%
    wilcox_test(reformulate("comparisontype3", j),p.adjust.method = "BH") %>%
    mutate(comparsion = "Trajectory 3 vs. non-Trajectory 3")
  
  stat.test <-rbind(stat.test,stat.test2)
  stat.test <-rbind(stat.test,stat.test3) 
  stat.test$p_scientific <- format(stat.test$p, scientific = TRUE, digits = 3)
  
  holddf<-rbind(holddf,stat.test)
  
  stat.test <- stat.test %>%
    add_significance("p") %>%
    add_xy_position(x = "visit",step.increase = 0.05,fun = "median_iqr",
                    dodge = 0.5) %>%
    mutate(color = case_when(comparsion == "Trajectory 1 vs. non-Trajectory 1" ~ "1",
                             comparsion == "Trajectory 2 vs. non-Trajectory 2" ~ "2",
                             comparsion == "Trajectory 3 vs. non-Trajectory 3" ~ "3"))
  
  newname <- str_replace(j,"_"," ")
  
  pic <- trajectory %>% 
    filter(feeding_m6==1)  %>%
    filter(str_detect(visit,"BM")) %>%
    ggplot(aes(x = as.numeric(visit), y = get(j),
               color = trajcluster))+ 
    geom_smooth(method = "loess", se = TRUE) +
    theme_classic() +
    xlab("Visit") +
    ylab(expression(log[10] ~ "Relative Abundance")) +
    scale_fill_manual(values = c("#e95280","#23b1a5","#e49b0f")) +
    scale_color_manual(values = c("#e95280","#23b1a5","#e49b0f"),
                       labels = c("1" = "Trajectory 1 vs. non-Trajectory 1", 
                                  "2" = "Trajectory 2 vs. non-Trajectory 2", 
                                  "3" = "Trajectory 3 vs. non-Trajectory 3")) +
    ggtitle(label = paste0(newname,"\n","exclusively breastfeeding"))+
    theme (plot.title = element_text (size = 14, face = "bold" ))+
    guides(color=guide_legend(title="Trajectory"))+
    scale_x_discrete(limits=c("BM1.5","BM6","BM12","BM24","BM36"))+
    theme(axis.text=element_text(size=20,face="bold", color = "black"), 
          axis.title=element_text(size=24,face="bold"), 
          legend.text = element_text(size=14,face="bold"), 
          legend.title = element_text(size=18,face="bold"), 
          plot.title = element_text(size=26, face="italic", hjust = 0.5), 
          axis.text.x = element_text(size=20,face="bold",angle = 315, hjust = -0.01))  +
    stat_pvalue_manual(
      stat.test,  
      y.position = -4.0, 
      size = 11, 
      label = "p.signif",
      tip.length = 0, 
      color = "color", 
      hide.ns = TRUE,
      vjust = -0.5, 
      fontface = "bold",  # Make the p-value bold
      hjust = 0.5,  # Adjust horizontally for better spacing if needed
      position = position_dodge(width = 2)
    )
  ggsave(paste0(j," exclusively breastfeeding.pdf"),width =25,height = 20,units = "cm")
}

##########################################
#plot de for mother
##########################################
setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/differential_abundance_results/differentiatingspecies_for_mother/")
motherdiff <- c("Bifidobacterium_pseudocatenulatum","Phocaeicola_coprocola",
                "Bacteroides_uniformis","Bacteroides_thetaiotaomicron")

for (j in motherdiff){
  
  
  stat.test <- trajectorymaternal %>%
    group_by(visit) %>%
    wilcox_test(reformulate("trajcluster", j ),p.adjust.method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(p.adj = round(p.adj,3)) 
  
  stat.test <- stat.test %>%
    add_xy_position(x = "visit",step.increase = 0.08,fun = "median_iqr",
                    dodge = 0.8) %>%
    mutate(color = case_when(group1 == "1" & group2 == "2" ~ "1",
                             group1 == "2" & group2 == "3" ~ "2",
                             group1 == "1" & group2 == "3" ~ "3")) 
  
  newname <- str_replace(j,"_"," ")
  
  pic <- trajectorymaternal %>% 
    dplyr::select(all_of(j),trajcluster, visit) %>%
    group_by(trajcluster, visit) %>%
    mutate(
      traj_median = median(get(j), na.rm=TRUE),
      traj_q1 = quantile(get(j), 0.25, na.rm=TRUE),
      traj_q3 = quantile(get(j), 0.75, na.rm=TRUE)
    ) %>%
    ggplot(aes(x = visit, y = get(j),
               color = trajcluster))+
    geom_boxplot(outlier.size = 0.2, color = "black",  # box border
                 lwd = 0.5,                            # line width
                 aes(fill = trajcluster)) +
    theme_classic() +
    xlab("Visit") +
    ylab("Relative abundance(log 10 transformation)") +
    scale_fill_manual(values = c("#e95280","#23b1a5","#e49b0f")) +
    scale_color_manual(values = c("#e95280","#23b1a5","#e49b0f")) +
    ggtitle(label = paste0(newname))+
    theme (plot.title = element_text (size = 14, face = "bold" ))+
    guides(color=guide_legend(title="Trajectory"))+
    theme(axis.text=element_text(size=20,face="bold", color = "black"), 
          axis.title=element_text(size=24,face="bold"), 
          legend.text = element_text(size=20,face="bold"), 
          legend.title = element_text(size=20,face="bold"), 
          plot.title = element_text(size=26, face="italic", hjust = 0.5), 
          axis.text.x = element_text(size=20,face="bold",angle = 315, hjust = -0.01))  +
    stat_pvalue_manual(
      stat.test,  
      y.position = -4.5, size = 6 ,label = "p.adj.signif", tip.length = 0, color = "color",hide.ns = T
    )
  ggsave(paste0(j,".pdf"),width =20,height = 20,units = "cm")
}
