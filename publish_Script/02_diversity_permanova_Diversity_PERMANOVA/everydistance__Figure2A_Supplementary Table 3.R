# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 02_diversity_permanova_Diversity_PERMANOVA
# Figure(s): Figure2A
# Table(s): Supplementary Table 3
# Purpose: Calculates microbiome distance/ordination summaries and PCoA figure outputs.
# Main input(s): newdata_phylo_NMrevision_0305.rds; dist_aitchison_pc_rda.rda; pcoa_res.rda; df_pcoa.rda; plus 1 more
# Main output(s): dist_aitchison_pc_rda.rda; species_scores.rda; pcoa_res.rda; df_pcoa.rda; plus 1 more
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(lme4)
library(vegan)
library(phyloseq)
library(zCompositions)   # cmultRepl
library(compositions)    # clr
library(rstatix)
library(ggpubr)
library(patchwork)
library(labelled)
library(ape)
library(ggrepel)

conflicted::conflicts_prefer(rstatix::filter)
conflicted::conflicts_prefer(compositions::dist)


data_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")
metadata <- as.data.frame(data_phylo@sam_data)
class(metadata) <- "data.frame"


X_pc <- as.matrix(data_phylo@otu_table) + 1e-6    # or +1, but smaller often better
X_clr_pc <- t(apply(X_pc, 1, clr))
dist_aitchison_pc <- dist(X_clr_pc, method="euclidean")


dist_aitchison_pc_rda <- as.data.frame(as.matrix(dist_aitchison_pc))
save(dist_aitchison_pc_rda,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/dist_aitchison_pc_rda.rda")


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/dist_aitchison_pc_rda.rda")
pcoa_res <- pcoa(dist_aitchison_pc_rda)   # PCoA on Aitchison dist (same as PCA on CLR)

# 1. Access the values data frame
val_df <- pcoa_res$values

fit <- envfit(pcoa_res$vectors[,1:2], X_clr_pc, permutations = 0)


species_scores <- as.data.frame(scores(fit, display = "vectors"))

save(species_scores,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/species_scores.rda")

# 2. Calculate % variance manually 
# (using the 'Eigenvalues' column)
var_explained <- (val_df$Eigenvalues / sum(val_df$Eigenvalues)) * 100

pc1_var <- var_explained[1] #21.19072
pc2_var <- var_explained[2] #5.907813

save(pcoa_res,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/pcoa_res.rda")


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/pcoa_res.rda")

df_pcoa <- data.frame(
  SampleID = rownames(pcoa_res[["vectors"]]),
  Axis1 = pcoa_res$vectors[,1],
  Axis2 = pcoa_res$vectors[,2]
)


save(df_pcoa,file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/df_pcoa.rda")


#####################################################
#start
#####################################################
conflicted::conflicts_prefer(compositions::dist)
data_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")
metadata <- as.data.frame(data_phylo@sam_data)
class(metadata) <- "data.frame"
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/df_pcoa.rda")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/species_scores.rda")


# add magnitude (strength of association)
species_scores$length <- sqrt(species_scores$Axis.1^2 + species_scores$Axis.2^2)

top_pc1_pos <- species_scores[order(-species_scores$Axis.1), ][1:5, ]
top_pc2_pos <- species_scores[order(-species_scores$Axis.2), ][1:5, ]

top_pc1_neg <- species_scores[order(species_scores$Axis.1), ][1:2, ]
top_pc2_neg <- species_scores[order(species_scores$Axis.2), ][0, ]

top_species <- bind_rows(
  top_pc1_pos,
  top_pc1_neg,
  top_pc2_pos,
  top_pc2_neg
)


arrow_scale <- 70


top_species <- top_species %>%
  rownames_to_column("Species") %>%
  mutate(
    Axis1 = Axis.1 * arrow_scale,
    Axis2 = Axis.2 * arrow_scale
  )



joinpcoa <- metadata %>%
  rownames_to_column("SampleID") %>%
  inner_join(df_pcoa,by = "SampleID")



joinpcoa$visit <- factor(joinpcoa$visit,
                         levels = c("MG1",
                                    "MG2",
                                    "MG3",
                                    "BM1.5",
                                    "BM6",
                                    "BM12",
                                    "BM24",
                                    "BM36")) 




Figure2a <- joinpcoa %>% 
  ggplot(aes( x = Axis1 ,y = Axis2, fill = visit,color = visit)) +
  geom_point(alpha = 0.5,size = 0.3) +
  theme_classic()+
  scale_fill_manual(values = c("#e1f6f4","#7db9b3","#166678",
                               "#8B8000","#ffd2a5","#ffa8b8",
                               "#d988bc","#66429b")) +
  scale_color_manual(values =c("#e1f6f4","#7db9b3","#166678",
                               "#8B8000","#ffd2a5","#ffa8b8",
                               "#d988bc","#66429b")) +
  xlab("PCoA axis 1 (%variance explained: 21.19%)") +
  ylab("PCoA axis 2 (%variance explained: 5.91%)") +
  stat_ellipse(linetype = 2,
               lwd = 1.0)+
  #coord_cartesian(ylim = c(-0.03,0.02)) +
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        plot.title = element_text(size=22, face="bold"),
        plot.subtitle=element_text(size=18,  color="black"))



Figure2a <- joinpcoa %>% 
  ggplot(aes(x = Axis1, y = Axis2, fill = visit, color = visit)) +
  geom_point(alpha = 0.5, size = 0.3) +
  

  geom_segment(data = top_species,
               aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.1, "cm")),
               color = "grey",
               linewidth = 0.2) +
  
  geom_text_repel(
    data = top_species,
    aes(x = Axis1, y = Axis2, label = Species),
    inherit.aes = FALSE,
    size = 6,
    fontface = "italic",
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  
  theme_classic() +
  scale_fill_manual(values = c("#e1f6f4","#7db9b3","#166678",
                               "#8B8000","#ffd2a5","#ffa8b8",
                               "#d988bc","#66429b")) +
  scale_color_manual(values = c("#e1f6f4","#7db9b3","#166678",
                                "#8B8000","#ffd2a5","#ffa8b8",
                                "#d988bc","#66429b")) +
  xlab("PCoA axis 1 (%variance explained: 21.19%)") +
  ylab("PCoA axis 2 (%variance explained: 5.91%)") +
  stat_ellipse(linetype = 2, lwd = 1.0) +
  theme(axis.text = element_text(size=16,face="bold", color="black"), 
        axis.title = element_text(size=16,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        plot.title = element_text(size=22, face="bold"),
        plot.subtitle = element_text(size=18, color="black"))



stat.test_1 <- joinpcoa %>%
  group_by(visit) %>%
  t_test(reformulate("class_1", "Axis1" ),detailed = TRUE,p.adjust.method ="BH") %>%
  add_significance("p.adj")   %>%
  filter(p.adj.signif != "ns")

stat.test_1 <- stat.test_1 %>%
  add_xy_position(x = "visit",step.increase = 0.07,fun = "median_iqr",
                  dodge = 0.8) %>%
  mutate(color = case_when(group1 == "1" & group2 == "2" ~ "1",
                           group1 == "2" & group2 == "3" ~ "2",
                           group1 == "1" & group2 == "3" ~ "3")) 




Figure2b_1 <- joinpcoa %>% 
  group_by(class_1,visit) %>%
  summarise(pcoa_mean = mean(Axis1),
            pcoa_sd = sd(Axis1)) %>%
  ggplot(aes(x = as.numeric(visit), y = pcoa_sd,color = as.factor(class_1), 
             fill = as.factor(class_1)))+
  geom_errorbar(mapping =  aes (ymax = pcoa_mean+pcoa_sd,
                                ymin = pcoa_mean-pcoa_sd),  
                width = 0.6,position = position_dodge(width = 0.9)) +
  geom_point(mapping = aes(y=pcoa_mean),
             position = position_dodge(width = 0.9),
             size = 4)+
  theme_classic()+
  xlab(paste0("")) +
  ylab("PCoA Axis1") +
  scale_fill_manual(values = c("#e95280","#23b1a5","#E49B0F","grey")) +
  scale_color_manual(values = c("#e95280","#23b1a5","#E49B0F","grey")) +
  theme(axis.text=element_text(size=16,face="bold", angle = 45, hjust = 1,color = "black"), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        legend.position="top",
        plot.title = element_text(size=22, face="bold"),
        plot.subtitle=element_text(size=18,  color="black"))+
 
  guides(fill="none")+
  guides(color=guide_legend(title="Trajectory",nrow=3))+
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("MG1","MG2","MG3","BM1.5","BM6","BM12","BM24","BM36"))+
  stat_pvalue_manual(
    stat.test_1, label = "p.adj.signif", tip.length = 0,
    color = "color",fill = "color",hide.ns = F, size= 6
  )




stat.test_2 <- joinpcoa %>%
  group_by(visit) %>%
  t_test(reformulate("class_1", "Axis2" ),detailed = TRUE,p.adjust.method ="BH") %>%
  add_significance("p.adj")   %>%
  filter(p.adj.signif != "ns")

stat.test_2  <- stat.test_2  %>%
  add_xy_position(x = "visit",step.increase = 0.08,fun = "median_iqr",
                  dodge = 0.6) %>%
  mutate(y.position=y.position+18)%>%
  mutate(color = case_when(group1 == "1" & group2 == "2" ~ "1",
                           group1 == "2" & group2 == "3" ~ "2",
                           group1 == "1" & group2 == "3" ~ "3")) 



Figure2b_2 <- joinpcoa %>% 
  group_by(class_1,visit) %>%
  summarise(pcoa_mean = mean(Axis1),
            pcoa_sd = sd(Axis1)) %>%
  ggplot(aes(x = as.numeric(visit), y = pcoa_sd,color = as.factor(class_1), 
             fill = as.factor(class_1)))+
  geom_errorbar(mapping =  aes (ymax = pcoa_mean+pcoa_sd,
                                ymin = pcoa_mean-pcoa_sd),  
                width = 0.5,position = position_dodge(width = 0.9)) +
  geom_point(mapping = aes(y=pcoa_mean),
             position = position_dodge(width = 0.9),
             size = 4)+
  theme_classic()+
  xlab(paste0("")) +
  ylab("PCoA Axis2") +
  scale_fill_manual(values = c("#e95280","#23b1a5","#E49B0F","grey")) +
  scale_color_manual(values = c("#e95280","#23b1a5","#E49B0F","grey"),
                                labels = c("1" = "Trajectory 1 vs. Trajectory 2", 
                                           "2" = "Trajectory 2 vs. Trajectory 3", 
                                           "3" = "Trajectory 1 vs. Trajectory 3")) +
  theme(axis.text=element_text(size=16,face="bold", angle = 45, hjust = 1,color = "black"), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text = element_text(size=18,face="bold"), 
        legend.title = element_text(size=18,face="bold"), 
        legend.position="top",
        plot.title = element_text(size=22, face="bold"),
        plot.subtitle=element_text(size=18,  color="black"))+
  
  guides(fill="none")+
  guides(color=guide_legend(title="Trajectory",nrow=3))+
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("MG1","MG2","MG3","BM1.5","BM6","BM12","BM24","BM36"))+
  stat_pvalue_manual(
    stat.test_2 , label = "p.adj.signif", tip.length = 0,
    color = "color",fill = "color",hide.ns = F, size= 6
  )


pic <- Figure2a/(Figure2b_1+Figure2b_2)

ggsave( plot = pic,filename = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure2A.pdf",width = 30,height =30,units = "cm")


stat.test <- as.data.frame(rbind(stat.test_1,stat.test_2)) %>%
  dplyr::select(1:17)
  
  
stat.test <- as.data.frame(lapply(stat.test, function(x) {
  if (is.list(x)) unlist(x) else x
}))
write.table(stat.test, file ="C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 3 pcoaAxiscompare.csv",
            sep=",",fileEncoding = "GBK",row.names = F,col.names = T,quote = F)
