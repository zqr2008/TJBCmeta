#!bin/Rscript

library(tidyverse)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)

stat.test <- df %>%
  group_by(Species) %>%
  wilcox_test(Mutation_Rate ~ Traj2) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(p.adj = round(p.adj,3))

stat.test <- stat.test %>%
  add_xy_position(x = "Traj2", step.increase = 0.3, dodge = 0.5)

fig <- df %>%
  ggplot(aes( x = Traj2 ,y = Mutation_Rate, group=Traj2, color = Traj2)) +
  geom_boxplot()+
  xlab("Trajectory") +
  ylab("Mutation Rate") +
  guides(color=guide_legend(title="Trajectory"))+
  scale_color_manual(name = "Trajectory", 
                     values = c("#e95280","#23b1a5","#ffdd7e"))+
  stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01, #vjust = 1,
    angle = 0, coord.flip=F, size = 4, hjust = 0.01
  )+
  theme_classic() +
  theme (plot.title = element_text (size = 14, face = "bold" ))+
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text = element_text(size=10,face="bold"), 
        legend.title = element_text(size=16,face="bold"), 
        axis.text.y =  element_text(face = "italic"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        axis.text.x = element_text(size=18,face="bold",angle = 45, hjust = 1))

pdf("Figure3E.pdf", 5, 4)
print(fig)
dev.off()