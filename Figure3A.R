library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)


# manually set species order
speciesorder <- c( "Klebsiella oxytoca",
                   "Enterococcus casseliflavus",
                   "Enterococcus faecium (t_SGB7967)",
                   "Enterococcus faecium (t_SGB7968)",
                   "Clostridium paraputrificum",
                   "Blautia hansenii (t_SGB4794)",
                   "Blautia hansenii (t_SGB4798)",
                   "Dysosmobacter welbionis",
                   "Blautia faecis",
                   "Faecalibacterium SGB15346",
                   "Agathobaculum butyriciproducens",
                   "Bacteroides thetaiotaomicron",
                   "Bacteroides fragilis (t_SGB1853)",
                   "Bacteroides fragilis (t_SGB1855)",
                   "Bifidobacterium bifidum",
                   "Roseburia faecis",
                   "Veillonella dispar",
                   "Veillonella atypica",
                   "Streptococcus parasanguinis",
                   "Ruminococcus gnavus",
                   "Haemophilus parainfluenzae (t_SGB9657)",
                   "Haemophilus parainfluenzae (t_SGB9663_group)",
                   "Phocaeicola vulgatus",
                   "Streptococcus salivarius (t_SGB8005)",
                   "Streptococcus salivarius (t_SGB8007_group)")

# mark species with more than one strains
DUP_spe <- c("Blautia hansenii","Haemophilus parainfluenzae","Streptococcus salivarius","Enterococcus faecium","Bacteroides fragilis")

# read Mean and SD of nGD calculated for strains with script <Figure3A_nGD_test.sh>
# results merged through rbind
merge_df <- read.table("Allmerge.MeanSD.tsv",header=T) %>%
  mutate(V1 = gsub("_"," ",gsub("^s__","",V1))) %>%
  mutate(Strain = gsub("__","_",Strain)) %>%
  mutate(X=case_when(V1 %in% DUP_spe ~ paste(V1," (",Strain,")",sep = ""),
                     .default = V1)) %>%
  dplyr::filter(X %in% speciesorder) %>%
  mutate(X = factor(X, levels=speciesorder))

# read wilcox test results of nGD comparison between groups calculated for strains with script <Figure3A_nGD_test.sh>
# results merged through rbind
stat.test <- read.table("Allmerge.stat.test.tsv",header = T, sep="\t") %>%
  mutate(across(c(10,12,15,16), as.numeric)) %>%
  mutate(xmin = case_when(xmin != 1 ~ xmin-0.2, .default = xmin)) %>%
  #mutate(xmax = case_when(xmax != 1 ~ xmax-0.4, .default = xmax)) %>%
  mutate(V1 = gsub("_"," ",gsub("^s__","",V1))) %>%
  mutate(Strain = gsub("__","_",Strain)) %>%
  mutate(X=case_when(V1 %in% DUP_spe ~ paste(V1," (",Strain,")",sep = ""),
                     .default = V1)) %>%
  dplyr::filter(X %in% speciesorder) %>%
  mutate(X = factor(X, levels=speciesorder)) %>%
  filter(!(group1=="None" | group2=="None"))

# modified x for plotting
idx=0
for (i in c(levels(factor(stat.test$X)))){
  stat.test[which(stat.test$X==i),"x"] <- stat.test[which(stat.test$X==i),"x"]+idx
  stat.test[which(stat.test$X==i),"xmin"] <- stat.test[which(stat.test$X==i),"xmin"]+idx
  stat.test[which(stat.test$X==i),"xmax"] <- stat.test[which(stat.test$X==i),"xmax"]+idx
  idx=idx+1
}

# show only comparison between inter-family and intra-family groups
# modify y.position and p.adj.signif for plotting
stat.test <- stat.test %>%
  filter(grepl("family",group1) & grepl("family",group2)) %>%
  mutate(y.position = 0.12) %>%
  mutate(p.adj.signif=gsub("\\*\\*\\*\\*","\\*\\*\\*",p.adj.signif))

# calculate ratio of inter-family and intra-family group
fig2 <- merge_df %>%
  filter(grepl("family",V12)) %>% 
  reshape2::dcast(X~V12,value.var = "Mean") %>%
  rename(inter=`Inter-family`,intra=`Intra-family`) %>%
  mutate(Y="ratio") %>%
  mutate(color=inter/intra) 

# plotting
fig <- merge_df %>%
  filter(V12!="None") %>%   # remove statistic for 'None' group
  merge(.,fig2,by="X") %>%
  ggplot(aes( x = X ,y = Mean, group=V12, color = V12)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.05, 
                position = position_dodge(1)) +
  geom_point(aes(color = V12), size=1, position = position_dodge(1)) +
  geom_tile(mapping=aes(x = X,y=0.14,fill=color),
            inherit.aes=F,show.legend = TRUE, 
            linewidth=0, width=0.5,height=0.008) +
  scale_fill_gradient2(name="Inter-family/Intra-family",
                       low = "blue", mid = "white", high = "orange",
                       midpoint=1,limits=c(0,6),breaks=c(0,1,6)) +
  scale_color_manual(values = c("#3ec1d3","#f38181","#d3d4d8","#b9bbdf")) +
  xlab(NULL) +
  ylab("nGD") +
  guides(color=guide_legend(title="Group"))+
  stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01, #vjust = 1,
    angle = 0, coord.flip=T, size = 4, hjust = 0.01
  )+
  scale_y_continuous(breaks = c(0,0.05,0.1)) +
  coord_flip() +
  theme_classic() +
  theme (plot.title = element_text (size = 14, face = "bold" ))+
  theme(axis.text=element_text(size=16,face="bold", color = "black"), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text = element_text(size=10,face="bold"), 
        legend.title = element_text(size=16,face="bold"), 
        axis.text.y =  element_text(face = "italic"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        axis.text.x = element_text(size=18,face="bold",angle = 45, hjust = 1))

# save plot
pdf("Fig3A.pdf", 12, 8)
print(fig)
dev.off()