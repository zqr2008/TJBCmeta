library(tidyverse)
library(dplyr)
library(ggplot2)

data <- readRDS("complete_phylo.rds")
sam_data <- data@sam_data

visit_level=c("MG1", "MG2","MG3", "BM1.5", "BM6", "BM12", "BM24", "BM36")
bar_df <- data.frame(sam_data[,c("visit","trajcluster")]) %>%
  group_by(trajcluster,visit) %>%
  summarise(Count=n()) %>%
  mutate(visit = factor(visit, levels=visit_level))

fig <- bar_df %>%
  ggplot(aes(visit, Count, fill=trajcluster)) +
  geom_bar(stat="identity",position="stack", color="grey", width=0.7,size=0.25)+
  geom_text(aes(label=Count),size=3, color = "#4d4545", 
            position=position_stack(vjust=0.8)) +
  scale_fill_manual(name = "Trajectory", values = c("#e95280","#23b1a5","#ffdd7e"),
                    labels = c("Trajectory 1", "Trajectory 2", "Trajectory 3"))+
  xlab("Visit") + 
  ylab("Sample Size") +
  theme_classic() +
  theme (plot.title = element_text (size = 14, face = "bold" ))+
  theme(axis.text=element_text(size=10, color = "black"), 
        axis.title=element_text(size=10), 
        legend.text = element_text(size=10), 
        legend.title = element_text(size=12), 
        axis.text.y =  element_text(),
        plot.title = element_text(size=20, hjust = 0.5),
        axis.text.x = element_text(size=10)
        )

pdf("FigureS1C.pdf", 6, 4)
print(fig)
dev.off()