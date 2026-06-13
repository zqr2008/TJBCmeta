# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 02_diversity_permanova_Diversity_PERMANOVA
# Figure(s): FigureS3
# Table(s): Supplementary Table 4; Supplementary Table 5
# Purpose: Calculates microbiome distance/ordination summaries and PCoA figure outputs.
# Main input(s): newdata_phylo_NMrevision_0305.rds; dist_aitchison_pc_rda.rda
# Main output(s): Supplementary Table 4 pcoa_effect_sizes_separeatevisit.csv; FigureS3.pdf; Supplementary Table 5 report_axis_meaning.csv
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
library(patchwork)



data_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")
metadata <- as.data.frame(data_phylo@sam_data)
class(metadata) <- "data.frame"

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/dist_aitchison_pc_rda.rda")

metadatasub <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,visit,class_1)

dist_aitchison_pc_withvisit <- dist_aitchison_pc_rda %>%
  rownames_to_column("sampleid") %>%
  left_join(metadatasub,by="sampleid") %>%
  column_to_rownames("sampleid")



# extract visit levels
visit_levels <- levels(factor(dist_aitchison_pc_withvisit$visit))

pcoa_list <- list()
pc1_var <- list()
pc2_var <- list()
for (v in visit_levels) {
  
  samples_v <- rownames(
    dist_aitchison_pc_withvisit[
      dist_aitchison_pc_withvisit$visit == v, ]
  )
  
  dist_sub <- as.dist(
    dist_aitchison_pc_withvisit[samples_v, samples_v]
  )
  
  pcoa_res <- pcoa(dist_sub)
  
  
  # 1. Access the values data frame
  val_df <- pcoa_res$values
  
  # 2. Calculate % variance manually 
  # (using the 'Eigenvalues' column)
  var_explained <- (val_df$Eigenvalues / sum(val_df$Eigenvalues)) * 100
  
  
  # To see the variance for the first two axes:
  pc1_var[[v]] <- var_explained[1]
  pc2_var[[v]] <- var_explained[2]
  
  pcoa_list[[v]] <- data.frame(
    sampleid = samples_v,
    PC1 = pcoa_res$vectors[,1],
    PC2 = pcoa_res$vectors[,2],
    visit = v
  )%>%
    
    left_join(metadatasub,by ="sampleid")
}



storetest<-data.frame()
plot_list <- list()


elementslist <- c("MG1","MG2","MG3","BM1.5","BM6","BM12","BM24","BM36")
for (elements in elementslist) {
  
  subject_number <- as.data.frame(
    table(pcoa_list[[elements]]$class_1)
  )
  
  textinfo <- paste0(
    "\n Trajectory 1, n=", subject_number[1,2],"\n",
    " Trajectory 2, n=", subject_number[2,2],"\n",
    " Trajectory 3, n=", subject_number[3,2]
  )
  
  pic <- pcoa_list[[elements]] %>%
    ggplot(aes(x = PC1, y = PC2,
               fill = as.factor(class_1),
               color = as.factor(class_1))) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_classic() +
    scale_fill_manual(name = "Trajectory",
                      values = c("#e95280","#23b1a5","#E49B0F","grey")) +
    scale_color_manual(name = "Trajectory",
                       values = c("#e95280","#23b1a5","#E49B0F","grey")) +
    xlab(paste0("PCoA axis 1 \n(%variance explained: ",round(pc1_var[[elements]],2),"%)")) +
    ylab(paste0("PCoA axis 2 \n(%variance explained: ",round(pc2_var[[elements]],2),"%)")) +
    stat_ellipse(type = "norm", linetype = 2) +
    theme(
      axis.text = element_text(size=16, face="bold", color="black"),
      axis.title = element_text(size=16, face="bold"),
      legend.text = element_text(size=18, face="bold"),
      legend.title = element_text(size=18, face="bold"),
      strip.text = element_text(size=14, face="bold")
    ) +
    ggtitle(paste0(elements, textinfo))
  
  plot_list[[elements]] <- pic
  
  # statistics
  stat.test1 <- pcoa_list[[elements]] %>%
    t_test(reformulate("class_1", "PC1"),
           detailed = TRUE,
           p.adjust.method ="BH") %>%
    add_significance("p.adj") %>%
    mutate(visit = elements)
  
  

  
  stat.test2 <- pcoa_list[[elements]] %>%
    t_test(reformulate("class_1", "PC2"),
           detailed = TRUE,
           p.adjust.method ="BH") %>%
    add_significance("p.adj") %>%
    mutate(visit = elements)
  
  storetest <- rbind(storetest, stat.test1, stat.test2)
}


pcoa_list[["MG1"]] %>%
  filter(class_1==1 | class_1==3) %>%
  cohen.d(PC1 ~ class_1, data=. )



write.table(storetest,file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 4 pcoa_effect_sizes_separeatevisit.csv",sep=",",fileEncoding="GBK",row.names = F)



final_plot <- wrap_plots(plot_list, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_plot

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS3.pdf",width = 40,height =25,units = "cm")




metadata <- metadata %>%
  rownames_to_column("sampleid")

cor <- stats::cor
storecor <- data.frame()
x=1

elementslist <- c("BM1.5","BM6","BM12","BM24","BM36","MG1","MG2","MG3")

for (elements in elementslist) {
  
     reexamineaxis <-   pcoa_list[[elements]] %>%
       left_join(metadata,by="sampleid")
     
    
 
     joinpcoa <- remove_labels(reexamineaxis)
     
     for (i in (1:ncol(joinpcoa))){

         
          specificvisit <- joinpcoa %>%
            dplyr::filter(visit==elements)

         
         if(is.numeric(joinpcoa[,i])){
           
           storecor[x,1] <- colnames(specificvisit)[i]
           storecor[x,2] <- cor(specificvisit[,i],specificvisit$PC1)
           storecor[x,3] <- cor(specificvisit[,i],specificvisit$PC2)
           storecor[x,4] <- cor.test(specificvisit[,i],specificvisit$PC1)[["p.value"]]
           storecor[x,5] <- cor.test(specificvisit[,i],specificvisit$PC2)[["p.value"]]
           
           storecor[x,6] <- elements
           
           x=x+1
           
         }
         
         if(is.factor(specificvisit[,i])){
           
           storecor[x,1] <- colnames(specificvisit)[i]
           storecor[x,2] <- cor(as.numeric(specificvisit[,i]),specificvisit$PC1)
           storecor[x,3] <- cor(as.numeric(specificvisit[,i]),specificvisit$PC2)
           storecor[x,4] <- cor.test(as.numeric(specificvisit[,i]),specificvisit$PC1)[["p.value"]]
           storecor[x,5] <- cor.test(as.numeric(specificvisit[,i]),specificvisit$PC2)[["p.value"]]
           storecor[x,6] <- elements
           x=x+1

       }
     }
     

  }


colnames(storecor)[2] <- "axis1_r"
colnames(storecor)[3] <- "axis2_r"
colnames(storecor)[4] <- "axis1_p"
colnames(storecor)[5] <- "axis2_p"


report <- storecor %>%
  filter( abs(axis1_r) >0.10 | abs(axis2_r) >0.10) %>%
  filter(str_detect(V1,"PC")==FALSE) %>%
  filter(V1=="weight_mo_first_check" | V1=="weight_mo_pre_preg" |
          V1=="UA.x" | V1=="gravida" | V1=="parity") 


write.table(report,file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Supplementary Table 5 report_axis_meaning.csv",sep=",",fileEncoding="GBK",row.names = F)
