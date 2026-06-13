# README — script.zip contents

Generated summary of each file in the archive. This README was auto-generated.

- Total files scanned: 48


---

## GT_GH_relationship_with_abundance.R

- Path: `/mnt/data/script_contents/GT_GH_relationship_with_abundance.R`

- Size: 5KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggpubr)
library(rstatix)
library(lme4)
library(lmerTest)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/spreadcazy.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
significant_results <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/statistics_results_final/offspring_output_species/significant_results.tsv")
metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,familyid, delivery_mode,feeding_type_28days)
```


---

## QCscript.R

- Path: `/mnt/data/script_contents/QCscript.R`

- Size: 10KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
library(tidyverse)
library(ggplot2)
library(readxl)
library(pheatmap)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/")
complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
meta_batch <- read_excel("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/forsubmit/meta_batch.xlsx")

#complete_phylo <- aggregate_taxa(complete_phylo, "Genus")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

expose <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)

expose$id_anon <- paste0("ID_", seq_len(nrow(expose)))
```


---

## Use_MG_predict_offspring_trajectory.R

- Path: `/mnt/data/script_contents/Use_MG_predict_offspring_trajectory.R`

- Size: 6KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(smotefamily)
library(ranger)
library(pROC)
library(dplyr)
library(caret)
library(sjmisc)
library(tibble)
library(tidyverse)

complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")



otudf <- as.data.frame(complete_phylo@otu_table)

trajectory <- as.data.frame(complete_phylo@sam_data)

class(trajectory) <- "data.frame"
class(otudf) <- "data.frame"
```


---

## abundance relationship with transmission.R

- Path: `/mnt/data/script_contents/abundance relationship with transmission.R`

- Size: 5KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggpubr)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/spreadcazy.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)

Bifidobacterium_pseudocatenulatum_list <- Bifidobacterium %>%
  filter(species == "Bifidobacterium pseudocatenulatum")
```


---

## bifidobacterium_ko.R

- Path: `/mnt/data/script_contents/bifidobacterium_ko.R`

- Size: 11KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
library(data.table)
library(tidyverse)
library(uwot)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(ggrepel)

my_colors <- brewer.pal(8, "Set2")  # 8 distinct colors


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/filtered_results/")
bifido <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/bifido.list", header=FALSE) %>%
  mutate(bifiname = str_split_fixed(V1,".fa",n=2)[,1])


complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
```


---

## caculateindex_new.R

- Path: `/mnt/data/script_contents/caculateindex_new.R`

- Size: 5KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(vegan)
library(phyloseq)
library(sjmisc)
library(picante)
library(ape)
library(igraph)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/")
load("unfiltercleantj.rda")
load("cleansampleidandid.rda")

filteredspecies <- cleantj
sampleidandid <-  cleansampleidandid
rotatesepcies <- filteredspecies %>%
  mutate(species = str_split_fixed(filteredspecies$clade_name,"s__", n = 2)[,2]) %>%
  dplyr::select(-c(clade_name,clade_taxid)) %>%
  column_to_rownames("species") %>%
  rotate_df() %>%
```


---

## cazy.R

- Path: `/mnt/data/script_contents/cazy.R`

- Size: 17KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggpubr)


###########################################################
#clustering heatmap for cazy 
###########################################################
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/spreadcazy_high_qual.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit,zwhz_1:zwhz_7)
```


---

## cazy_PV.R

- Path: `/mnt/data/script_contents/cazy_PV.R`

- Size: 4KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Phocaeicola_vulgatus.rda")

complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")




metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)

```


---

## cazy_bifidum_comparsions.R

- Path: `/mnt/data/script_contents/cazy_bifidum_comparsions.R`

- Size: 16KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium_same.rda")

complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")




metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)
```


---

## compareBPsize.R

- Path: `/mnt/data/script_contents/compareBPsize.R`

- Size: 11KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(viridis)   # better colors
library(ggpubr)
library(purrr)
library(patchwork)

combined_strain_comparison <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/new_result/combined_strain_comparison.tsv") %>%
  mutate(tlevel = str_split_fixed(Strain,"t__",n=2)[,2])

detailed_strain_comparison <- read.csv("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/new_result/detailed_strain_comparison.csv")


# 1. Total Original vs Trimmed (line plot with labels)
df_long_total <- detailed_strain_comparison %>%
  dplyr::select(SampleID, Total_Original, Total_Trimmed) %>%
  pivot_longer(-SampleID, names_to="Type", values_to="Count")

```


---

## compareBPsize_ourdata.R

- Path: `/mnt/data/script_contents/compareBPsize_ourdata.R`

- Size: 12KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(viridis)   # better colors
library(ggpubr)
library(purrr)
library(patchwork)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/new_result/")

final_10_benchmark <- read_excel("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/new_result/final 10 benchmark.xlsx") %>%
  mutate(SampleID = paste0(样品名称,".mp4"))

finalsample<-c("MYF0006847.mp4","MYF0054133.mp4","MYF0071181.mp4",
               "MYF0071417.mp4","MYF0071670.mp4","MYF0046187.mp4",
               "MYF0071466.mp4","MYF0073097.mp4")
#delte host rate too high, and resequence error exsit
ourdata_compare150bp100bp.all.t.summary <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/new_result/ourdata_compare150bp100bp.all.t.summary.tsv") %>%
  filter(SampleID %in% finalsample) 
```


---

## draw_separeatevisit.R

- Path: `/mnt/data/script_contents/draw_separeatevisit.R`

- Size: 19KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(Maaslin2)
library(microViz)
library(tidyverse)
library(microbiome)
library(sjmisc)
library(tibble)
library(MetBrewer)
library(labelled)

##########################################
#handle data
##########################################


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")


data_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

```


---

## everydistance.R

- Path: `/mnt/data/script_contents/everydistance.R`

- Size: 9KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
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

data_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
metadata <- as.data.frame(data_phylo@sam_data)
class(metadata) <- "data.frame"


X_pc <- as.matrix(data_phylo@otu_table) + 1e-6    # or +1, but smaller often better
X_clr_pc <- t(apply(X_pc, 1, clr))
```


---

## feeding_organize.R

- Path: `/mnt/data/script_contents/feeding_organize.R`

- Size: 7KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(phyloseq)
library(dplyr)
library(randomForest)
library(ggplot2)
library(microbiome)
library(microViz)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(haven)

m6 <-read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/m6.sav")
y1 <- read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/y1.sav")

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/who_m6.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/who_y1.rda")

complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

```


---

## figure1_bubblepermonva_new.R

- Path: `/mnt/data/script_contents/figure1_bubblepermonva_new.R`

- Size: 15KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(viridis)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/batchscript")



######################################################
#PREPROCESS CODE
######################################################
BM_hold <- list()
for (kk in c("BM1\\.5.*\\.rda","^BM6.*\\.rda$","^BM12.*\\.rda$","^BM24.*\\.rda$","^BM36.*\\.rda$")){
  
  myFiles <- list.files(pattern=kk)
  
  for (i in myFiles){
    load(i)  
```


---

## generateR_permonva.R

- Path: `/mnt/data/script_contents/generateR_permonva.R`

- Size: 2KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/batchscript")


writeLines('
rm(list=ls())
library(vegan)
library(ggplot2)
library(tidyverse)
library(phyloseq)

setwd("/storeData/mchri/zhaiqiangrong/zhaiqiangrong/fecalmeta/secondbatchTJBC/permonva")

load("dist_aitchison_pc_df.rda")

bray <- dist_aitchison_pc_df
data_phylo <- readRDS("ps_filtered.rds")

```


---

## heatmaptraj.R

- Path: `/mnt/data/script_contents/heatmaptraj.R`

- Size: 11KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
library(data.table)
library(ggplot2)
library(dplyr)
library(forcats)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(rstatix)
library(sjmisc)
library(tibble)
library(MetBrewer)
library(ggpubr)


df <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/offspring_interaction_significant_results.tsv")

# Clean column names if necessary
df <- df %>%
  mutate(feature = gsub("_", " ", feature),
         metadata = gsub("_", " ", metadata))
```


---

## high_quality_GH_GT.R

- Path: `/mnt/data/script_contents/high_quality_GH_GT.R`

- Size: 16KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium_same.rda")

complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")




metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit)
```


---

## instrain_releaserate.R

- Path: `/mnt/data/script_contents/instrain_releaserate.R`

- Size: 13KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
library(readr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rstatix)


complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

#complete_phylo <- aggregate_taxa(complete_phylo, "Genus")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,visit) %>%
  mutate(Samp1 = sampleid)

```


---

## instrain_transmission.R

- Path: `/mnt/data/script_contents/instrain_transmission.R`

- Size: 13KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
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
```


---

## microbiotaage.R

- Path: `/mnt/data/script_contents/microbiotaage.R`

- Size: 7KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(phyloseq)
library(dplyr)
library(randomForest)
library(ggplot2)
library(microbiome)
library(microViz)
library(tidyverse)
library(ranger)
library(ggpubr)
library(rstatix)
library(readxl)
library(lubridate)


sampletime <- read_excel("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/0723导出母子粪便样本_数据位置_20251017更新采样登记时间.xlsx") %>%
  dplyr::rename(samplingtime = 登记时间) %>%
  dplyr::select(sampleid,samplingtime)


```


---

## new_FigureS1.R

- Path: `/mnt/data/script_contents/new_FigureS1.R`

- Size: 2KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
library(tidyverse)
library(dplyr)
library(ggplot2)

data <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
sam_data <- as.data.frame(data@sam_data)
class(sam_data) <- "data.frame"

visit_level=c("MG1", "MG2","MG3", "BM1.5", "BM6", "BM12", "BM24", "BM36")
bar_df <- sam_data %>%
  dplyr::select(class_1, visit) %>%
  group_by(class_1, visit) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(visit) %>%
  mutate(Total_visit = sum(Count),
         Percent = Count / Total_visit * 100) %>%
  mutate(visit = factor(visit, levels = visit_level))


forstatsitic <- bar_df %>%
```


---

## new_maslin2.R

- Path: `/mnt/data/script_contents/new_maslin2.R`

- Size: 4KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(Maaslin2)
library(microViz)
library(tidyverse)
library(microbiome)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")


data_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

data_phylo <- aggregate_taxa(data_phylo, "Genus")



supple <- phyloseq_validate(data_phylo, remove_undetected = TRUE)
supple <- tax_fix(supple)

mother_ps <- supple %>%
  ps_filter(
```


---

## new_maslin2_motherheatmap.R

- Path: `/mnt/data/script_contents/new_maslin2_motherheatmap.R`

- Size: 6KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")
significant_results1 <- read.delim("statistics_results_final/mother_output_species_comparsion1/significant_results.tsv")
significant_results2 <- read.delim("statistics_results_final/mother_output_species_comparsion2/significant_results.tsv")
significant_results3 <- read.delim("statistics_results_final/mother_output_species_comparsion3/significant_results.tsv")

significant_results_select <- rbind(significant_results1,significant_results2)
significant_results_select <- rbind(significant_results_select,significant_results3) %>%
  filter(str_detect(metadata,"comparsion")) %>%
  filter(qval<0.1)


write.table(significant_results_select,file = "mother_speices_significant_results_select.csv",
            fileEncoding = "GBK",sep = ",",quote = F,col.names = T,row.names = F)
```


---

## new_maslin2_motherheatmap_three_comparsion.R

- Path: `/mnt/data/script_contents/new_maslin2_motherheatmap_three_comparsion.R`

- Size: 5KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")
significant_results <- read.delim("statistics_results_final/mother_output_species/significant_results.tsv")


significant_results_BMI <- significant_results %>%
  filter(metadata == "BMI_mo") %>%
  filter(qval < 0.1) %>%
  mutate(trajectory = paste0(metadata,value)) %>%
  pivot_wider(id_cols = feature,
              names_from = value,
              values_from = coef ) %>%
  mutate(`BMI status` = case_when(BMI_mo < 0 ~ "Negatively associated with mother's BMI",
                                  BMI_mo > 0 ~ "Positively associated with mother's BMI"))
```


---

## new_maslin2_offspringheatmap.R

- Path: `/mnt/data/script_contents/new_maslin2_offspringheatmap.R`

- Size: 5KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")

significant_results1 <- read.delim("statistics_results_final/offspring_output_species_comparsion1/significant_results.tsv")
significant_results2 <- read.delim("statistics_results_final/offspring_output_species_comparsion2/significant_results.tsv")
significant_results3 <- read.delim("statistics_results_final/offspring_output_species_comparsion3/significant_results.tsv")


significant_results_select <- rbind(significant_results1,significant_results2)
significant_results_select <- rbind(significant_results_select,significant_results3) %>%
  filter(str_detect(metadata,"comparsion")) %>%
  filter(qval<0.1)

write.table(significant_results_select,file = "offspring_speices_significant_results_select.csv",
```


---

## new_maslin2_offspringheatmap_three_comparsion.R

- Path: `/mnt/data/script_contents/new_maslin2_offspringheatmap_three_comparsion.R`

- Size: 4KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")

significant_results <- read.delim("statistics_results_final/offspring_output_species/significant_results.tsv")





significant_results_select <- significant_results %>%
  filter(metadata == "trajcluster") %>%
  filter(qval < 0.1)  %>%
  mutate(trajectory =  case_when(value == "traj3" ~ "Trajectory 3 vs. Trajectory 2",
                                 value == "traj1" ~ "Trajectory 1 vs. Trajectory 2"))
```


---

## newheatmapde_offspring.R

- Path: `/mnt/data/script_contents/newheatmapde_offspring.R`

- Size: 5KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")


significant_results <- read.delim("all_results_withoutinteraction_offspring.tsv")




significant_results_select <- significant_results %>%
  filter(metadata == "trajcluster") %>%
  filter(qval_individual <0.1 ) %>%
  filter(is.na(error)==TRUE)  %>%
  mutate(Trajectory = case_when(value =="traj1" ~ "Trajectory 1 vs. Trajectory 2",
```


---

## para.R

- Path: `/mnt/data/script_contents/para.R`

- Size: 4KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium_same.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Parabacteroides_distasonis.rda")

Parabacteroides_distasonis <- Parabacteroides_distasonis %>%
  dplyr::select(genome,pair_visit)

cazy_annotation_genome.parabac <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/8.cazy_annotation_genome.parabac.tsv", header=FALSE)

complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")



```


---

## pcoa1.5.R

- Path: `/mnt/data/script_contents/pcoa1.5.R`

- Size: 4KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(microViz)
library(tidyverse)
library(microbiome)
library(cluster)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision")


data_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")


supple <- phyloseq_validate(data_phylo, remove_undetected = TRUE)
supple <- tax_fix(supple)



offspring_1.5 <- supple %>%
  ps_filter(
    str_detect(visit,"BM6") | str_detect(visit,"MG3") 
```


---

## permonvafortransmission.R

- Path: `/mnt/data/script_contents/permonvafortransmission.R`

- Size: 8KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
# REQUIRED LIBRARIES

library(ggplot2)
library(tidyverse)
library(vegan)
library(umap)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/who_m6.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/allvisittransmission.rda")

#d28 <- read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/d28.sav")
#d42 <- read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/d42.sav")
#y1 <- read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/y1.sav")
#m6 <- read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/m6.sav")



# Original column names
#orig_names <- colnames(m6)
```


---

## phenotypefortable1_new.R

- Path: `/mnt/data/script_contents/phenotypefortable1_new.R`

- Size: 29KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(gtsummary)
library(tableone)
library(labelled)
library(rstatix)
library(writexl)
library(officer)
library(flextable)
library(readxl)
library(scales)
setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/")





complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

complete_phylo_pheno <- as.data.frame(complete_phylo@sam_data)
```


---

## pick_samples_for_revision.R

- Path: `/mnt/data/script_contents/pick_samples_for_revision.R`

- Size: 692B

- Detected type: R script (.R)

- Preview (first 30 lines):

```
library(tidyverse)
library(phyloseq)
library(microViz)

`%nin%` = Negate(`%in%`)

newdata_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/newdata_phylo.rds")

forrevision <- phyloseq_validate(newdata_phylo, remove_undetected = TRUE)



forrevision <- forrevision %>%
  ps_filter(
    selected == "Yes"
  ) %>%
  ps_filter(
    is.na(class_1) == FALSE
  )

```


---

## processDDS_at1year.R

- Path: `/mnt/data/script_contents/processDDS_at1year.R`

- Size: 6KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())

library(dplyr)
library(stringr)
library(tidyr)
library(haven)
library(Hmisc)


y1 <- read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/y1.sav")
#REMOVE less quailifed qustioniarie

new_names <- sapply(y1, function(x) label(x))



label_map <- c(
  IDchild = "IDchild",
  VAR1 = "VAR1",
  
```


---

## processDDS_at6months.R

- Path: `/mnt/data/script_contents/processDDS_at6months.R`

- Size: 3KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())

library(dplyr)
library(stringr)
library(tidyr)

m6 <-read_sav("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/喂养数据250617/m6.sav")

cols_to_fill <- c("rice_cereal", "egg_yolk","egg_white","vegetable_puree",
                  "fruit_puree","liver_puree","fish_puree","minced_meat","mashed_rice")

#REMOVE less quailifed qustioniarie
m6 <- m6   %>%
  filter(rowMeans(is.na(across(all_of(cols_to_fill)))) <= 0.4)

# 1️⃣ Define your structured columns

# 2️⃣ Recode structured columns: 1 = intake, 0 = no intake
m6 <- m6 %>%
  mutate(across(all_of(cols_to_fill), ~ ifelse(. >= 2, 1, 0), .names = "intake_{col}"))
```


---

## releaserate.R

- Path: `/mnt/data/script_contents/releaserate.R`

- Size: 13KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rstatix)


complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

#complete_phylo <- aggregate_taxa(complete_phylo, "Genus")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid")

speciestable  <- as.data.frame(complete_phylo@otu_table)
class(speciestable) <- "data.frame"
```


---

## releaserate_genus.R

- Path: `/mnt/data/script_contents/releaserate_genus.R`

- Size: 7KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(sjmisc)

complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

complete_phylo <- aggregate_taxa(complete_phylo, "Genus")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid")

speciestable  <- as.data.frame(complete_phylo@otu_table)
class(speciestable) <- "data.frame"
```


---

## reorganizemeta_new.R

- Path: `/mnt/data/script_contents/reorganizemeta_new.R`

- Size: 4KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(readxl)
library(haven)
library(memisc)
library(labelled)


`%ni%` <- Negate(`%in%`)

setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta")


tjmetarevision <- read_delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/tjmeta_merged_1add2add3add4add5.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)



cleantj <- tjmetarevision 
```


---

## retention.R

- Path: `/mnt/data/script_contents/retention.R`

- Size: 15KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)


complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

#complete_phylo <- aggregate_taxa(complete_phylo, "Genus")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"

metadata <- metadata %>%
  rownames_to_column("sampleid")

speciestable  <- as.data.frame(complete_phylo@otu_table)
class(speciestable) <- "data.frame"

```


---

## revision_everydistance.R

- Path: `/mnt/data/script_contents/revision_everydistance.R`

- Size: 6KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```

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



data_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
metadata <- as.data.frame(data_phylo@sam_data)
```


---

## somaticgrowth_with_GT_GH.R

- Path: `/mnt/data/script_contents/somaticgrowth_with_GT_GH.R`

- Size: 6KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggpubr)
library(broom)

###########################################################
#clustering heatmap for cazy 
###########################################################
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/spreadcazy.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"
metadata <- metadata %>%
  rownames_to_column("sampleid") %>%
  dplyr::select(sampleid,class_1,zwhz_1:zwhz_7)
```


---

## spearate_specific_transmission_new.R

- Path: `/mnt/data/script_contents/spearate_specific_transmission_new.R`

- Size: 18KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggh4x)
library(patchwork)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/")

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/allvisittransmission.rda")

DEsig.transmission_event<-allvisittransmission %>%
  mutate(familyid= str_remove(FamilyID1,"F")) 





subcomplete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
```


---

## three_transmit_status.R

- Path: `/mnt/data/script_contents/three_transmit_status.R`

- Size: 5KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggh4x)
library(patchwork)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/")

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/allvisittransmission.rda")

DEsig.transmission_event<-allvisittransmission %>%
  mutate(familyid= str_remove(FamilyID1,"F")) 





subcomplete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
```


---

## top_transmission.R

- Path: `/mnt/data/script_contents/top_transmission.R`

- Size: 6KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
rm(list = ls())
# REQUIRED LIBRARIES
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidytext) 
library(haven)
library(labelled)
library(vegan)

`%ni%` <- function(x, table) {
  !(x %in% table)
}


transmission <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/metadata/transmission.info2")


```


---

## transmission_genus.R

- Path: `/mnt/data/script_contents/transmission_genus.R`

- Size: 11KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggh4x)
library(patchwork)
library(microbiome)
library(sjmisc)



setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/")

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/allvisittransmission.rda")

DEsig.transmission_event<-allvisittransmission %>%
  mutate(familyid= str_remove(FamilyID1,"F")) 



```


---

## transmission_rate_separeatevisit.R

- Path: `/mnt/data/script_contents/transmission_rate_separeatevisit.R`

- Size: 3KB

- Detected type: R script (.R)

- Preview (first 30 lines):

```
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/instrain_transmission.rda")

library(ggrepel)


instrain_transmission <- instrain_transmission %>%
  dplyr::rename(Visit = V2,
                `Next visit` = V3,
                Trajectory = V4,
                `Transmitted number` = V5,
                `Transmissibility` = V6,
                `Not transmitted number` = V7) 


instrain_transmission$Visit <- factor(instrain_transmission$Visit,
                                      levels = c("BM1.5","BM6","BM12","BM24","BM36"))

instrain_transmission %>%
  filter(V1 == "Bifidobacterium pseudocatenulatum") %>%
  ggplot(aes(x = Visit, y = Transmissibility,
```


---

## transmisson_and_abundance_overall.R

- Path: `/mnt/data/script_contents/transmisson_and_abundance_overall.R`

- Size: 9KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(readr)
library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggpubr)
library(rstatix)
library(lme4)
library(lmerTest)

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/spreadcazy.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/Bifidobacterium.rda")
complete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
significant_results <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/statistics_results_final/offspring_output_species/significant_results.tsv")

significant_results_maternal <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/statistics_results_final/mother_output_species/significant_results.tsv")



```


---

## various_rate_new.R

- Path: `/mnt/data/script_contents/various_rate_new.R`

- Size: 9KB

- Detected type: R script (.R)

- Inferred actions / keywords: loop for

- Preview (first 30 lines):

```
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(viridis)
library(rstatix)
library(ggpubr)
library(data.table)

###############################################
#otu table
###############################################
  
#setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/")

subcomplete_phylo <-  readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")
otudata <- as.data.frame(subcomplete_phylo@otu_table)
class(otudata) <- "data.frame"


metadata <- as.data.frame(subcomplete_phylo@sam_data)
```


---
