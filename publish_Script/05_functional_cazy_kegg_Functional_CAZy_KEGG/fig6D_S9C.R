


#========= Script for Fig 6D ===============


df_mag <- read.table("BP_cazy.tsv", header = F, sep="\t") %>%
  filter(grepl("Bifidobacterium pseudocatenulatum",V9)) %>%
  filter(V10=="HIGH") %>%
  filter(V11=="Transmit") %>%
  filter(V8==3) %>%
  select(V1,V7) %>%
  rename(V6=V7)

df_GHCOUNT <- rbind(df_ex,df_mag) %>%
  reshape2::dcast(V1~V6,fun.aggregate = length) %>%
  mutate(sum_GH=rowSums(across(starts_with("GH")))) %>%
  mutate(sum_GH = rowSums(across(starts_with("GH")))) %>%
  mutate(`GH Count` = case_when(sum_GH >= 50 ~ "High GH gene count (>=50)",
                                sum_GH < 50 ~ "Low GH gene count (<50)")) %>%
  select(V1, `GH Count`)

df <- rbind(df_ex,df_mag) %>%
  dplyr::filter(V6 %in% c("GH1","GH172","GH2_10","GH104","GH121","GH127","GH27", "GH3+3.2.1.21","GH3",
                          "GH30_2","GH31_15","GH32","GH36","GH42","GH43_29","GH8"
  )
  | str_detect(V6,"^GH13_.*|^GH30_.*|^GH31_.*|^GH43_.*")
  | str_detect(V6, "^GH29[_]*|^GH95[_]*|^GH141[_]*|^GH33[_]*|^GH2_.*|^GH42[_]*|^GH20[_]*|
                             ^GH136[_]*|^GH112[_]*|^GH18[_]*|^GH109[_]*")
  ) %>%
  reshape2::dcast(V1~V6,fun.aggregate = length) %>%
  select(where(~sum(. != 0, na.rm = TRUE) > 0))
bin_matrix <- df %>% column_to_rownames("V1") %>% as.matrix()

gene_dist <- vegdist(bin_matrix, method = "bray") 
nj_tree <- nj(gene_dist)

tree <- nj_tree
plot_tree <- tree
plot_tree$edge.length <- pmax(plot_tree$edge.length, 0)
ultra_tree <- chronos(plot_tree)
hclust_tree <- as.hclust(ultra_tree)

annotation_row <- df %>%
  merge(.,df_GHCOUNT,by="V1") %>%
  select(V1,`GH Count`) %>%
  mutate(samp = gsub("\\..*","",V1)) %>%
  mutate(samp = case_when(grepl("MYF",V1) ~ gsub("^D","",samp), .default = samp)) %>%
  mutate(samp = case_when(grepl("MYF",V1) ~ gsub("-1","",samp), .default = samp)) %>%
  merge(., metadata[,c("feeding_type_28days","BMI_mo","UA.x","TG",
                       #"preterm",
                       "delivery_mode","class_1")],
        by.x="samp",by.y=0,all.x=T) %>%
  select(-samp) %>%
  column_to_rownames("V1") %>%
  mutate(feeding_type_28days = case_when(is.na(feeding_type_28days) ~ "Missing", 
                                         .default = paste0(toupper(
                                           substring(feeding_type_28days, 1, 1)), 
                                           substring(feeding_type_28days, 2)))) %>%
  mutate(delivery_mode = case_when(delivery_mode=="Caesarean section" ~ "C-section", 
                                   .default = delivery_mode)) %>%
  mutate(class_1 = paste("Trajectory ", class_1, sep="")) %>%
  mutate(across(c(1,2,6:7),as.factor)) %>%
  mutate(across(c(3:5),as.numeric)) %>%
  rename(Trajectory = class_1, `Pregnancy BMI` = BMI_mo, UA = UA.x,
         `Delivery mode` = delivery_mode, `Feeding type` = feeding_type_28days)
annotation_row <- annotation_row[row.names(bin_matrix),]
annotation_row <- annotation_row[,c("GH Count","TG","UA","Pregnancy BMI",
                                    "Feeding type","Delivery mode","Trajectory")]

UA_palette <- colorRampPalette(c("#e9f1f6", "#B72289"))(100)
TG_palette <- colorRampPalette(c("#e9f1f6", "#D55E00"))(100)
BMI_mo_palette <- colorRampPalette(c("#e9f1f6","#B22222"))(100)
annotation_color <- list(
  Trajectory = c("Trajectory 1" = "#e95280", 
                 "Trajectory 2" = "#23b1a5", 
                 "Trajectory 3"="#e49b0f"),
  `Delivery mode` = c("C-section" = "#B0B0B0", "Vaginal delivery" = "#16317d"),
  `Feeding type` = c("Exclusive breastfeeding" = "#006400",
                     "Exclusive formula fed" = "#911f27",
                     "Mixed feeding" = "#60c6db",
                     "Missing" = "#e9f1f6"),
  UA = UA_palette,
  TG = TG_palette,
  `Pregnancy BMI` = BMI_mo_palette,
  `GH Count` = c("High GH gene count (>=50)" = "#fdaa89", "Low GH gene count (<50)" = "#aab8d7")
)

col.matrix <- as.data.frame(t(bin_matrix))
col_dist = dist(col.matrix)
hclust_1 <- hclust(col_dist)
col_cluster <- hclust_1 %>% 
  as.dendrogram %>% 
  ladderize(right = TRUE) %>% 
  as.hclust

bin_matrix     <- bin_matrix[hclust_tree$labels, , drop = FALSE]
annotation_row <- annotation_row[hclust_tree$labels, , drop = FALSE]
pheatmap ::pheatmap(t(bin_matrix), annotation_col = annotation_row, cluster_col = hclust_tree,
                    treeheight_col = 50, show_colnames =F, annotation_colors = annotation_color,
                    border_color = "white",cutree_cols =10, cluster_rows = col_cluster, cutree_rows=3,
                    fontsize = 12)


#========= Script for sup Fig 9C =============== 

df_ex_COUNT <- read.table("validation.cazy.tsv", 
                          header = F, sep="\t") %>%
  filter(V7==3) %>%
  select(V1,V6) %>%
  reshape2::dcast(V1~V6,fun.aggregate = length) %>%
  mutate(sum_GH=rowSums(across(starts_with("GH")))) %>%
  mutate(sum_GH = rowSums(across(starts_with("GH")))) %>%
  mutate(`GH Count` = case_when(sum_GH >= 50 ~ "High GH gene count (>=50)",
                                sum_GH < 50 ~ "Low GH gene count (<50)")) %>%
  select(V1, sum_GH, `GH Count`)

df_ex <- read.table("validation.cazy.tsv", 
                    header = F, sep="\t") %>%
  filter(V7==3) %>%
  merge(.,strain_rename_key,by=1) %>%
  select(strain_internal,V6) %>%
  dplyr::filter(V6 %in% c("GH1","GH172","GH2_10","GH104","GH121","GH127","GH27", "GH3+3.2.1.21","GH3",
                          "GH30_2","GH31_15","GH32","GH36","GH42","GH43_29","GH8"
  )
  | str_detect(V6,"^GH13_.*|^GH30_.*|^GH31_.*|^GH43_.*")
  | str_detect(V6, "^GH29[_]*|^GH95[_]*|^GH141[_]*|^GH33[_]*|^GH2_.*|^GH42[_]*|^GH20[_]*|
                             ^GH136[_]*|^GH112[_]*|^GH18[_]*|^GH109[_]*")
  ) %>%
  reshape2::dcast(strain_internal~V6,fun.aggregate = length) %>%
  select(where(~sum(. != 0, na.rm = TRUE) > 0))
bin_matrix <- df_ex %>% column_to_rownames("strain_internal") %>% as.matrix()

gene_dist <- vegdist(bin_matrix, method = "bray") 
nj_tree <- nj(gene_dist)

growth <- read.csv(file = "growth.csv") %>%
  select(-c(sd_OD600_blank_subtracted)) %>%
  reshape2::dcast(strainname_plot+GH_count+genecoutcate~condition_formal) %>%
  merge(.,strain_rename_key, by.x=1, by.y=2,all=T)

annotation_row <- growth %>%
  mutate(genecoutcate=case_when(genecoutcate == "Low gene counts (<50)" ~ "Low GH gene count (<50)",
                                genecoutcate == "High gene counts (>=50)" ~ "High GH gene count (>=50)")) %>%
  mutate(across(c(3),as.factor)) %>%
  mutate(across(c(2,4:10),as.numeric)) %>%
  rename(`GH Count` = GH_count) %>%
  #rbind(.,growth_external) %>%
  #column_to_rownames("strain") %>%
  #select(-c(strainname_plot,`GH Count`)) %>% 
  column_to_rownames("strainname_plot") %>%
  select(-c(strain,`GH Count`)) %>% 
  rename(`GH Count` = genecoutcate)

annotation_row <- annotation_row[row.names(bin_matrix),]
annotation_row <- annotation_row[,c("GH Count","HMO Mix","LNT","LNnT",
                                    "6'-SL","3'-FL","2'-FL","MRS")]

palette1 <- colorRampPalette(c("#e9f1f6", "#dc3023"))(100)
palette2 <- colorRampPalette(c("#e9f1f6", "#ff8936"))(100)
palette3 <- colorRampPalette(c("#e9f1f6","#ffc64b"))(100)
palette4 <- colorRampPalette(c("#e9f1f6","#0aa344"))(100)
palette5 <- colorRampPalette(c("#e9f1f6","#48c0a3"))(100)
palette6 <- colorRampPalette(c("#e9f1f6","#177cb0"))(100)
palette7 <- colorRampPalette(c("#e9f1f6","#801dae"))(100)
annotation_color <- list(
  #`GH Count` = UA_palette,
  `MRS` = palette1,
  `2'-FL` = palette2,
  `3'-FL` = palette3,
  `6'-SL` = palette4,
  `LNnT` = palette5,
  `LNT` = palette6,
  `HMO Mix` = palette7,
  `GH Count` = c("High GH gene count (>=50)" = "#fdaa89", "Low GH gene count (<50)" = "#aab8d7")
)

tree <- nj_tree
plot_tree <- tree
plot_tree$edge.length <- pmax(plot_tree$edge.length, 0)
ultra_tree <- chronos(plot_tree)
hclust_tree <- as.hclust(ultra_tree)

col.matrix <- as.data.frame(t(bin_matrix))
col_dist = dist(col.matrix)
hclust_1 <- hclust(col_dist)
col_cluster <- hclust_1 %>% 
  as.dendrogram %>% 
  ladderize(right = TRUE) %>% 
  as.hclust

bin_matrix     <- bin_matrix[hclust_tree$labels, , drop = FALSE]
annotation_row <- annotation_row[hclust_tree$labels, , drop = FALSE]
pheatmap ::pheatmap(bin_matrix, annotation_row = annotation_row, cluster_rows = hclust_tree,
                    treeheight_row = 100, show_rownames =F,annotation_colors = annotation_color,
                    border_color = "white",cutree_rows =8, cluster_cols = col_cluster, cutree_cols=3,
                    fontsize = 12,angle_col = "90")




