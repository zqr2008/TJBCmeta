#!bin/bash
# for strain <Ex_SGB> of <species>
t=Ex_SGB
species=$(grep ${t} SE_species.list | cut -f 1 )

# add group label
# Intra-individual(mother): a pair of fecal samples from the same mother subject
# Intra-individual(offspring): a pair of fecal samples from the same infant subject
# Intra-family: a pair of fecal samples from the same family; one is mother's, another is offspring's
# Intra-family: a pair of fecal samples from different families; one is mother's, another is offspring's
# None: a pair of fecal samples does not belong to any group above 
awk -v OFS="\t" 'NR==FNR{s[$1]=$3"\t"$5; f[$1]=$4}NR!=FNR{print $1,s[$1],f[$1],$2,s[$2],f[$2],$3}' ./1w8.metadata.forR ./output/${t}/${t}_nGD.tsv | awk '{if($2==$6 && $3=="mother"){print "'${species}'\t'${t}'\t"$0"\tIntra-individual(mother)"} else if($2==$6 && $3=="offspring"){print "'${species}'\t'${t}'\t"$0"\tIntra-individual(offspring)"}else if($2!=$6 && $4==$8 && ($3"-"$7=="mother-offspring" || $3"-"$7=="offspring-mother")){print "'${species}'\t'${t}'\t"$0"\tIntra-family"}else if ($4!=$8 && ($3"-"$7=="mother-offspring" || $3"-"$7=="offspring-mother")){print "'${species}'\t'${t}'\t"$0"\tInter-family"}else{print "'${species}'\t'${t}'\t"$0"\tNone"}}' > ./output/${t}/${t}_nGD_labeled.tsv && \

# generate the R script for determine group differences
cat > ./output/${t}/${t}.R <<EOF

library(tidyverse)
library(dplyr)
library(rstatix)

df <- read.table("./output/${t}/${t}_nGD_labeled.tsv",header = F, sep = "\t")
df <- df %>%
  mutate(across(c(V12), as.factor)) %>%	#V12 is the group label
  mutate(across(c(V11), as.numeric))	#V11 is the nGD

# determine group differences
stat.test <- df %>%
  group_by(V1) %>%			#V1 is the species name
  wilcox_test(V11 ~ V12) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(p.adj = round(p.adj,3))

stat.test <- stat.test %>%
  add_xy_position(x = "V1", step.increase = 0.3, dodge = 0.5)
stat.test <- stat.test %>% mutate(across(c(1:ncol(stat.test)), as.character))
write.table(stat.test,file="./output/${t}/${t}_stat.test.tsv", quote = F, row.names = F, sep = "\t")

# calculate group Mean and SD
MeanSD <- df %>%
  group_by(V1, V12) %>%
  summarise(Mean=mean(V11),SD=sd(V11))
write.table(MeanSD, file="./output/${t}/${t}_MeanSD.tsv", quote = F, row.names = F, sep = "\t")

EOF

# run R script
Rscript ./output/${t}/${t}.R
