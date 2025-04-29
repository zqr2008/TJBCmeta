
rm(list=ls())
library(vegan)
library(ggplot2)
library(tidyverse)
library(phyloseq)

setwd("/hwfsxx2/ST_HN/P17Z10200N0306/zhaiqiangrong/tjmetasecondround/")

load("bray.rda")

data_phylo <- readRDS("complete_phylo.rds")


sepvisit <- subset_samples(data_phylo,visit == "MG1")


subsample <- as.data.frame(sepvisit@otu_table)


phenotype <- sepvisit@sam_data[,xxxxxx]

covert.species <- bray[colnames(bray) %in% rownames(subsample),
                         rownames(bray) %in% rownames(subsample)]


#create dataframe to hold result
store_result = matrix(rep(0,ncol(phenotype)*3),ncol = 3)
store_result[,1] = colnames(phenotype)


#loop the permonva
for(i in 1:nrow(store_result)){
  if(sum(is.na(phenotype[,i]))==nrow(phenotype)){
    next
  }
  
  m = covert.species[!is.na(phenotype[,i]),]
  n = phenotype[,i][!is.na(phenotype[,i])]
  n = as.matrix(n)
  
  if(dim(unique(n))[1]==1){
    next
  }
  
  a = adonis(m~n,permutations = 999,method = "bray")
  
  store_result[i,2] = a$aov.tab[1,5]
  store_result[i,3] = a$aov.tab[1,6]
  }
  
           
  visitname <- store_result
           
  save(visitname,file = "storestorestore.rda")
  
