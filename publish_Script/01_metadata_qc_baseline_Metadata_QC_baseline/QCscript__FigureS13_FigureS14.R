# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 01_metadata_qc_baseline_Metadata_QC_baseline
# Figure(s): FigureS13; FigureS14
# Table(s): No direct submitted table matched from filename
# Purpose: Generates quality-control summaries and QC figure/table outputs.
# Main input(s): 7.genome_QC_GTDB.update.addFID.addReRunQC.short.reNewQClass.include.addFdrepMAG.profile
# Main output(s): FigureS13 QC_genomes1.pdf; FigureS14 QC_genomes2.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
library(tidyverse)
library(ggplot2)
library(readxl)
library(pheatmap)



genome_QC_GTDB <-read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/genome_qc_inputs/7.genome_QC_GTDB.update.addFID.addReRunQC.short.reNewQClass.include.addFdrepMAG.profile")

genome_QC_GTDB$Quality_class <- genome_QC_GTDB$Quality_class_F
genome_QC_GTDB$Completeness <- genome_QC_GTDB$Completeness
genome_QC_GTDB$Contamination<- genome_QC_GTDB$Contamination


df <- genome_QC_GTDB

#######################################################
# 1. p1 = Stacked histogram of Completeness
#######################################################

library(patchwork)
library(ggplot2)

# Define the bins once → ensures perfect alignment
binwidth <- 2   # adjust if you want finer or coarser bins
xmin <- floor(min(df$Completeness, na.rm = TRUE))
xmax <- ceiling(max(df$Completeness, na.rm = TRUE))

p1 <- ggplot(df, aes(x = Completeness, fill = Quality_class)) +
  geom_histogram(binwidth = binwidth, position = "stack", color = "black") +
  scale_x_continuous(limits = c(xmin, xmax)) +
  theme_classic() +
  labs(
    title = "Distribution of Completeness (%)",
    x = "Completeness (%)",
    y = "Count"
  )+
  scale_fill_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))

#######################################################
# 2. p2 = Scatter plot with EXACT SAME x-axis scale
#######################################################

p2 <- ggplot(df, aes(x = Completeness, y = Contamination)) +
  geom_point(alpha = 0.7, size = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax)) +   # identical x-axis
  theme_classic() +
  labs(
    x = "Completeness (%)",
    y = "Contamination (%)"
  )

p1/p2
ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS13 QC_genomes1.pdf",width = 7, height = 7)

df$Quality_class <- factor(df$Quality_class, levels = c("LOW", "MEDIUM", "HIGH"))

## Plot 1: N50 by Quality_class
pic3 <- ggplot(df, aes(x = Quality_class, y = N50, color = Quality_class, fill = Quality_class)) +
  geom_boxplot(width = 0.2, outlier.colour = "black", outlier.size = 0.1,
               alpha =  0.2) +
  theme_classic() +
  labs(title = "N50 by Quality Class",
       x = "Quality Class",
       y = "N50")+
  scale_fill_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))+
  scale_color_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))

## Plot 1: N50 by Quality_class
pic4 <- ggplot(df, aes(x = Quality_class, y = Strain.heterogeneity, color = Quality_class, fill = Quality_class)) +
  geom_boxplot(width = 0.2, outlier.colour = "black", outlier.size = 0.1,
               alpha =  0.2) +
  theme_classic() +
  labs(title = "Strain heterogeneity by Quality Class",
       x = "Quality Class",
       y = "Strain.heterogeneity(%)")+
  scale_fill_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))+
  scale_color_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))

pic3/pic4
ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS14 QC_genomes2.pdf",width = 7, height = 7)

  

pic5 <- ggplot(df, aes(x = Quality_class, y =Total.length, color = Quality_class, fill = Quality_class)) +
  geom_boxplot(width = 0.2, outlier.colour = "black", outlier.size = 0.1,
               alpha =  0.2) +
  theme_classic() +

  scale_fill_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))+
  scale_color_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))



pic6 <- ggplot(df, aes(x = Quality_class, y = X..contigs.....5000.bp., color = Quality_class, fill = Quality_class)) +
  geom_boxplot(width = 0.2, outlier.colour = "black", outlier.size = 0.1,
               alpha =  0.2) +
  theme_classic() +
  labs(title = "contigs count (>=5000bp)",
       x = "Quality Class",
       y = "contigs count (>=5000bp)")+
  scale_fill_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))+
  scale_color_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))

pic7 <- ggplot(df, aes(x = Quality_class, y = X..contigs.....10000.bp., color = Quality_class, fill = Quality_class)) +
  geom_boxplot(width = 0.2, outlier.colour = "black", outlier.size = 0.1,
               alpha =  0.2) +
  theme_classic() +
  labs(title = "contigs count (>=10000bp)",
       x = "Quality Class",
       y = "contigs count (>=10000bp)")+
  scale_fill_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))+
  scale_color_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))


pic8 <- ggplot(df, aes(x = Quality_class, y = X..contigs.....50000.bp., color = Quality_class, fill = Quality_class)) +
  geom_boxplot(width = 0.2, outlier.colour = "black", outlier.size = 0.1,
               alpha =  0.2) +
  theme_classic() +
  labs(title = "contigs count (>=50000bp)",
       x = "Quality Class",
       y = "contigs count (>=50000bp)")+
  scale_fill_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))+
  scale_color_manual(values = c(
    "HIGH" = "#1A9850",      # blue
    "MEDIUM" =  "#FDBF6F",    # grey
    "LOW" =  "#D73027"       # red
  ))


(pic5+pic6)/(pic7+pic8)


