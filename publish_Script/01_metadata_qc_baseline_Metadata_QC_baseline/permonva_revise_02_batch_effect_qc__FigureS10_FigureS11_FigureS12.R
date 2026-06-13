# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 01_metadata_qc_baseline_Metadata_QC_baseline
# Figure(s): FigureS10; FigureS11; FigureS12
# Table(s): No direct submitted table matched from filename
# Purpose: Generates quality-control summaries and QC figure/table outputs.
# Main input(s): newdata_phylo_NMrevision_0305.rds; df_pcoa.rda
# Main output(s): FigureS12.pdf; FigureS10.pdf; FigureS11.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
# ---- Source Rmd chunk: setup (permonva_revise.Rmd lines 16-25) ----
rm(list = ls())

knitr::opts_chunk$set(echo = TRUE)
setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision")

set.seed(1000)

# ---- Source Rmd chunk: loading pkgs (permonva_revise.Rmd lines 28-38) ----
library(tidyverse)
library(readxl)
library(phyloseq)
library(sjmisc)
library(lubridate)
library(labelled)
library(conflicted)
conflicts_prefer(base::setdiff)

# ---- Source Rmd chunk: batch effect plot (permonva_revise.Rmd lines 121-260) ----
cols10 <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # grey
  "#000000", # black
  "#66C2A5"  # extra (teal)
)


complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"


load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/pcoa_distance_objects/df_pcoa.rda")


# ensure SampleID is rownames in metadata
metadata$SampleID <- rownames(metadata)
qc_vars <- grep("^QC__", colnames(metadata), value = TRUE)


# merge
df_plot <- df_pcoa %>%
  left_join(metadata, by = "SampleID")




df_long <- df_plot %>%
  select(SampleID, Axis1, Axis2, all_of(qc_vars)) %>%
  pivot_longer(
    cols = all_of(qc_vars),
    names_to = "QC_metric",
    values_to = "QC_value"
  )


df_long <- df_long %>%
  left_join(metadata,
            by = "SampleID")
# sanity check
stopifnot(nrow(df_plot) == nrow(df_pcoa))


qc_vars <- grep("^QC__", colnames(df_plot), value = TRUE)

df_long <- df_long %>%
  group_by(QC_metric) %>%
  mutate(QC_group = ntile(QC_value, 3)) %>%  # tertiles
  ungroup()



ggplot(df_long, aes(x = Axis1, y = Axis2, color = factor(QC_group))) +
  geom_point(size = 0.5,alpha=0.3) +
  facet_wrap(~ QC_metric, scales = "free") +
  theme_bw() +
  labs(color = "QC tertile")

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS12.pdf",width = 30, height = 20, units = "cm")


ggplot(df_long, aes(x = QC_value)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~ QC_metric, scales = "free") +
  theme_bw() +
  labs(
    x = "QC value",
    y = "Density",
    title = "Distribution of QC metrics"
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS10.pdf",width = 30, height = 20, units = "cm")



ggplot(df_long, aes(x = QC_value, fill = trajectory)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ QC_metric, scales = "free") +
  theme_bw() +
  labs(fill = "Trajectory")+
      scale_fill_manual(values = cols10)

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS11.pdf",width = 30, height = 20, units = "cm")

