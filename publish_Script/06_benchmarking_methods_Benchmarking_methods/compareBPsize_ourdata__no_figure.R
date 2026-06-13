# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 06_benchmarking_methods_Benchmarking_methods
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Benchmarks pipeline/read-size effects and exports benchmarking plots.
# Main input(s): final 10 benchmark.xlsx; ourdata_compare150bp100bp.all.t.summary.tsv; ourdata_compare150bp100bp.all.diff.only_share.summary; tjmetamerged100bp_15_20251117.txt; plus 1 more
# Main output(s): ourdata_Benchmarking1.pdf; ourdata_Benchmarking2.pdf; ourdata_Benchmarking3.pdf; ourdata_Benchmarking4.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(viridis)   # better colors
library(ggpubr)
library(purrr)
library(patchwork)

setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/")

final_10_benchmark <- read_excel("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/final 10 benchmark.xlsx") %>%
  mutate(SampleID = paste0(样品名称,".mp4"))

finalsample<-c("MYF0006847.mp4","MYF0054133.mp4","MYF0071181.mp4",
               "MYF0071417.mp4","MYF0071670.mp4","MYF0046187.mp4",
               "MYF0071466.mp4","MYF0073097.mp4")
#delte host rate too high, and resequence error exsit
ourdata_compare150bp100bp.all.t.summary <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/ourdata_compare150bp100bp.all.t.summary.tsv") %>%
  filter(SampleID %in% finalsample) 
  




# 1. Total Original vs Trimmed (line plot with labels)
df_long_total <- ourdata_compare150bp100bp.all.t.summary %>%
  dplyr::select(SampleID, Total_150bp, Tqotal_100bp) %>%
  pivot_longer(-SampleID, names_to="Type", values_to="Count")

p1 <- ggplot(df_long_total, aes(x=SampleID, y=Count, group=Type, color=Type)) +
  geom_line(size=1) +
  geom_point(size=2) +
  geom_text(data = subset(df_long_total, Type=="Total_150bp"),
            aes(label=Count), vjust=-1, size=3) +
  geom_text(data = subset(df_long_total, Type=="Tqotal_100bp"),
            aes(label=Count), vjust=1.5, size=3) +
  scale_color_viridis(discrete=TRUE, option="H") +
  theme_classic(base_size=13) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="top") +
  labs(title="Total Strains 150bp and 100bp",
       x="Sample ID", y="Total Count", color="")
# 2. Unique Counts (bar plot with labels)
df_long_unique <- ourdata_compare150bp100bp.all.t.summary %>%
  dplyr::select(SampleID, Uniq_to_150bp, Uniq_to_100bp) %>%
  pivot_longer(-SampleID, names_to="Type", values_to="Count")

p2 <- ggplot(df_long_unique, aes(x=SampleID, y=Count, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  geom_text(aes(label=Count), 
            position=position_dodge(width=0.8), 
            vjust=-0.3, size=3) +
  scale_fill_viridis(discrete=TRUE, option="H") +
  theme_classic(base_size=13) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="top") +
  labs(title="Unique Strain Counts", x="Sample ID", y="Count", fill="")




## ---- Trimmed percentages ----
plot_data_trimmed <- ourdata_compare150bp100bp.all.t.summary %>%
  dplyr::rename(Total = Tqotal_100bp,
                Unique = Uniq_to_100bp,
                Shared = Share) %>%
  dplyr::select(SampleID, Total, Unique, Shared) %>%
  pivot_longer(cols = c("Shared", "Unique"),
               names_to = "Category", values_to = "Count") %>%
  mutate(Percent = ifelse(Total > 0, Count / Total * 100, 0))

cols_trimmed <- viridis(2, option = "H")
names(cols_trimmed) <- c("Shared", "Unique")

p_trimmed <- ggplot(plot_data_trimmed, aes(x = SampleID, y = Percent, fill = Category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(Percent >= 3, paste0(round(Percent,1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 4,color="white") +
  scale_fill_manual(values = cols_trimmed,
                    labels = c("Shared" = "Shared strains", "Unique" = "Unique (100bp)")) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Shared vs Unique (100bp)",
       x = "Sample ID", y = "Percent of Total 100bp", fill = "")


## ---- Original percentages ----
plot_data_original <- ourdata_compare150bp100bp.all.t.summary %>%
  dplyr::rename(Total = Total_150bp,
                Unique = Uniq_to_150bp,
                Shared = Share) %>%
  dplyr::select(SampleID, Total, Unique, Shared) %>%
  pivot_longer(cols = c("Shared", "Unique"),
               names_to = "Category", values_to = "Count") %>%
  mutate(Percent = ifelse(Total > 0, Count / Total * 100, 0))

cols_original <- c("Shared" = viridis(1, option = "H"),
                   "Unique" = "firebrick")

p_original <- ggplot(plot_data_original, aes(x = SampleID, y = Percent, fill = Category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(Percent >= 3, paste0(round(Percent,1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 4,color="white") +
  scale_fill_manual(values = cols_original,
                    labels = c("Shared" = "Shared strains", "Unique" = "Unique 150bp")) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Shared vs Unique (150bp)",
       x = "Sample ID", y = "Percent of Total 150bp", fill = "")



# Combine into 2x2
final_plot <- (p1 | p2) / (p_trimmed |p_original)
final_plot

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/ourdata_Benchmarking1.pdf",width =30,height = 30,units = "cm")


combined_strain_comparison <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/ourdata_compare150bp100bp.all.diff.only_share.summary") %>%
  mutate(tlevel = str_split_fixed(taxonomy,"t__",n=2)[,2])



p <- combined_strain_comparison %>%
  filter(Diff != 0) %>%
  ggplot(aes(x = SampleID, y = Diff)) +
  # Create boxplot with outlier points hidden (we'll add them separately)
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
  # Add jittered points with transparency
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "steelblue") +
  # Add mean points
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  # Add horizontal line at y=0 for reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  # Label outliers using IQR method
  geom_label_repel(
    data = . %>%
      group_by(SampleID) %>%
      mutate(
        Q1 = quantile(Diff, 0.25),
        Q3 = quantile(Diff, 0.75),
        IQR = Q3 - Q1,
        is_outlier = Diff < (Q1 - 1.5 * IQR) | Diff > (Q3 + 1.5 * IQR)
      ) %>%
      filter(is_outlier),
    aes(label = tlevel),
    size = 2.5,
    box.padding = 0.5,
    max.overlaps = 20,
    segment.color = "gray",
    min.segment.length = 0.2
  ) +
  # Customize theme for academic appearance
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  # Add labels and title
  labs(
    title = "Distribution of Strain Abundance Differences Between 150bp and 100bp Samples",
    x = "Sample ID",
    y = "Difference in Relative Abundance (100bp result - 150bp result)"
    #caption = "Outliers defined as values beyond 1.5×IQR from the quartiles\nRed diamonds represent mean values"
  ) 
# Optional: add faceting if you have many samples
# facet_wrap(~ SampleID, scales = "free_x", ncol = 4) +
# Adjust y-axis limits if needed
# coord_cartesian(ylim = c(-2, 2))

# Display the plot
print(p)

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/ourdata_Benchmarking2.pdf",width =28,height = 20,units = "cm")


tjmetamerged100bp_15_20251117 <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/tjmetamerged100bp_15_20251117.txt") %>%
  filter(str_detect(clade_name,"t__")) %>%
  mutate(tlevel = str_split_fixed(clade_name,"t__",n=2)[,2]) %>%
  dplyr::select( any_of(finalsample),tlevel) 


tjmetamerged150bp_15_20251029 <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/tjmetamerged150bp_15_20251029.txt")%>%
  filter(str_detect(clade_name,"t__")) %>%
  mutate(tlevel = str_split_fixed(clade_name,"t__",n=2)[,2])%>%
  dplyr::select( any_of(finalsample),tlevel) 


df100 <- tjmetamerged100bp_15_20251117
df150 <- tjmetamerged150bp_15_20251029

# sample columns (all except tlevel)
sample_cols <- setdiff(colnames(df100), "tlevel")

# reshape to long format
long100 <- df100 %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "SampleID",
               values_to = "Count100")

long150 <- df150 %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "SampleID",
               values_to = "Count150")

# join two long tables by tlevel + sampleID
merged <- long100 %>%
  left_join(long150, by = c("tlevel", "SampleID"))

# filter only cases where:
#   100bp != 0 AND 150bp == 0
filtered <- merged %>%
  filter(Count100 != 0 & Count150 == 0)

# add total count across samples for each tlevel
result <- filtered %>%
  group_by(tlevel) %>%
  mutate(total_count = sum(Count100)) %>%
  ungroup() %>%
  select(t__SGB = tlevel,
         total_count,
         SampleID,
         Count = Count100)



# Plot
ggplot(result, aes(x = t__SGB, y = Count, fill = SampleID)) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "H") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top") +
  labs(x = "Strain (t__SGB)", 
       y = "Sum of the Relative Abundance across all samples", 
       fill = "Sample ID",
       title = "Unique 100bp result Strain Relative Abundance Across Samples")

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/ourdata_Benchmarking3.pdf",width =25,height = 18,units = "cm")



# Reshape data
df_long <- combined_strain_comparison %>%
  dplyr::select(tlevel, AB_150bp, AB_100bp) %>%
  tidyr::pivot_longer(cols = c(AB_150bp, AB_100bp),
                      names_to = "Type", values_to = "Abundance") %>%
  dplyr::mutate(Type = factor(Type, levels=c("AB_150bp","AB_100bp")))

# Count observations per tlevel and filter >=3
top_candidates <- df_long %>%
  dplyr::group_by(tlevel) %>%
  dplyr::filter(n() >= 3) %>%     # keep only tlevels with >=3 obs
  dplyr::summarise(median_diff = median(Abundance[Type=="Trimmed"] - Abundance[Type=="Original"], na.rm=TRUE)) %>%
  dplyr::arrange(desc(median_diff)) %>%
  dplyr::slice_head(n=10) %>%
  dplyr::pull(tlevel)

# Filter top 10 tlevels
df_top10 <- df_long %>%
  dplyr::filter(tlevel %in% top_candidates) %>%
  dplyr::mutate(tlevel = factor(tlevel, levels = top_candidates))

# Calculate Wilcoxon p-values per tlevel
p_values <- ggpubr::compare_means(
  Abundance ~ Type,
  data = df_top10,
  method = "wilcox.test",
  group.by = "tlevel"
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(y.position = max(df_top10$Abundance[df_top10$tlevel == tlevel]) * 1.05) %>%
  dplyr::ungroup()

# Plot
ggplot(df_top10, aes(x=Type, y=Abundance, fill=Type)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) +
  geom_jitter(width=0.2, size=1, alpha=0.6) +
  geom_text(data=p_values, aes(x=1.5, y=y.position+2, label=paste0("P =", p.format)), inherit.aes = FALSE) +
  scale_fill_viridis(discrete=TRUE, option="H") +
  facet_wrap(~tlevel, scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(x="", y="Abundance",
       title="150bp vs 100bp Abundance for Top 10 Different Strains (>=3 Observations)")


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/ourdata_Benchmarking4.pdf",width =25,height = 18,units = "cm")
