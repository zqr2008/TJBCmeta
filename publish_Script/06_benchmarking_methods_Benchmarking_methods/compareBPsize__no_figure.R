# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 06_benchmarking_methods_Benchmarking_methods
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Benchmarks pipeline/read-size effects and exports benchmarking plots.
# Main input(s): combined_strain_comparison.tsv; detailed_strain_comparison.csv
# Main output(s): Benchmarking1.pdf; Benchmarking2.pdf; Benchmarking3.pdf; Benchmarking4.pdf
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

combined_strain_comparison <- read.delim("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/combined_strain_comparison.tsv") %>%
  mutate(tlevel = str_split_fixed(Strain,"t__",n=2)[,2])

detailed_strain_comparison <- read.csv("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/detailed_strain_comparison.csv")


# 1. Total Original vs Trimmed (line plot with labels)
df_long_total <- detailed_strain_comparison %>%
  dplyr::select(SampleID, Total_Original, Total_Trimmed) %>%
  pivot_longer(-SampleID, names_to="Type", values_to="Count")

p1 <- ggplot(df_long_total, aes(x=SampleID, y=Count, group=Type, color=Type)) +
  geom_line(size=1) +
  geom_point(size=2) +
  geom_text(data = subset(df_long_total, Type=="Total_Original"),
            aes(label=Count), vjust=-1, size=3) +
  geom_text(data = subset(df_long_total, Type=="Total_Trimmed"),
            aes(label=Count), vjust=1.5, size=3) +
  scale_color_viridis(discrete=TRUE, option="H") +
  theme_classic(base_size=13) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="top") +
  labs(title="Total Strains Before and After Trimming",
       x="Sample ID", y="Total Count", color="")
# 2. Unique Counts (bar plot with labels)
df_long_unique <- detailed_strain_comparison %>%
  dplyr::select(SampleID, Unique_Original_Count, Unique_Trimmed_Count) %>%
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
plot_data_trimmed <- detailed_strain_comparison %>%
  dplyr::rename(Total = Total_Trimmed,
         Unique = Unique_Trimmed_Count,
         Shared = Shared_Strains) %>%
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
                    labels = c("Shared" = "Shared strains", "Unique" = "Unique (Trimmed)")) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Shared vs Unique (Trimmed)",
       x = "Sample ID", y = "Percent of Total Trimmed", fill = "")


## ---- Original percentages ----
plot_data_original <- detailed_strain_comparison %>%
  dplyr::rename(Total = Total_Original,
         Unique = Unique_Original_Count,
         Shared = Shared_Strains) %>%
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
                    labels = c("Shared" = "Shared strains", "Unique" = "Unique (Original)")) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Shared vs Unique (Original)",
       x = "Sample ID", y = "Percent of Total Original", fill = "")



# Combine into 2x2
final_plot <- (p1 | p2) / (p_trimmed |p_original)
final_plot
ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/Benchmarking1.pdf",width =30,height = 30,units = "cm")






# Directory with files
file_dir <- "C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/"

# List all ERR*_unique_trimmed.txt files
file_list <- list.files(file_dir, pattern="^ERR.*_unique_trimmed\\.txt$", full.names = TRUE)

# Function to read a single file
read_unique_trimmed <- function(file) {
  df <- read.delim(file, header=FALSE, stringsAsFactors = FALSE)
  
  df <- df %>%
    dplyr::mutate(
      SampleID = gsub("_unique_trimmed\\.txt$", "", basename(file)),
      # Extract t__SGB from V1
      t__SGB = sub(".*(t__SGB\\d+).*", "\\1", V1),
      # Use V3 as numeric y-axis value
      Count = as.numeric(V3)
    ) %>%
    dplyr::select(SampleID, t__SGB, Count)
  
  return(df)
}

# Load all files
all_data <- purrr::map_dfr(file_list, read_unique_trimmed)

# Optional: order SGB by overall Count
all_data <- all_data %>%
  dplyr::group_by(t__SGB) %>%
  dplyr::summarise(total_count = sum(Count, na.rm = TRUE)) %>%
  dplyr::arrange(desc(total_count)) %>%
  dplyr::mutate(t__SGB = factor(t__SGB, levels = t__SGB)) %>%
  dplyr::right_join(all_data, by = "t__SGB")

# Summarize total Count per strain for ordering
strain_order <- all_data %>%
  group_by(t__SGB) %>%
  summarise(TotalCount = sum(Count)) %>%
  arrange(desc(TotalCount)) %>%
  pull(t__SGB)

# Convert t__SGB to factor with levels ordered by total count
all_data <- all_data %>%
  mutate(t__SGB = factor(t__SGB, levels = strain_order))

# Plot
ggplot(all_data, aes(x = t__SGB, y = Count, fill = SampleID)) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "H") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top") +
  labs(x = "Strain (t__SGB)", 
       y = "Relative Abundance of Unique Trimmed Strain", 
       fill = "Sample ID",
       title = "Unique Trimmed Strain Relative Abundance Across Samples")


ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/Benchmarking2.pdf",width =25,height = 18,units = "cm")




p <- combined_strain_comparison %>%
  filter(Difference != 0) %>%
  ggplot(aes(x = SampleID, y = Difference)) +
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
        Q1 = quantile(Difference, 0.25),
        Q3 = quantile(Difference, 0.75),
        IQR = Q3 - Q1,
        is_outlier = Difference < (Q1 - 1.5 * IQR) | Difference > (Q3 + 1.5 * IQR)
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
    title = "Distribution of Strain Abundance Differences Between Original and Trimmed Samples",
    x = "Sample ID",
    y = "Difference in Relative Abundance (Trimmed - Original)"
    #caption = "Outliers defined as values beyond 1.5×IQR from the quartiles\nRed diamonds represent mean values"
  ) 
  # Optional: add faceting if you have many samples
  # facet_wrap(~ SampleID, scales = "free_x", ncol = 4) +
  # Adjust y-axis limits if needed
  # coord_cartesian(ylim = c(-2, 2))
  
  # Display the plot
  print(p)
  
  
ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/Benchmarking3.pdf",width =28,height = 20,units = "cm")
  
  

  # Reshape data
df_long <- combined_strain_comparison %>%
    dplyr::select(tlevel, Original, Trimmed) %>%
    tidyr::pivot_longer(cols = c(Original, Trimmed),
                        names_to = "Type", values_to = "Abundance") %>%
    dplyr::mutate(Type = factor(Type, levels=c("Original","Trimmed")))
  
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
         title="Original vs Trimmed Abundance for Top 10 Different Strains (>=3 Observations)")

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/04_analysis_results/benchmarking_results/new_result/Benchmarking4.pdf",width =30,height = 25,units = "cm")
  
  