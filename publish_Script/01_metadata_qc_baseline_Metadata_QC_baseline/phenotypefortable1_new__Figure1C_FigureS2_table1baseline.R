# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 01_metadata_qc_baseline_Metadata_QC_baseline
# Figure(s): Figure1C; FigureS2
# Table(s): table1baseline
# Purpose: Generates quality-control summaries and QC figure/table outputs.
# Main input(s): newdata_phylo_NMrevision_0305.rds
# Main output(s): Figure1C.pdf; FigureS2 whzscore_sparate.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
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

conflicted::conflicts_prefer(dplyr::first)

setwd("C:/Users/zqr20/Documents/tjmeta/BIG_revision/")

complete_phylo <-  readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

complete_phylo_pheno <- as.data.frame(complete_phylo@sam_data)
class(complete_phylo_pheno) <- "data.frame"


complete_phylo_pheno <- remove_labels(complete_phylo_pheno)


showtable1 <- c("class_1","Age_mo","BMI_mo","TG",
                "LDLC","HDLC","height_mo","height_fa","edu_mo","smoking_mo",
                "drinking_mo_18topreg","GDM_group2023","SBP_first_check",
                "DBP_first_check","delivery_mode_kid","preterm","gestational_week_delivery",
                "gender_kid","birth_weight","birth_length",
                "SGA","LGA","macrosomia","feeding_type_28days","antibiotic_usage_180day")


fortable1 <- complete_phylo_pheno %>%
  distinct(familyid,.keep_all = T) %>%
  dplyr::select(any_of(showtable1)) %>%
  mutate(gestational_week_delivery = as.numeric(gestational_week_delivery)) %>%
  mutate(class_1= paste0("Trajectory ", class_1))



#################################################################
#overall description, only for writing, please do not delete
#################################################################
tableaddinfor2 <- complete_phylo_pheno %>%
  distinct(familyid, .keep_all = TRUE) %>%
  group_by(sjid_kid) %>%
  summarise(
    trajcluster = first(class_1),
    obsesityby36month = case_when(
      any(c_across(zwhz_1:zwhz_7) > 2, na.rm = TRUE) ~ "overweight",
      any(c_across(zwhz_1:zwhz_7) < -2, na.rm = TRUE) ~ "lean",
      TRUE ~ "normal"
    ),
    .groups = "drop"
  ) %>%
  tbl_summary(
    by = trajcluster,
    include = obsesityby36month
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct"))


ft <- as_flex_table(tableaddinfor2)
save_as_docx("Trajectory healthy" = ft, path = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/07_manuscript_and_submission/manuscript_versions/table1STrajectory_healthy.docx")

plot_df <- complete_phylo_pheno %>%
  distinct(familyid, .keep_all = TRUE) %>%
  group_by(sjid_kid) %>%
  summarise(
    trajcluster = first(class_1),
    obesityby36month = case_when(
      any(c_across(zwhz_1:zwhz_7) > 2, na.rm = TRUE) ~ "Overweight",
      any(c_across(zwhz_1:zwhz_7) < -2, na.rm = TRUE) ~ "Underweight",
      TRUE ~ "Normal"
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(trajcluster), !is.na(obesityby36month)) %>%
  count(trajcluster, obesityby36month, name = "n") %>%
  group_by(trajcluster) %>%
  mutate(
    total_n = sum(n),
    Percent = 100 * n / total_n,
    label = paste0(n, "\n(", round(Percent, 1), "%)")
  ) %>%
  ungroup() %>%
  mutate(
    Trajectory = paste0("Trajectory ", trajcluster),
    obesityby36month = factor(
      obesityby36month,
      levels = c("Underweight", "Normal", "Overweight")
    )
  )

matched_colors <- c(
  "Underweight" = "#C51B7D",  # strong magenta
  "Normal"      = "#1B9E77",  # green-teal
  "Overweight"  = "#D95F02"   # burnt orange
)

ggplot(plot_df, aes(x = Trajectory, y = Percent, fill = obesityby36month)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = matched_colors) +
  labs(
    x = "Trajectory",
    y = "Proportion (%)",
    fill = "Obesity Category"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure1C.pdf",width = 4, height = 6)



tableaddinfor2separate <- complete_phylo_pheno %>%
  distinct(familyid, .keep_all = TRUE) %>%
  group_by(sjid_kid) %>%
  summarise(
    trajcluster = first(class_1),
    across(
      zwhz_1:zwhz_7,
      ~ case_when(
        .x > 2 ~ "overweight",
        .x < -2 ~ "lean",
        TRUE ~ "normal"
      ),
      .names = "cat_{col}"   # will generate cat_zwhz_1 … cat_zwhz_7
    ),
    .groups = "drop"
  ) %>%
  tbl_summary(
    by = trajcluster,
    include = starts_with("cat_")
  ) %>%
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct"))






# Step 1: classify each zwhz_* into categories
zvars <- paste0("zwhz_", 1:7)

complete_phylo_pheno_long <- complete_phylo_pheno %>%
  distinct(familyid, .keep_all = TRUE) %>%
  dplyr::select(sjid_kid, class_1, all_of(zvars)) %>%
  pivot_longer(cols = all_of(zvars), names_to = "timepoint", values_to = "zwhz") %>%
  mutate(
    obesity_cat = case_when(
      zwhz > 2  ~ "Overweight",
      zwhz < -2 ~ "Underweight",
      !is.na(zwhz) ~ "Normal",
      TRUE ~ NA_character_
    )
  )

# Step 2: summarize each timepoint by trajcluster
tableaddinfor2 <- complete_phylo_pheno_long %>%
  group_by(class_1, timepoint) %>%
  summarise(
    n = n(),
    overweight = sum(obesity_cat == "Overweight", na.rm = TRUE),
    lean = sum(obesity_cat == "Underweight", na.rm = TRUE),
    normal = sum(obesity_cat == "Normal", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Overweight_prop = overweight / n,
    Underweight_prop = lean / n,
    Normal_prop = normal / n
  )


matched_colors <- c(
  "Underweight" = "#C51B7D",  # strong magenta
  "Normal"      = "#1B9E77",  # green-teal
  "Overweight"  = "#D95F02"   # burnt orange
)
plot_df <- tableaddinfor2 %>%
  dplyr::select(class_1, timepoint, Underweight_prop, Normal_prop, Overweight_prop) %>%
  pivot_longer(cols = ends_with("_prop"),
               names_to = "category",
               values_to = "prop") %>%
  mutate(category = factor(category,
                           levels = c("Underweight_prop","Normal_prop","Overweight_prop"),
                           labels = c("Underweight","Normal","Overweight")))



ggplot(plot_df, aes(x = timepoint, y = prop, fill = category, group = category)) +
  geom_col(position = "fill") +
  
  # proportion labels next to the bar
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_fill(vjust = 0.05), # slightly above the bar
            size = 4) +
  # connecting lines across timepoints
  geom_line(aes(y = prop, color = category),
            position = position_fill(),
            size = 0.7,
            show.legend = FALSE) +
  scale_fill_manual(values = matched_colors)+
  
  facet_wrap(~ class_1, ncol = 3) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Measurement",
    y = "Proportion",
    fill = "Phenotype"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/FigureS2 whzscore_sparate.pdf",width = 14, height = 8)






table1 <- fortable1 %>%
  mutate(
    across(where(is.character), ~ stringr::str_to_sentence(.x)),
    across(where(is.factor), ~ factor(stringr::str_to_sentence(as.character(.x))))
  ) %>%
  tbl_summary(
    by = class_1,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ c(2, 2),
      all_categorical() ~ c(0, 2)
    ),
    label = list(
      Age_mo ~ "Age, years",
      BMI_mo ~ "Prepregnancy BMI, kg/m\u00B2",
      height_mo ~ "Maternal height, cm",
      height_fa ~ "Paternal height, cm",
      edu_mo ~ "Education level",
      smoking_mo ~ "Current or former smoking",
      drinking_mo_18topreg ~ "Drinking during pregnancy",
      GDM_group2023 ~ "GDM",
      SBP_first_check ~ "SBP at first prenatal check-up, mmHg",
      DBP_first_check ~ "DBP at first prenatal check-up, mmHg",
      delivery_mode_kid ~ "Delivery mode",
      preterm ~ "Preterm",
      gestational_week_delivery ~ "Week of gestation",
      gender_kid ~ "Gender of child",
      birth_weight ~ "Birth weight of child, g",
      birth_length ~ "Birth length of child, cm",
      macrosomia ~ "Macrosomia",
      feeding_type_28days ~ "Feeding type",
      antibiotic_usage_180day ~ "Antibiotic exposure during the first 180 days of life"
    )
  ) %>%
  add_n() %>%
  add_p(
    test = list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "chisq.test.no.correct"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_header(p.value ~ "***P* value**") %>%
  bold_labels() %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "gestational_week_delivery" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "gestational_week_delivery" & row_type == "label",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "gestational_week_delivery" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "birth_weight" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "birth_weight" & row_type == "label",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "birth_weight" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "delivery_mode_kid" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "delivery_mode_kid" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "preterm" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "preterm" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "SGA" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "SGA" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "LGA" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "LGA" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "macrosomia" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "macrosomia" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "height_mo" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "height_mo" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_table_styling(
    column = stat_3,
    rows = variable == "TG" & row_type == "label",
    footnote = "Significantly different between Trajectory 1 and Trajectory 3"
  ) %>%
  
  modify_footnote(
    update = label ~ paste(
      "Continuous variables were compared using the Kruskal–Wallis rank-sum test for the overall group difference.",
      "Categorical variables were compared using Pearson's Chi-squared test for the overall group difference.",
      "Pairwise comparisons were subsequently performed between two trajectories.",
      "BMI, body mass index;",
      "GDM, gestational diabetes mellitus;",
      "PGDM, Pregestational Diabetes Mellitus;",
      "SBP, systolic blood pressure;",
      "DBP, diastolic blood pressure;",
      "SGA, small for gestational age;",
      "LGA, large for gestational age."
    )
  )
ft <- as_flex_table(table1)

ft <- ft %>%
  italic(j = "p.value", part = "body") %>%
  bold(j = "p.value", part = "body") %>%
  italic(j = "p.value", part = "header") %>%
  bold(j = "p.value", part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

doc <- officer::read_docx()

doc <- officer::body_add_par(
  doc,
  value = "Table 1. Baseline characteristics",
  style = "Normal"
)

doc <- flextable::body_add_flextable(doc, value = ft)

print(
  doc,
  target = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/table1baseline.docx"
)



#############################################################
#below are tests used for annotated for above
#############################################################

pairwise.t.test(fortable1[["gestational_week_delivery"]], fortable1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortable1[["birth_weight"]], fortable1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortable1[["height_mo"]], fortable1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortable1[["TG"]], fortable1[["trajcluster"]], p.adj = "BH")



xtab <- table(fortable1$trajcluster,fortable1$delivery_mode_kid)
pairwise_prop_test(xtab,p.adjust.method = "BH")

xtab <- table(fortable1$trajcluster,fortable1$preterm)
pairwise_prop_test(xtab,p.adjust.method = "BH")


xtab <- table(fortable1$trajcluster,fortable1$SGA)
pairwise_prop_test(xtab,p.adjust.method = "BH")

xtab <- table(fortable1$trajcluster,fortable1$LGA)
pairwise_prop_test(xtab,p.adjust.method = "BH")

xtab <- table(fortable1$trajcluster,fortable1$macrosomia)
pairwise_prop_test(xtab,p.adjust.method = "BH")

#exceeding two columns, can not compute
#xtab <- table(fortable1$trajcluster,fortable1$feeding_type_28days)
#pairwise_prop_test(xtab)

xtab <- table(fortable1$trajcluster,fortable1$antibiotic_usage_180day)
pairwise_prop_test(xtab,p.adjust.method = "BH")


####################################################################
#for reviewer 3's comments
####################################################################
fortable1 %>%ggplot(aes(x= trajcluster, y=gestational_week_delivery,fill=trajcluster))+
  geom_boxplot() +
  theme_classic()+
  scale_fill_manual(values = c("Trajectory 1" = "#e95280", "Trajectory 2" = "#23b1a5", "Trajectory 3" = "#E49B0F")) +
  scale_color_manual(values = c("Trajectory 1" = "#e95280", "Trajectory 2" = "#23b1a5", "Trajectory 3" = "#E49B0F")) +
  ylab("Week of gestation")+
  xlab("Trajectory")+
  theme(
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  labs(fill = "Trajectory")+
  theme (plot.title = element_text (size = 16, face = "bold" ))



fortable1 %>%
  ggplot(aes(x = gestational_week_delivery, fill = trajcluster)) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Trajectory 1" = "#e95280",
    "Trajectory 2" = "#23b1a5",
    "Trajectory 3" = "#E49B0F"
  )) +
  xlab("Week of gestation") +
  ylab("Density") +
  theme(
    strip.text = element_text(size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(fill = "Trajectory")

fortableS1 <- complete_phylo_pheno %>%
  distinct(familyid,.keep_all = T) %>%
  dplyr::select(trajcluster,ends_with("_1"),ends_with("_2"),
                ends_with("_3"),ends_with("_4"),ends_with("_5")) %>%
  dplyr::select(-c(starts_with("zBMI")))




tables1 <- fortableS1 %>%
  tbl_summary(
    by = trajcluster, # split table by group
    missing = "no" # don't list missing data separately
    ) |> 
  add_n() |> # add column with total number of non-missing observations
  add_p(test = list(all_categorical() ~ "chisq.test.no.correct")) |> # test for a difference between groups
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels() %>%
  modify_caption("**Table S1. Trajectories Characteristics**") |>
  modify_table_styling(
    column = stat_1,
    rows = variable == "weight_1",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "weight_1",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "weight_1",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "height_1",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "height_1",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "height_1",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zhaz_1",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zhaz_1",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zhaz_1",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwaz_1",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwaz_1",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwaz_1",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%

  modify_table_styling(
    column = stat_1,
    rows = variable == "zwhz_1",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwhz_1",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwhz_1",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  )   %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "weight_2",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "weight_2",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "weight_2",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "height_2",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "height_2",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "height_2",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zhaz_2",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zhaz_2",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zhaz_2",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  )%>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwaz_2",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwaz_2",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwaz_2",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwhz_2",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwhz_2",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwhz_2",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "weight_3",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "weight_3",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "weight_3",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "height_3",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "height_3",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "height_3",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zhaz_3",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zhaz_3",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zhaz_3",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwhz_3",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwhz_3",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwhz_3",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwaz_3",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwaz_3",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwaz_3",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "weight_4",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "weight_4",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "weight_4",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "height_4",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "height_4",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "height_4",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zhaz_4",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zhaz_4",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zhaz_4",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  )%>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwaz_4",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwaz_4",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwaz_4",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwhz_4",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwhz_4",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwhz_4",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "weight_5",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "weight_5",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "weight_5",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "height_5",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "height_5",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "height_5",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zhaz_5",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zhaz_5",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zhaz_5",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  )%>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwaz_5",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwaz_5",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwaz_5",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
    column = stat_1,
    rows = variable == "zwaz_5",
    footnote = "Significantly different between Trajectory 1 and Trajectory 2"
  ) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwaz_5",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwaz_5",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  modify_table_styling(
  column = stat_1,
  rows = variable == "zwhz_5",
  footnote = "Significantly different between Trajectory 1 and Trajectory 2"
) %>%
  modify_table_styling(
    column = stat_2,
    rows = variable == "zwhz_5",
    footnote = "Significantly different between Trajectory 2 and Trajectory 3"
  ) %>%
  modify_table_styling(
    column = stat_3,
    rows = variable == "zwhz_5",
    footnote = "Significantly different between Trajectory 3 and Trajectory 1"
  ) %>%
  as_gt() |> 
  gt::gtsave(filename = "tableS1baseline.docx")



pairwise.t.test(fortableS1[["delivery_mode_kid"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["height_1"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zhaz_1"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zwaz_1"]], fortableS1[["trajcluster"]], p.adj = "BH")


pairwise.t.test(fortableS1[["weight_2"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["height_2"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zhaz_2"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zwaz_2"]], fortableS1[["trajcluster"]], p.adj = "BH")



pairwise.t.test(fortableS1[["weight_3"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["height_3"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zhaz_3"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zwaz_3"]], fortableS1[["trajcluster"]], p.adj = "BH")


pairwise.t.test(fortableS1[["weight_4"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["height_4"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zhaz_4"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zwaz_4"]], fortableS1[["trajcluster"]], p.adj = "BH")




pairwise.t.test(fortableS1[["weight_5"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["height_5"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zhaz_5"]], fortableS1[["trajcluster"]], p.adj = "BH")
pairwise.t.test(fortableS1[["zwaz_5"]], fortableS1[["trajcluster"]], p.adj = "BH")

