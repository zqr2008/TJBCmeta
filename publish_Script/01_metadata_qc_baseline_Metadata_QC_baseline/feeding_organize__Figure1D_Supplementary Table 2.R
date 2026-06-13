# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 01_metadata_qc_baseline_Metadata_QC_baseline
# Figure(s): Figure1D
# Table(s): Supplementary Table 2
# Purpose: Generates quality-control summaries and QC figure/table outputs.
# Main input(s): who_m6.rda; who_y1.rda; newdata_phylo_NMrevision_0305.rds
# Main output(s): Figure1D.pdf; Supplementary Table 2 multivariate association analyses adjusting for potential confounders.csv
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
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
library(broom)


m6 <-read_sav("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/metadata/metadata/喂养数据250617/m6.sav")
y1 <- read_sav("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/metadata/metadata/喂养数据250617/y1.sav")

load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/growth_feeding_objects/who_m6.rda")
load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/growth_feeding_objects/who_y1.rda")

complete_phylo <- readRDS("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds")

#complete_phylo <- aggregate_taxa(complete_phylo, "Genus")

metadata <- as.data.frame(complete_phylo@sam_data)
class(metadata) <- "data.frame"


who_m6$IDchild <- as.character(who_m6$IDchild )


overlapping6month <- metadata %>%
  dplyr::rename(IDchild = sjid_kid) %>%
  left_join(who_m6,by ="IDchild")



overlapping1y <- metadata %>%
  dplyr::rename(IDchild = sjid_kid) %>%
  inner_join(who_y1,by ="IDchild")


overlapping6month$class_1 <- factor(overlapping6month$class_1)




m6$IDchild <- as.character(m6$IDchild)

smallmetadata <- metadata  %>%
  dplyr::rename(IDchild = sjid_kid) %>%
  dplyr::select(IDchild,class_1,antibiotic_usage_180day,feeding_type_28days,delivery_mode,gestational_week_delivery,TG,height_mo,class_1) 

studypheno6month <- smallmetadata %>%
  inner_join(m6,by ="IDchild")
  



df <- studypheno6month %>%
  mutate(class_1 = factor(class_1)) %>%
  mutate(feeding_types_6_months= case_when(feeding_m6 ==1~"exclusive_breastfeeding",
                                           feeding_m6 ==2~"exclusive_formula",
                                           feeding_m6 ==3~"mixed_feeding")) %>%
  dplyr::select(IDchild, antibiotic_usage_180day,feeding_type_28days ,feeding_types_6_months,delivery_mode,TG,height_mo,gestational_week_delivery,class_1) %>%
  distinct(IDchild,.keep_all = T) %>%
  mutate(class_1 = factor(class_1, levels = c(1,2,3))) %>%
  mutate(antibiotic_usage_180day=as.factor(antibiotic_usage_180day)) %>%
  mutate(delivery_mode = factor(delivery_mode,levels=c("Vaginal delivery","Caesarean section")))




df <- df %>%
  mutate(comparisontype1 = case_when(class_1==1~1,
                                     TRUE~0)) %>%
  mutate(comparisontype2 = case_when(class_1==2~1,
                                     TRUE~0)) %>%
  mutate(comparisontype3 = case_when(class_1==3~1,
                                     TRUE~0))
  
  
df$delivery_mode <- relevel(df$delivery_mode, ref = "Vaginal delivery" )

# Fit multinomial model
fit <- glm(comparisontype1 ~ antibiotic_usage_180day+feeding_type_28days + feeding_types_6_months +
                  delivery_mode + gestational_week_delivery+TG+height_mo,
           family = binomial(link = "logit"),
                data = df)

fit2 <- glm(comparisontype2 ~ antibiotic_usage_180day+feeding_type_28days + feeding_types_6_months +
             delivery_mode + gestational_week_delivery+TG+height_mo,
           family = binomial(link = "logit"),
           data = df)


fit3 <- glm(comparisontype3 ~ antibiotic_usage_180day+feeding_type_28days + feeding_types_6_months +
              delivery_mode + gestational_week_delivery+TG+height_mo,
            family = binomial(link = "logit"),
            data = df)



# Get tidy coefficients
tidy_fit <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(comparsion = "Trajectory 1 vs. non-Trajectory 1")

tidy_fit2 <- broom::tidy(fit2, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(comparsion = "Trajectory 2 vs. non-Trajectory 2")


tidy_fit3 <- broom::tidy(fit3, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(comparsion = "Trajectory 3 vs. non-Trajectory 3")


tidy_all <- rbind(tidy_fit,tidy_fit2) 
tidy_all <- rbind(tidy_all,tidy_fit3)


# Dodge width
dodge_width <- 0.6

# Add a small offset to place p-values slightly to the right of CI
pval_offset <- 1.05

tidy_all %>%
  filter(term!="(Intercept)") %>%
  filter(p.value<0.05) %>%
  ggplot(aes(x = estimate, y = term, color = comparsion)) +
  # Points
  geom_point(position = position_dodge(width = dodge_width), size = 3) +
  # Error bars
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = dodge_width),
                 height = 0.2) +
  # P-values aligned with error bars (dodged)
  geom_text(
    aes(
      x = conf.high * pval_offset,
      label = ifelse(
        p.value < 0.001,
        "italic(P) < 0.001",
        paste0("italic(P) == ", sprintf("%.3f", p.value))
      )
    ),
    parse = TRUE,
    position = position_dodge(width = dodge_width),
    hjust = 0,
    size = 3
  )+
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = c("Trajectory 1 vs. non-Trajectory 1" = "#e95280",
                                "Trajectory 2 vs. non-Trajectory 2" = "#23b1a5",
                                "Trajectory 3 vs. non-Trajectory 3" = "#E49B0F")) +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.3)))+
  labs(
    x = "Odds Ratio 95% CI",
    y = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )+
  scale_y_discrete(labels = c(
    "height_mo"="Maternal height",
    "antibiotic_usage_180day1"="Antibiotic use within 180 days",
    "feeding_type_28daysexclusive formula fed" = "Exclusive formula feeding at 28 days",
    "feeding_type_28daysmixed feeding" = "Mixed feeding at 28 days",
    "feeding_types_6_monthsexclusive_formula" = "Exclusive formula feeding at 6 months",
    "feeding_types_6_monthsmixed_feeding" = "Mixed feeding at 6 months",
    "delivery_modeCaesarean section" = "C-section",
    "gestational_week_delivery" = "Gestational week at delivery"
  ))+
  guides(
    color = guide_legend(nrow = 3),
    fill = guide_legend(nrow = 3)
  )

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure1D.pdf",width =20,height = 10,units = "cm")


write.table(tidy_all, file = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision//Supplementary Table 2 multivariate association analyses adjusting for potential confounders.csv",
            fileEncoding = "GBK",row.names = F,col.names = T,
            quote = F,sep = ",")
