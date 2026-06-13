# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): 20260608第三次反馈结果_24株菌.xlsx; 20260529strain1_4.xlsx; 20260601strain5_8.xlsx; firstbatch.rda
# Main output(s): allgrow.rda; allgrowth.pdf
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END

rm(list = ls())
library(readxl)


X20260608第三次反馈结果_24株菌 <- read_excel("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/secondbatch/20260608第三次反馈结果_24株菌.xlsx", 
                                    skip = 20)

X20260529strain1_4 <- read_excel("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/secondbatch/20260529strain1_4.xlsx") %>%
  dplyr::select(16:28) %>%
  filter(is.na(培养基...16)==FALSE)


X20260601strain5_8 <- read_excel("C:/Users/zqr20/Documents/tjmeta/BIG_revision/01_raw_inputs/strain_inputs/secondbatch/20260601strain5_8.xlsx") %>%
dplyr::select(16:28)%>%
  filter(is.na(培养基...16)==FALSE)

library(tidyverse)

make_allstrain <- function(df) {
  
  condition_col <- names(df)[1]
  od_cols <- names(df)[-1]
  
  if (length(od_cols) %% 3 != 0) {
    stop("The number of OD600 columns is not divisible by 3. Please check the input dataframe.")
  }
  
  # Build strict column map:
  # every 3 columns = one strain
  col_map <- tibble(
    raw_col = od_cols,
    col_index = seq_along(od_cols),
    strain_block = ceiling(col_index / 3),
    rep_id = paste0("rep", ((col_index - 1) %% 3) + 1)
  )
  
  # Extract strain name only from the first column of each 3-column block
  strain_names <- col_map %>%
    filter(rep_id == "rep1") %>%
    mutate(
      strainname = str_extract(raw_col, "^Strain\\d+")
    ) %>%
    select(strain_block, strainname)
  
  col_map <- col_map %>%
    left_join(strain_names, by = "strain_block") %>%
    select(raw_col, strainname, rep_id, col_index, strain_block)
  
  if (any(is.na(col_map$strainname))) {
    print(col_map, n = Inf)
    stop("Some strain names are NA. Please check whether the first column of each 3-column block starts with 'Strain'.")
  }
  
  rep_check <- col_map %>%
    count(strainname)
  
  if (any(rep_check$n != 3)) {
    print(rep_check)
    stop("Some strains do not have exactly 3 replicate columns.")
  }
  
  allstrain_long <- df %>%
    rename(condition = all_of(condition_col)) %>%
    pivot_longer(
      cols = all_of(od_cols),
      names_to = "raw_col",
      values_to = "OD600"
    ) %>%
    left_join(col_map, by = "raw_col") %>%
    mutate(
      OD600 = as.numeric(OD600)
    ) %>%
    select(
      condition,
      strainname,
      rep_id,
      OD600
    )
  
  allstrain <- allstrain_long %>%
    pivot_wider(
      id_cols = c(condition, strainname),
      names_from = rep_id,
      values_from = OD600
    ) %>%
    rowwise() %>%
    mutate(
      mean_OD600_blank_subtracted = mean(c(rep1, rep2, rep3), na.rm = TRUE),
      sd_OD600_blank_subtracted = sd(c(rep1, rep2, rep3), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(
      condition,
      strainname,
      rep1,
      rep2,
      rep3,
      mean_OD600_blank_subtracted,
      sd_OD600_blank_subtracted
    )
  
  return(allstrain)
}



allstrain1_4 <- make_allstrain(X20260529strain1_4)

allstrain5_8 <- make_allstrain(X20260601strain5_8)

allstrainother <- make_allstrain(X20260608第三次反馈结果_24株菌)

allstrain <- bind_rows(
  allstrain1_4,
  allstrain5_8,
  allstrainother
)


allstrain <- allstrain %>%
  mutate(
    condition_formal = recode(
      condition,
      "2'-FL" = "2'-FL",
      "3'-FL" = "3'-FL",
      "LNT" = "LNT",
      "LNnT" = "LNnT",
      "6'-SL" = "6'-SL",
      "HMOS" = "HMO Mix",
      "无葡萄糖MRS" = "sMRS",
      "MRS" = "MRS"
    )
  ) 



# Order strains by 3-FL OD600 value, highest to lowest
strain_order <- allstrain %>%
  filter(condition_formal == "3'-FL") %>%
  arrange(desc(mean_OD600_blank_subtracted)) %>%
  pull(strainname)

# Optional: keep HMO/media order fixed on x-axis
condition_order <- c("2'-FL", "3'-FL", "LNT", "LNnT", "6'-SL", "HMO Mix", "sMRS", "MRS")

allstrain_plot <- allstrain %>%
  mutate(
    strainname = factor(strainname, levels = strain_order),
    condition_formal = factor(condition_formal, levels = condition_order)
  )




load("C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/strain_growth_objects/firstbatch.rda")

allstrain_plot <- rbind(allstrain,firstbatch)


strain_order <- allstrain_plot %>%
  filter(condition_formal == "3'-FL") %>%
  arrange(desc(mean_OD600_blank_subtracted)) %>%
  pull(strainname)

allstrain_plot <- allstrain_plot %>%
  mutate(
    strainname = factor(strainname, levels = strain_order),
    condition_formal = factor(condition_formal, levels = condition_order)
  ) %>%
  filter(condition_formal !="sMRS")


ggplot(allstrain_plot,
       aes(
         x = condition_formal,
         y = mean_OD600_blank_subtracted,
         fill = strainname
       )) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(
    aes(
      ymin = mean_OD600_blank_subtracted - sd_OD600_blank_subtracted,
      ymax = mean_OD600_blank_subtracted + sd_OD600_blank_subtracted
    ),
    width = 0.2
  ) +
  facet_wrap(
    ~ strainname,
    nrow = 1,
    scales = "fixed"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = expression("24h OD"[600] ~ "(blank-subtracted)")
  )

save(allstrain_plot,file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/strain_growth_objects/allgrow.rda")

ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/allgrowth.pdf",width = 40, height = 8)
