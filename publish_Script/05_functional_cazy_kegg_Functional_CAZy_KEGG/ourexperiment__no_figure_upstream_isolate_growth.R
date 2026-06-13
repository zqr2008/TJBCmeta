# >>> PUBLISH_ANNOTATION_START
# Generated: 2026-06-13
# Category: 05_functional_cazy_kegg_Functional_CAZy_KEGG
# Figure(s): no direct figure
# Table(s): No direct submitted table matched from filename
# Purpose: Processes CAZy functional profiles and Bifidobacterium functional figures/tables.
# Main input(s): none detected
# Main output(s): firstbatch.rda
# Static warning(s): 0 potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv
# Scope note: annotation added to this published copy only; original script body below is unchanged.
# <<< PUBLISH_ANNOTATION_END
rm(list = ls())
library(xml2)
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)

pzfx <- r"(C:\Users\zqr20\Documents\tjmeta\BIG_revision\01_raw_inputs\strain_inputs\假小链原始数据\假小链26.3.27.pzfx)"
doc <- read_xml(pzfx)
ns <- xml_ns(doc)

tables <- xml_find_all(doc, ".//d1:Table", ns)


allstrain <- data.frame()
for (tb in tables) {
  strain <- xml_text(xml_find_first(tb, "./d1:Title", ns), trim = TRUE)
  if (strain == "") next
  
  ycols <- xml_find_all(tb, "./d1:YColumn", ns)
  
  out <- lapply(ycols, function(yc) {
    condition <- xml_text(xml_find_first(yc, "./d1:Title", ns), trim = TRUE)
    vals <- as.numeric(xml_text(xml_find_all(yc, ".//d1:d", ns)))
    
    data.frame(
      condition = condition,
      rep1 = vals[1],
      rep2 = vals[2],
      rep3 = vals[3],
      mean_OD600_blank_subtracted = mean(vals, na.rm = TRUE),
      sd_OD600_blank_subtracted = sd(vals),
      stringsAsFactors = FALSE
    )
  }) |> bind_rows() %>%
    mutate(strainname = strain)
    
  allstrain <- rbind(allstrain,out)
  safe <- gsub("[^A-Za-z0-9_.-]", "_", strain)

}



allstrain <- allstrain %>%
  mutate(
    condition_formal = recode(
      condition,
      "2-FL" = "2'-FL",
      "3-FL" = "3'-FL",
      "LNT" = "LNT",
      "LNnT" = "LNnT",
      "6-SL" = "6'-SL",
      "HMOs(混合糖）" = "HMO Mix",
      "sMRS" = "sMRS",
      "MRS" = "MRS"
    )
  ) #%>%
  #filter(strainname != "ZHB-PYG-1") %>%
  #filter(strainname != "TQ-M-4")


strain_colors <- c(
  "GBWPFCBAD1"  = "#D94E7A",  # muted rose
  "CBXQCBAD1"   = "#2AA89C",  # teal
  "Y0808-CBA-1" = "#4C72B0",  # soft blue
  "YCRPC1-N-11" = "#55A868",  # muted green
  "YCRPC9-N-11" = "#C44E52",  # brick red
  "YF2-M136-4"  = "#8172B2"
  
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
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = expression("24h OD"[600] ~ "(blank-subtracted)")
   )

firstbatch <- allstrain_plot 
save(firstbatch, file="C:/Users/zqr20/Documents/tjmeta/BIG_revision/02_processed_objects/strain_growth_objects/firstbatch.rda")


#ggsave("C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision/Figure6B OD.pdf",width = 14, height = 8)
