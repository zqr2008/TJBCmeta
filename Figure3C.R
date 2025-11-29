rm(list = ls())
library(smotefamily)
library(ranger)
library(pROC)
library(dplyr)
library(caret)
library(sjmisc)
library(tibble)
library(tidyverse)

complete_phylo <- readRDS("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/ps_filtered.rds")



otudf <- as.data.frame(complete_phylo@otu_table)

trajectory <- as.data.frame(complete_phylo@sam_data)

class(trajectory) <- "data.frame"
class(otudf) <- "data.frame"

otudf <- otudf %>%
  rownames_to_column("sampleid")

trajectory <- trajectory %>%
  rownames_to_column("sampleid") %>%
  left_join(otudf,by = "sampleid")   %>%
  mutate(trajcluster= as.factor(class_1))


# Identify numeric columns that are also in otudf
num_cols <- intersect(
  colnames(otudf),
  trajectory %>% 
    dplyr::select(where(is.numeric)) %>% colnames()
)

# Apply log10 transformation with pseudo-count
trajectory <- trajectory %>%
  mutate(across(all_of(num_cols), ~ log10(. + 1e-6)))



trajectorymaternal <- trajectory %>%
  filter(str_detect(visit,"MG"))

trajectory <- trajectory %>% 
  filter(str_detect(visit,"BM"))




# -------------------------------
# Function to train RF + SMOTE + get ROC
# -------------------------------
# -------------------------------
# Function to train RF + SMOTE + get ROC
# -------------------------------
get_roc_data <- function(df, species_cols, outcome_col, outcome_value, train_frac = 0.8, ntree = 1000) {
  
  # 1️⃣ Make binary outcome
  df <- df %>%
    mutate(temp_class = as.factor(ifelse(.data[[outcome_col]] == outcome_value, 1, 0)))
  
  # 2️⃣ Prepare model data
  data_model <- df %>%
    dplyr::select(all_of(species_cols), temp_class) %>%
    dplyr::select(-sampleid)
  
  # 3️⃣ Split train/test
  trainIndex <- createDataPartition(data_model$temp_class, p = train_frac, list = FALSE)
  train <- data_model[trainIndex, ]
  test  <- data_model[-trainIndex, ]
  
  train$temp_class <- as.factor(train$temp_class)
  
  # 4️⃣ Apply SMOTE
  train_smote <- SMOTE(train[, !(names(train) %in% "temp_class")], train$temp_class, K = 5, dup_size = 5)
  train_balanced <- train_smote$data
  colnames(train_balanced)[ncol(train_balanced)] <- "temp_class"
  train_balanced$temp_class <- as.factor(train_balanced$temp_class)
  
  # 5️⃣ Train random forest
  rf <- ranger(
    temp_class ~ ., 
    data = train_balanced,
    importance = "impurity",
    probability = TRUE,
    num.trees = ntree,
    num.threads = parallel::detectCores() - 1
  )

  # 6️⃣ Predict probabilities on test
  test_num <- test %>% dplyr::select(-temp_class)
  pred <- predict(rf, data = test_num)
  prob <- pred$predictions[, "1"]
  
  # 7️⃣ Compute ROC and AUC
  roc_obj <- roc(test$temp_class, prob)
  auc_value <- auc(roc_obj)
  
  # 7️⃣1️⃣ Compute 95% CI for AUC
  auc_ci <- ci.auc(roc_obj)  # default 95% CI
  ci_lower <- auc_ci[1]
  ci_upper <- auc_ci[3]
  
  # 8️⃣ Extract ROC data for ggplot
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    outcome = paste0("Trajectory ", outcome_value, 
                     " (AUC=", round(auc_value,3),
                     ", 95% CI [", round(ci_lower,3), "-", round(ci_upper,3), "])")
  )
  
  # 9️⃣ Return list with ROC, AUC, and CI
  list(
    roc_df = roc_df,
    auc = auc_value,
    auc_ci_lower = ci_lower,
    auc_ci_upper = ci_upper,
    rf_model = rf,
    test_data = test
  )
}


# -------------------------------
# Species columns
# -------------------------------
species_cols <- colnames(otudf)

# -------------------------------
# Run for class_1 == 3
# -------------------------------
roc3 <- get_roc_data(trajectorymaternal, species_cols, outcome_col = "class_1", outcome_value = 3)

# -------------------------------
# Run for class_1 == 1
# -------------------------------
roc1 <- get_roc_data(trajectorymaternal, species_cols, outcome_col = "class_1", outcome_value = 1)

# -------------------------------
# Combine ROC data
# -------------------------------
roc_combined <- rbind(roc3$roc_df, roc1$roc_df)

# -------------------------------
# Plot both ROC curves
# -------------------------------
ggplot(roc_combined, aes(x = fpr, y = tpr, color = outcome)) +
  geom_line(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("#e95280","#E49B0F")) +
  labs(
    title = "Maternal microbiome species abundance predict trajectory",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Outcome"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )+
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )


importance_df <- data.frame(
  species = names(roc3$rf_model$variable.importance),
  importance = roc3$rf_model$variable.importance
) %>%
  arrange(desc(importance))


top20 <- importance_df %>% slice_max(order_by = importance, n = 20)

ggplot(top20, aes(x = reorder(species, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 20 Predictive Species for Trajectory 3",
    x = "Species",
    y = "Variable Importance"
  ) +
  theme_classic(base_size = 13)




importance_df <- data.frame(
  species = names(roc1$rf_model$variable.importance),
  importance = roc1$rf_model$variable.importance
) %>%
  arrange(desc(importance))


top20 <- importance_df %>% slice_max(order_by = importance, n = 20)

ggplot(top20, aes(x = reorder(species, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 20 Predictive Species for Trajectory 1",
    x = "Species",
    y = "Variable Importance"
  ) +
  theme_classic(base_size = 13)

