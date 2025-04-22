library(readr)
library(dplyr)
library(ggplot2)
library(cutpointr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  args[2] = "out.txt"
}

metadata <- args[1]
nGD <- args[2]
outdir <- args[3]

# Reading the metadata
md <- read_tsv(file = metadata, col_names = T, show_col_types = F)
# reading the tsv table of pairwise distances
nGD <- read_tsv(file = nGD, col_names = F, show_col_types = F)

# Adding the metadata to the table
nGD <- left_join(nGD %>% select(sampleID_1 = X1, everything()),
                 md %>% select(sampleID_1 = sampleID,
                               subjectID_1 = subjectID,
                               family1 = family,
                               days_from_first_collection_1 = days_from_first_collection,
                               Dataset_1 = Dataset))
nGD <- left_join(nGD %>% select(sampleID_2 = X2, everything()),
                 md %>% select(sampleID_2 = sampleID,
                               subjectID_2 = subjectID,
                               family2 = family,
                               days_from_first_collection_2 = days_from_first_collection,
                               Dataset_2 = Dataset))

# Computing time difference between sample (important for longitudinal samples)
nGD$time_day_diff <- abs(nGD$days_from_first_collection_1 - nGD$days_from_first_collection_2)

# Annotating pairs of samples. Are they related? Are they from the same individual?
nGD$same_individual <- ifelse(nGD$subjectID_1 == nGD$subjectID_2, "same_individual", "different_individual")
nGD$related <- ifelse(nGD$family1 == nGD$family2, "related", "unrelated")
nGD$same_dataset <- ifelse(nGD$Dataset_1 == nGD$Dataset_2, "same_dtset", "different_dtset")


# Keeping only the training data
nGD_training <- nGD 

nGD_training <- rbind(nGD_training %>% filter(same_individual == "same_individual") %>%
                 filter(time_day_diff <= 180) %>%
                 group_by(subjectID_1) %>% arrange(subjectID_1, time_day_diff) %>% slice_head(n = 1) %>% ungroup(),
             nGD_training %>% filter(same_individual != "same_individual",  related == "unrelated"
                                     ) %>%
                 group_by(subjectID_1, subjectID_2) %>% slice_head(n = 1) %>% ungroup())


table(nGD_training$same_individual)

res_youden <- cutpointr(data = nGD_training, x = X3, class = same_individual, method = maximize_metric, metric = youden)
quantile_5pc <- nGD_training %>% filter(related == "unrelated") %>% pull(X3) %>% quantile(0.05)
quantile_3pc <- nGD_training %>% filter(related == "unrelated") %>% pull(X3) %>% quantile(0.03)
fig <- ggplot(data = nGD_training) +
    geom_density(aes(x = X3, fill = same_individual), alpha = 0.5 ) +
    geom_vline(aes(xintercept = res_youden$optimal_cutpoint, color = "youden")) +
    geom_vline(aes(xintercept = quantile_5pc, color = "quantile_5pc"), linetype = "dotted", lwd = 1) +
    geom_vline(aes(xintercept = quantile_3pc, color = "quantile_3pc"), linetype = "dotdash", lwd = 1) +
    theme_minimal() + xlab("StrainPhlAn nGD") + ylab("") +
    ggtitle("") +
    theme(legend.title =  element_blank(), legend.position = "bottom") +
    scale_color_manual(name = "statistics", values = c(youden = "blue", quantile_5pc = "red", quantile_3pc = "green"))
pdf(paste(outdir,"youden.pdf",sep="/"),width = 8,height = 4)
print(fig)
dev.off()

out <- data.frame("prevalence"=table(nGD_training$same_individual)["same_individual"],
                  "youden_index"=res_youden$optimal_cutpoint,
                  "quantile_5pc"=quantile_5pc, "quantile_3pc"=quantile_3pc)
write.table(out,paste(outdir,"youden_result.tsv", sep="/"), row.names = FALSE, quote = FALSE, sep="\t")

threshold=NA
if(table(nGD_training$same_individual)["same_individual"]>=50){
        if(res_youden$optimal_cutpoint <= quantile_5pc){
                res=prop.table(table(nGD_training %>% filter(related == "unrelated") %>% pull(X3) >= res_youden$optimal_cutpoint))
                print(res)
                print(table(nGD_training %>% filter(related == "unrelated") %>% pull(X3) >= res_youden$optimal_cutpoint))
                threshold=res["FALSE"]
        }else{
                threshold=0.05
        }
}else{
        threshold=0.03
}

fileConn<-file(paste(outdir,"threshold.tsv", sep="/"))
writeLines(as.character(threshold), fileConn)
close(fileConn)

