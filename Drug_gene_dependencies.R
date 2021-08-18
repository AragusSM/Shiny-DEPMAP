
library(ggplot2)
library(data.table)
library(ggpubr)
install.packages(c("GGally", "network"))
theme_set(theme_pubr())
## load dependency data ##
setwd("C:/Users/ningz/OneDrive/Desktop/My files/DrugResponse2021/crispr_dep_map")

# Achilles+Sanger dependency data
CRISPR_dep <- read.delim(file = "CRISPR_gene_effect.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
#CHRONOS version
CRISPR_dep <- read.delim(file = "CRISPR_gene_effect_Chronos.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# cell-line info
sample_info_dep <- read.delim(file = "sample_info.csv", stringsAsFactors = FALSE, header = TRUE, sep=",", row.names = 1)

## load drug response data (primary set) ##
setwd("C:/Users/ningz/OneDrive/Desktop/My files/DrugResponse2021/PRISM")

# drug response data
drug_resp_primary <- read.delim(file = "primary-screen-replicate-collapsed-logfold-change.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# cell-line info
sample_info_primary <- read.delim(file = "primary-screen-cell-line-info.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# drug info
treatment_info_primary <- read.delim(file = "primary-screen-replicate-collapsed-treatment-info.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

length(intersect(row.names(CRISPR_dep),row.names(drug_resp_primary)))



#load functions
setwd("C:/Users/ningz/OneDrive/Desktop/My files/DrugResponse2021")
source("NameChanger.R")
source("DepCorrGenerator.R")
source("correlation_query.R")

update_gene <- cell_name(CRISPR_dep, sample_info_dep, name_col = 2)
update_gene <- gene_name(update_gene)


update_drug <- cell_name(drug_resp_primary, sample_info_dep, name_col = 2)
update_drug <- drug_name(update_drug, treatment_info_primary, name_col = 2)

correlations_table <- data.table(read.csv("gene_drug_dependencies_chronos.csv"))
correlations_table <- correlations[,-1]
correlations_table <- correlations(update_gene, update_drug)
write.csv(correlations_table, file = "gene_drug_dependencies_chronos.csv")

summary(correlations)

corlist <- target(correlations_table, corr = 0.3)
corlist1 <- target(correlations_table2, corr = 0.5)
corlist1 <- target(correlations_table, corr = 0.5, drugname = "idasanutlin")
corlist2 <- target(correlations_table, corr = 0.2, gene = "BRAF")


#make sure correlation function is correct,  then summarize data and query drug targets, then use shiny

write.csv(corlist, file = "correlations_0.3.csv")
write.csv(corlist1, file = "correlations_0.4.csv")

v1 <- as.numeric(as.matrix(correlations[,2:length(correlations)]))
hist(v1)
v2 <- as.numeric(as.matrix(correlations[,2:10]))
ggplot()+aes(v2)+geom_histogram()


summary(v1)

#############################################################################

source("correlation_plots.R")


dp <- drug_corplot(corlist1, "0.4")
dp

gp <- gene_corplot(corlist1, "0.4")
gp

dp1 <- drug_corplot(corlist, "0.3")
dp1

gp1 <- gene_corplot(corlist, "0.3")
gp1


target_correlations <- targetgene(correlations_table, treatment_info_primary)
target_correlations1 <- remove_dupes(target_correlations)

target_correlations3 <- remove_dupes(target_correlations2)

ggplot(target_correlations, aes(correlation)) + geom_histogram() + xlim(-0.7,0.7)

# excluded <- correlations_table
# columns <- as.numeric(target_correlations$column)
# rows <- as.numeric(target_correlations$row)
# for(i in 10001:nrow(target_correlations)){
#   column <- columns[i]
#   row <- rows[i]
#   excluded[row, column] <- NA
# }
# 
# s1 <- sample(excluded, nrow(target_correlations), replace = TRUE)

#delete target pairs from the main correlation table
excluded <- difference(correlations_table, target_correlations)
excluded1 <- difference(correlations_table2, target_correlations2)

#use the excluded table and target correlation as desired
sample_plot <- cor_sample(excluded1, target_correlations2)
sample_plot$plot_env$sample

#which(is.na(sample_plot$plot_env$sample))

p1 <- quantile(sample_plot$plot_env$sample, probs = c(.025, .975))
p4 <- quantile(sample_plot$plot_env$sample, probs = c(.05, .95))

#quantile(sample_plot$plot_env$sample, probs = c(.0025, .9975))

sample_plot1 <- cor_sample(excluded1, target_correlations2)
sample_plot1$plot_env$sample

p2 <- quantile(sample_plot1$plot_env$sample, probs = c(.025, .975))
p5 <- quantile(sample_plot1$plot_env$sample, probs = c(.05, .95))

sample_plot2 <- cor_sample(excluded1, target_correlations2)
sample_plot2$plot_env$sample

p3 <- quantile(sample_plot2$plot_env$sample, probs = c(.025, .975))
p6 <- quantile(sample_plot2$plot_env$sample, probs = c(.05, .95))

lower_cutoff <- mean(c(p1[1], p2[1], p3[1]))
upper_cutoff <- mean(c(p1[2], p2[2], p3[2]))

lower_cutoff2 <- mean(c(p4[1], p5[1], p6[1]))
upper_cutoff2 <- mean(c(p4[2], p5[2], p6[2]))

combined <- combine_plots2(target_correlations, sample_plot)
combined

combined2 <- combine_plots3(target_correlations, sample_plot)
combined2

summary(target_correlations$correlation)
summary(sample_plot$plot_env$sample)

dat <- data.frame(correlation = c(as.numeric(sample_plot$plot_env$sample), 
                                  as.numeric(sample_plot1$plot_env$sample),
                                  as.numeric(sample_plot2$plot_env$sample)),
                  Data = rep(c("S1", "S2", "S3"),each = length(sample_plot$plot_env$sample)))

ggplot(dat, aes(x=correlation, color=Data, fill = Data)) + 
  geom_density(alpha = 0.4, size=2) +
  theme_gray() + scale_fill_brewer(palette = "Pastel2") + 
  scale_alpha(range = c(0.4, 0.8)) +
  labs(title = "sample distributions", x = "Correlation", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))

############################################################

#generate target pairs above cutoff from the target_correlations table
valuable_targets <- read.csv("valuable_targets_5_percent.csv")
valuable_targets <- valuable_targets[, -1]
valuable_targets <- subset(target_correlations1, as.numeric(correlation) < lower_cutoff | as.numeric(correlation) > upper_cutoff)
valuable_targets <- valuable_targets[order(-correlation)]
write.csv(valuable_targets, file = "valuable_targets_5_percent.csv")

valuable_targets2 <- read.csv("valuable_targets_10_percent.csv")
valuable_targets2 <- valuable_targets2[, -1]
valuable_targets2 <- subset(target_correlations1, as.numeric(correlation) < lower_cutoff2 | as.numeric(correlation) > upper_cutoff2)
valuable_targets2 <- valuable_targets2[order(-correlation)]
write.csv(valuable_targets2, file = "valuable_targets_10_percent.csv")

valuable_targets3 <- read.csv("valuable_targets_5_percent_spearman.csv")
valuable_targets3 <- valuable_targets3[, -1]
valuable_targets3 <- subset(target_correlations3, as.numeric(correlation) < lower_cutoff | as.numeric(correlation) > upper_cutoff)
valuable_targets3 <- valuable_targets3[order(-correlation)]
write.csv(valuable_targets3, file = "valuable_targets_5_percent_spearman.csv")

valuable_targets4 <- read.csv("valuable_targets_10_percent_spearman.csv")
valuable_targets4 <- valuable_targets4[, -1]
valuable_targets4 <- subset(target_correlations3, as.numeric(correlation) < lower_cutoff2 | as.numeric(correlation) > upper_cutoff2)
valuable_targets4 <- valuable_targets4[order(-correlation)]
write.csv(valuable_targets4, file = "valuable_targets_10_percent_spearman.csv")

true_corlist <- target2(correlations_table, corr1 = -0.1368361, corr2 = 0.1448784)
true_corlist

#legacy network code
# adj_matrix <- adjacency(valuable_targets)
# net1 <- dgnet(adj_matrix, valuable_targets)
# net1
