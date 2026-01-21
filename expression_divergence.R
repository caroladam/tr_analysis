# ---
# title: "Gene expression divergence"
# date: "2025-12-02"
# author: Adam, Carolina L.
# ---

# Running limma with recommendations from Law et al. (2018)

library(BiocManager)
BiocManager::install("limma")
BiocManager::install("biomaRt")

library(biomaRt)
library(limma)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrastr)
library(rstatix)
library(ggpubr)

expr <- read.table("", header=T, sep="\t")

# keep only human and chimp data
expr <- expr %>% select(starts_with("hsa"), starts_with("ptr"))

# Get ID colums
id_cols <- c("hsa", "ptr")
expr_ids <- expr %>% select(all_of(id_cols))

# get sample values
expr_samples <- expr %>% select(-all_of(id_cols))

# Tidy table to get one expression value (RPKM) per row 
expr_long <- expr_samples %>%
  mutate(hsa_id = expr_ids$hsa,
         ptr_id = expr_ids$ptr) %>%
  pivot_longer(
    cols = -c(hsa_id, ptr_id),
    names_to = "sample",
    values_to = "expr") %>%
  separate(sample, into = c("species", "tissue", "sex", "indiv"), sep = "_") %>%
  mutate(sample_id = paste(species, tissue, sex, indiv, sep="_")) %>%
  mutate(ind_id = paste(species, sex, indiv, sep="_"))

# write long format table
write.table(expr_long, "rpkm_human_and_chimp_long_table.txt", col.names=T, row.names=F, sep="\t", quote=F)

# See if expression values are normalized per individual per tissue
ggplot(expr_long, aes(x = sample_id, y = expr + 1, fill = species)) +
  geom_boxplot(outlier.shape = 1, width=0.1, position = position_dodge(width = 0.8), show.legend = FALSE) + 
  geom_violin(trim=F, show.legend=F, alpha=0.5) +
  facet_wrap(~ tissue, scales="free") +
  scale_fill_manual(values = c("hsa" = "#5C608A", "ptr" = "#8A865D")) +
  scale_x_discrete(name = "sample") +
  scale_y_continuous(trans="log1p", name = "Expression (log1p))") +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(colour = "black", size = 6, angle=90, vjust=0.5),
    axis.text.y = element_text(colour = "black", size = 6))

# Filter: Only test differential expression for genes that are detectably expressed in all tissues in both species.
expr_long_filtered <- expr_long %>%
  group_by(species, tissue, hsa_id) %>%
  mutate(expressed = any(expr >= 1)) %>% # in at least one sample across tissues
  group_by(hsa_id) %>%
  filter(all(expressed)) %>% # appears in both species
  ungroup() %>%
  dplyr::select(-expressed)

# Get average gene expression per spp per tissue
avg_expr <- expr_long_filtered %>%
  group_by(hsa_id, tissue, species) %>%
  summarize(mean_expr = mean(expr), .groups = "drop") %>%
  pivot_wider(
    names_from = species,
    values_from = mean_expr,
    names_prefix = "avg_")

###################################################################
################# Running limma per tissue ########################
###################################################################

# start here with long format table
# tissues: br, cb, ht, kd, lv, ts

tissue_code <- "br"
expr_tissue <- expr_long_filtered %>% filter(tissue == tissue_code)

# Put the table back into wide format to have one gene per row
expr_wide <- expr_tissue %>%
  dplyr::select(hsa_id, sample_id, expr) %>%
  pivot_wider(names_from = sample_id, values_from = expr)

# Expression matrix to prep for running limma
mat <- expr_wide %>% 
  dplyr::select(-hsa_id) %>%
  as.matrix()

# I am adding hsa gene names as row names because it is more comparable with our data
rownames(mat) <- expr_wide$hsa_id

# log2 transform plus offset - recommendation by authors
mat_log <- log2(mat + 1)

# create design matrix
sample_info <- tibble(
  sample_id = colnames(mat_log),
  species   = ifelse(grepl("^hsa", colnames(mat_log)), "human", "chimp"))

design <- model.matrix(~ 0 + species, data = sample_info)
colnames(design) <- c("chimp", "human")

# create contrast
# this tells limma to take the average expression in one species and subtract the average from the other
contr_matrix <- makeContrasts(human_vs_chimp = human - chimp, levels = design)

# Run limma
fit <- lmFit(mat_log, design) # fir a linear model per gene per species
fit2 <- contrasts.fit(fit, contrasts = contr_matrix) # apply the contrast to the linear model
fit2 <- eBayes(fit2) # performs empirical Bayes moderation of the standard errors - uses info from all genes to stabilize variance per gene

plotSA(fit2) # you should not see mean-variance dependence - variance should not be dependent on mean expression values

# access significance using adjusted p-value
de <- decideTests(fit2)
summary(de) # just number of genes

# Convert DE results (decideTests) to tibbles
de_df <- as.data.frame(de) %>%
  tibble::rownames_to_column("hsa_id") %>%
  dplyr::rename(de = human_vs_chimp)

# for more strict definition
# TREAT tests whether the log2 fold change is significantly greater than the threshold (lfc) 
fit2_treat <- treat(fit2, lfc=1) # log-FC > 1 (equal to a 2-fold difference)

de_treat <- decideTests(fit2_treat)
summary(de_treat)

# Convert DE treat results to tibble
de_treat_df <- as.data.frame(de_treat) %>%
  tibble::rownames_to_column("hsa_id") %>%
  dplyr::rename(de_treat = human_vs_chimp)

# get you limma results 
tt <- topTable(fit2, coef = 1, number = Inf) %>%
  tibble::rownames_to_column("hsa_id") %>%
  dplyr::mutate(tissue = tissue_code) 

# get your limma results with "treat"
tt_treat <- topTable(fit2_treat, coef = 1, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("hsa_id") %>%
  dplyr::rename(
    logFC_treat = logFC,
    AveExpr_treat = AveExpr,
    t_treat = t,
    P.Value_treat = P.Value,
    adj.P.Val_treat = adj.P.Val)

#Output
# logFC - log2 fold change (the way contrasts are set up, positive values are > humans and negative values are > chimps)
# AveExpr - average log2 expression across all samples in the tissue
# t - moderated t-statistic for the contrast (hsa x ptr)
# B - bayes results. Higher values = more likely DE
    # B = logP(gene DE)/ 1-P(gene DE)

# Merge into a final table with all values
tt_final <- tt %>%
  left_join(tt_treat, by = "hsa_id") %>%
  left_join(de_df, by = "hsa_id") %>%
  left_join(de_treat_df, by = "hsa_id") %>%
  dplyr::select(-logFC_treat, -AveExpr_treat)

write.table(tt_final, "limma_results_brain.txt", col.names=T, row.names=F, sep="\t", quote=F)
# You can loop through each tissue to get the full table results
