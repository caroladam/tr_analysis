# ---
# Title: get_tr_ratios.R
# Date: 2025
# Author: Adam, Carolina L.
# Purpose: Calculate HKA-like divergence-diversity TR ratio
# Output:
  # A tab-separated file with columns:
  # tr_id, mean_homo_len, mean_pantro_len, var_homo, var_pantro, var_between, ratio_homo, ratio_pantro
# ---

required_packages <- c("dplyr", "data.table", "tidyr", "doParallel", "foreach")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
lapply(required_packages, install_if_missing)

library(dplyr)
library(data.table)
library(tidyr)
library(doParallel)
library(foreach)

# Load input file
df <- fread("human_pantro_long_format.txt", header = FALSE,
                 col.names = c("species", "tr", "motif", "allele_length", "copy_c", "constancy_score"), stringAsFactors = F)

# Remove TRs with mean length < 11bp
total_trs_before_mean <- n_distinct(df$tr)

valid_mean <- df %>%
  group_by(tr, species) %>%
  summarise(Mean = mean(allele_length, na.rm = TRUE), .groups = "drop") %>%
  filter(Mean >= 11)

df <- df %>%
  semi_join(valid_mean, by = c("tr", "species"))

trs_after_mean_filter <- n_distinct(df$tr)
cat("TRs removed by mean length filter:", total_trs_before_mean - trs_after_mean_filter, "\n")
cat("TRs remaining:", trs_after_mean_filter, "\n\n")

# Remove invariant loci
total_trs_before_var <- n_distinct(df$tr)

valid_trs <- df %>%
  group_by(tr) %>%
  summarise(Variance = var(allele_length, na.rm = TRUE), .groups = "drop") %>%
  filter(Variance > 0)

df <- df %>%
  filter(tr %in% valid_trs$tr)

trs_after_var_filter <- n_distinct(df$tr)
cat("TRs removed by zero variance filter:", total_trs_before_var - trs_after_var_filter, "\n")
cat("TRs remaining:", trs_after_var_filter, "\n\n")

# Remove outlier loci function
remove_top_n <- function(x, n = 3) {
  if (length(x) <= n) return(x)
  med <- median(x, na.rm = TRUE)
  # Order by absolute distance from the median
  idx <- order(abs(x - med), decreasing = TRUE)  
  # Drop the n most divergent
  x[-idx[1:n]]}

# Variance functions
safe_var <- function(x) {
  v <- var(x, na.rm = TRUE)
  ifelse(is.na(v), NA, v)}

safe_ratio <- function(var_within, var_between, epsilon) {
  denom <- ifelse(is.na(var_within) || var_within < epsilon, epsilon, var_within)
  return(var_between / denom)}

# Start loop
start_time <- Sys.time()

output_file <- "homo_pantro_trv_ratios.csv"

# Create output file with header if not exists
if (!file.exists(output_file)) {
  write.table(data.frame(
    tr = character(),
    mean_homo = numeric(),
    mean_pantro = numeric(),
    var_homo = numeric(),
    var_pantro = numeric(),
    var_between = numeric(),
    ratio_homo = numeric(),
    ratio_pantro = numeric()), output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)}

# Set up parallel backend
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

tr_groups <- split(df, sort(df$tr))
lock <- tempfile()

foreach(group = tr_groups, .combine = rbind, .packages = c("dplyr")) %dopar% {
  tr_id <- unique(group$tr)
  
  epsilon <- 1

  # Remove top 3 allele outliers per species
  homo_alleles <- remove_top_n(group$allele_length[group$species == "homo"], n = 3)
  pantro_alleles <- remove_top_n(group$allele_length[group$species == "pantro"], n = 3)

  # Skip TRs with zero variance after outlier removal
  all_alleles <- c(homo_alleles, pantro_alleles)
  if(length(all_alleles) <= 1 || var(all_alleles) == 0) return(NULL)

  var_homo <- safe_var(homo_alleles)
  var_pantro <- safe_var(pantro_alleles)

  mean_homo <- mean(homo_alleles)
  mean_pantro <- mean(pantro_alleles)

  n_h <- length(homo_alleles)
  n_p <- length(pantro_alleles)
  mean_overall <- (mean_homo * n_h + mean_pantro * n_p) / (n_h + n_p)

  var_between <- (n_h * (mean_homo - mean_overall)^2 +
                  n_p * (mean_pantro - mean_overall)^2)

  # Compute ratios
  ratio_homo <- safe_ratio(var_homo, var_between, epsilon)
  ratio_pantro <- safe_ratio(var_pantro, var_between, epsilon)

  res <- data.frame(
    tr = tr_id,
    mean_homo = mean_homo,
    mean_pantro = mean_pantro,
    var_homo = var_homo,
    var_pantro = var_pantro,
    var_between = var_between,
    ratio_homo = ratio_homo,
    ratio_pantro = ratio_pantro)

  while (file.exists(lock)) Sys.sleep(0.01)
  file.create(lock)
  write.table(res, file = output_file, sep = "\t", row.names = FALSE,
              quote = FALSE, col.names = FALSE, append = TRUE)
  file.remove(lock)
}

stopCluster(cl)
end_time <- Sys.time()
cat("Elapsed time:", end_time - start_time, "\n")
