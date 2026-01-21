required_packages <- c("dplyr", "tidyr", "doParallel", "foreach")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

# Load the input file
df <- read.table("human_chimp_exp_sumstats.txt", header = FALSE,
                 col.names = c("tr", "ind", "motif_len", "allele_length", "motif", "copy_n", "species"),
                 stringsAsFactors = FALSE)

# Convert copy_n to numeric, remove NA or invalid rows
df$copy_n <- as.numeric(df$copy_n)
df <- df[!is.na(df$copy_n), ]
df <- df[!is.na(df$species), ]

# Filter to TRs with both species represented
species_count <- df %>%
  group_by(tr) %>%
  summarise(Species_Count = n_distinct(species), .groups = "drop")

df <- df %>%
  filter(tr %in% species_count$tr[species_count$Species_Count >= 2])

# Remove invariant loci across all samples
valid_trs <- df %>%
  group_by(tr) %>%
  summarise(Variance = var(copy_n, na.rm = TRUE), .groups = "drop") %>%
  filter(Variance > 0)

df <- df %>%
  filter(tr %in% valid_trs$tr)

# Set pseudocount to avoid division by zero
epsilon <- 0.01

# Safe variance function
safe_var <- function(x) {
  v <- var(x, na.rm = TRUE)
  if (is.na(v)) return(NA)
  if (v == 0) return(0)
  return(v)
}

# Safe ratio function
safe_ratio <- function(var_within, var_between) {
  if (is.na(var_within) || is.na(var_between)) return(NA)
  if (var_within == 0) {
    return(var_between / epsilon)
  } else {
    return(var_between / var_within)
  }
}

start_time <- Sys.time()

output_file <- "tr_variance_ratios.csv"

# Create output file with header if not exists
if (!file.exists(output_file)) {
  write.table(data.frame(
    tr = character(),
    mean_homo = numeric(),
    mean_chimp = numeric(),
    var_homo = numeric(),
    var_chimp = numeric(),
    var_between = numeric(),
    ratio_homo = numeric(),
    ratio_chimp = numeric()
  ), output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}

# Set up parallel backend
n_cores <- min(6, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

tr_groups <- split(df, df$tr)

# Avoid corrupting the file
lock <- tempfile()

foreach(group = tr_groups, .combine = rbind, .packages = c("dplyr")) %dopar% {
  if (nrow(group) == 0) return(NULL)

  tr_id <- unique(group$tr)

  # Require both species to be present
  if (sum(group$species == "homo") == 0 || sum(group$species == "chimp") == 0) return(NULL)

  mean_homo <- mean(group$copy_n[group$species == "homo"], na.rm = TRUE)
  mean_chimp <- mean(group$copy_n[group$species == "chimp"], na.rm = TRUE)

  var_homo <- safe_var(group$copy_n[group$species == "homo"])
  var_chimp <- safe_var(group$copy_n[group$species == "chimp"])

  n_h <- sum(group$species == "homo")
  n_p <- sum(group$species == "chimp")

  if ((n_h + n_p) == 0) return(NULL)

  mean_overall <- (mean_homo * n_h + mean_chimp * n_p) / (n_h + n_p)

  var_between <- (n_h * (mean_homo - mean_overall)^2 +
                  n_p * (mean_chimp - mean_overall)^2)

  ratio_homo <- safe_ratio(var_homo, var_between)
  ratio_chimp <- safe_ratio(var_chimp, var_between)

  res <- data.frame(
    tr = tr_id,
    mean_homo = mean_homo,
    mean_chimp = mean_chimp,
    var_homo = var_homo,
    var_chimp = var_chimp,
    var_between = var_between,
    ratio_homo = ratio_homo,
    ratio_chimp = ratio_chimp
  )

  # Write to file with lock
  while (file.exists(lock)) Sys.sleep(0.01)
  file.create(lock)
  write.table(res, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
  file.remove(lock)

  return(NULL)  # No need to return anything for combine
}

stopCluster(cl)

end_time <- Sys.time()
cat("Elapsed time:", end_time - start_time, "\n")
