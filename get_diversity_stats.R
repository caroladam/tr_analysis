# Load required libraries
library(poppr)
library(vcfR)
library(adegenet)
library(hierfstat)
library(tidyr)
library(dplyr)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if required arguments are provided
if (length(args) < 3) {
  stop("Usage: Rscript script_name.R <vcf_file> <popinfo_file> <output_file>")
}

vcf_file <- args[1]
popinfo_file <- args[2]
output_file <- args[3]

# Read VCF file
vcf <- read.vcfR(vcf_file)

# Read population information file
popinfo <- read.csv(popinfo_file, header = TRUE, sep = "\t")

# Convert VCF to genind object and assign population information
df <- vcfR2genind(vcf, ind.names = popinfo$Individual, pop = popinfo$Individual)

# Calculate basic statistics
basic_stats <- basic.stats(df)

ho_pop <- basic_stats$Ho
ho_pop_df <- as.data.frame(ho_pop)
ho_pop_df$locus <- rownames(ho_pop_df)

# Gather the data into long format
ho_pop_long <- pivot_longer(ho_pop_df, cols = -locus, names_to = "Individual", values_to = "Heterozygosity")

# Merge with population info
ho_pop_long <- merge(ho_pop_long, popinfo, by = "Individual", all.x=T)

# Write heterozygosity values to file
write.table(ho_pop_long, output_file, sep="\t", row.names=F, col.names=F)

# Get mean values per individual maintaining population label
mean_het <- ho_pop_long %>%
  group_by(Individual) %>%  # Group by Individual and Population
  summarize(mean_het = mean(Heterozygosity, na.rm = TRUE), .groups = 'drop')  # Calculate mean heterozygosity

# Write mean heterozygosity per individual to a file
mean_het_file <- paste0("mean_het_", basename(vcf_file), ".txt")
write.table(mean_het, mean_het_file, sep="\t", row.names=F, col.names=F)

cat("Mean heterozygosity for individuals saved to:", mean_het_file, "\n")
