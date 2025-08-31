library(dplyr)
library(ggplot2)
library(ggsignif)

df <- read.table("human_chimp_mean_het_per_locus.txt")

df <- df %>%
  rename(
    species        = V1,
    TR_id          = V2,
    feature        = V3,
    heterozygosity = V4,
    pathogenicity  = V5,
    motif_length   = V6,
    motif_seq      = V7)

# Define groups
df <- df %>%
  mutate(group = paste(species, pathogenicity, sep = "_"))

# Wilcoxon tests
test1 <- wilcox.test(heterozygosity ~ pathogenicity, data = df)
test2 <- wilcox.test(heterozygosity ~ species, data = filter(df, pathogenicity == "pathogenic"))
test3 <- wilcox.test(heterozygosity ~ species, data = filter(df, pathogenicity == "non_pathogenic"))

# Annotation values for plot
max_y <- max(df$heterozygosity, na.rm = TRUE)

pvals <- data.frame(
  group1 = c("non_pathogenic", "chimp_pathogenic", "chimp_non_pathogenic"),
  group2 = c("pathogenic", "homo_pathogenic", "homo_non_pathogenic"),
  y.position = c(max_y * 1.05, max_y * 1.15, max_y * 1.25),
  p.adj = c(test1$p.value, test2$p.value, test3$p.value),
  p.label = c(
    scales::pvalue(test1$p.value),
    scales::pvalue(test2$p.value),
    scales::pvalue(test3$p.value)))

ggplot(df, aes(x = factor(pathogenicity), y = heterozygosity, fill = species, colour = species)) + 
  geom_boxplot(alpha = 0.8, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               position = position_dodge(width = 0.75)) +
  
  # pathogenic vs non-pathogenic
  geom_signif(
    comparisons = list(c("non_pathogenic", "pathogenic")),
    map_signif_level = TRUE,
    y_position = max_y * 1.05,
    tip_length = 0.03,
    textsize = 6) +
  # chimp vs human within pathogenic
  geom_signif(
    annotations = ifelse(test2$p.value < 0.001, "***", 
                  ifelse(test2$p.value < 0.01, "**", 
                  ifelse(test2$p.value < 0.05, "*", "ns"))),
    y_position = max_y * 1.15,
    xmin = 2 - 0.3, xmax = 2 + 0.3,
    tip_length = 0.03,
    textsize = 6) +
  # chimp vs human within non-pathogenic
  geom_signif(
    annotations = ifelse(test3$p.value < 0.001, "***", 
                  ifelse(test3$p.value < 0.01, "**", 
                  ifelse(test3$p.value < 0.05, "*", "ns"))),
    y_position = max_y * 1.25,
    xmin = 1 - 0.3, xmax = 1 + 0.3,
    tip_length = 0.03,
    textsize = 6) +
  scale_x_discrete(name = "Pathogenic Status") +
  scale_y_continuous(name = "Expected Heterozygosity") +
  scale_fill_manual(values = mycolors) +
  scale_colour_manual(values = mycolors) +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 21),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20))
