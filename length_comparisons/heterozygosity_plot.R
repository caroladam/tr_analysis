required_packages <- c("ggplot2", "dplyr", "tidyr", "ggsignif", "stringr")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(required_packages, install_if_missing))

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(ggsignif)

df <- read.table("homo_pantro_df_heteroz.txt", header= F)

df <- df %>%
  mutate(gc_content = str_count(toupper(V7), "[GC]") / str_length(V7) * 100)

df <- df %>%
  mutate(V9 = case_when(
    V6 %in% 1:6 ~ as.character(V6),
    V6 >= 7 & V6 <= 20 ~ "7-20",
    V6 >= 21 & V6 <= 30 ~ "21-30", 
    V6 >= 31 & V6 <= 40 ~ "31-40",
    V6 >= 41 & V6 <= 50 ~ "41-50",
    V6 >= 51 & V6 <= 100 ~ "51-100",
    V6 > 100 ~ ">100")) %>%
  mutate(V9 = factor(V9, levels = c(
    as.character(1:6), "7-20", "21-30", "31-40", "41-50", "51-100", "101-200", ">100")))

df <- df %>%
  mutate(V10 = interaction(V1, V5, sep = "_"))

df <- df %>%
  rename(
    spp = V1,
    tr_id = V2,
    Ho = V3,
    He = V4,
    patogenic = V5,
    motif_size = V6,
    motif_seq = V7,
    feature = V8,
    GC = gc_content,
    motif_cat = V9,
    spp_pat_status = V10) %>%
  mutate(
    feature = as.factor(feature),
    spp = as.factor(spp),
    motif_size = as.numeric(motif_size),
    GC = as.numeric(GC))

mycolors=c('homo'='#5C608A','chimp'='#8A865D')

p <- df %>%
  ggplot(aes(x = factor(motif_cat), y = He, fill = spp, colour = spp)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, show.legend = FALSE, position = position_dodge(width = 0.75)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1, show.legend = FALSE, position = position_dodge(width = 0.75)) +
  #facet_wrap(~ feature, scale="free_y") +
  scale_x_discrete(name = "Motif length (bp)") +
  scale_y_continuous(name = "Genetic diversity") +
  scale_fill_manual(values = mycolors) +
  scale_colour_manual(values = mycolors) +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 8, family="Helvetica"),
    strip.background = element_blank(),
    text = element_text(size = 7),
    axis.title.x = element_text(size = 8, family="Helvetica"),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_text(colour = "black", size = 8))

ggsave("homo_pantro_heteroz_motifs.pdf", plot = p, device = cairo_pdf, width = 6, height = 6)

# Plot expected heterozygosity relative to the pathogenic status of TR

# Pathogenic vs Non-pathogenic (all species)
test1 <- wilcox.test(He ~ patogenic, data = df)

# chimp vs human within pathogenic
test2 <- wilcox.test(He ~ spp, data = filter(df, patogenic == "patogenic"))

# chimp vs human within non-pathogenic
test3 <- wilcox.test(He ~ spp, data = filter(df, patogenic == "non_patogenic"))

# coordinates for plot
max_y <- max(df$He, na.rm = TRUE)

pvals <- data.frame(
  group1 = c("non_patogenic", "chimp_patogenic", "chimp_non_patogenic"),
  group2 = c("patogenic", "homo_patogenic", "homo_non_patogenic"),
  y.position = c(max_y * 1.05, max_y * 1.15, max_y * 1.25),
  p.adj = c(test1$p.value, test2$p.value, test3$p.value),
  p.label = c(
    scales::pvalue(test1$p.value),
    scales::pvalue(test2$p.value),
    scales::pvalue(test3$p.value)))

p <- ggplot(df, aes(x = factor(patogenic), y = He, fill = spp, colour= spp)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, position = position_dodge(width = 0.75), show.legend=F) +
  #geom_point_rast(data = outliers_df, aes(x = factor(patogenic), y = He, colour = spp, fill = spp, alpha = 0.01), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75), shape = 16, size=1, alpha = 0.01, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, position = position_dodge(width = 0.75), show.legend=F) +
  # pathogenic vs non-pathogenic (across x-axis groups)
  geom_signif(
    comparisons = list(c("non_patogenic", "patogenic")),
    map_signif_level = TRUE, y_position = max_y * 1.05, size = 0.4, tip_length = 0.03, textsize = 4, family="Helvetica", color = "black") +
  # chimp vs human within pathogenic
  geom_signif(
    annotations = ifelse(test2$p.value < 0.001, "***", ifelse(test2$p.value < 0.01, "**", ifelse(test2$p.value < 0.05, "*", "ns"))),
    y_position = max_y * 1.15, xmin = 2 - 0.3, xmax = 2 + 0.3, size = 0.4, tip_length = 0.03, textsize = 4, family="Helvetica", color = "black") +
  # chimp vs human within non-pathogenic
  geom_signif(
    annotations = ifelse(test3$p.value < 0.001, "***", ifelse(test3$p.value < 0.01, "**", ifelse(test3$p.value < 0.05, "*", "ns"))),
    y_position = max_y * 1.25, xmin = 1 - 0.3, xmax = 1 + 0.3, size = 0.4, tip_length = 0.03, textsize = 4, family="Helvetica", color = "black") +
  theme(legend.position = "none") +
  guides(fill = "none", colour = "none") +
  scale_x_discrete(labels = c("non_patogenic" = "Non-pathogenic", "patogenic" = "Pathogenic")) +
  scale_y_continuous(name = "Genetic diversity") +
  scale_fill_manual(values = mycolors) +
  scale_colour_manual(values = mycolors) +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(family="Helvetica"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 8, family="Helvetica"),
    axis.text.x = element_text(colour = "black", size = 8, family="Helvetica"),
    axis.text.y = element_text(colour = "black", size = 8, family="Helvetica"))

ggsave("het_patogenic_trs.pdf", plot = p, device = cairo_pdf, width = 3, height = 6)
