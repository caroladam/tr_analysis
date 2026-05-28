required_packages <- c("ggplot2", "dplyr", "tidyr", "viridis", "stringr", "ggrastr")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(tidyr)
library(stringr)
library(ggrastr)

df<-read.table("human_chimp_avg_len.bed")

color_mapping <- c("cds" = "darkgoldenrod1", "promoter" = "dodgerblue", "5utr" = "firebrick", "3utr" = "coral1", "intron"= "mediumseagreen", "intergenic" = "darkviolet")

df <- df %>% arrange(factor(feature, levels = c("intergenic", "intron", "promoter", "cds", "3utr", "5utr")))

# Main plot
p <- ggplot(df, aes(x = homo_len, y = pantro_len, colour = feature)) + 
  geom_abline(linewidth=2, slope = 1, intercept = 0, color = "black") +
  #scale_x_continuous(breaks=c(10,50,100,500,1000,2000,5000,10000), trans="log1p", name = "Human mean TR allele length") +
  #scale_y_continuous(breaks=c(10,50,100,500,1000,2000,5000,10000), trans="log1p", name = "Chimpanzee mean TR allele length") +
  scale_x_continuous(name = "Human mean TR allele length") +
  scale_y_continuous(name = "Chimpanzee mean TR allele length") +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"))

# Plot categories with CI lines
plot_category_with_ci <- function(category, color) {
  cat_data <- df[df$feature == category, ]
  list(
    rasterize(geom_linerange(data = cat_data, aes(ymin = pantro_p5, ymax = pantro_p95), color = color, alpha = 0.5, show.legend = FALSE), dpi=500),
    rasterize(geom_linerange(data = cat_data, aes(xmin = homo_p5, xmax = homo_p95), color = color, alpha = 0.5, show.legend = FALSE), dpi=500),
    ggrastr::geom_point_rast(data = cat_data, size = 2, alpha = 0.4, color = color, shape = 16, raster.dpi = 500, show.legend = FALSE),
    geom_smooth(data = cat_data, method = "lm", color = color, linewidth = 1, alpha = 0.8, se = FALSE))}

# Add each category layer
p <- p + 
    plot_category_with_ci("intergenic", color_mapping["intergenic"]) +
    plot_category_with_ci("intron", color_mapping["intron"]) +
    plot_category_with_ci("promoter", color_mapping["promoter"]) +
    plot_category_with_ci("cds", color_mapping["cds"]) +
    plot_category_with_ci("3utr", color_mapping["3utr"]) +
    plot_category_with_ci("5utr", color_mapping["5utr"]) 

ggsave("human_chimp_avg_len.pdf", plot = p, device = cairo_pdf, width = 6, height = 6)

## Heatmap
p <- ggplot(df, aes(x = homo_len, y = pantro_len)) + 
    geom_hex(bins = 50, aes(fill = after_stat(count), alpha = after_stat(count)), show.legend=F) +
    scale_fill_viridis_c(option = "F", direction=-1, name = "Density (log)", trans="log10") +
    scale_alpha_continuous(range = c(0.8, 1), guide = "none") +
    geom_abline(slope = 1, intercept = 0, color = "red", alpha = 1) +
    scale_x_continuous(breaks=c(10,25,50,100,250,500,1000,2500,5000,7000), trans = "log1p", name = "Human TR length") +
    scale_y_continuous(breaks=c(10,25,50,100,250,500,1000,2500,5000,5000), trans = "log1p", name = "Chimpanzee TR length") +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"))

ggsave("human_chimp_avg_len_heatmap.pdf", plot = p, device = cairo_pdf, width = 6, height = 6)

# Plot within-species variance x mean allele len divergence between species
df <- df %>%
  mutate(
    slope = homo_len - pantro_len,
    direction = case_when(
      slope > 0  ~ "human_exp",
      slope < 0  ~ "chimp_exp",
      TRUE       ~ "zero"))

df<-df%>%mutate(homo_central_90 = homo_p95 - homo_p5)
df<-df%>%mutate(chimp_central_90 = pantro_p95 - pantro_p5)

df <- df %>%
  mutate(
    mean_length = (homo_len + pantro_len) / 2,
    relative_diff = abs(homo_len - pantro_len) / mean_length)

df<-df%>%mutate(relative_var = homo_var/homo_len)

# highlight trait-associated TRs
trait <- df[df$id %in% c("MUC1", "ACAN", "TCHH", "TENT5A", "MACF1", "RRBP1", "SPDYE3", "KRTAP5-1", "KRTAP5-5", "NACA", "PDZD7", "PHETA1", "PHGR1", "PPP1R15A", "RERE", "TNFRSF10C", "GIGYF2", "ZNF470", "IRF5", "CCDC40", "UBC", "CLEC4M", "MUC21", "MUC22", "ABCD3", "AFF2", "AR", "ARX1", "ARX2", "ATN1", "ATXN1", "ATXN10", "ATXN2", "ATXN3", "ATXN7", "ATXN8OS", "BEAN1","C9ORF72", "CACNA1A", "CBL", "CNBP", "COMP", "CSTB", "DAB1", "DIP2B", "DMD", "EIF4A3", "FGF14", "FMR1", "FOXL2", "FXN", "GIPC1", "GLS", "HOXA131", "HOXA132", "HOXA133", "HOXD13", "HTT", "JPH3", "LRP12", "MARCHF6", "NIPA1", "NOP56", "NOTCH2NLC", "NUTM2B", "PABPN1", "PHOX2B", "PPP2R2B", "PRDM12", "PRNP", "RAPGEF2", "RILPL1", "RUNX2", "SAMD12", "SOX3", "TBP", "TBX1", "TCF4", "THAP11", "TNRC6A", "XYLT1", "YEATS2", "ZFHX3", "ZIC2","ZIC3"), ]

# remove LPA 
p <- ggplot(df%>%filter(id!="LPA"), aes(x = slope, y = homo_central_90)) +
    geom_hex(bins = 60) +
    scale_fill_viridis(option = "G", direction = -1, name = "Density (log10)", trans="log10") +
    geom_smooth(method = "lm", color = "red", linewidth = 1, se = T, level = 0.99, show.legend=F) +
    geom_vline(xintercept = 0, linetype="dashed") + 
    theme_bw() +
    scale_x_continuous(name = "Mean allele length divergence") +
    scale_y_continuous(name = "90% IPR of human allele lengths") +
    coord_cartesian(ylim = c(0, NA)) +
    theme(
        text = element_text(size=5, family = "Helvetica", colour="black"),
        axis.text = element_text(size=5, family = "Helvetica", colour="black"),
        axis.title = element_text(size=5, family = "Helvetica", colour="black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.96, 0.2),
        legend.key.size = unit(0.4, "cm"),
        axis.line = element_line(linewidth = 1, colour = "black"))

p <- p +
  geom_point(
    data = trait, aes(x = slope, y = homo_central_90), color = "red", size = 3, shape = 8, stroke = 0.5) +
    ggrepel::geom_text_repel(data = trait, aes(x = slope, y = homo_central_90, label = id), segment.color = "black", min.segment.length = 0, segment.size= 0.3, color = "black", size = 2, fontface = "plain")

ggsave("human_var_x_divergence_heatmap.pdf", plot = p, device = cairo_pdf, width = 6, height = 6)

df %>%
  summarize(
    cor = cor(slope, homo_central_90, use = "complete.obs"),
    p = cor.test(slope, homo_central_90)$p.value) %>%
    mutate(p_adj = p.adjust(p, method = "fdr"))
