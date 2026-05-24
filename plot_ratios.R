required_packages <- c("ggplot2", "hexbin", "dplyr", "ggrepel", "patchwork", "viridis")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

library(ggplot2)
library(hexbin)
library(dplyr)
library(ggrepel)
library(patchwork)
library(viridis)

# input
input_main_file <- "all_categories_intersect.txt"

# outliers
highlight_trids <- c("TR1334497", "TR140813")
highlight_points <- df[df$TRID %in% highlight_trids, ]

# Main hex plot
main <- ggplot(df, aes(x = ratio_homo, y = ratio_pantro)) + 
  geom_hex(bins = 50, aes(fill = after_stat(count), alpha = after_stat(count)), show.legend = F) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Density (log)", trans = "log10") +
  scale_alpha_continuous(range = c(0.8, 1), guide = "none") +
  scale_x_continuous(name = "Human ratio", trans = "log1p", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 1e7, 2e8)) +
  scale_y_continuous(name = "Chimpanzee ratio", trans = "log1p", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 1e7, 2e8, 2e9)) +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(colour = "black", size = 16, angle = 45, vjust = 0.5),
    axis.text.y = element_text(colour = "black", size = 16)) 

# plot corner boxes with limits of the top categories
bounds_high_both <- highlight_data %>%
  filter(category == "high_both") %>%
  summarise(x_min = min(ratio_homo), x_max = max(ratio_homo), y_min = min(ratio_pantro), y_max = max(ratio_pantro))
bounds_low_both <- highlight_data %>%
  filter(category == "low_both") %>%
  summarise(x_min = min(ratio_homo), x_max = max(ratio_homo), y_min = min(ratio_pantro), y_max = max(ratio_pantro))

# add minimum variance
bounds_low_both <- bounds_low_both %>% mutate(x_min = pmax(x_min, 1), y_min = pmax(y_min, 1))

main <- main +
  # high–high
  geom_segment(
    data = bounds_high_both, aes(x = x_min, xend = x_max, y = y_min, yend = y_min), linetype = "dashed") +
  geom_segment(
    data = bounds_high_both, aes(x = x_min, xend = x_min, y = y_min, yend = y_max), linetype = "dashed") +
  # low–low
  geom_segment(
    data = bounds_low_both, aes(x = x_min, xend = x_max, y = y_max, yend = y_max), linetype = "dashed") +
  geom_segment(
    data = bounds_low_both, aes(x = x_max, xend = x_max, y = y_min, yend = y_max), linetype = "dashed")

# Marginal histogram of human ratio (top)
xhist <- ggplot(df, aes(x = ratio_homo)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "gray70", color = "black") +
  scale_x_continuous(trans = "log1p") +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(colour = "black", size = 10))

# HMarginal histogram of chimp ratio (right)
yhist <- ggplot(df, aes(x = ratio_pantro)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "gray70", color = "black") +
  scale_x_continuous(trans = "log1p") +
  coord_flip() +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(colour = "black", size = 10))

# Add points with outliers
main <- main +
  geom_point(
    data = highlight_points,
    aes(x = ratio_homo, y = ratio_pantro),
    color = "red", size = 1, shape = 8, stroke = 1) +
  ggrepel::geom_text_repel(
    data = highlight_points,
    aes(x = ratio_homo, y = ratio_pantro, label = V14),
    color = "black", size = 4, fontface = "bold")

final_plot <- (xhist + plot_spacer() + main + yhist) + 
  plot_layout(widths = c(4, 1), heights = c(1, 4))
