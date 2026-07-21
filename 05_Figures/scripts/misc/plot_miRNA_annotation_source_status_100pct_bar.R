# ============================================================
# 100% stacked bar plot for miRNA annotation category by source
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", script_args[grep(file_arg, script_args)])
if (length(script_path) == 0) {
  script_path <- "SoymiR-atlas/05_Figures/scripts/misc/plot_miRNA_annotation_source_status_100pct_bar.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
module_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

input_file <- file.path(
  module_dir,
  "input/plotting_data/misc/2814_annotation_hierarchical_counts_source_annotation_category.tsv"
)
out_dir <- file.path(
  module_dir,
  "results/intermediate_figures/misc"
)
final_dir <- file.path(module_dir, "results/final_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(input_file, show_col_types = FALSE)

category_order <- c(
  "reference_matched",
  "reference_locus_variant",
  "known_family_new_locus",
  "known_family_new_locus_variant",
  "novel_family_new_locus",
  "novel_family_new_locus_variant",
  "unannotated_opposite_arm_product",
  "unannotated_opposite_arm_variant"
)

legend_order <- c(
  "reference_matched",
  "reference_locus_variant",
  "known_family_new_locus",
  "known_family_new_locus_variant",
  "novel_family_new_locus",
  "novel_family_new_locus_variant",
  "unannotated_opposite_arm_product",
  "unannotated_opposite_arm_variant"
)

category_labels <- c(
  "reference_matched" = "Reference matched",
  "reference_locus_variant" = "Reference locus variant",
  "unannotated_opposite_arm_product" = "Unannotated opposite-arm",
  "unannotated_opposite_arm_variant" = "Unannotated opposite-arm variant",
  "known_family_new_locus" = "Known-family new locus",
  "known_family_new_locus_variant" = "Known-family new locus variant",
  "novel_family_new_locus" = "Novel-family new locus",
  "novel_family_new_locus_variant" = "Novel-family new locus variant"
)

source_labels <- c(
  "miRbase" = "miRBase",
  "pmiren" = "PmiREN",
  "soymir" = "SoymiR"
)

category_colors <- c(
  "reference_matched" = "#4D4D4D",
  "reference_locus_variant" = "#0072B2",
  "known_family_new_locus" = "#56B4E9",
  "known_family_new_locus_variant" = "#009E73",
  "novel_family_new_locus" = "#E69F00",
  "novel_family_new_locus_variant" = "#F0E442",
  "unannotated_opposite_arm_product" = "#D55E00",
  "unannotated_opposite_arm_variant" = "#CC79A7"
)

df_plot <- df %>%
  mutate(
    Source = factor(Source, levels = c("soymir", "pmiren", "miRbase")),
    annotation_category = factor(annotation_category, levels = category_order)
  ) %>%
  arrange(Source, annotation_category) %>%
  group_by(Source) %>%
  mutate(
    Total = sum(Count),
    Proportion = Count / Total,
    label_pos = cumsum(Proportion) - Proportion / 2,
    label_pos = ifelse(Source == "soymir" & annotation_category == "unannotated_opposite_arm_product", 0.005, label_pos),
    label = ifelse(
      Count < 6 |
        (Source == "soymir" & annotation_category == "unannotated_opposite_arm_product") |
        (Source == "miRbase" & annotation_category %in% c("novel_family_new_locus", "novel_family_new_locus_variant")),
      "",
      as.character(Count)
    ),
    label_color = ifelse(annotation_category == "reference_matched", "white", "black")
  ) %>%
  ungroup()

df_total <- df_plot %>%
  distinct(Source, Total)

p <- ggplot(df_plot, aes(x = Source, y = Proportion, fill = annotation_category)) +
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    aes(y = label_pos, label = label, color = label_color),
    size = 4.0,
    show.legend = FALSE
  ) +
  geom_text(
    data = df_total,
    aes(x = Source, y = 1.035, label = paste0("n = ", Total)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.2
  ) +
  coord_flip(ylim = c(0, 1), clip = "off") +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_color_identity() +
  scale_x_discrete(labels = source_labels) +
  scale_fill_manual(
    values = category_colors,
    breaks = legend_order,
    labels = category_labels[legend_order],
    name = NULL
  ) +
  labs(x = NULL, y = "Proportion within source") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, color = "black", margin = margin(r = 10)),
    axis.text = element_text(size = 13, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.line.y = element_blank(),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.55, "cm"),
    legend.spacing.x = unit(1.0, "cm"),
    legend.margin = margin(t = 5),
    plot.margin = margin(15, 80, 15, 35)
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = FALSE,
      override.aes = list(color = NA)
    )
  )

ggsave(
  filename = file.path(out_dir, "miRNA_annotation_source_annotation_category_100pct_bar.pdf"),
  plot = p,
  width = 9,
  height = 5.8,
  units = "in",
  device = pdf
)

ggsave(
  filename = file.path(out_dir, "miRNA_annotation_source_annotation_category_100pct_bar.png"),
  plot = p,
  width = 9,
  height = 5.8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(final_dir, "Figure_2.png"),
  plot = p,
  width = 9,
  height = 5.8,
  units = "in",
  dpi = 600
)
