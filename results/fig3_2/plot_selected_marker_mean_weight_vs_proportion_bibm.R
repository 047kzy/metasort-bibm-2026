suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

input_path <- "results/fig3_2/version_2/metasort_weight_rerun_sqrt_tol0p005_avg2p5/target_celltype_weight_trends/target_marker_weight_by_mixture.csv"
plot_dir <- "results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots"

font_family <- "Times New Roman"

selected_targets <- tibble::tribble(
  ~dataset, ~tissue, ~cell_type, ~plot_title, ~panel_color, ~file_stem,
  "CELlxGENE", "Basal_Zone_Of_Heart", "vein_endothelial_cell",
    "Heart: vein endothelial cell", "#E45756",
    "heart_vein_endothelial_marker_mean_weight_vs_proportion_nonzero_bibm",
  "Tabula-Sapiens", "Skin", "CD1c_positive_myeloid_dendritic_cell",
    "Skin: CD1c+ myeloid DC", "#E45756",
    "skin_cd1c_myeloid_marker_mean_weight_vs_proportion_nonzero_bibm"
)

plot_data <- read_csv(input_path, show_col_types = FALSE) %>%
  inner_join(
    selected_targets,
    by = c("dataset", "tissue", "cell_type")
  ) %>%
  filter(truth_proportion > 0)

if (nrow(plot_data) == 0) {
  stop("No non-zero selected tissue/cell-type rows were found in target_marker_weight_by_mixture.csv")
}

summary_data <- plot_data %>%
  group_by(dataset, tissue, cell_type, plot_title, file_stem) %>%
  summarise(
    n = n(),
    mean_truth_proportion = mean(truth_proportion, na.rm = TRUE),
    mean_marker_weight = mean(target_marker_mean_weight, na.rm = TRUE),
    pearson_r = cor(truth_proportion, target_marker_mean_weight, method = "pearson"),
    .groups = "drop"
  ) %>%
  mutate(
    stat_label = sprintf(
      "Pearson r = %.2f",
      pearson_r
    )
  )

label_data <- summary_data %>%
  mutate(
    x_pos = min(plot_data$truth_proportion, na.rm = TRUE) +
      0.05 * diff(range(plot_data$truth_proportion, na.rm = TRUE)),
    y_pos = max(plot_data$target_marker_mean_weight, na.rm = TRUE) -
      0.05 * diff(range(plot_data$target_marker_mean_weight, na.rm = TRUE))
  )

write_csv(
  plot_data,
  file.path(plot_dir, "selected_marker_mean_weight_vs_proportion_nonzero_data.csv")
)
write_csv(
  summary_data,
  file.path(plot_dir, "selected_marker_mean_weight_vs_proportion_nonzero_summary.csv")
)

make_plot <- function(target) {
  one_data <- plot_data %>%
    filter(file_stem == target$file_stem)
  one_summary <- summary_data %>%
    filter(file_stem == target$file_stem)

  p <- ggplot(
    one_data,
    aes(x = truth_proportion, y = target_marker_mean_weight)
  ) +
    geom_point(size = 2.0, alpha = 0.76, stroke = 0, color = target$panel_color) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.75,
      formula = y ~ x,
      color = target$panel_color
    ) +
    scale_x_continuous(
      labels = number_format(accuracy = 0.01),
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.08, 0.10))
    ) +
    labs(
      title = target$plot_title,
      subtitle = one_summary$stat_label,
      x = "True cell-type proportion",
      y = "Target marker mean weight"
    ) +
    theme_classic(base_family = font_family, base_size = 13.5) +
    theme(
      plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5, margin = margin(b = 7), lineheight = 0.95),
      plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 8), color = "black"),
      axis.title = element_text(size = 13.5),
      axis.text = element_text(size = 12, color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      panel.grid.major = element_line(linewidth = 0.24, color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 12, 8, 12)
    )

  pdf_path <- file.path(plot_dir, paste0(target$file_stem, ".pdf"))
  png_path <- file.path(plot_dir, paste0(target$file_stem, ".png"))

  ggsave(
    pdf_path,
    p,
    device = cairo_pdf,
    width = 3.55,
    height = 2.75,
    units = "in",
    family = font_family
  )

  ggsave(
    png_path,
    p,
    width = 3.55,
    height = 2.75,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  cat("Wrote:\n", pdf_path, "\n", png_path, "\n", sep = "")
}

for (row_index in seq_len(nrow(selected_targets))) {
  make_plot(selected_targets[row_index, ])
}

cat("\nSummary:\n")
print(summary_data)
