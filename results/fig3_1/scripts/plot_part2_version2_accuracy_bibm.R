suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

version_dir <- "results/fig3_1/version_2"
plot_dir <- file.path(version_dir, "plots")
metrics_path <- file.path(version_dir, "all_methods_version2_metrics.tsv")
plot_font <- "Times New Roman"

metrics <- read.delim(metrics_path, check.names = FALSE, stringsAsFactors = FALSE)

method_levels <- c(
  "NNLS",
  "DWLS",
  "LinDeconSeq",
  "BayesPrism",
  "HiDecon",
  "HiDecon_hierarchy",
  "MuSic",
  "MuSic_hierarchy",
  "MetaSort_jk_v2",
  "MetaSort_jk_hierarchy_v2"
)
method_labels <- c(
  "NNLS",
  "DWLS",
  "LinDeconSeq",
  "BayesPrism",
  "HiDecon",
  "HiDecon-H",
  "MuSiC",
  "MuSiC-H",
  "MetaSort",
  "MetaSort-H"
)
names(method_labels) <- method_levels

tissue_labels <- c(
  Salivary_Gland = "Salivary gland",
  Thymus = "Thymus"
)

metrics <- metrics[metrics$method %in% method_levels, ]
metrics$method <- factor(metrics$method, levels = method_levels)
metrics$method_label <- factor(
  method_labels[as.character(metrics$method)],
  levels = unname(method_labels)
)
metrics$tissue_label <- factor(
  tissue_labels[metrics$tissue],
  levels = unname(tissue_labels)
)

fill_values <- c(
  "NNLS" = "#D9D9D9",
  "DWLS" = "#D9D9D9",
  "LinDeconSeq" = "#D9D9D9",
  "BayesPrism" = "#D9D9D9",
  "HiDecon" = "#D9D9D9",
  "HiDecon-H" = "#C7C7C7",
  "MuSiC" = "#D9D9D9",
  "MuSiC-H" = "#C7C7C7",
  "MetaSort" = "#4C78A8",
  "MetaSort-H" = "#E45756"
)

base_theme <- theme_classic(base_family = plot_font, base_size = 9) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, margin = margin(r = 4)),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.tag = element_text(size = 10, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(4, 4, 4, 4)
  )

plot_one_tissue <- function(tissue_name, title_text) {
  ggplot(
    metrics[metrics$tissue_label == tissue_name, ],
    aes(x = method_label, y = Accuracy, fill = method_label)
  ) +
    geom_boxplot(
      width = 0.56,
      outlier.shape = NA,
      linewidth = 0.35,
      color = "black",
      alpha = 0.95
    ) +
    geom_point(
      position = position_jitter(width = 0.09, height = 0, seed = 1),
      size = 0.55,
      alpha = 0.45,
      color = "grey20",
      stroke = 0
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_y_continuous(
      limits = c(0.35, 0.85),
      breaks = seq(0.4, 0.8, by = 0.1),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(title = title_text, y = "Accuracy (1 - L1 / 2)") +
    base_theme
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, file_stem, width, height) {
  pdf_path <- file.path(plot_dir, paste0(file_stem, ".pdf"))
  png_path <- file.path(plot_dir, paste0(file_stem, ".png"))

  ggsave(
    pdf_path,
    plot,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in",
    family = plot_font
  )

  ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  c(pdf_path, png_path)
}

p_salivary <- plot_one_tissue("Salivary gland", "Salivary gland")
p_thymus <- plot_one_tissue("Thymus", "Thymus")
p_two_panel <- p_salivary + p_thymus + plot_layout(ncol = 2, widths = c(1, 1))

written <- c(
  save_plot(
    p_salivary,
    "all_methods_version2_accuracy_boxplot_bibm_salivary_gland",
    width = 3.7,
    height = 3.2
  ),
  save_plot(
    p_thymus,
    "all_methods_version2_accuracy_boxplot_bibm_thymus",
    width = 3.7,
    height = 3.2
  ),
  save_plot(
    p_two_panel,
    "all_methods_version2_accuracy_boxplot_bibm_two_panel",
    width = 7.4,
    height = 3.2
  )
)

cat("Wrote:\n", paste(written, collapse = "\n"), "\n", sep = "")
