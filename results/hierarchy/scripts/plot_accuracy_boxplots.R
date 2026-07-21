suppressPackageStartupMessages({
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
data_path <- file.path(project_dir, "data", "outputs", "all_methods_version2_metrics.tsv")
plot_dir <- file.path(project_dir, "plots")
plot_font <- "Times New Roman"

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
  NNLS = "NNLS",
  DWLS = "DWLS",
  LinDeconSeq = "LinDeconSeq",
  BayesPrism = "BayesPrism",
  HiDecon = "HiDecon",
  HiDecon_hierarchy = "HiDecon-H",
  MuSic = "MuSiC",
  MuSic_hierarchy = "MuSiC-H",
  MetaSort_jk_v2 = "MetaSort",
  MetaSort_jk_hierarchy_v2 = "MetaSort-H"
)
tissue_labels <- c(
  Salivary_Gland = "Salivary gland",
  Thymus = "Thymus"
)

metrics <- read.delim(data_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("tissue", "method", "mixture", "L1", "Accuracy", "MAE", "RMSE")
if (!identical(names(metrics), required_columns)) {
  stop("Unexpected columns in ", data_path)
}
if (!identical(sort(unique(metrics$tissue)), sort(names(tissue_labels)))) {
  stop("The input does not contain exactly the expected tissues.")
}
if (!identical(sort(unique(metrics$method)), sort(method_levels))) {
  stop("The input does not contain exactly the expected methods.")
}
group_sizes <- table(metrics$tissue, metrics$method)
if (any(group_sizes != 30L)) {
  stop("Each tissue-method combination must contain exactly 30 mixtures.")
}
numeric_columns <- c("L1", "Accuracy", "MAE", "RMSE")
if (any(!is.finite(as.matrix(metrics[numeric_columns])))) {
  stop("Metric columns must contain only finite numeric values.")
}
if (!isTRUE(all.equal(metrics$Accuracy, 1 - metrics$L1 / 2, tolerance = 1e-12))) {
  stop("Accuracy values are inconsistent with 1 - L1 / 2.")
}

metrics$method <- factor(metrics$method, levels = method_levels)
metrics$method_label <- factor(
  unname(method_labels[as.character(metrics$method)]),
  levels = unname(method_labels)
)
metrics$tissue_label <- factor(
  unname(tissue_labels[metrics$tissue]),
  levels = unname(tissue_labels)
)

fill_values <- c(
  NNLS = "#D9D9D9",
  DWLS = "#D9D9D9",
  LinDeconSeq = "#D9D9D9",
  BayesPrism = "#D9D9D9",
  HiDecon = "#D9D9D9",
  `HiDecon-H` = "#C7C7C7",
  MuSiC = "#D9D9D9",
  `MuSiC-H` = "#C7C7C7",
  MetaSort = "#4C78A8",
  `MetaSort-H` = "#E45756"
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
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 4)
    ),
    plot.margin = margin(4, 4, 4, 4)
  )

plot_one_tissue <- function(tissue_name) {
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
    labs(title = tissue_name, y = "Accuracy (1 - L1 / 2)") +
    base_theme
}

save_plot <- function(tissue_key) {
  tissue_name <- unname(tissue_labels[[tissue_key]])
  file_stem <- paste0(
    "all_methods_version2_accuracy_boxplot_bibm_",
    tolower(tissue_key)
  )
  figure <- plot_one_tissue(tissue_name)

  ggsave(
    file.path(plot_dir, paste0(file_stem, ".pdf")),
    figure,
    device = cairo_pdf,
    width = 3.7,
    height = 3.2,
    units = "in",
    family = plot_font
  )
  ggsave(
    file.path(plot_dir, paste0(file_stem, ".png")),
    figure,
    width = 3.7,
    height = 3.2,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(names(tissue_labels), save_plot))
cat("Wrote four figure files to ", plot_dir, "\n", sep = "")
