suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Run this file with Rscript")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
project_dir <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(project_dir, "data")
plot_dir <- file.path(project_dir, "figures")

font_family <- "Times New Roman"
base_font_size <- 9
horizontal_width_in <- 7.0
horizontal_height_in <- 2.35
single_column_width_in <- 3.55
single_column_height_in <- 3.75

detail_path <- file.path(data_dir, "rare_celltype_truth_prediction.csv")

detail <- read_csv(detail_path, show_col_types = FALSE)

method_levels <- c("MetaSort_jk", "DWLS", "MuSic", "LinDeconSeq", "bayesprism", "NNLS", "HiDecon")

safe_cor <- function(x, y) {
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  cor(x, y, method = "pearson")
}

cors <- detail %>%
  filter(method %in% method_levels) %>%
  group_by(dataset, tissue, cell_type, method) %>%
  summarise(
    n = n(),
    pearson_r = safe_cor(truth, prediction),
    mae = mean(abs(prediction - truth)),
    .groups = "drop"
  )

method_labels <- c(
  MetaSort_jk = "MetaSort",
  DWLS = "DWLS",
  MuSic = "MuSiC",
  LinDeconSeq = "LinDeconSeq",
  bayesprism = "BayesPrism",
  NNLS = "NNLS",
  HiDecon = "HiDecon"
)

method_colors <- c(
  MetaSort_jk = "#E45756",
  DWLS = "#4C78A8",
  MuSic = "#54A24B",
  LinDeconSeq = "#B279A2",
  bayesprism = "#F2CF5B",
  NNLS = "#72B7B2",
  HiDecon = "#8C8C8C"
)

format_name <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

format_dataset_label <- function(x) {
  ifelse(x == "CELlxGENE", "CELLxGENE", x)
}

format_r <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.2f", x))
}

format_mae <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.3f", x))
}

wrap_title <- function(x) {
  paste(strwrap(x, width = 80), collapse = "\n")
}

make_plot <- function(
  dataset_name,
  tissue_name,
  cell_type_name,
  file_stem,
  output_suffix = "_bibm",
  width_in = horizontal_width_in,
  height_in = horizontal_height_in,
  facet_nrow = 1,
  facet_ncol = NULL,
  point_size = 0.9,
  title_size = 10,
  axis_title_size = 8.8,
  axis_text_size = 7.2,
  strip_size = 7.2,
  panel_spacing_x = 0.55,
  panel_spacing_y = 0.62,
  title_wrap_width = 80
) {
  plot_data <- detail %>%
    filter(
      dataset == dataset_name,
      tissue == tissue_name,
      cell_type == cell_type_name,
      method %in% method_levels
    ) %>%
    mutate(
      method = factor(method, levels = method_levels)
    )

  stat_data <- cors %>%
    filter(
      dataset == dataset_name,
      tissue == tissue_name,
      cell_type == cell_type_name,
      method %in% method_levels
    ) %>%
    mutate(
      method = factor(method, levels = method_levels),
      facet_label = sprintf(
        "%s\nr=%s\nMAE=%s",
        method_labels[as.character(method)],
        format_r(pearson_r),
        format_mae(mae)
      )
    ) %>%
    arrange(method)

  label_map <- setNames(stat_data$facet_label, as.character(stat_data$method))
  plot_data$facet_label <- factor(
    label_map[as.character(plot_data$method)],
    levels = stat_data$facet_label
  )

  max_value <- max(plot_data$truth, plot_data$prediction, na.rm = TRUE)
  axis_upper <- ceiling((max_value * 1.05) / 0.02) * 0.02
  axis_upper <- max(axis_upper, 0.08)
  axis_limit <- axis_upper * 1.04
  breaks <- seq(0, axis_upper, length.out = 3)
  x_label_fn <- function(x) {
    ifelse(
      abs(x) < sqrt(.Machine$double.eps),
      "",
      scales::number_format(accuracy = 0.01)(x)
    )
  }

  title_text <- sprintf(
    "%s, %s: %s",
    format_dataset_label(dataset_name),
    format_name(tissue_name),
    format_name(cell_type_name)
  )
  title_text <- paste(strwrap(title_text, width = title_wrap_width), collapse = "\n")

  facet_layer <- if (is.null(facet_ncol)) {
    facet_wrap(~ facet_label, nrow = facet_nrow, axes = "margins", axis.labels = "margins")
  } else {
    facet_wrap(~ facet_label, ncol = facet_ncol, axes = "margins", axis.labels = "margins")
  }

  p <- ggplot(plot_data, aes(x = truth, y = prediction, color = method)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = 0.35,
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_point(size = point_size, alpha = 0.68, stroke = 0) +
    facet_layer +
    scale_color_manual(values = method_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, axis_limit),
      breaks = breaks,
      labels = x_label_fn,
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, axis_limit),
      breaks = breaks,
      labels = scales::number_format(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.01))
    ) +
    coord_equal(clip = "off") +
    labs(
      title = title_text,
      x = "True proportion",
      y = "Predicted proportion"
    ) +
    theme_classic(base_family = font_family, base_size = base_font_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5, margin = margin(b = 5), lineheight = 0.95),
      axis.title = element_text(size = axis_title_size),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      axis.text = element_text(size = axis_text_size, color = "black"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = strip_size, face = "bold", margin = margin(b = 3), lineheight = 0.95),
      panel.grid.major = element_line(linewidth = 0.22, color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(panel_spacing_x, "lines"),
      panel.spacing.y = unit(panel_spacing_y, "lines"),
      plot.margin = margin(4, 5, 4, 5)
    )

  pdf_path <- file.path(plot_dir, paste0(file_stem, output_suffix, ".pdf"))
  png_path <- file.path(plot_dir, paste0(file_stem, output_suffix, ".png"))

  ggsave(
    pdf_path,
    p,
    device = cairo_pdf,
    width = width_in,
    height = height_in,
    units = "in",
    family = font_family
  )

  ggsave(
    png_path,
    p,
    width = width_in,
    height = height_in,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  cat("Wrote:\n", pdf_path, "\n", png_path, "\n", sep = "")
}

make_plot(
  dataset_name = "CELlxGENE",
  tissue_name = "Basal_Zone_Of_Heart",
  cell_type_name = "vein_endothelial_cell",
  file_stem = "CELlxGENE__Basal_Zone_Of_Heart__vein_endothelial_cell"
)

make_plot(
  dataset_name = "Tabula-Sapiens",
  tissue_name = "Skin",
  cell_type_name = "CD1c_positive_myeloid_dendritic_cell",
  file_stem = "Tabula-Sapiens__Skin__CD1c_positive_myeloid_dendritic_cell"
)
