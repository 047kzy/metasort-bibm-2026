suppressPackageStartupMessages({
  library(ggplot2)
})

plot_font <- "Times New Roman"
base_dir <- "results/fig3_1"
version_dir <- file.path(base_dir, "version_2")
plot_dir <- file.path(version_dir, "plots")
tissues <- c("Salivary_Gland", "Thymus")
n_mixtures <- 30L

method_files <- c(
  NNLS = "NNLS.txt",
  DWLS = "DWLS.txt",
  LinDeconSeq = "LinDeconSeq.txt",
  BayesPrism = "bayesprism.txt",
  HiDecon = "HiDecon.txt",
  HiDecon_hierarchy = "HiDecon_hierarchy.txt",
  MuSic = "MuSic.txt",
  MuSic_hierarchy = "MuSic_hierarchy.txt",
  MetaSort_jk_v2 = "MetaSort_jk_avg2p5_tol0005.txt",
  MetaSort_jk_hierarchy_v2 = "MetaSort_jk_hierarchy_avg2p5_tol0005.txt"
)

read_matrix <- function(path) {
  as.matrix(read.table(path, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, comment.char = ""))
}

clean_names <- function(x) {
  x <- gsub('^"|"$', "", x)
  sub("^Bulk_([0-9]+)$", "Mixture\\1", x)
}

normalize_rows <- function(x) {
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  sums <- rowSums(x)
  ok <- is.finite(sums) & sums > 0
  x[ok, ] <- x[ok, , drop = FALSE] / sums[ok]
  x
}

calc_metrics <- function(truth, pred, tissue, method) {
  rownames(truth) <- clean_names(rownames(truth))
  colnames(truth) <- clean_names(colnames(truth))
  rownames(pred) <- clean_names(rownames(pred))
  colnames(pred) <- clean_names(colnames(pred))

  common_mixtures <- intersect(colnames(truth), rownames(pred))
  common_mixtures <- common_mixtures[seq_len(min(length(common_mixtures), n_mixtures))]
  common_cell_types <- intersect(rownames(truth), colnames(pred))
  if (length(common_mixtures) == 0 || length(common_cell_types) == 0) {
    stop("No overlap for ", tissue, " / ", method)
  }

  pred_sub <- normalize_rows(pred[common_mixtures, common_cell_types, drop = FALSE])
  truth_sub <- t(truth[common_cell_types, common_mixtures, drop = FALSE])

  rows <- lapply(common_mixtures, function(mix) {
    diff <- pred_sub[mix, ] - truth_sub[mix, ]
    l1 <- sum(abs(diff))
    data.frame(
      tissue = tissue,
      method = method,
      mixture = mix,
      L1 = l1,
      Accuracy = 1 - l1 / 2,
      MAE = mean(abs(diff)),
      RMSE = sqrt(mean(diff^2)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

result_path_for <- function(tissue, method, file_name) {
  if (grepl("_v2$", method)) {
    return(file.path(version_dir, tissue, "result", file_name))
  }
  first30_path <- file.path(base_dir, tissue, "result", paste0(method, "_first30.txt"))
  if (file.exists(first30_path)) {
    return(first30_path)
  }
  file.path(base_dir, tissue, "result", file_name)
}

all_metrics <- list()
for (tissue in tissues) {
  truth <- read_matrix(file.path(base_dir, tissue, "bulkRatio.txt"))
  for (method in names(method_files)) {
    result_path <- result_path_for(tissue, method, method_files[[method]])
    if (!file.exists(result_path)) {
      warning("Skipping missing result: ", result_path)
      next
    }
    pred <- read_matrix(result_path)
    all_metrics[[length(all_metrics) + 1L]] <- calc_metrics(truth, pred, tissue, method)
  }
}

metrics <- do.call(rbind, all_metrics)
metrics$method <- factor(metrics$method, levels = names(method_files))
metrics$tissue <- factor(metrics$tissue, levels = tissues)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  metrics,
  file.path(version_dir, "all_methods_version2_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- aggregate(
  cbind(Accuracy, L1, MAE, RMSE) ~ tissue + method,
  data = metrics,
  FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x))
)
summary_flat <- data.frame(
  tissue = summary$tissue,
  method = summary$method,
  mean_accuracy = summary$Accuracy[, "mean"],
  median_accuracy = summary$Accuracy[, "median"],
  sd_accuracy = summary$Accuracy[, "sd"],
  mean_l1 = summary$L1[, "mean"],
  median_l1 = summary$L1[, "median"],
  sd_l1 = summary$L1[, "sd"],
  mean_mae = summary$MAE[, "mean"],
  median_mae = summary$MAE[, "median"],
  sd_mae = summary$MAE[, "sd"],
  mean_rmse = summary$RMSE[, "mean"],
  median_rmse = summary$RMSE[, "median"],
  sd_rmse = summary$RMSE[, "sd"],
  stringsAsFactors = FALSE
)
write.table(
  summary_flat,
  file.path(version_dir, "all_methods_version2_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_metric <- function(metric, x_label, file_stem, higher_is_better = TRUE) {
  p <- ggplot(metrics, aes(x = .data[[metric]], y = method, fill = method)) +
    geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(height = 0.16, size = 0.8, alpha = 0.45, color = "gray20") +
    facet_wrap(~ tissue, ncol = 1, scales = "free_x") +
    labs(x = x_label, y = NULL) +
    theme_bw(base_family = plot_font, base_size = 11) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray92", color = "gray65"),
      panel.grid.major.y = element_blank()
    )
  if (higher_is_better) {
    p <- p + coord_cartesian(xlim = c(max(0, min(metrics[[metric]]) - 0.02), 1))
  }
  ggsave(file.path(plot_dir, paste0(file_stem, ".png")), p, width = 9.5, height = 8.5, dpi = 300)
  ggsave(
    file.path(plot_dir, paste0(file_stem, ".pdf")),
    p,
    device = cairo_pdf,
    width = 9.5,
    height = 8.5,
    family = plot_font
  )
}

plot_metric("Accuracy", "Accuracy (1 - L1 / 2)", "all_methods_version2_accuracy_boxplot", TRUE)
plot_metric("MAE", "Mean absolute error", "all_methods_version2_mae_boxplot", FALSE)

cat("Wrote version_2 all-method metrics and boxplots to ", version_dir, "\n", sep = "")
