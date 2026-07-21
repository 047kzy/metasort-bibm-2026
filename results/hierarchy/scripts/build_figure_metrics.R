script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_dir <- file.path(project_dir, "data", "inputs")
results_dir <- file.path(project_dir, "results")
output_dir <- file.path(project_dir, "data", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "all_methods_version2_metrics.tsv")
summary_path <- file.path(output_dir, "all_methods_version2_summary.tsv")
tissues <- c("Salivary_Gland", "Thymus")
n_mixtures <- 30L

method_files <- c(
  NNLS = "NNLS.tsv",
  DWLS = "DWLS.tsv",
  LinDeconSeq = "LinDeconSeq.tsv",
  BayesPrism = "bayesprism.tsv",
  HiDecon = "HiDecon.tsv",
  HiDecon_hierarchy = "HiDecon_hierarchy.tsv",
  MuSic = "MuSic.tsv",
  MuSic_hierarchy = "MuSic_hierarchy.tsv",
  MetaSort_jk_v2 = "MetaSort_jk_avg2p5_tol0005.tsv",
  MetaSort_jk_hierarchy_v2 = "MetaSort_jk_hierarchy_avg2p5_tol0005.tsv"
)

read_matrix <- function(path) {
  as.matrix(read.table(path, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, quote = "", comment.char = ""))
}
clean_names <- function(x) sub("^Bulk_([0-9]+)$", "Mixture\\1", gsub('^"|"$', "", x))
normalize_rows <- function(x) {
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  totals <- rowSums(x)
  valid <- is.finite(totals) & totals > 0
  x[valid, ] <- x[valid, , drop = FALSE] / totals[valid]
  x
}

all_metrics <- list()
for (tissue in tissues) {
  truth <- read_matrix(file.path(input_dir, tissue, "bulkRatio.tsv"))
  rownames(truth) <- clean_names(rownames(truth))
  colnames(truth) <- clean_names(colnames(truth))
  for (method in names(method_files)) {
    pred <- read_matrix(file.path(results_dir, tissue, method_files[[method]]))
    rownames(pred) <- clean_names(rownames(pred))
    colnames(pred) <- clean_names(colnames(pred))
    mixtures <- intersect(colnames(truth), rownames(pred))
    mixtures <- mixtures[seq_len(min(length(mixtures), n_mixtures))]
    cell_types <- intersect(rownames(truth), colnames(pred))
    if (length(mixtures) == 0L || length(cell_types) == 0L) {
      stop("No matching mixtures or cell types for ", tissue, " / ", method)
    }
    pred_sub <- normalize_rows(pred[mixtures, cell_types, drop = FALSE])
    truth_sub <- t(truth[cell_types, mixtures, drop = FALSE])
    rows <- lapply(mixtures, function(mixture) {
      difference <- pred_sub[mixture, ] - truth_sub[mixture, ]
      l1 <- sum(abs(difference))
      data.frame(
        tissue = tissue,
        method = method,
        mixture = mixture,
        L1 = l1,
        Accuracy = 1 - l1 / 2,
        MAE = mean(abs(difference)),
        RMSE = sqrt(mean(difference^2)),
        stringsAsFactors = FALSE
      )
    })
    all_metrics[[length(all_metrics) + 1L]] <- do.call(rbind, rows)
  }
}

metrics <- do.call(rbind, all_metrics)
write.table(metrics, output_path, sep = "\t", quote = FALSE, row.names = FALSE)

summary_input <- metrics
summary_input$tissue <- factor(summary_input$tissue, levels = tissues)
summary_input$method <- factor(summary_input$method, levels = names(method_files))
summary <- aggregate(
  cbind(Accuracy, L1, MAE, RMSE) ~ tissue + method,
  data = summary_input,
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
write.table(summary_flat, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(
  "Wrote ", nrow(metrics), " metric rows and ", nrow(summary_flat),
  " summary rows.\n", sep = ""
)
