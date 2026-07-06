suppressPackageStartupMessages(library(LinDeconSeq))

base_dir <- "results/fig3_1"
out_dir <- file.path(base_dir, "hierarchy_reference_level_analysis")
tissues <- c("Salivary_Gland", "Thymus")
levels <- c("coarse", "leaf")
n_mixtures <- 30L
method <- "LinDeconSeq"

read_matrix <- function(path) {
  read.table(path, header = TRUE, row.names = 1, sep = "\t",
             check.names = FALSE, quote = "", comment.char = "")
}

normalize_nonnegative <- function(x) {
  x <- pmax(as.numeric(x), 0)
  total <- sum(x)
  if (!is.finite(total) || total <= 0) return(rep(1 / length(x), length(x)))
  x / total
}

normalize_rows <- function(mat) {
  mat <- as.matrix(mat)
  rn <- rownames(mat)
  cn <- colnames(mat)
  mat[!is.finite(mat)] <- 0
  out <- t(apply(mat, 1, normalize_nonnegative))
  rownames(out) <- rn
  colnames(out) <- cn
  out
}

score_prediction <- function(pred, truth) {
  pred <- normalize_nonnegative(pred)
  truth <- normalize_nonnegative(truth)
  diff <- pred - truth
  l1 <- sum(abs(diff))
  c(L1 = l1, Accuracy = 1 - l1 / 2, MAE = mean(abs(diff)), RMSE = sqrt(mean(diff^2)))
}

score_table <- function(pred, truth, tissue, level) {
  rows <- lapply(seq_len(nrow(pred)), function(i) {
    s <- score_prediction(pred[i, colnames(pred)], truth[colnames(pred), rownames(pred)[i]])
    data.frame(
      tissue = tissue,
      level = level,
      method = method,
      mixture = rownames(pred)[i],
      mixture_index = i,
      n_reference_celltypes = ncol(pred),
      L1 = s[["L1"]],
      Accuracy = s[["Accuracy"]],
      MAE = s[["MAE"]],
      RMSE = s[["RMSE"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

load_mapping <- function(tissue, level) {
  mapping_df <- read.table(file.path(out_dir, "hierarchy_reference_level_node_mapping.tsv"),
                           header = TRUE, sep = "\t", quote = "",
                           comment.char = "", check.names = FALSE,
                           stringsAsFactors = FALSE)
  sub <- mapping_df[mapping_df$tissue == tissue & mapping_df$level == level, , drop = FALSE]
  mapping <- strsplit(sub$leaf_celltypes, "\\|", fixed = FALSE)
  names(mapping) <- sub$node
  mapping
}

aggregate_leaf_prediction <- function(leaf_pred, mapping) {
  mat <- sapply(mapping, function(leaves) rowSums(leaf_pred[, leaves, drop = FALSE]))
  as.data.frame(mat, check.names = FALSE)
}

run_lindeconseq <- function(signature, bulk) {
  common <- intersect(rownames(bulk), rownames(signature))
  fractions <- deconSeq(
    as.matrix(bulk[common, , drop = FALSE]),
    as.matrix(signature[common, , drop = FALSE]),
    verbose = FALSE
  )
  fractions <- as.matrix(fractions)
  if (nrow(fractions) != ncol(bulk) && ncol(fractions) == ncol(bulk)) {
    fractions <- t(fractions)
  }
  rownames(fractions) <- colnames(bulk)
  colnames(fractions) <- colnames(signature)
  normalize_rows(fractions)
}

all_metrics <- list()
all_consistency <- list()
all_logs <- list()

for (tissue in tissues) {
  message("[LOAD] ", tissue)
  case_dir <- file.path(base_dir, tissue)
  tissue_out <- file.path(out_dir, tissue)
  bulk <- read_matrix(file.path(case_dir, "bulk.txt"))
  bulk <- bulk[, seq_len(min(n_mixtures, ncol(bulk))), drop = FALSE]
  level_pred <- list()

  for (level in levels) {
    pred_path <- file.path(tissue_out, paste0(tissue, "_", level, "_", method, "_prediction.tsv"))
    truth <- read_matrix(file.path(tissue_out, paste0(tissue, "_", level, "_truth.tsv")))
    start <- Sys.time()
    status <- "OK"
    error <- ""
    pred <- NULL
    tryCatch({
      if (file.exists(pred_path)) {
        pred <- as.matrix(read_matrix(pred_path))
        if (nrow(pred) != ncol(bulk)) stop("existing prediction has wrong row count")
        message("[LOAD] ", tissue, " / ", level, " / ", method)
      } else {
        signature <- read_matrix(file.path(tissue_out, paste0(tissue, "_", level, "_reference_signature.tsv")))
        message("[RUN] ", tissue, " / ", level, " / ", method)
        pred <- run_lindeconseq(signature, bulk)
      }
      pred <- pred[, rownames(truth), drop = FALSE]
      rownames(pred) <- colnames(bulk)
      write.table(pred, pred_path, sep = "\t", quote = FALSE, col.names = NA)
      level_pred[[level]] <- pred
      all_metrics[[length(all_metrics) + 1L]] <- score_table(pred, truth, tissue, level)
    }, error = function(e) {
      status <<- "ERROR"
      error <<- conditionMessage(e)
      message("[ERROR] ", tissue, " / ", level, " / ", method, ": ", error)
    })
    all_logs[[length(all_logs) + 1L]] <- data.frame(
      tissue = tissue,
      level = level,
      method = method,
      status = status,
      rows = if (is.null(pred)) NA_integer_ else nrow(pred),
      cols = if (is.null(pred)) NA_integer_ else ncol(pred),
      seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 3),
      error = error,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(level_pred[["leaf"]])) {
    for (level in c("coarse")) {
      if (is.null(level_pred[[level]])) next
      mapping <- load_mapping(tissue, level)
      leaf_agg <- aggregate_leaf_prediction(level_pred[["leaf"]], mapping)
      direct_pred <- level_pred[[level]]
      rows <- lapply(seq_len(nrow(direct_pred)), function(i) {
        s <- score_prediction(direct_pred[i, ], leaf_agg[rownames(direct_pred)[i], colnames(direct_pred)])
        data.frame(
          tissue = tissue,
          level = level,
          method = method,
          comparison = "direct_ref_level_vs_leaf_aggregated",
          mixture = rownames(direct_pred)[i],
          mixture_index = i,
          L1 = s[["L1"]],
          Consistency = s[["Accuracy"]],
          MAE = s[["MAE"]],
          RMSE = s[["RMSE"]],
          stringsAsFactors = FALSE
        )
      })
      all_consistency[[length(all_consistency) + 1L]] <- do.call(rbind, rows)
    }
  }
}

metrics_df <- do.call(rbind, all_metrics)
consistency_df <- do.call(rbind, all_consistency)
log_df <- do.call(rbind, all_logs)

write.table(metrics_df, file.path(out_dir, "hierarchy_reference_level_lindeconseq_metrics_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(consistency_df, file.path(out_dir, "hierarchy_reference_level_lindeconseq_consistency_vs_leaf_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(log_df, file.path(out_dir, "hierarchy_reference_level_lindeconseq_run_log.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

metric_summary <- aggregate(
  cbind(Accuracy, MAE, RMSE, L1) ~ tissue + level + method + n_reference_celltypes,
  data = metrics_df,
  FUN = function(x) c(mean = mean(x), median = stats::median(x), sd = stats::sd(x))
)
metric_summary <- do.call(data.frame, metric_summary)
write.table(metric_summary, file.path(out_dir, "hierarchy_reference_level_lindeconseq_metric_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

consistency_summary <- aggregate(
  cbind(Consistency, MAE, RMSE, L1) ~ tissue + level + method,
  data = consistency_df,
  FUN = function(x) c(mean = mean(x), median = stats::median(x), sd = stats::sd(x))
)
consistency_summary <- do.call(data.frame, consistency_summary)
write.table(consistency_summary, file.path(out_dir, "hierarchy_reference_level_lindeconseq_consistency_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

print(metric_summary)
print(consistency_summary)
