parse_cli_args <- function(argv) {
  parsed <- list(
    data_root = "data/benchmark_inputs",
    max_mixtures = NA_integer_
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--") || i == length(argv)) {
      stop("Arguments must be provided as --key value pairs")
    }
    value <- argv[[i + 1L]]
    name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    parsed[[name]] <- value
    i <- i + 2L
  }
  if (!is.null(parsed$max_mixtures) && !is.na(parsed$max_mixtures)) {
    parsed$max_mixtures <- as.integer(parsed$max_mixtures)
  }
  parsed
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$output_dir) || is.null(args$dataset) || is.null(args$tissue) || is.null(args$method)) {
  stop("--output-dir, --dataset, --tissue, and --method are required")
}

read_matrix <- function(path) {
  read.table(path, header = TRUE, row.names = 1, sep = "\t",
             check.names = FALSE, quote = "", comment.char = "")
}

normalize_nonnegative <- function(x) {
  x <- pmax(as.numeric(x), 0)
  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    return(rep(1 / length(x), length(x)))
  }
  x / total
}

normalize_rows <- function(mat) {
  mat <- as.matrix(mat)
  rn <- rownames(mat)
  cn <- colnames(mat)
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
  pearson <- NA_real_
  if (length(pred) > 1 && sd(pred) > 0 && sd(truth) > 0) {
    pearson <- suppressWarnings(cor(pred, truth))
  }
  c(
    l1 = l1,
    accuracy = 1 - l1 / 2,
    mae = mean(abs(diff)),
    rmse = sqrt(mean(diff^2)),
    pearson = pearson
  )
}

write_prediction <- function(pred, out_dir, method) {
  prediction_dir <- file.path(out_dir, "r_baseline_predictions")
  dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(pred, file.path(prediction_dir, paste0(method, ".txt")),
              sep = "\t", quote = FALSE, col.names = NA)
}

run_lindeconseq <- function(signature, bulk) {
  if (!requireNamespace("LinDeconSeq", quietly = TRUE)) {
    stop("LinDeconSeq is not available in this R environment")
  }
  suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(MASS)
  })
  common <- intersect(rownames(signature), rownames(bulk))
  fractions <- LinDeconSeq::deconSeq(
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

run_music <- function(sc, labels, subjects, bulk, cell_types) {
  if (!requireNamespace("MuSiC", quietly = TRUE) ||
      !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("MuSiC/SingleCellExperiment/S4Vectors are not available in this R environment")
  }
  suppressPackageStartupMessages({
    library(MuSiC)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
  })
  common <- intersect(rownames(sc), rownames(bulk))
  sc <- as.matrix(sc[common, , drop = FALSE])
  bulk <- as.matrix(bulk[common, , drop = FALSE])
  keep <- labels %in% cell_types
  sc <- sc[, keep, drop = FALSE]
  labels <- labels[keep]
  subjects <- subjects[keep]

  if (length(unique(subjects)) < 2) {
    subjects <- paste0("sample", ((seq_along(labels) - 1) %% 3) + 1)
  }

  meta <- data.frame(
    cell_type = labels,
    sampleID = subjects,
    row.names = colnames(sc),
    stringsAsFactors = FALSE
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = sc),
    colData = S4Vectors::DataFrame(meta)
  )
  est <- MuSiC::music_prop(
    bulk.mtx = bulk,
    sc.sce = sce,
    clusters = "cell_type",
    samples = "sampleID",
    select.ct = cell_types,
    verbose = FALSE
  )
  pred <- as.matrix(est[["Est.prop.weighted"]])
  pred <- pred[, cell_types, drop = FALSE]
  normalize_rows(pred)
}

sampled_dir <- file.path(args$output_dir, "sampled_cases", args$dataset, args$tissue)
original_dir <- file.path(args$data_root, args$dataset, args$tissue)
bulk <- read_matrix(file.path(sampled_dir, "sampled_bulk.txt"))
truth <- read_matrix(file.path(sampled_dir, "sampled_bulkRatio.txt"))
signature <- read_matrix(file.path(sampled_dir, "selected_cellTypeExpr.txt"))

if (!is.na(args$max_mixtures)) {
  keep_mixtures <- head(colnames(bulk), args$max_mixtures)
  bulk <- bulk[, keep_mixtures, drop = FALSE]
  truth <- truth[, keep_mixtures, drop = FALSE]
}

cell_types <- intersect(colnames(signature), rownames(truth))
signature <- signature[, cell_types, drop = FALSE]
truth <- truth[cell_types, colnames(bulk), drop = FALSE]
original_bulk_names <- colnames(bulk)
if (ncol(bulk) == 1L) {
  bulk <- cbind(bulk, bulk[, 1, drop = FALSE])
  colnames(bulk) <- c(original_bulk_names, paste0(original_bulk_names, "__dup"))
}

start <- Sys.time()
method <- args$method
if (method == "LinDeconSeq") {
  pred <- run_lindeconseq(signature, bulk)
} else if (method == "MuSiC") {
  sc <- read_matrix(file.path(original_dir, "singleCellExpr.txt"))
  selected_genes <- rownames(signature)
  sc <- sc[intersect(selected_genes, rownames(sc)), , drop = FALSE]
  labels_df <- read.table(file.path(original_dir, "singleCellLabels.txt"), header = TRUE, sep = "\t",
                          check.names = FALSE, quote = "", comment.char = "",
                          stringsAsFactors = FALSE)
  labels <- as.character(labels_df[[1]])
  subjects_path <- file.path(original_dir, "singleCellSubjects.txt")
  if (file.exists(subjects_path)) {
    subjects_df <- read.table(subjects_path, header = TRUE, sep = "\t",
                              check.names = FALSE, quote = "", comment.char = "",
                              stringsAsFactors = FALSE)
    subjects <- as.character(subjects_df[[1]])
  } else {
    subjects <- rep("sample1", length(labels))
  }
  pred <- run_music(sc, labels, subjects, bulk, cell_types)
} else {
  stop("Unsupported method: ", method)
}
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

pred <- pred[original_bulk_names, cell_types, drop = FALSE]
bulk_for_metrics <- bulk[, original_bulk_names, drop = FALSE]
write_prediction(pred, sampled_dir, method)

metric_rows <- lapply(seq_len(nrow(pred)), function(i) {
  mixture <- rownames(pred)[i]
  scores <- score_prediction(pred[i, cell_types], truth[cell_types, mixture])
  data.frame(
    dataset = args$dataset,
    tissue = args$tissue,
    case = paste(args$dataset, args$tissue, sep = "/"),
    mixture = mixture,
    mixture_index = i,
    method = method,
    n_genes = nrow(signature),
    n_cell_types = length(cell_types),
    l1 = unname(scores[["l1"]]),
    accuracy = unname(scores[["accuracy"]]),
    mae = unname(scores[["mae"]]),
    rmse = unname(scores[["rmse"]]),
    pearson = unname(scores[["pearson"]]),
    runtime_sec = elapsed / nrow(pred),
    stringsAsFactors = FALSE
  )
})
metrics <- do.call(rbind, metric_rows)

out_file <- file.path(sampled_dir, paste0("r_baseline_", method, "_metrics_long.csv"))
write.csv(metrics, out_file, row.names = FALSE)

summary <- aggregate(
  cbind(l1, accuracy, mae, rmse, pearson, runtime_sec) ~ dataset + tissue + case + method,
  data = metrics,
  FUN = mean,
  na.rm = TRUE
)
summary$mixtures <- nrow(metrics)
summary_file <- file.path(sampled_dir, paste0("r_baseline_", method, "_case_summary.csv"))
write.csv(summary, summary_file, row.names = FALSE)

cat("Wrote", out_file, "\n")
cat("Wrote", summary_file, "\n")
