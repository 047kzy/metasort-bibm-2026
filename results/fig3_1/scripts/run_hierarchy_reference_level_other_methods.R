suppressPackageStartupMessages({
  library(jsonlite)
  library(nnls)
  library(quadprog)
  library(Biobase)
  library(MuSiC)
  library(SingleCellExperiment)
  library(HiDecon)
  library(BayesPrism)
})

base_dir <- "results/fig3_1"
out_dir <- file.path(base_dir, "hierarchy_reference_level_analysis")
tissues <- c("Salivary_Gland", "Thymus")
levels <- c("coarse", "leaf")
n_mixtures <- 30L

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

score_prediction <- function(pred, truth) {
  pred <- normalize_nonnegative(pred)
  truth <- normalize_nonnegative(truth)
  diff <- pred - truth
  l1 <- sum(abs(diff))
  c(L1 = l1, Accuracy = 1 - l1 / 2, MAE = mean(abs(diff)), RMSE = sqrt(mean(diff^2)))
}

strip_level_prefix <- function(x) {
  sub("^L[0-9]+: ", "", x)
}

clean_node_name <- function(x) {
  gsub(" / ", "__", strip_level_prefix(x), fixed = TRUE)
}

leaf_names <- function(node) {
  children <- node$children
  if (is.null(children) || length(children) == 0) return(node$name)
  unlist(lapply(children, leaf_names), use.names = FALSE)
}

nodes_at_depth <- function(root, depth) {
  current <- root$children
  if (depth == 1L) return(current)
  for (i in seq_len(depth - 1L)) {
    next_nodes <- list()
    for (node in current) {
      children <- node$children
      if (is.null(children) || length(children) == 0) {
        next_nodes[[length(next_nodes) + 1L]] <- node
      } else {
        next_nodes <- c(next_nodes, children)
      }
    }
    current <- next_nodes
  }
  current
}

hierarchy_mappings <- function(tissue, fine_cell_types) {
  tree <- fromJSON(file.path(base_dir, "hierarchies", paste0(tissue, "_hierarchy.json")),
                   simplifyVector = FALSE)
  mappings <- list()
  mapping <- list()
  used <- character()
  for (node in tree$children) {
    name <- clean_node_name(node$name)
    if (name %in% used) name <- paste0(name, "_", length(used) + 1L)
    used <- c(used, name)
    leaves <- intersect(leaf_names(node), fine_cell_types)
    if (length(leaves) > 0) mapping[[name]] <- leaves
  }
  mappings[["coarse"]] <- mapping
  mappings[["leaf"]] <- as.list(fine_cell_types)
  names(mappings[["leaf"]]) <- fine_cell_types
  mappings
}

group_labels <- function(labels, mapping) {
  label_to_group <- unlist(lapply(names(mapping), function(group) {
    stats::setNames(rep(group, length(mapping[[group]])), mapping[[group]])
  }))
  unname(label_to_group[labels])
}

aggregate_truth <- function(truth, mapping) {
  mat <- sapply(mapping, function(leaves) colSums(truth[leaves, , drop = FALSE]))
  t(mat)
}

aggregate_leaf_prediction <- function(leaf_pred, mapping) {
  mat <- sapply(mapping, function(leaves) rowSums(leaf_pred[, leaves, drop = FALSE]))
  as.data.frame(mat, check.names = FALSE)
}

run_nnls <- function(signature, bulk) {
  common <- intersect(rownames(signature), rownames(bulk))
  sig <- as.matrix(signature[common, , drop = FALSE])
  b <- as.matrix(bulk[common, , drop = FALSE])
  out <- matrix(NA_real_, nrow = ncol(b), ncol = ncol(sig),
                dimnames = list(colnames(b), colnames(sig)))
  for (i in seq_len(ncol(b))) {
    fit <- nnls::nnls(sig, b[, i])
    out[i, ] <- normalize_nonnegative(coef(fit))
  }
  out
}

solve_qp_nonnegative <- function(S, B) {
  D <- crossprod(S) + diag(1e-8, ncol(S))
  d <- crossprod(S, B)
  A <- diag(ncol(S))
  sol <- tryCatch(
    solve.QP(D, d, A, rep(0, ncol(S)))$solution,
    error = function(e) rep(0, ncol(S))
  )
  normalize_nonnegative(sol)
}

solve_dwls_one <- function(S, B) {
  solution <- solve_qp_nonnegative(S, B)
  denom <- as.vector(S %*% solution)
  denom[denom <= 0 | !is.finite(denom)] <- min(denom[denom > 0 & is.finite(denom)], na.rm = TRUE)
  ws <- (1 / denom)^2
  ws <- ws / min(ws[is.finite(ws)])
  ws[!is.finite(ws)] <- max(ws[is.finite(ws)])
  cap <- stats::quantile(ws, 0.95, na.rm = TRUE)
  ws <- pmin(ws, cap)
  W <- diag(as.numeric(ws))
  D <- t(S) %*% W %*% S + diag(1e-8, ncol(S))
  d <- t(S) %*% W %*% B
  sol <- tryCatch(
    solve.QP(D / norm(D, "2"), d / norm(D, "2"), diag(ncol(S)), rep(0, ncol(S)))$solution,
    error = function(e) solution
  )
  normalize_nonnegative(sol)
}

run_dwls <- function(signature, bulk) {
  common <- intersect(rownames(signature), rownames(bulk))
  sig <- as.matrix(signature[common, , drop = FALSE])
  b <- as.matrix(bulk[common, , drop = FALSE])
  out <- matrix(NA_real_, nrow = ncol(b), ncol = ncol(sig),
                dimnames = list(colnames(b), colnames(sig)))
  for (i in seq_len(ncol(b))) {
    out[i, ] <- solve_dwls_one(sig, b[, i])
  }
  out
}

run_music <- function(sc, labels, bulk) {
  common <- intersect(rownames(sc), rownames(bulk))
  sc <- as.matrix(sc[common, , drop = FALSE])
  bulk <- as.matrix(bulk[common, , drop = FALSE])
  sce <- SingleCellExperiment(assays = list(counts = sc))
  set.seed(42)
  meta <- data.frame(
    cell_type = labels,
    sampleID = paste0("sample", sample.int(3, ncol(sce), replace = TRUE)),
    row.names = colnames(sce),
    stringsAsFactors = FALSE
  )
  colData(sce) <- S4Vectors::DataFrame(meta)
  est <- music_prop(
    bulk.mtx = bulk,
    sc.sce = sce,
    clusters = "cell_type",
    samples = "sampleID",
    select.ct = unique(labels),
    verbose = FALSE
  )
  normalize_rows(est[["Est.prop.weighted"]])
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

run_hidecon_flat <- function(sc, labels, bulk) {
  common <- intersect(rownames(sc), rownames(bulk))
  sc <- as.matrix(sc[common, , drop = FALSE])
  bulk <- as.matrix(bulk[common, , drop = FALSE])
  type_order <- unique(labels)
  k <- length(type_order)
  B <- list(matrix(1, nrow = 1, ncol = k, dimnames = list("root", type_order)),
            diag(k))
  colnames(B[[2]]) <- type_order
  rownames(B[[2]]) <- type_order
  select_res <- select_HiDecon(
    bulk = bulk,
    ref = sc,
    B = B,
    cell_type = labels,
    type_order = type_order,
    normalized = FALSE
  )
  normalize_rows(select_res$res)
}

run_bayesprism <- function(sc, labels, bulk) {
  common <- intersect(rownames(sc), rownames(bulk))
  sc <- as.matrix(sc[common, , drop = FALSE])
  bulk <- t(as.matrix(bulk[common, , drop = FALSE]))
  keep_genes <- rowSums(sc) > 0 & colSums(bulk) > 0
  sc <- t(sc[keep_genes, , drop = FALSE])
  bulk <- bulk[, keep_genes, drop = FALSE]
  stat <- BayesPrism::get.exp.stat(
    sc.dat = sc,
    cell.type.labels = labels,
    cell.state.labels = labels,
    pseudo.count = 0.1,
    cell.count.cutoff = 1,
    n.cores = 1
  )
  markers <- BayesPrism::select.marker(sc.dat = sc, stat = stat, pval.max = 0.05, lfc.min = 0.1)
  prism <- BayesPrism::new.prism(
    reference = markers,
    mixture = bulk,
    input.type = "count.matrix",
    cell.type.labels = labels,
    cell.state.labels = labels,
    key = NULL,
    outlier.cut = 1,
    outlier.fraction = 0.2
  )
  bp <- BayesPrism::run.prism(prism = prism, n.cores = 1)
  normalize_rows(BayesPrism::get.fraction(bp = bp, which.theta = "final", state.or.type = "type"))
}

write_prediction <- function(pred, tissue, level, method) {
  tissue_out <- file.path(out_dir, tissue)
  dir.create(tissue_out, recursive = TRUE, showWarnings = FALSE)
  write.table(pred, file.path(tissue_out, paste0(tissue, "_", level, "_", method, "_prediction.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
}

score_table <- function(pred, truth, tissue, level, method) {
  rows <- lapply(seq_len(nrow(pred)), function(i) {
    s <- score_prediction(pred[i, colnames(pred)], truth[colnames(pred), rownames(pred)[i]])
    data.frame(tissue = tissue, level = level, method = method,
               mixture = rownames(pred)[i], mixture_index = i,
               n_reference_celltypes = ncol(pred),
               L1 = s[["L1"]], Accuracy = s[["Accuracy"]],
               MAE = s[["MAE"]], RMSE = s[["RMSE"]],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

partition_key <- function(mapping) {
  paste(vapply(mapping, paste, collapse = "|", FUN.VALUE = character(1)), collapse = ";;")
}

all_metrics <- list()
all_consistency <- list()
all_logs <- list()

for (tissue in tissues) {
  message("[LOAD] ", tissue)
  case_dir <- file.path(base_dir, tissue)
  tissue_out <- file.path(out_dir, tissue)
  dir.create(tissue_out, recursive = TRUE, showWarnings = FALSE)
  bulk <- read_matrix(file.path(case_dir, "bulk.txt"))
  bulk <- bulk[, seq_len(min(n_mixtures, ncol(bulk))), drop = FALSE]
  truth_leaf <- read_matrix(file.path(case_dir, "bulkRatio.txt"))
  sc <- read_matrix(file.path(case_dir, "singleCellExpr.txt"))
  labels <- read.table(file.path(case_dir, "singleCellLabels.txt"), header = TRUE, sep = "\t",
                       check.names = FALSE, quote = "", comment.char = "",
                       stringsAsFactors = FALSE)[[1]]
  mappings <- hierarchy_mappings(tissue, rownames(truth_leaf))
  pred_cache <- list()
  map_cache <- list()
  level_pred_cache <- list()

  for (level in levels) {
    mapping <- mappings[[level]]
    key <- partition_key(mapping)
    truth <- aggregate_truth(truth_leaf, mapping)
    ref_path <- file.path(tissue_out, paste0(tissue, "_", level, "_reference_signature.tsv"))
    signature <- read_matrix(ref_path)
    grouped <- group_labels(labels, mapping)
    keep <- !is.na(grouped)
    sc_level <- sc[, keep, drop = FALSE]
    labels_level <- grouped[keep]

    methods <- list(
      NNLS = function() run_nnls(signature, bulk),
      DWLS = function() run_dwls(signature, bulk),
      MuSic = function() run_music(sc_level, labels_level, bulk),
      HiDecon = function() run_hidecon_flat(sc_level, labels_level, bulk),
      BayesPrism = function() run_bayesprism(sc_level, labels_level, bulk)
    )

    for (method in names(methods)) {
      pred_path <- file.path(tissue_out, paste0(tissue, "_", level, "_", method, "_prediction.tsv"))
      start <- Sys.time()
      status <- "OK"
      error <- ""
      pred <- NULL
      tryCatch({
        if (file.exists(pred_path)) {
          pred <- as.matrix(read_matrix(pred_path))
          if (nrow(pred) != ncol(bulk)) stop("existing prediction has wrong row count")
          message("[LOAD] ", tissue, " / ", level, " / ", method)
        } else if (!is.null(pred_cache[[method]]) && !is.null(pred_cache[[method]][[key]])) {
          cached <- pred_cache[[method]][[key]]
          source_mapping <- map_cache[[method]][[key]]
          source_by_leaves <- setNames(names(source_mapping), vapply(source_mapping, paste, collapse = "|", FUN.VALUE = character(1)))
          pred <- sapply(mapping, function(leaves) cached[, source_by_leaves[[paste(leaves, collapse = "|")]]])
          pred <- as.matrix(pred)
          rownames(pred) <- rownames(cached)
          message("[REUSE] ", tissue, " / ", level, " / ", method, " identical partition")
        } else {
          message("[RUN] ", tissue, " / ", level, " / ", method)
          pred <- methods[[method]]()
        }
        pred <- normalize_rows(pred)
        pred <- pred[, names(mapping), drop = FALSE]
        rownames(pred) <- colnames(bulk)
        write_prediction(pred, tissue, level, method)
        all_metrics[[length(all_metrics) + 1L]] <- score_table(pred, truth, tissue, level, method)
        if (is.null(pred_cache[[method]])) pred_cache[[method]] <- list()
        if (is.null(map_cache[[method]])) map_cache[[method]] <- list()
        if (is.null(level_pred_cache[[method]])) level_pred_cache[[method]] <- list()
        pred_cache[[method]][[key]] <- pred
        map_cache[[method]][[key]] <- mapping
        level_pred_cache[[method]][[level]] <- pred
      }, error = function(e) {
        status <<- "ERROR"
        error <<- conditionMessage(e)
        message("[ERROR] ", tissue, " / ", level, " / ", method, ": ", error)
      })
      all_logs[[length(all_logs) + 1L]] <- data.frame(
        tissue = tissue, level = level, method = method, status = status,
        rows = if (is.null(pred)) NA_integer_ else nrow(pred),
        cols = if (is.null(pred)) NA_integer_ else ncol(pred),
        seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 3),
        error = error,
        stringsAsFactors = FALSE
      )
    }
  }

  for (method in names(level_pred_cache)) {
    leaf_pred <- level_pred_cache[[method]][["leaf"]]
    if (is.null(leaf_pred)) next
    for (level in c("coarse")) {
      direct_pred <- level_pred_cache[[method]][[level]]
      if (is.null(direct_pred)) next
      leaf_agg <- aggregate_leaf_prediction(leaf_pred, mappings[[level]])
      rows <- lapply(seq_len(nrow(direct_pred)), function(i) {
        s <- score_prediction(direct_pred[i, ], leaf_agg[rownames(direct_pred)[i], colnames(direct_pred)])
        data.frame(tissue = tissue, level = level, method = method,
                   comparison = "direct_ref_level_vs_leaf_aggregated",
                   mixture = rownames(direct_pred)[i], mixture_index = i,
                   L1 = s[["L1"]], Consistency = s[["Accuracy"]],
                   MAE = s[["MAE"]], RMSE = s[["RMSE"]],
                   stringsAsFactors = FALSE)
      })
      all_consistency[[length(all_consistency) + 1L]] <- do.call(rbind, rows)
    }
  }
}

metrics_df <- do.call(rbind, all_metrics)
consistency_df <- do.call(rbind, all_consistency)
log_df <- do.call(rbind, all_logs)

write.table(metrics_df, file.path(out_dir, "hierarchy_reference_level_other_methods_metrics_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(consistency_df, file.path(out_dir, "hierarchy_reference_level_other_methods_consistency_vs_leaf_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(log_df, file.path(out_dir, "hierarchy_reference_level_other_methods_run_log.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

metric_summary <- aggregate(
  cbind(Accuracy, MAE, RMSE, L1) ~ tissue + level + method + n_reference_celltypes,
  data = metrics_df,
  FUN = function(x) c(mean = mean(x), median = stats::median(x), sd = stats::sd(x))
)
metric_summary <- do.call(data.frame, metric_summary)
write.table(metric_summary, file.path(out_dir, "hierarchy_reference_level_other_methods_metric_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

consistency_summary <- aggregate(
  cbind(Consistency, MAE, RMSE, L1) ~ tissue + level + method,
  data = consistency_df,
  FUN = function(x) c(mean = mean(x), median = stats::median(x), sd = stats::sd(x))
)
consistency_summary <- do.call(data.frame, consistency_summary)
write.table(consistency_summary, file.path(out_dir, "hierarchy_reference_level_other_methods_consistency_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

print(metric_summary)
print(consistency_summary)
