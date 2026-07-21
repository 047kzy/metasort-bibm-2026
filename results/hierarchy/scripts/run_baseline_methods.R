script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_dir <- file.path(project_dir, "data", "inputs")
results_dir <- file.path(project_dir, "results")
tissues <- c("Salivary_Gland", "Thymus")
n_mixtures <- 30L

read_matrix <- function(path) {
  read.table(path, header = TRUE, row.names = 1, sep = "\t",
             check.names = FALSE, quote = "", comment.char = "")
}

read_labels <- function(path) {
  x <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
  as.character(x[[1]])
}

normalize_nonnegative <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  total <- sum(x)
  if (!is.finite(total) || total <= 0) return(rep(0, length(x)))
  x / total
}

normalize_rows <- function(x) {
  x <- as.matrix(x)
  rn <- rownames(x)
  cn <- colnames(x)
  out <- t(apply(x, 1, normalize_nonnegative))
  rownames(out) <- rn
  colnames(out) <- cn
  out
}

require_method_package <- function(package, method) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(method, " requires the R package '", package, "'.")
  }
}

load_case <- function(tissue, single_cell = FALSE) {
  case_dir <- file.path(input_dir, tissue)
  bulk <- read_matrix(file.path(case_dir, "bulk.tsv.gz"))
  bulk <- bulk[, seq_len(min(n_mixtures, ncol(bulk))), drop = FALSE]
  signature <- read_matrix(file.path(case_dir, "signature.tsv.gz"))
  out <- list(case_dir = case_dir, bulk = bulk, signature = signature)
  if (single_cell) {
    out$sc <- read_matrix(file.path(case_dir, "singleCellExpr.tsv.gz"))
    out$labels <- read_labels(file.path(case_dir, "singleCellLabels.tsv.gz"))
    if (length(out$labels) != ncol(out$sc)) {
      stop(tissue, ": label count does not match single-cell matrix columns.")
    }
  }
  out
}

write_result <- function(mat, tissue, method) {
  tissue_dir <- file.path(results_dir, tissue)
  dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(mat, file.path(tissue_dir, paste0(method, ".tsv")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

run_nnls <- function(tissue) {
  require_method_package("nnls", "NNLS")
  dat <- load_case(tissue)
  genes <- intersect(rownames(dat$signature), rownames(dat$bulk))
  sig <- as.matrix(dat$signature[genes, , drop = FALSE])
  bulk <- as.matrix(dat$bulk[genes, , drop = FALSE])
  res <- t(vapply(colnames(bulk), function(mixture) {
    fit <- nnls::nnls(sig, bulk[, mixture])
    normalize_nonnegative(fit$x)
  }, numeric(ncol(sig))))
  rownames(res) <- colnames(bulk)
  colnames(res) <- colnames(sig)
  write_result(res, tissue, "NNLS")
  res
}

run_lindeconseq <- function(tissue) {
  require_method_package("LinDeconSeq", "LinDeconSeq")
  dat <- load_case(tissue)
  genes <- intersect(rownames(dat$bulk), rownames(dat$signature))
  res <- LinDeconSeq::deconSeq(
    as.matrix(dat$bulk[genes, , drop = FALSE]),
    as.matrix(dat$signature[genes, , drop = FALSE]),
    verbose = FALSE
  )
  res <- as.matrix(res)
  if (nrow(res) != ncol(dat$bulk) && ncol(res) == ncol(dat$bulk)) res <- t(res)
  rownames(res) <- colnames(dat$bulk)
  colnames(res) <- colnames(dat$signature)
  res <- normalize_rows(res)
  write_result(res, tissue, "LinDeconSeq")
  res
}

load_dwls_functions <- function() {
  dwls_path <- Sys.getenv("DWLS_R", unset = "")
  if (!nzchar(dwls_path) || !file.exists(dwls_path)) {
    stop(
      "DWLS requires DWLS.R from lxy-simulation-and-deconvolution. ",
      "Set DWLS_R to its absolute path; see README.md."
    )
  }
  lines <- readLines(dwls_path, warn = FALSE)
  if (length(lines) < 428L) stop("DWLS.R is shorter than the pinned source version.")
  eval(parse(text = lines[seq_len(428L)]), envir = .GlobalEnv)
}

run_dwls <- function(tissue) {
  if (!exists("solveDampenedWLS", mode = "function")) load_dwls_functions()
  dat <- load_case(tissue)
  signature <- as.matrix(dat$signature)
  bulk <- as.matrix(dat$bulk)
  res <- t(vapply(colnames(bulk), function(mixture) {
    bulk_single <- setNames(as.numeric(bulk[, mixture]), rownames(bulk))
    trimmed <- trimData(signature, bulk_single)
    normalize_nonnegative(solveDampenedWLS(trimmed$sig, trimmed$bulk))
  }, numeric(ncol(signature))))
  rownames(res) <- colnames(bulk)
  colnames(res) <- colnames(signature)
  write_result(res, tissue, "DWLS")
  res
}

run_hidecon <- function(tissue) {
  require_method_package("HiDecon", "HiDecon")
  dat <- load_case(tissue, single_cell = TRUE)
  bulk <- as.matrix(dat$bulk)
  sc <- as.matrix(dat$sc)
  storage.mode(bulk) <- "numeric"
  storage.mode(sc) <- "numeric"
  type_order <- unique(dat$labels)
  hierarchy <- list(NULL, diag(length(type_order)))
  selected <- HiDecon::select_HiDecon(
    bulk = bulk,
    ref = sc,
    B = hierarchy,
    cell_type = dat$labels,
    type_order = type_order,
    normalized = FALSE
  )
  write_result(selected$res, tissue, "HiDecon")
  selected$res
}

run_music <- function(tissue) {
  for (package in c(
    "MuSiC", "SingleCellExperiment", "SummarizedExperiment", "Biobase", "S4Vectors"
  )) {
    require_method_package(package, "MuSiC")
  }
  dat <- load_case(tissue, single_cell = TRUE)
  bulk <- as.matrix(dat$bulk)
  sc <- as.matrix(dat$sc)
  storage.mode(bulk) <- "numeric"
  storage.mode(sc) <- "numeric"
  sample_data <- data.frame(sample = colnames(bulk), row.names = colnames(bulk))
  bulk_set <- Biobase::ExpressionSet(
    bulk,
    phenoData = methods::new("AnnotatedDataFrame", data = sample_data)
  )
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = sc))
  set.seed(42)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(data.frame(
    cell_type = dat$labels,
    sampleID = paste0("sample", sample.int(3, ncol(sce), replace = TRUE)),
    row.names = colnames(sce),
    stringsAsFactors = FALSE
  ))
  marker_df <- read.delim(
    file.path(dat$case_dir, "markers.tsv.gz"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  markers <- if (all(c("cell_type", "gene") %in% names(marker_df))) {
    split(marker_df$gene, marker_df$cell_type)
  } else {
    unique(marker_df[[1]])
  }
  estimate <- MuSiC::music_prop(
    bulk.mtx = Biobase::exprs(bulk_set),
    sc.sce = sce,
    markers = markers,
    clusters = "cell_type",
    samples = "sampleID",
    select.ct = unique(dat$labels)
  )
  res <- estimate[["Est.prop.weighted"]]
  write_result(res, tissue, "MuSic")
  res
}

run_bayesprism <- function(tissue) {
  require_method_package("BayesPrism", "BayesPrism")
  dat <- load_case(tissue, single_cell = TRUE)
  bulk <- t(as.matrix(dat$bulk))
  sc <- t(as.matrix(dat$sc))
  storage.mode(bulk) <- "numeric"
  storage.mode(sc) <- "numeric"
  n_cores <- as.integer(Sys.getenv("BAYESPRISM_CORES", unset = "4"))
  stat <- BayesPrism::get.exp.stat(
    sc.dat = sc,
    cell.type.labels = dat$labels,
    cell.state.labels = dat$labels,
    cell.count.cutoff = 1,
    n.cores = n_cores
  )
  markers <- BayesPrism::select.marker(
    sc.dat = sc,
    stat = stat,
    pval.max = 0.05,
    lfc.min = 0.1
  )
  prism <- BayesPrism::new.prism(
    reference = markers,
    mixture = bulk,
    input.type = "count.matrix",
    cell.type.labels = dat$labels,
    cell.state.labels = dat$labels,
    key = NULL,
    outlier.cut = 1,
    outlier.fraction = 0.2
  )
  fit <- BayesPrism::run.prism(prism = prism, n.cores = n_cores)
  res <- BayesPrism::get.fraction(fit, which.theta = "final", state.or.type = "type")
  write_result(res, tissue, "bayesprism")
  res
}

methods <- list(
  NNLS = run_nnls,
  DWLS = run_dwls,
  LinDeconSeq = run_lindeconseq,
  BayesPrism = run_bayesprism,
  HiDecon = run_hidecon,
  MuSic = run_music
)
requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) requested <- names(methods)
unknown <- setdiff(requested, names(methods))
if (length(unknown) > 0L) stop("Unknown methods: ", paste(unknown, collapse = ", "))

log_rows <- list()
for (tissue in tissues) {
  for (method in requested) {
    message("[RUN] ", tissue, " / ", method)
    start <- Sys.time()
    status <- "OK"
    error <- ""
    out_dim <- c(NA_integer_, NA_integer_)
    tryCatch({
      res <- methods[[method]](tissue)
      out_dim <- dim(res)
    }, error = function(e) {
      status <<- "ERROR"
      error <<- conditionMessage(e)
      message("[ERROR] ", tissue, " / ", method, ": ", error)
    })
    log_rows[[length(log_rows) + 1L]] <- data.frame(
      tissue = tissue,
      method = method,
      status = status,
      rows = out_dim[1],
      cols = out_dim[2],
      seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 3),
      error = error,
      stringsAsFactors = FALSE
    )
  }
}

log_df <- do.call(rbind, log_rows)
write.table(log_df, file.path(results_dir, "baseline_methods_run_log.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(log_df)
