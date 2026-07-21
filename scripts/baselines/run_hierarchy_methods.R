suppressPackageStartupMessages({
  library(Biobase)
  library(HiDecon)
  library(MuSiC)
  library(SingleCellExperiment)
  library(jsonlite)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
base_dir <- normalizePath(file.path(repo_dir, "results", "hierarchy"), mustWork = TRUE)
input_dir <- file.path(base_dir, "data", "inputs")
result_dir <- file.path(base_dir, "results")
hierarchy_dir <- file.path(base_dir, "data", "hierarchies")
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

strip_level_prefix <- function(x) {
  sub("^L[0-9]+: ", "", x)
}

strip_l2_prefix <- function(x) {
  x <- strip_level_prefix(x)
  sub("^.* / ", "", x)
}

load_hierarchy <- function(path, type_order) {
  tree <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  l1_nodes <- tree$children
  l1_names <- vapply(l1_nodes, function(node) strip_level_prefix(node$name), character(1))

  leaf_parent_l1 <- character()

  get_leaf_names <- function(node) {
    children <- node$children
    if (is.null(children) || length(children) == 0) {
      return(node$name)
    }
    unlist(lapply(children, get_leaf_names), use.names = FALSE)
  }

  for (i in seq_along(l1_nodes)) {
    l1 <- l1_names[[i]]
    leaves <- get_leaf_names(l1_nodes[[i]])
    for (leaf_name in leaves) {
      leaf_parent_l1[[leaf_name]] <- l1
    }
  }

  missing_leaves <- setdiff(type_order, names(leaf_parent_l1))
  extra_leaves <- setdiff(names(leaf_parent_l1), type_order)
  if (length(missing_leaves) > 0 || length(extra_leaves) > 0) {
    stop(
      "Hierarchy leaves do not match observed cell types. Missing: ",
      paste(missing_leaves, collapse = ", "),
      " Extra: ",
      paste(extra_leaves, collapse = ", ")
    )
  }

  b1 <- matrix(1, nrow = 1, ncol = length(l1_names),
               dimnames = list("root", l1_names))
  b2 <- matrix(0, nrow = length(l1_names), ncol = length(type_order),
               dimnames = list(l1_names, type_order))
  for (leaf in type_order) {
    b2[leaf_parent_l1[[leaf]], leaf] <- 1
  }

  clusters_type <- split(type_order, vapply(type_order, function(x) leaf_parent_l1[[x]], character(1)))
  clusters_type <- clusters_type[l1_names]
  clusters_type <- clusters_type[!vapply(clusters_type, is.null, logical(1))]

  list(
    B = list(b1, b2),
    leaf_to_l1 = unlist(leaf_parent_l1),
    clusters_type = clusters_type
  )
}

read_case <- function(tissue) {
  case_dir <- file.path(input_dir, tissue)
  bulk <- read_matrix(file.path(case_dir, "bulk.tsv.gz"))
  keep <- seq_len(min(n_mixtures, ncol(bulk)))
  bulk <- bulk[, keep, drop = FALSE]

  sc <- read_matrix(file.path(case_dir, "singleCellExpr.tsv.gz"))
  labels <- read_labels(file.path(case_dir, "singleCellLabels.tsv.gz"))
  if (length(labels) != ncol(sc)) {
    stop(tissue, ": label count does not match singleCellExpr columns")
  }

  bulk_mat <- as.matrix(bulk)
  sc_mat <- as.matrix(sc)
  storage.mode(bulk_mat) <- "numeric"
  storage.mode(sc_mat) <- "numeric"

  list(
    case_dir = case_dir,
    result_dir = file.path(result_dir, tissue),
    bulk = bulk_mat,
    sc = sc_mat,
    labels = labels
  )
}

write_result <- function(mat, out_dir, method) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(mat, file.path(out_dir, paste0(method, "_hierarchy.tsv")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

run_hidecon_hierarchy <- function(dat, hierarchy) {
  select_res <- select_HiDecon(
    bulk = dat$bulk,
    ref = dat$sc,
    B = hierarchy$B,
    cell_type = dat$labels,
    type_order = unique(dat$labels),
    normalized = FALSE
  )
  write_result(select_res$res, dat$result_dir, "HiDecon")
  select_res$res
}

load_group_markers <- function(case_dir, clusters_type, genes) {
  marker_file <- file.path(case_dir, "markers.tsv.gz")
  if (!file.exists(marker_file)) {
    stop("Missing marker file: ", marker_file)
  }
  marker_df <- read.delim(marker_file, sep = "\t", header = TRUE,
                          check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("cell_type", "gene") %in% colnames(marker_df))) {
    stop("Expected marker file columns: cell_type, gene")
  }
  out <- lapply(clusters_type, function(leaves) {
    unique(intersect(marker_df$gene[marker_df$cell_type %in% leaves], genes))
  })
  empty <- names(out)[vapply(out, length, integer(1)) == 0]
  if (length(empty) > 0) {
    stop("No marker genes found for hierarchy groups: ", paste(empty, collapse = ", "))
  }
  out
}

run_music_hierarchy <- function(dat, hierarchy) {
  sce <- SingleCellExperiment(assays = list(counts = dat$sc))
  set.seed(42)
  cell_groups <- hierarchy$leaf_to_l1[dat$labels]
  if (anyNA(cell_groups)) {
    stop("Some cell labels are missing from hierarchy group map")
  }
  colData(sce) <- S4Vectors::DataFrame(data.frame(
    cell_type = dat$labels,
    hierarchy_group = as.character(cell_groups),
    sampleID = paste0("sample", sample.int(3, ncol(sce), replace = TRUE)),
    row.names = colnames(sce),
    stringsAsFactors = FALSE
  ))

  sample_label <- data.frame(sample = colnames(dat$bulk), row.names = colnames(dat$bulk))
  bulk_set <- ExpressionSet(dat$bulk, phenoData = new("AnnotatedDataFrame", data = sample_label))
  genes <- intersect(rownames(dat$bulk), rownames(dat$sc))
  group_markers <- load_group_markers(dat$case_dir, hierarchy$clusters_type, genes)

  est <- music_prop.cluster(
    bulk.mtx = exprs(bulk_set),
    sc.sce = sce,
    group.markers = group_markers,
    groups = "hierarchy_group",
    clusters = "cell_type",
    samples = "sampleID",
    clusters.type = hierarchy$clusters_type,
    verbose = TRUE
  )
  res <- est[["Est.prop.weighted.cluster"]]
  write_result(res, dat$result_dir, "MuSic")
  res
}

log_rows <- list()
for (tissue in tissues) {
  message("[LOAD] ", tissue)
  dat <- read_case(tissue)
  hierarchy_path <- file.path(hierarchy_dir, paste0(tissue, "_hierarchy.json"))
  hierarchy <- load_hierarchy(hierarchy_path, unique(dat$labels))

  for (method in c("HiDecon", "MuSic")) {
    message("[RUN] ", tissue, " / ", method, "_hierarchy")
    start <- Sys.time()
    status <- "OK"
    error <- ""
    out_dim <- c(NA_integer_, NA_integer_)
    tryCatch({
      res <- if (method == "HiDecon") {
        run_hidecon_hierarchy(dat, hierarchy)
      } else {
        run_music_hierarchy(dat, hierarchy)
      }
      out_dim <- dim(res)
    }, error = function(e) {
      status <<- "ERROR"
      error <<- conditionMessage(e)
      message("[ERROR] ", tissue, " / ", method, "_hierarchy: ", error)
    })
    log_rows[[length(log_rows) + 1L]] <- data.frame(
      tissue = tissue,
      method = paste0(method, "_hierarchy"),
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
write.table(log_df, file.path(result_dir, "hierarchy_methods_run_log.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(log_df)
