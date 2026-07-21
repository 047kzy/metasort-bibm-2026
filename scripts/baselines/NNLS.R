conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
if (dir.exists(conda_lib)) .libPaths(conda_lib)
suppressPackageStartupMessages(library(nnls))

root_path <- Sys.getenv(
  "ROOT_PATH",
  unset = ""
)
if (!nzchar(root_path)) stop("Set ROOT_PATH to a dataset directory containing tissue folders.")

# 初始化计数器
total_exp <- 0
processed_count <- 0
skipped_count <- 0
failed_count <- 0

# 获取该路径下的所有文件夹
folders <- list.dirs(root_path, full.names = TRUE, recursive = FALSE)
total_exp <- length(folders)

read_gene_matrix <- function(file_path, label) {
  if (!file.exists(file_path)) {
    stop(label, " 文件不存在: ", file_path)
  }

  df <- read.delim(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  message("[", label, "] dim(gene x sample/cell) = ", paste(dim(df), collapse = " x "))

  mtx <- as.matrix(df)
  storage.mode(mtx) <- "numeric"
  message("[", label, "] rownames(head) = ", paste(head(rownames(mtx)), collapse = ", "))
  message("[", label, "] colnames(head) = ", paste(head(colnames(mtx)), collapse = ", "))

  mtx
}

build_signature_matrix <- function(sc_mat, cell_types, markers = NULL) {
  if (length(cell_types) != ncol(sc_mat)) {
    stop("单细胞表达矩阵列数与标签长度不一致。")
  }

  valid_ct <- !is.na(cell_types) & nzchar(as.character(cell_types))
  if (!all(valid_ct)) {
    sc_mat <- sc_mat[, valid_ct, drop = FALSE]
    cell_types <- cell_types[valid_ct]
  }

  genes <- rownames(sc_mat)
  if (is.null(genes)) {
    stop("singleCellExpr 需要包含基因行名。")
  }

  if (!is.null(markers)) {
    marker_genes <- unique(unlist(markers, use.names = FALSE))
    marker_genes <- marker_genes[!is.na(marker_genes) & nzchar(marker_genes)]
    genes <- intersect(genes, marker_genes)
    message("[markers] after intersection with sc genes = ", length(genes))
  }

  if (length(genes) < 2) {
    stop("可用于构建签名矩阵的基因过少。")
  }

  ct_levels <- unique(as.character(cell_types))
  signature <- sapply(ct_levels, function(ct) {
    idx <- which(as.character(cell_types) == ct)
    rowMeans(sc_mat[genes, idx, drop = FALSE])
  })

  signature <- as.matrix(signature)
  rownames(signature) <- genes
  colnames(signature) <- ct_levels

  signature
}

solve_nnls_bulk <- function(signature_mat, bulk_vec) {
  common_genes <- intersect(rownames(signature_mat), names(bulk_vec))
  if (length(common_genes) < 2) {
    stop("bulk 与 signature 的交集基因过少，无法进行 NNLS。")
  }

  sig <- signature_mat[common_genes, , drop = FALSE]
  bulk_sub <- as.numeric(bulk_vec[common_genes])

  fit <- nnls::nnls(as.matrix(sig), bulk_sub)
  coef <- fit$x
  names(coef) <- colnames(sig)
  coef[coef < 0] <- 0

  total <- sum(coef)
  if (!is.finite(total) || total <= 0) {
    stop("NNLS 解的总和无效，可能是输入矩阵退化或全零。")
  }

  coef / total
}

# 遍历每个文件夹
for (folder in folders) {
  result_folder <- file.path(folder, "result")
  result_file <- file.path(result_folder, "NNLS.txt")

  if (file.exists(result_file)) {
    message(paste("[Skipped]", folder, "- NNLS结果已存在"))
    skipped_count <- skipped_count + 1
    next
  }

  message(paste("\n[Processing]", folder))

  tryCatch({
    message("---- DEBUG START ----")
    message("[folder] ", folder)

    # ========== 1) bulk (TSV, gene x mixture) ==========
    bulk_file <- file.path(folder, "bulk.txt")
    message("[bulk_file] ", bulk_file)
    bulk.df <- read.delim(bulk_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    message("[bulk.df] dim(gene x sample) = ", paste(dim(bulk.df), collapse = " x "))
    message("[bulk.df] rownames(head) = ", paste(head(rownames(bulk.df)), collapse = ", "))
    message("[bulk.df] colnames(head) = ", paste(head(colnames(bulk.df)), collapse = ", "))

    bulk.mtx <- as.matrix(bulk.df)
    storage.mode(bulk.mtx) <- "numeric"
    message("[bulk.mtx] dim(gene x sample) = ", paste(dim(bulk.mtx), collapse = " x "))

    # ========== 2) single cell counts (TSV, gene x cell) ==========
    sc_dat_file <- file.path(folder, "singleCellExpr.txt")
    message("[sc_dat_file] ", sc_dat_file)
    sc.counts <- read_gene_matrix(sc_dat_file, "sc.counts")

    # ========== 3) labels (single column, order matches sc columns) ==========
    sc_ct_file <- file.path(folder, "singleCellLabels.txt")
    message("[sc_ct_file] ", sc_ct_file)

    lab.df <- read.delim(sc_ct_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    message("[lab.df] dim = ", paste(dim(lab.df), collapse = " x "))
    message("[lab.df] colnames = ", paste(colnames(lab.df), collapse = ", "))
    print(head(lab.df, 5))

    if (nrow(lab.df) != ncol(sc.counts)) {
      stop("labels 行数与 singleCellExpr 列数不一致：请检查数据是否配套。")
    }
    cell_type_vec <- lab.df[[1]]

    # ========== 4) markers (TSV) ==========
    marker_file <- file.path(folder, "markers.txt")
    markers <- NULL
    if (file.exists(marker_file)) {
      marker.df <- read.delim(marker_file, sep = "\t", header = TRUE,
                              check.names = FALSE, stringsAsFactors = FALSE)
      message("[marker.df] dim = ", paste(dim(marker.df), collapse = " x "))
      message("[marker.df] colnames = ", paste(colnames(marker.df), collapse = ", "))

      if (all(c("cell_type", "gene") %in% colnames(marker.df))) {
        markers <- split(marker.df$gene, marker.df$cell_type)
      } else {
        markers <- unique(marker.df[[1]])
      }
      message("[markers] loaded genes = ", length(unique(unlist(markers, use.names = FALSE))))
    } else {
      message("[markers] not found, use all genes")
    }

    # ========== 5) build signature matrix ==========
    message("[signature] building mean expression per cell type ...")
    signature_mat <- build_signature_matrix(sc.counts, cell_type_vec, markers)
    message("[signature] dim(gene x cell_type) = ", paste(dim(signature_mat), collapse = " x "))
    message("[signature] colnames(head) = ", paste(head(colnames(signature_mat)), collapse = ", "))

    common_genes <- intersect(rownames(signature_mat), rownames(bulk.mtx))
    message("[common_genes] signature ∩ bulk = ", length(common_genes))
    if (length(common_genes) < 2) {
      stop("signature 与 bulk 的交集基因过少，无法进行反卷积。")
    }

    # ========== 6) NNLS ==========
    message("[NNLS] running deconvolution sample by sample ...")
    result_mat <- matrix(NA_real_, nrow = ncol(bulk.mtx), ncol = ncol(signature_mat))
    rownames(result_mat) <- colnames(bulk.mtx)
    colnames(result_mat) <- colnames(signature_mat)

    for (sample_name in colnames(bulk.mtx)) {
      message("[NNLS] sample = ", sample_name)
      coef <- solve_nnls_bulk(signature_mat, bulk.mtx[, sample_name])
      result_mat[sample_name, names(coef)] <- coef
      message("[NNLS] sample sum = ", signif(sum(coef), 6))
    }

    message("[result_mat] dim = ", paste(dim(result_mat), collapse = " x "))
    print(head(result_mat, 3))

    # ========== 7) save ==========
    if (!dir.exists(result_folder)) dir.create(result_folder, recursive = TRUE)
    write.table(result_mat, result_file, sep = "\t", row.names = TRUE, col.names = TRUE)

    processed_count <- processed_count + 1
    message("[Completed] ", folder)
    message("---- DEBUG END ----")

  }, error = function(e) {
    failed_count <<- failed_count + 1
    message("---- DEBUG ERROR ----")
    message("处理失败: ", folder)
    message("错误信息: ", e$message)
  })
}

# 输出最终统计
message("\n======== 运行统计 ========")
message(paste("总实验数量:", total_exp))
message(paste("成功处理:", processed_count))
message(paste("跳过已存在:", skipped_count))
message(paste("失败数量:", failed_count))
message(paste("未完成数量:", total_exp - processed_count - skipped_count))
