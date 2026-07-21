conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
if (dir.exists(conda_lib)) .libPaths(conda_lib)
suppressPackageStartupMessages(library(LinDeconSeq))
suppressPackageStartupMessages(library(dplyr))


root_dir <- Sys.getenv("ROOT_PATH", unset = "")
if (!nzchar(root_dir)) stop("Set ROOT_PATH to a dataset directory containing tissue folders.")
required_files <- c("bulk.txt", "signature.txt", "bulkRatio.txt")
tissue_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
tissue_dirs <- tissue_dirs[sapply(tissue_dirs, function(d) all(file.exists(file.path(d, required_files))))]
tissue_names <- basename(tissue_dirs)
cat("发现可运行组织数:", length(tissue_dirs), "\n")
print(tissue_names)

run_one_tissue <- function(tissue_dir, verbose = FALSE) {
  tissue <- basename(tissue_dir)
  bulk <- read.table(file.path(tissue_dir, "bulk.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  sig <- read.table(file.path(tissue_dir, "signature.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

  common_genes <- intersect(rownames(bulk), rownames(sig))
  bulk_mat <- as.matrix(bulk[common_genes, , drop = FALSE])
  sig_mat <- as.matrix(sig[common_genes, , drop = FALSE])

  fractions <- deconSeq(bulk_mat, sig_mat, verbose = verbose)

  out_dir <- file.path(tissue_dir, "result")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- file.path(out_dir, "LinDeconSeq.txt")
  write.table(fractions, file = out_file, sep = "\t", quote = FALSE, col.names = NA)

  data.frame(
    Tissue = tissue,
    GenesUsed = nrow(sig_mat),
    Samples = nrow(fractions),
    CellTypes = ncol(fractions),
    OutputFile = out_file,
    Status = "OK",
    Error = "",
    stringsAsFactors = FALSE
  )
}

run_results <- lapply(tissue_dirs, function(td) {
  tissue <- basename(td)
  out_file <- file.path(td, "result", "LinDeconSeq.txt")

  # 参考 MuSic：结果已存在则跳过
  if (file.exists(out_file)) {
    cat("\n[SKIP]", tissue, "- LinDeconSeq结果已存在\n")
    return(data.frame(
      Tissue = tissue,
      GenesUsed = NA_integer_,
      Samples = NA_integer_,
      CellTypes = NA_integer_,
      OutputFile = out_file,
      Status = "SKIPPED",
      Error = "",
      stringsAsFactors = FALSE
    ))
  }

  cat("\n[RUN]", tissue, "\n")
  tryCatch(
    run_one_tissue(td, verbose = FALSE),
    error = function(e) data.frame(
      Tissue = tissue,
      GenesUsed = NA_integer_,
      Samples = NA_integer_,
      CellTypes = NA_integer_,
      OutputFile = "",
      Status = "ERROR",
      Error = as.character(e),
      stringsAsFactors = FALSE
    )
  )
})
manifest <- do.call(rbind, run_results)
print(manifest)
cat("运行完成。", "\n")

bench_out_dir <- file.path(root_dir, "result")
if (!dir.exists(bench_out_dir)) dir.create(bench_out_dir, recursive = TRUE)
manifest_file <- file.path(bench_out_dir, "LinDeconSeq.txt")
write.table(manifest, file = manifest_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("已保存运行清单到:", manifest_file, "\n")

sessionInfo()
