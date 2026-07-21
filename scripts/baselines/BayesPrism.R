conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
if (dir.exists(conda_lib)) .libPaths(conda_lib)
suppressPackageStartupMessages(library(BayesPrism))

root_path <- Sys.getenv(
  "BAYESPRISM_ROOT_PATH",
  unset = ""
)
if (!nzchar(root_path)) stop("Set BAYESPRISM_ROOT_PATH to a dataset directory containing tissue folders.")

# 初始化计数器
total_exp <- 0
processed_count <- 0
skipped_count <- 0

# 获取该路径下的所有文件夹
folders <- list.dirs(root_path, full.names = TRUE, recursive = FALSE)
required_files <- c("bulk.txt", "singleCellExpr.txt", "singleCellLabels.txt")
folders <- folders[vapply(
  folders,
  function(folder) all(file.exists(file.path(folder, required_files))),
  logical(1)
)]
total_exp <- length(folders)

# 遍历每个文件夹
for (folder in folders) {
  # 检查结果文件是否存在
  result_file <- file.path(folder, "result", "bayesprism.txt")
  
  # 如果结果文件已存在则跳过
  if (file.exists(result_file)) {
    message(paste("[Skipped]", folder, "- BayesPrism结果已存在"))
    skipped_count <- skipped_count + 1
    next  # 跳过当前循环
  }
  
  # 如果不存在则继续处理
  message(paste("\n[Processing]", folder))
  
  # 数据加载和处理流程
  tryCatch({
    message("---- DEBUG START ----")
    message("[folder] ", folder)

    # ========== 1) bulk ==========
    bulk_file <- file.path(folder, "bulk.txt")
    message("[bulk_file] ", bulk_file)

    bk.df <- read.delim(bulk_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    message("[bk.df(raw)] dim = ", paste(dim(bk.df), collapse = " x "))

    # 转置（bulk）
    bk.dat <- t(as.matrix(bk.df))
    storage.mode(bk.dat) <- "numeric"
    message("[bk.dat] dim(after transpose) = ", paste(dim(bk.dat), collapse = " x "))
    message("[bk.dat] rownames(head) = ", paste(head(rownames(bk.dat)), collapse = ", "))
    message("[bk.dat] colnames(head) = ", paste(head(colnames(bk.dat)), collapse = ", "))

    # ========== 2) single cell expr ==========
    sc_dat_file <- file.path(folder, "singleCellExpr.txt")
    message("[sc_dat_file] ", sc_dat_file)

    sc.df <- read.delim(sc_dat_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    message("[sc.df(raw)] dim = ", paste(dim(sc.df), collapse = " x "))

    # 转置（single cell）
    sc.dat <- t(as.matrix(sc.df))
    storage.mode(sc.dat) <- "numeric"
    message("[sc.dat] dim(after transpose) = ", paste(dim(sc.dat), collapse = " x "))
    message("[sc.dat] rownames(head) = ", paste(head(rownames(sc.dat)), collapse = ", "))
    message("[sc.dat] colnames(head) = ", paste(head(colnames(sc.dat)), collapse = ", "))
    
    # ========== 3) labels ==========
    sc_ct_file <- file.path(folder, "singleCellLabels.txt")
    message("[sc_ct_file] ", sc_ct_file)

    # labels 只有一列，第一行是列名；顺序与 sc.dat 的行顺序一致
    ct.df <- read.delim(sc_ct_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    message("[labels(raw)] dim = ", paste(dim(ct.df), collapse = " x "))
    message("[labels(raw)] colnames = ", paste(colnames(ct.df), collapse = ", "))

    if (nrow(ct.df) != nrow(sc.dat)) {
      stop("labels 行数与 singleCellExpr 转置后的行数不一致：请检查数据是否配套。")
    }
    cell.type.labels <- ct.df[[1]]
    message("[labels] n = ", length(cell.type.labels))
    message("[labels] table(top):")
    print(sort(table(cell.type.labels), decreasing = TRUE)[1:min(10, length(unique(cell.type.labels)))])

    # ========== 4) BayesPrism ==========
    message("[BayesPrism] get.exp.stat ...")
    diff.exp.stat <- get.exp.stat(
      sc.dat = sc.dat,
      cell.type.labels = cell.type.labels,
      cell.state.labels = cell.type.labels,
      cell.count.cutoff = 1,
      n.cores = 24
    )
    message("[BayesPrism] get.exp.stat done")

    message("[BayesPrism] select.marker ...")
    sc.dat.filtered.pc.sig <- select.marker(
      sc.dat = sc.dat,
      stat = diff.exp.stat,
      pval.max = 0.05,
      lfc.min = 0.1
    )
    message("[select.marker] filtered dim = ", paste(dim(sc.dat.filtered.pc.sig), collapse = " x "))
    message("[select.marker] marker genes(head) = ", paste(head(rownames(sc.dat.filtered.pc.sig)), collapse = ", "))

    message("[BayesPrism] new.prism ...")
    myPrism <- new.prism(
      reference = sc.dat.filtered.pc.sig,
      mixture = bk.dat,
      input.type = "count.matrix",
      cell.type.labels = cell.type.labels,
      cell.state.labels = cell.type.labels,
      key = NULL,
      outlier.cut = 1,
      outlier.fraction = 0.2
    )
    message("[BayesPrism] new.prism done")

    message("[BayesPrism] run.prism ...")
    bp.res <- run.prism(prism = myPrism, n.cores = 24)
    message("[BayesPrism] run.prism done")

    theta <- get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")
    message("[theta] dim = ", paste(dim(theta), collapse = " x "))
    print(head(theta, 3))

    # ========== 5) save ==========
    result_folder <- file.path(folder, "result")
    if (!dir.exists(result_folder)) dir.create(result_folder, recursive = TRUE)
    message("[save] ", result_file)
    write.table(theta, result_file, sep = "\t", row.names = TRUE, col.names = TRUE)

    processed_count <- processed_count + 1
    message(paste("[Completed]", folder))
    message("---- DEBUG END ----")

    }, error = function(e) {
    message("---- DEBUG ERROR ----")
    message(paste("处理失败:", folder))
    message(paste("错误信息:", e$message))
    })
}

# 输出最终统计
message("\n======== 运行统计 ========")
message(paste("总实验数量:", total_exp))
message(paste("成功处理:", processed_count))
message(paste("跳过已存在:", skipped_count))
message(paste("失败数量:", total_exp - processed_count - skipped_count))
