conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
if (dir.exists(conda_lib)) .libPaths(conda_lib)
suppressPackageStartupMessages(library(HiDecon))

## 无层次化结构
root_path <- Sys.getenv(
  "ROOT_PATH",
  unset = ""
)
if (!nzchar(root_path)) stop("Set ROOT_PATH to a dataset directory containing tissue folders.")
b_root_path <- Sys.getenv("HIDECON_B_ROOT", unset = "")
overwrite <- identical(tolower(Sys.getenv("HIDECON_OVERWRITE", unset = "0")), "1")
allow_non_level <- identical(tolower(Sys.getenv("HIDECON_ALLOW_NON_LEVEL", unset = "0")), "1")

hidecon_label_aliases <- list(
  Anterior_Cingulate_Cortex = list(
    level2 = c(
      cerebral_cortex_endothelial_cell = "Vascular / mesenchymal lineage",
      corticothalamic_projecting_glutamatergic_cortical_neuron = "Glutamatergic_cortical_neurons",
      l5_extratelencephalic_projecting_glutamatergic_cortical_neuron = "Glutamatergic_cortical_neurons",
      l6b_glutamatergic_cortical_neuron = "Glutamatergic_cortical_neurons",
      near_projecting_glutamatergic_cortical_neuron = "Glutamatergic_cortical_neurons",
      vascular_leptomeningeal_cell = "Vascular / mesenchymal lineage"
    ),
    level3 = c(
      Glutamatergic_interneurons = "Glutamatergic_cortical_neurons",
      Vascular_Mesenchymal_lineage = "Vascular / mesenchymal lineage"
    )
  ),
  Bladder = list(
    level2 = c(Lymphoid_lineage = "Lymphocyte")
  ),
  Lymph_Node = list(
    level2 = c(T_lineage = "T_cell")
  ),
  Skin = list(
    level2 = c(Lymphoid_lineage = "T_cell")
  )
)

apply_hidecon_aliases <- function(labels, tissue_name, level_name) {
  aliases <- hidecon_label_aliases[[tissue_name]][[level_name]]
  if (is.null(aliases)) {
    return(list(labels = labels, reverse = character(0)))
  }

  recoded <- labels
  for (from in names(aliases)) {
    recoded[recoded == from] <- aliases[[from]]
  }

  alias_targets <- unname(aliases)
  one_to_one_targets <- alias_targets[!duplicated(alias_targets) & !duplicated(alias_targets, fromLast = TRUE)]
  reverse <- names(aliases)[match(one_to_one_targets, aliases)]
  names(reverse) <- one_to_one_targets

  list(labels = recoded, reverse = reverse)
}

# 初始化计数器
total_exp <- 0
processed_count <- 0
skipped_count <- 0

# 获取该路径下的所有文件夹
folders <- list.dirs(root_path, full.names = TRUE, recursive = FALSE)
required_files <- c("bulk.txt", "singleCellExpr.txt", "singleCellLabels.txt")
folders <- folders[
  (allow_non_level | grepl("^level[0-9]+$", basename(folders))) &
    vapply(
      folders,
      function(folder) all(file.exists(file.path(folder, required_files))),
      logical(1)
    )
]
total_exp <- length(folders)

# 遍历每个文件夹
for (folder in folders) {
  # 检查结果文件是否存在
  result_file <- file.path(folder, "result", "HiDecon.txt")
  
  # 如果结果文件已存在则跳过
  if (file.exists(result_file) && !overwrite) {
    message(paste("[Skipped]", folder, "HiDecon结果已存在"))
    skipped_count <- skipped_count + 1
    next  # 跳过当前循环
  }
  
  # 如果不存在则继续处理
  message(paste("\n[Processing]", folder))
  
  # 数据加载和处理流程
  tryCatch({
    bulk_file <- file.path(folder, "bulk.txt")
    bk.dat <- read.delim(bulk_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    bulk.matrix <- as.matrix(bk.dat)  # gene x mixture
    storage.mode(bulk.matrix) <- "numeric"

    sc_dat_file <- file.path(folder, "singleCellExpr.txt")
    sc.df <- read.delim(sc_dat_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    sc_df <- as.matrix(sc.df)  # gene x cell
    storage.mode(sc_df) <- "numeric"
   
    sc_ct_file <- file.path(folder, "singleCellLabels.txt")
    lab.df <- read.delim(sc_ct_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(lab.df) != ncol(sc_df)) {
      stop("labels 行数与 singleCellExpr 列数不一致：请检查数据是否配套。")
    }
    cell.type.labels <- lab.df[[1]]
    order.type <- unique(cell.type.labels)

    if (nzchar(b_root_path)) {
      tissue_name <- basename(normalizePath(root_path, mustWork = FALSE))
      level_name <- basename(folder)
      b_file <- file.path(
        b_root_path,
        tissue_name,
        level_name,
        "HiDecon_B",
        paste0(tissue_name, "_", level_name, "_HiDecon_B.rds")
      )
      if (!file.exists(b_file)) {
        stop("HiDecon B文件不存在: ", b_file)
      }
      b_obj <- readRDS(b_file)
      B <- b_obj$B
      order.type <- b_obj$type_order
      alias_result <- list(reverse = character(0))
      if (!setequal(unique(cell.type.labels), order.type)) {
        alias_result <- apply_hidecon_aliases(cell.type.labels, tissue_name, level_name)
        cell.type.labels <- alias_result$labels
      }
      message("[HiDecon B] ", b_file)
    } else {
      num_cell_types <- length(unique(cell.type.labels))
      B <- list(
        NULL,  # 占位符（满足函数参数结构要求）
        diag(num_cell_types)  # 创建单位矩阵（表示独立细胞类型）
      )
    }

    select.res <- select_HiDecon(bulk = bulk.matrix, ref = sc_df, B = B, cell_type = cell.type.labels,
                                  type_order = order.type, normalized = F)
        
    # 创建结果目录
    result_folder <- file.path(folder, "result")
    if (!dir.exists(result_folder)) dir.create(result_folder)
    
    # 保存结果
    if (exists("alias_result") && length(alias_result$reverse) > 0) {
      colnames(select.res$res) <- ifelse(
        colnames(select.res$res) %in% names(alias_result$reverse),
        alias_result$reverse[colnames(select.res$res)],
        colnames(select.res$res)
      )
    }
    write.table(select.res$res, result_file, sep = "\t", row.names=TRUE, col.names=TRUE)
    
    processed_count <- processed_count + 1
    message(paste("[Completed]", folder))
    
  }, error = function(e) {
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
