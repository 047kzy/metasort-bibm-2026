configure_music_libs <- function() {
    conda_prefix <- Sys.getenv("CONDA_PREFIX")
    if (!nzchar(conda_prefix)) {
        return(invisible(NULL))
    }

    conda_r_lib <- file.path(conda_prefix, "lib", "R", "library")
    if (dir.exists(conda_r_lib)) {
        Sys.setenv(R_LIBS_USER = conda_r_lib)
        .libPaths(unique(c(conda_r_lib, .Library)))
    }
}

load_required_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(
            "Missing required R package '", pkg, "'. ",
            "Run this script with the rebuilt conda environment, e.g. ",
            "`conda run -n music Rscript <script>`.",
            call. = FALSE
        )
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

configure_music_libs()
invisible(lapply(
    c("MuSiC", "SingleCellExperiment", "Biobase", "S4Vectors"),
    load_required_package
))

root_path <- Sys.getenv(
    "ROOT_PATH",
    unset = ""
)
if (!nzchar(root_path)) stop("Set ROOT_PATH to a dataset directory containing tissue folders.")

# 初始化计数器
total_exp <- 0
processed_count <- 0   
skipped_count <- 0

# 获取该路径下的实验文件夹
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
  result_file <- file.path(folder, "result", "MuSic.txt")
  
  # 如果结果文件已存在则跳过
  if (file.exists(result_file)) {
    message(paste("[Skipped]", folder, "- MuSic结果已存在"))
    skipped_count <- skipped_count + 1
    next  # 跳过当前循环
  }
  
  # 如果不存在则继续处理
  message(paste("\n[Processing]", folder))
  
  # 数据加载和处理流程
  tryCatch({
    message("---- DEBUG START ----")
    message("[folder] ", folder)

    # ========== 1) bulk (TSV, gene x mixture) ==========
    bulk_file <- file.path(folder, "bulk.txt")
    message("[bulk_file] ", bulk_file)

    bulk.df <- read.delim(bulk_file, sep = "\t", header = TRUE,
                            row.names = 1, check.names = FALSE)
    message("[bulk.df] dim(gene x sample) = ", paste(dim(bulk.df), collapse=" x "))
    message("[bulk.df] rownames(head) = ", paste(head(rownames(bulk.df)), collapse=", "))
    message("[bulk.df] colnames(head) = ", paste(head(colnames(bulk.df)), collapse=", "))

    # 已知格式：gene x mixture
    bulk.mtx <- as.matrix(bulk.df)
    storage.mode(bulk.mtx) <- "numeric"
    message("[bulk.mtx] dim(gene x sample) = ", paste(dim(bulk.mtx), collapse=" x "))
    message("[bulk.mtx] colnames(samples head) = ", paste(head(colnames(bulk.mtx)), collapse=", "))
    message("[bulk.mtx] rownames(genes head) = ", paste(head(rownames(bulk.mtx)), collapse=", "))

    sample.label <- data.frame(sample = colnames(bulk.mtx), row.names = colnames(bulk.mtx))
    sl <- new("AnnotatedDataFrame", data = sample.label)
    bk.set <- ExpressionSet(bulk.mtx, phenoData = sl)

    # ========== 2) single cell counts (TSV, gene x cell) ==========
    sc_dat_file <- file.path(folder, "singleCellExpr.txt")
    message("[sc_dat_file] ", sc_dat_file)

    sc.df <- read.delim(sc_dat_file, sep = "\t", header = TRUE,
                        row.names = 1, check.names = FALSE)
    message("[sc.df] dim(gene x cell) = ", paste(dim(sc.df), collapse=" x "))

    sc.counts <- as.matrix(sc.df)  # gene x cell
    storage.mode(sc.counts) <- "numeric"
    message("[sc.counts] dim(gene x cell) = ", paste(dim(sc.counts), collapse=" x "))

    sce <- SingleCellExperiment(assays = list(counts = sc.counts))
    message("[sce] nrow(gene) = ", nrow(sce), " ; ncol(cell) = ", ncol(sce))

    # ========== 3) labels (single column, order matches sc columns) ==========
    sc_ct_file <- file.path(folder, "singleCellLabels.txt")
    message("[sc_ct_file] ", sc_ct_file)

    lab.df <- read.delim(sc_ct_file, sep = "\t", header = TRUE,
                        check.names = FALSE, stringsAsFactors = FALSE)
    message("[lab.df] dim = ", paste(dim(lab.df), collapse=" x "))
    message("[lab.df] colnames = ", paste(colnames(lab.df), collapse=", "))
    print(head(lab.df, 5))

    # labels 文件为单列，且与 singleCellExpr 列顺序一致
    if (nrow(lab.df) != ncol(sce)) {
        stop("labels 行数与 singleCellExpr 列数不一致：请检查数据是否配套。")
    }
    cell_type_vec <- lab.df[[1]]

    # ========== 4) subjects (random) ==========
    set.seed(42)  # 保证可复现；想每次不同就删掉这行
    # n_pseudo_samples <- max(2, min(5, ncol(bulk.mtx)))  # 伪样本数：至少2个，最多10个（可按需调整）
    sample_id_vec <- paste0("sample", sample.int(3, ncol(sce), replace = TRUE))

    # ========== 5) markers (TSV) ==========
    marker_file <- file.path(folder, "markers.txt")
    markers <- NULL
    if (file.exists(marker_file)) {
        marker.df <- read.delim(marker_file, sep = "\t", header = TRUE,
                                check.names = FALSE, stringsAsFactors = FALSE)
        message("[marker.df] dim = ", paste(dim(marker.df), collapse=" x "))
        message("[marker.df] colnames = ", paste(colnames(marker.df), collapse=", "))
        if (all(c("cell_type", "gene") %in% colnames(marker.df))) {
            markers <- split(marker.df$gene, marker.df$cell_type)
        } else {
            # 兼容只有一列基因的情况
            markers <- unique(marker.df[[1]])
        }
        message("[markers] loaded genes = ", length(unique(unlist(markers))))
    } else {
        message("[markers] not found, use all genes")
    }

    meta <- data.frame(
        cell_type = cell_type_vec,
        sampleID = sample_id_vec,
        row.names = colnames(sce),
        stringsAsFactors = FALSE
    )

    message("[meta] nrow = ", nrow(meta), " ; sce ncol = ", ncol(sce))
    message("[cell_type] table(top):")
    print(sort(table(meta$cell_type), decreasing = TRUE)[1:min(10, length(unique(meta$cell_type)))])

    colData(sce) <- S4Vectors::DataFrame(meta)

    # ========== 6) MuSiC ==========
    message("[MuSiC] running music_prop ...")
    Est.prop.test <- music_prop(
        bulk.mtx = exprs(bk.set),   # gene x bulk
        sc.sce = sce,              # gene x cell
        markers = markers,
        clusters = "cell_type",
        samples = "sampleID",
        select.ct = unique(meta$cell_type)
    )

    message("[MuSiC] names = ", paste(names(Est.prop.test), collapse=", "))
    result_music <- Est.prop.test[["Est.prop.weighted"]]
    message("[result_music] dim = ", paste(dim(result_music), collapse=" x "))
    print(head(result_music, 3))

    # ========== 7) save ==========
    result_folder <- file.path(folder, "result")
    if (!dir.exists(result_folder)) dir.create(result_folder, recursive = TRUE)
    result_file <- file.path(result_folder, "MuSic.txt")
    write.table(result_music, result_file, sep = "\t", row.names = TRUE, col.names = TRUE)

    processed_count <- processed_count + 1
    message("[Completed] ", folder)
    message("---- DEBUG END ----")

    }, error = function(e) {
    message("处理失败: ", folder)
    message("错误信息: ", e$message)
    })
}
# 输出最终统计
message("\n======== 运行统计 ========")
message(paste("总实验数量:", total_exp))
message(paste("成功处理:", processed_count))
message(paste("跳过已存在:", skipped_count))
message(paste("失败数量:", total_exp - processed_count - skipped_count))
