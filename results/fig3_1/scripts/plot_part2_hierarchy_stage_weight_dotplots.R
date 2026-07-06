suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)][1])
base_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()
analysis_dir <- file.path(base_dir, "metasort_original_weight_de_analysis")
hierarchy_dir <- file.path(base_dir, "hierarchies")
plot_dir <- file.path(analysis_dir, "plots")
top_n <- 10L

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

node_leaves <- function(node) {
  children <- node$children
  if (is.null(children) || length(children) == 0) {
    return(node$name)
  }
  unique(unlist(lapply(children, node_leaves), use.names = FALSE))
}

collect_node_map <- function(node, acc = list()) {
  acc[[node$name]] <- node_leaves(node)
  children <- node$children
  if (!is.null(children) && length(children) > 0) {
    for (child in children) {
      acc <- collect_node_map(child, acc)
    }
  }
  acc
}

clean_node_label <- function(x) {
  x <- sub("^L[0-9]+: ", "", x)
  x <- sub("^([^/]+)/\\s*\\1$", "\\1", x)
  x
}

scale_by_gene <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      return(rep(0, length(x)))
    }
    (x - mean(x, na.rm = TRUE)) / s
  }))
  pmax(pmin(scaled, 2), -2)
}

make_stage_dot_data <- function(expr, labels, genes, node_names, node_map, gene_order) {
  avg_expr <- sapply(node_names, function(node_name) {
    leaves <- node_map[[node_name]]
    if (is.null(leaves)) {
      leaves <- node_name
    }
    cells <- labels %in% leaves
    rowMeans(expr[genes, cells, drop = FALSE], na.rm = TRUE)
  })
  pct_expr <- sapply(node_names, function(node_name) {
    leaves <- node_map[[node_name]]
    if (is.null(leaves)) {
      leaves <- node_name
    }
    cells <- labels %in% leaves
    rowMeans(expr[genes, cells, drop = FALSE] > 0, na.rm = TRUE) * 100
  })
  avg_expr_scaled <- scale_by_gene(log1p(avg_expr))

  dot_data <- expand.grid(
    gene = rownames(avg_expr),
    node = node_names,
    stringsAsFactors = FALSE
  )
  dot_data$node_label <- clean_node_label(dot_data$node)
  dot_data$avg_log1p_expr <- as.vector(log1p(avg_expr))
  dot_data$avg_log1p_expr_scaled <- as.vector(avg_expr_scaled)
  dot_data$pct_expr <- as.vector(pct_expr)
  dot_data$gene <- factor(dot_data$gene, levels = rev(gene_order))
  dot_data$node_label <- factor(dot_data$node_label, levels = clean_node_label(node_names))
  dot_data
}

weights_long <- read.delim(
  gzfile(file.path(analysis_dir, "metasort_gene_weights_long.tsv.gz")),
  check.names = FALSE
)
hier_weights <- weights_long[weights_long$method == "MetaSort_jk_hierarchy", ]
hier_weights$solve_index <- as.integer(hier_weights$solve_index)

top_records <- list()

for (tissue_name in sort(unique(hier_weights$tissue))) {
  expr_path <- file.path(base_dir, tissue_name, "singleCellExpr.txt")
  labels_path <- file.path(base_dir, tissue_name, "singleCellLabels.txt")
  hierarchy_path <- file.path(hierarchy_dir, sprintf("%s_hierarchy.json", tissue_name))

  expr <- read.delim(expr_path, row.names = 1, check.names = FALSE)
  labels <- read.delim(labels_path, check.names = FALSE)[[1]]
  if (length(labels) != ncol(expr)) {
    stop(sprintf("%s labels do not match expression columns", tissue_name))
  }

  hierarchy <- fromJSON(hierarchy_path, simplifyVector = FALSE)
  node_map <- collect_node_map(hierarchy)
  leaf_names <- unique(labels)
  for (leaf in leaf_names) {
    node_map[[leaf]] <- leaf
  }

  tissue_weights <- hier_weights[hier_weights$tissue == tissue_name, ]
  stage_keys <- unique(tissue_weights[, c("solve_index", "visible_nodes")])
  stage_keys <- stage_keys[order(stage_keys$solve_index), ]

  for (stage_idx in seq_len(nrow(stage_keys))) {
    solve_index <- stage_keys$solve_index[stage_idx]
    visible_nodes <- stage_keys$visible_nodes[stage_idx]
    node_names <- strsplit(visible_nodes, "\\|", fixed = FALSE)[[1]]

    stage_weights <- tissue_weights[
      tissue_weights$solve_index == solve_index &
        tissue_weights$visible_nodes == visible_nodes &
        tissue_weights$gene %in% rownames(expr),
    ]
    gene_summary <- aggregate(
      meta_weight ~ tissue + method + solve_index + visible_nodes + gene,
      stage_weights,
      function(x) c(
        mean_meta_weight = mean(x),
        median_meta_weight = median(x),
        max_meta_weight = max(x),
        fraction_weight_gt_baseline = mean(x > 10)
      )
    )
    gene_summary <- do.call(data.frame, gene_summary)
    names(gene_summary) <- sub("^meta_weight\\.", "", names(gene_summary))
    gene_summary <- gene_summary[
      order(
        -gene_summary$mean_meta_weight,
        -gene_summary$fraction_weight_gt_baseline,
        gene_summary$gene
      ),
    ]
    selected <- head(gene_summary, top_n)
    if (nrow(selected) == 0) {
      next
    }
    selected$stage_label <- sprintf(
      "Stage %d: %s",
      solve_index,
      paste(clean_node_label(node_names), collapse = " vs ")
    )
    top_records[[paste(tissue_name, solve_index, sep = "__")]] <- selected

    gene_order <- selected$gene
    dot_data <- make_stage_dot_data(expr, labels, selected$gene, node_names, node_map, gene_order)

    p <- ggplot(dot_data, aes(x = node_label, y = gene)) +
      geom_point(aes(size = pct_expr, color = avg_log1p_expr_scaled), alpha = 0.9) +
      scale_size_area(max_size = 6.2, limits = c(0, 100), name = "Pct expressed") +
      scale_color_gradient2(
        low = "#2c7bb6",
        mid = "grey95",
        high = "#d7191c",
        midpoint = 0,
        limits = c(-2, 2),
        name = "Scaled log1p\nmean expr"
      ) +
      labs(
        title = sprintf(
          "%s hierarchy stage %d: top %d high-weight genes",
          tissue_name,
          solve_index,
          top_n
        ),
        subtitle = paste(clean_node_label(node_names), collapse = " vs "),
        x = NULL,
        y = NULL
      ) +
      theme_bw(base_size = 8) +
      theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6.5),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6)
      )

    width <- max(6.2, 0.8 * length(node_names) + 3.8)
    height <- max(4.6, 0.34 * nrow(selected) + 2.2)
    stem <- file.path(
      plot_dir,
      sprintf("%s__hierarchy_stage%d__top%d_high_weight_gene_dotplot", tissue_name, solve_index, top_n)
    )
    ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 260, limitsize = FALSE)
    ggsave(paste0(stem, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  }
}

top_table <- do.call(rbind, top_records)
top_table <- top_table[order(top_table$tissue, top_table$solve_index, -top_table$mean_meta_weight), ]
write.table(
  top_table,
  file.path(analysis_dir, sprintf("metasort_hierarchy_stage_top%d_high_weight_genes.tsv", top_n)),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(
  aggregate(
    gene ~ tissue + solve_index + visible_nodes,
    top_table,
    function(x) length(unique(x))
  )
)
cat(sprintf("Wrote hierarchy-stage top%d dotplots to %s\n", top_n, plot_dir))
