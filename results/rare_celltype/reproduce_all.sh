#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Rscript "$project_dir/code/plot_rare_celltype_bibm.R"
Rscript "$project_dir/code/plot_selected_marker_mean_weight_vs_proportion_bibm.R"

echo "All retained rare-cell-type figures and result tables were regenerated."
