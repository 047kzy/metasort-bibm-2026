# Data dictionary

## Figure-source tables

### `data/outputs/all_methods_version2_metrics.tsv`

One row per tissue, deconvolution method, and mixture (600 rows total).

| Column | Meaning |
|---|---|
| `tissue` | Dataset: `Salivary_Gland` or `Thymus`. |
| `method` | Internal method identifier used by the scripts. |
| `mixture` | Mixture identifier. |
| `L1` | L1 distance between predicted and true cell-type proportions. |
| `Accuracy` | Accuracy calculated as `1 - L1 / 2`. |
| `MAE` | Mean absolute error across cell types. |
| `RMSE` | Root mean squared error across cell types. |

### `data/outputs/all_methods_version2_summary.tsv`

One row per tissue and method (20 rows total). Columns beginning with `mean_`,
`median_`, and `sd_` give the corresponding statistic across 30 mixtures.

## Method inputs

Each tissue under `data/inputs/` contains:

| File | Rows × columns | Meaning |
|---|---|---|
| `bulk.tsv.gz` | genes × mixtures | Bulk expression matrix; method runners use the first 30 mixtures. |
| `bulkRatio.tsv` | cell types × mixtures | Ground-truth cell-type proportions. |
| `signature.tsv.gz` | genes × cell types | Signature/reference matrix used by NNLS, DWLS, LinDeconSeq, and MetaSort. |
| `singleCellExpr.tsv.gz` | genes × cells | Single-cell reference expression used by BayesPrism, HiDecon, and MuSiC. |
| `singleCellLabels.tsv.gz` | cells × 1 | Cell-type label for each single-cell reference column. |
| `markers.tsv.gz` | marker records | Marker genes used by MuSiC and hierarchy-aware MuSiC. |

Files under `data/hierarchies/` define the cell-type trees used by HiDecon-H,
MuSiC-H, and MetaSort-H.

## Prediction results

Each file under `results/<tissue>/` is a 30-mixture prediction matrix. Rows are
mixtures and columns are cell types. The ten result files correspond to the ten
methods shown in the figures. `scripts/build_figure_metrics.R` intersects and
orders mixture and cell-type names before normalizing prediction rows and
calculating metrics.

## Figure outputs

Files under `plots/` are the manuscript figures. PDF files are vector outputs;
PNG files are 600-dpi raster outputs at 2220 × 1920 pixels.
