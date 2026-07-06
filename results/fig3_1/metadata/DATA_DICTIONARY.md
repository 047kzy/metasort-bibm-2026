# Data dictionary

## `results/version_2/all_methods_version2_metrics.tsv`

Per-mixture evaluation metrics for each tissue and method.

| Column | Meaning |
|---|---|
| `tissue` | Tissue dataset: `Salivary_Gland` or `Thymus`. |
| `method` | Deconvolution method name. |
| `mixture` | Mixture identifier. |
| `L1` | L1 distance between predicted and true cell-type proportions. |
| `Accuracy` | Accuracy score, calculated as `1 - L1 / 2`. |
| `MAE` | Mean absolute error across cell types. |
| `RMSE` | Root mean squared error across cell types. |

## `results/version_2/all_methods_version2_summary.tsv`

Method-level summary statistics grouped by tissue.

| Column pattern | Meaning |
|---|---|
| `mean_*` | Mean metric value across mixtures. |
| `median_*` | Median metric value across mixtures. |
| `sd_*` | Standard deviation across mixtures. |

## Per-tissue prediction files

Located under:

- `results/version_2/Salivary_Gland/result/`
- `results/version_2/Thymus/result/`

Files ending in `.txt` contain predicted cell-type proportions. Rows correspond
to mixtures and columns correspond to cell types or hierarchy nodes.

## Gene-weight files

Files ending in `_gene_weights.tsv.gz` contain per-gene MetaSort weight outputs.
Files ending in `_weight_summary.csv` contain summarized gene-weight statistics.

## Processed input files

| File | Meaning |
|---|---|
| `bulk.txt` | Bulk mixture expression matrix. |
| `bulkRatio.txt` | Ground-truth mixture cell-type proportions. |
| `cellTypeExpr.txt` | Cell-type-level expression reference. |
| `signature.txt` | Signature matrix used for deconvolution. |
| `markers.txt` | Marker-gene table. |
| `sigGenes.txt` | Signature gene list. |
| `nCellTypes.txt` | Number of retained cell types. |

## Figures

PDF files are vector outputs for publication assembly. PNG files are raster
previews or high-resolution bitmap exports.
