# Hierarchy-aware deconvolution benchmark

This module contains the Salivary Gland and Thymus benchmark used to compare
ten baseline and hierarchy-aware deconvolution methods across 30 mixtures per
tissue. Accuracy is defined as `1 - L1 / 2`.

## Contents

- `data/inputs/`: compressed bulk, truth, signature, and single-cell inputs.
- `data/hierarchies/`: cell-type trees used by hierarchy-aware methods.
- `data/outputs/`: retained mixture-level metrics and method summaries.
- `results/`: prediction matrices for the ten compared methods.
- `config/`: retained MetaSort parameter snapshots.
- `plots/`: manuscript PNG and PDF figures.
- `scripts/`: method runners, metric construction, and plotting code.
- `docs/DATA_DICTIONARY.md`: table and matrix definitions.

## Reproduce metrics and figures

From the repository root:

```bash
Rscript results/hierarchy/scripts/build_figure_metrics.R
Rscript results/hierarchy/scripts/plot_accuracy_boxplots.R
```

The figure workflow requires R and `ggplot2`. The committed figures use Times
New Roman and Cairo PDF output; a different installed font may produce small
layout differences.

## Optional method reruns

```bash
Rscript results/hierarchy/scripts/run_baseline_methods.R
Rscript scripts/baselines/run_hierarchy_methods.R
METASORT_JK_PATH=. python results/hierarchy/scripts/run_metasort.py
METASORT_JK_PATH=. python results/hierarchy/scripts/run_metasort_hierarchy.py
```

The R runners require the corresponding external method packages. Set
`DWLS_R=scripts/baselines/DWLS.R` for the shared DWLS implementation. The other
historical non-MetaSort runners are also archived under `scripts/baselines/`.
The hierarchy-aware HiDecon and MuSiC calls are in
`scripts/baselines/run_hierarchy_methods.R`.
The MetaSort runners use the root `metasort` package when `METASORT_JK_PATH=.`
is set. Rerun method outputs before rebuilding the metrics and figures.
