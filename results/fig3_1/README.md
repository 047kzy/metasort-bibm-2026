# Figure 3A-B: Hierarchy-Aware Deconvolution

This directory contains the Salivary Gland and Thymus hierarchy-aware benchmark
results used for Figure 3A-B in the manuscript.

## Contents

- `data/processed_inputs/`: processed input matrices for Salivary Gland and
  Thymus, including bulk mixtures, ground-truth proportions, signature matrices,
  marker genes, and cell-type counts.
- `data/hierarchies/`: tissue-specific hierarchy JSON files.
- `results/version_2/`: retained result tables, per-tissue predictions, metrics,
  run logs, and MetaSort configurations.
- `figures/version_2/`: figure outputs for the hierarchy benchmark.
- `scripts/`: MetaSort, baseline, and plotting wrappers used to produce the
  retained outputs.
- `metadata/`: manifest, checksums, file-selection notes, and data dictionary.

## Key Manuscript Files

- `results/version_2/all_methods_version2_metrics.tsv`
- `results/version_2/all_methods_version2_summary.tsv`
- `results/version_2/metasort_jk_avg2p5_tol0005_summary.tsv`
- `results/version_2/metasort_jk_hierarchy_avg2p5_tol0005_summary.tsv`
- `results/version_2/metasort_hierarchy_stage_top10_high_weight_genes.tsv`
- `figures/version_2/all_methods_version2_accuracy_boxplot_bibm_salivary_gland.pdf`
- `figures/version_2/all_methods_version2_accuracy_boxplot_bibm_thymus.pdf`
- `figures/version_2/all_methods_version2_accuracy_boxplot_bibm_two_panel.pdf`

## Scripts

```bash
python scripts/run_part2_metasort_jk_avg2p5_tol0005.py
python scripts/run_part2_metasort_jk_hierarchy_avg2p5_tol0005.py
Rscript scripts/plot_part2_version2_accuracy_bibm.R
Rscript scripts/plot_part2_hierarchy_stage_weight_dotplots.R
```

Baseline wrappers are retained for transparency, but external R/Python packages
must be installed separately.
