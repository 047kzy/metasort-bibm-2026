# Figure 3C-F: Rare Cell Types and Marker Weights

This directory contains rare-cell-type recovery and marker-weight trend results
used for Figure 3C-F.

## Contents

- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5/`:
  30-tissue MetaSort and baseline benchmark summaries, including rare-cell-type
  analysis tables.
- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/`:
  selected rare-cell figures and source data for manuscript panels.
- `version_2/metasort_weight_rerun_sqrt_tol0p005_avg2p5/`: rerun outputs used
  to extract MetaSort gene weights and marker-weight trends.
- `plot_rare_celltype_bibm.R`: plotting script for rare-cell recovery panels.
- `plot_selected_marker_mean_weight_vs_proportion_bibm.R`: plotting script for
  marker-weight trend panels.

## Key Manuscript Files

- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5/rare_celltype_analysis/summary.md`
- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5/rare_celltype_analysis/tissue_celltype_method_summary.csv`
- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/CELlxGENE__Basal_Zone_Of_Heart__vein_endothelial_cell_bibm.png`
- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/Tabula-Sapiens__Skin__CD1c_positive_myeloid_dendritic_cell_bibm.png`
- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/selected_marker_mean_weight_vs_proportion_data.csv`
- `version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/skin_marker_mean_weight_vs_proportion_data.csv`
- `version_2/metasort_weight_rerun_sqrt_tol0p005_avg2p5/target_celltype_weight_trends/target_marker_weight_by_mixture.csv`
- `version_2/metasort_weight_rerun_sqrt_tol0p005_avg2p5/target_celltype_weight_trends/target_marker_weight_trend_summary.csv`

## Scripts

```bash
Rscript plot_rare_celltype_bibm.R
Rscript plot_selected_marker_mean_weight_vs_proportion_bibm.R
```
