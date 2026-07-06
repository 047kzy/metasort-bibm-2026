# Manuscript Results Map

This file maps the Results section of the manuscript to the retained code,
processed data, result tables, and figure source files in this anonymous
artifact.

## Datasets and Setup

Full processed inputs for the 30-tissue benchmark are provided separately from
this anonymous code/results repository. The repository retains the result
tables, figure source data, and selected processed inputs already embedded in
the figure-specific result folders.

The MetaSort implementation and reusable scripts are in:

- `code/metasort/`
- `code/scripts/`
- `code/hierarchies/`

## Figure 2: Tongue Curvature and Iterative Optimization

Directory:

- `results/fig2/`

Manuscript claims covered:

- oracle curvature-guided update in the Tongue mixture;
- comparison with simplex-constrained least squares;
- iterative MetaSort trajectory from initialization to convergence;
- final comparison with MuSiC, DWLS, LinDeconSeq, and NNLS.

Key files:

- `results/fig2/publication_figures/bibm_true_init_metasort_vs_simplex_zscore_l1_summary.csv`
- `results/fig2/publication_figures/bibm_true_init_metasort_vs_simplex_zscore_proportions.csv`
- `results/fig2/publication_figures/bibm_metasort_tongue_case_methods.pdf`
- `results/fig2/publication_figures/bibm_metasort_tongue_case_iteration_0_13.csv`
- `results/fig2/publication_figures/bibm_metasort_tongue_case_combined_revised_values.csv`
- `results/fig2/publication_figures/bibm_metasort_tongue_case_combined_revised_with_true_init_l1.pdf`

Regeneration script:

- `results/fig2/merge_tongue_case_combined_with_true_init_l1.py`

## Figure 3A-B: Hierarchical Cell-Type Composition

Directory:

- `results/fig3_1/`

Manuscript claims covered:

- Salivary Gland benchmark with direct and hierarchical MetaSort;
- Thymus benchmark with direct and hierarchical MetaSort;
- comparison with BayesPrism, DWLS, MuSiC, LinDeconSeq, and other retained
  baselines;
- hierarchy-stage high-weight gene summaries.

Key files:

- `results/fig3_1/results/version_2/all_methods_version2_metrics.tsv`
- `results/fig3_1/results/version_2/all_methods_version2_summary.tsv`
- `results/fig3_1/results/version_2/metasort_jk_avg2p5_tol0005_summary.tsv`
- `results/fig3_1/results/version_2/metasort_jk_hierarchy_avg2p5_tol0005_summary.tsv`
- `results/fig3_1/results/version_2/metasort_hierarchy_stage_top10_high_weight_genes.tsv`
- `results/fig3_1/figures/version_2/all_methods_version2_accuracy_boxplot_bibm_salivary_gland.pdf`
- `results/fig3_1/figures/version_2/all_methods_version2_accuracy_boxplot_bibm_thymus.pdf`
- `results/fig3_1/figures/version_2/all_methods_version2_accuracy_boxplot_bibm_two_panel.pdf`

Regeneration scripts:

- `results/fig3_1/scripts/run_part2_metasort_jk_avg2p5_tol0005.py`
- `results/fig3_1/scripts/run_part2_metasort_jk_hierarchy_avg2p5_tol0005.py`
- `results/fig3_1/scripts/plot_part2_version2_accuracy_bibm.R`
- `results/fig3_1/scripts/plot_part2_hierarchy_stage_weight_dotplots.R`

## Figure 3C-F: Rare Cell Types and Marker Weights

Directory:

- `results/fig3_2/`

Manuscript claims covered:

- rare vein endothelial cell recovery in CELlxGENE Basal Zone of Heart;
- rare CD1c-positive myeloid dendritic cell recovery in Tabula Sapiens Skin;
- marker-weight trends as a function of true cell-type proportion;
- correlation and MAE summaries for selected rare-cell examples.

Key files:

- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5/rare_celltype_analysis/summary.md`
- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5/rare_celltype_analysis/tissue_celltype_method_summary.csv`
- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/CELlxGENE__Basal_Zone_Of_Heart__vein_endothelial_cell_bibm.png`
- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/Tabula-Sapiens__Skin__CD1c_positive_myeloid_dendritic_cell_bibm.png`
- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/selected_marker_mean_weight_vs_proportion_data.csv`
- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/selected_marker_mean_weight_vs_proportion_bibm.pdf`
- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/skin_marker_mean_weight_vs_proportion_data.csv`
- `results/fig3_2/version_2/nar30_all100_part2_nohier_config_sqrt_tol0p005_avg2p5_rare_celltype_plots/skin_marker_mean_weight_vs_proportion_bibm.pdf`

Regeneration scripts:

- `results/fig3_2/plot_rare_celltype_bibm.R`
- `results/fig3_2/plot_selected_marker_mean_weight_vs_proportion_bibm.R`

## Figure 4: Donor-Mismatched Reference Recovery

Directory:

- `results/fig4/`

This directory contains the `bulk_TSP14__reference_TSP10` experiment only, as
used in the manuscript.

Manuscript claims covered:

- pseudo-bulk mixtures from donor TSP14;
- initial reference from donor TSP10;
- recovered bulk-domain reference profiles;
- comparison with the true TSP14 reference and the original TSP10 reference.

Key files:

- `results/fig4/input_reference_donor_reference.tsv`
- `results/fig4/true_bulk_donor_reference.tsv`
- `results/fig4/recovered_bulk_domain_reference.tsv`
- `results/fig4/gene_recovery_metrics.tsv`
- `results/fig4/gene_metrics.tsv`
- `results/fig4/recovered_reference_nnls_iteration_metrics.tsv`
- `results/fig4/initial_metasort_proportions.tsv`
- `results/fig4/anchor_genes.txt`
- `results/fig4/metrics.json`

## Notes on Baselines

Baseline result files are retained where used by the manuscript figures and
summary tables. External baseline methods require their original package
environments; the artifact retains wrapper scripts, parameters, input/output
formats, and summarized outputs rather than forcing all third-party tools into
one environment.
