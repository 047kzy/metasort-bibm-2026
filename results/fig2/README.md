# Figure 2: Tongue Curvature and Iterative Optimization

This directory contains the Tongue Mixture case used for Figure 2.

## Manuscript Panels

- Figure 2A: oracle curvature-guided optimization initialized at the true
  composition, compared with simplex-constrained least squares on z-score
  normalized data.
- Figure 2B: final L1 error comparison among MetaSort, MuSiC, DWLS,
  LinDeconSeq, and NNLS.
- Figure 2C: iterative MetaSort trajectory from least-squares initialization to
  convergence after 13 iterations.

## Final Figure

- `publication_figures/bibm_metasort_tongue_case_combined_revised_with_true_init_l1.pdf`
- `publication_figures/bibm_metasort_tongue_case_combined_revised_with_true_init_l1.png`
- `publication_figures/bibm_metasort_tongue_case_combined_revised_with_true_init_l1.svg`
- `publication_figures/bibm_metasort_tongue_case_combined_revised_with_true_init_l1.tiff`

## Source Data

- `publication_figures/bibm_true_init_metasort_vs_simplex_zscore_l1_summary.csv`
- `publication_figures/bibm_true_init_metasort_vs_simplex_zscore_proportions.csv`
- `publication_figures/bibm_true_init_metasort_vs_simplex_zscore_curvature_log_eig_summary.csv`
- `publication_figures/bibm_true_init_metasort_vs_simplex_zscore_curvature_log_eig_values.csv`
- `publication_figures/bibm_metasort_tongue_case_combined_revised_values.csv`
- `publication_figures/bibm_metasort_tongue_case_iteration_0_13.csv`

## Regeneration

```bash
python merge_tongue_case_combined_with_true_init_l1.py
```

The script regenerates the final combined figure from the retained source-data
CSV files.
