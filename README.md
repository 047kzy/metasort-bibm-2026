# MetaSort Anonymous Review Artifact

This repository contains the anonymous code, retained result tables, figure
source data, and reproducibility wrappers for the BIBM submission. It is
organized to follow the Results section of the manuscript.

## Contents

- `code/`: MetaSort source code, benchmark wrappers, plotting scripts, hierarchy
  files, and Python requirements.
- `results/fig2/`: Tongue mixture curvature and iterative optimization case.
- `results/fig3_1/`: Salivary Gland and Thymus hierarchy-aware benchmark.
- `results/fig3_2/`: rare-cell-type and marker-weight analyses for Figure 3.
- `results/fig4/`: `bulk_TSP14__reference_TSP10` expression-profile and
  reference-recovery analysis results.
- `MANUSCRIPT_RESULTS_MAP.md`: statement-to-file map for the Results section.

## Minimal Python Setup

From the repository root:

```bash
python -m pip install -r code/requirements.txt
```

The core implementation is in `code/metasort/`. Baseline wrappers for external
methods are retained under `code/scripts/` and `results/*/scripts/`, but those
methods require their own R/Python package environments.

## Reproducibility Notes

The repository keeps the exact wrapper scripts, parameters, summary tables, and
figure outputs used for the manuscript. Full benchmark input data will be
provided separately. Paths such as `<H5AD_DATA_ROOT>` and `<CONDA_ROOT>` are
anonymous placeholders for optional raw h5ad files and local software
environments that are not required for inspecting the retained result files.

For review, the directly inspectable material is:

- final and intermediate result tables in `results/fig*/`;
- figure source data CSV/TSV files in each figure directory;
- plotting scripts used to regenerate retained figures;
- MetaSort source code and benchmark wrappers in `code/`.

## Directory Mapping

The current artifact uses figure-based names:

- `fig2`: curvature-guided optimization in the Tongue mixture case;
- `fig3_1`: hierarchy-aware deconvolution for Salivary Gland and Thymus;
- `fig3_2`: rare-cell-type recovery and marker-weight trends;
- `fig4`: `bulk_TSP14__reference_TSP10` expression-profile and reference-recovery
  analysis.

No original Git history, personal remote URL, notebook execution output, local
cache files, or author-identifying account metadata is included in this review
artifact.
