# Rare-cell-type and marker-weight analysis

This module contains the retained data, statistics, plotting code, historical
experiment drivers, and four manuscript figures used for the rare-cell-type
benchmark and marker-weight analysis.

## Contents

- `data/inputs/`: complete non-result inputs for the retained Heart and Skin
  tissues, including bulk, truth, signature, single-cell, subject, and marker
  files.
- `data/*.csv`: retained truth/prediction values and marker weights for the two
  cases.
- `results/`: marker-weight summary table.
- `code/experiment/`: MetaSort execution, rare-cell-type aggregation, and
  marker-weight extraction scripts.
- `code/plot_*.R`: retained-result plotting scripts.
- `figures/`: four figures in PNG and PDF formats.
- `reproduce_all.sh`: ordered entry point for all plots and result tables.
- `MANIFEST.sha256`: integrity manifest for retained files.

The shared NNLS, DWLS, HiDecon, LinDeconSeq, MuSiC, and BayesPrism execution
scripts are stored once under `../../scripts/baselines/`.

## Regenerate retained figures

From the repository root:

```bash
./results/rare_celltype/reproduce_all.sh
```

R dependencies are `dplyr`, `ggplot2`, `readr`, and `scales`. The R figures use
Times New Roman and Cairo PDF output. Regeneration overwrites matching files in
`figures/` and `results/`.

This command starts from the retained truth, prediction, and marker-weight
tables. It does not rerun the upstream deconvolution methods.

## Archived experiment drivers

- `code/experiment/run_nar30_first10_compare.py` runs MetaSort, reads the six
  baseline prediction files, and writes benchmark summaries. Despite its
  historical filename, this experiment used `--n-mixtures 100`.
- `code/experiment/analyze_nar_rare_celltypes.py` combines truth and all seven
  methods' predictions into the rare-cell-type analysis tables.
- `code/experiment/rerun_nar_version2_weights.py` reruns selected MetaSort cases
  while retaining the per-gene weights used by the marker-weight analysis.

The MetaSort algorithm itself is imported from the repository-level
`metasort/` package. The two retained tissue inputs are available under
`data/inputs/CELlxGENE/Basal_Zone_Of_Heart/` and
`data/inputs/Tabula-Sapiens/Skin/`. The archived execution scripts were not run
during release assembly.

Verify the retained artifact before regeneration with:

```bash
cd results/rare_celltype && sha256sum -c MANIFEST.sha256
```
