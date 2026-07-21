# Bulk-reference recovery experiment

This module contains the retained data and figures comparing the recovered
TSP14 bulk-domain reference with the TSP10 input reference.

## Contents

- `data/`: input, true, and recovered reference matrices; cell-type RMSE; and
  representative-gene selections.
- `figures/`: two component figures and the combined manuscript figure in PNG
  and PDF formats.
- `scripts/generate_figure.py`: figure-only reproduction from retained data.
- `scripts/run_skin_anchor_reference_recovery_test.py`: optional upstream
  MetaSort experiment wrapper.

Expression matrices use `log1p(CP10000)`, with genes in rows and cell types in
columns.

## Reproduce the figure

From the repository root, after installing the root `requirements.txt`:

```bash
python results/reference_recovery/scripts/generate_figure.py
```

Set `BIBM_FONT_DIR=/path/to/fonts` to use local Times New Roman files. The
script otherwise falls back to an installed serif font.

## Optional experiment rerun

The benchmark signature is retained under
`results/rare_celltype/data/inputs/Tabula-Sapiens/Skin/signature.txt`; the raw
Skin h5ad is not redistributed. With that raw input available, run:

```bash
PYTHONPATH=. python results/reference_recovery/scripts/run_skin_anchor_reference_recovery_test.py \
  --h5ad /path/to/Skin.h5ad \
  --output-dir results/reference_recovery/rerun_outputs
```

The retained files in `data/` are derived outputs and are not overwritten by
the command above.
