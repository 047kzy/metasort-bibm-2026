# Selected Tongue oracle-initialization experiment

This module retains the selected Tabula Sapiens Tongue `Mixture10` case and the
three-panel figure comparing a true-initialized MetaSort step, final method
errors, and the MetaSort iteration trace.

## Contents

- `data/selected_case/`: expression, bulk, truth, and selected-gene inputs.
- `data/figure_source/`: the three tables read by the final plotting script.
- `results/case_metrics/`: method-level metrics for the selected case.
- `results/method_predictions/`: standardized per-cell-type predictions.
- `results/metasort_run/` and `results/true_init_run/`: retained run traces.
- `results/final_figure/`: final PDF and PNG outputs.
- `scripts/`: export, validation, derivation, and plotting utilities.

## Reproduce and validate

From the repository root:

```bash
python -m pip install -r results/oracle_experiment/requirements.txt
python results/oracle_experiment/scripts/export_method_predictions.py
python results/oracle_experiment/scripts/validate_retained_data.py
python results/oracle_experiment/scripts/make_final_figure.py
```

Add `--with-tiff` to the last command only when a 600-dpi TIFF is required.

The optional upstream derivation below uses the root MetaSort package and
refreshes the true-initialization diagnostics and Panel A source table:

```bash
PYTHONPATH=. python results/oracle_experiment/scripts/derive_true_init_panel_data.py
```

Before regenerating files, the retained artifact can be verified with:

```bash
cd results/oracle_experiment && sha256sum -c SHA256SUMS
```

Exploratory searches and unrelated cases are intentionally excluded.
