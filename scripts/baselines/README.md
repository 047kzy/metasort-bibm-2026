# Baseline deconvolution runners

This directory centralizes the execution scripts for the seven
non-MetaSort methods used by the retained benchmarks:

- `NNLS.R`
- `DWLS.R`
- `HiDecon.R`
- `LinDeconSeq.R`
- `MuSic.R`
- `BayesPrism.R`
- `SimplexLS.py`
- `run_hierarchy_methods.R` (hierarchy-aware HiDecon and MuSiC)

The scripts expect tissue directories containing the original benchmark file
layout (`bulk.txt`, `signature.txt`, `singleCellExpr.txt`,
`singleCellLabels.txt`, and optional marker or subject files). Set `ROOT_PATH`
to the dataset directory for all scripts except BayesPrism, which uses
`BAYESPRISM_ROOT_PATH`. Each method writes its prediction matrix under the
tissue's `result/` directory.

`SimplexLS.py` is a thin command-line wrapper around the repository's existing
`solve_simplex_constrained_ls` implementation. It supports both the raw-input
variant and the normalized joint-z-score variant used in the oracle experiment;
the solver itself is not duplicated here.

From the repository root, the two variants are selected explicitly:

```bash
PYTHONPATH=. python scripts/baselines/SimplexLS.py --data-root /path/to/data --mode raw
PYTHONPATH=. python scripts/baselines/SimplexLS.py --data-root /path/to/data --mode zscore
```

They write `SimplexLS_raw.txt` and `SimplexLS_zscore.txt`, respectively.

`run_hierarchy_methods.R` reconstructs the hierarchy matrices from the retained
JSON trees and runs the hierarchy-aware HiDecon and MuSiC variants for the
Salivary Gland and Thymus benchmark. From the repository root:

```bash
Rscript scripts/baselines/run_hierarchy_methods.R
```

These files preserve the experiment-layer method calls. They were not executed
while assembling this release. Required R packages must be installed
from their respective upstream projects.
