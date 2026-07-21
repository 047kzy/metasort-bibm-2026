# MetaSort Anonymous Review Artifact

This repository contains the anonymous MetaSort implementation and the
retained data, scripts, and figures used for the BIBM submission.

## Repository layout

```text
.
├── metasort/                  # Core Python package
├── requirements.txt           # Core Python dependencies
├── scripts/baselines/         # Shared non-MetaSort method runners
└── results/
    ├── hierarchy/             # Hierarchy-aware deconvolution benchmark
    ├── rare_celltype/         # Rare-cell-type and marker-weight figures
    ├── oracle_experiment/     # Selected Tongue oracle-initialization case
    └── reference_recovery/    # Bulk-reference recovery experiment
```

Each result module has a short README describing its retained inputs, outputs,
dependencies, and optional upstream reruns. Paths and commands below are
relative to the repository root.

## Setup

Install the Python dependencies used by the core package and most plotting
workflows:

```bash
python -m pip install -r requirements.txt
```

The hierarchy and rare-cell-type modules also use R. Their module READMEs list
the required R packages and external benchmark-method dependencies.

## Reproduce the retained figures

Hierarchy benchmark:

```bash
Rscript results/hierarchy/scripts/build_figure_metrics.R
Rscript results/hierarchy/scripts/plot_accuracy_boxplots.R
```

Rare-cell-type analysis:

```bash
./results/rare_celltype/reproduce_all.sh
```

Selected Tongue oracle-initialization case:

```bash
python -m pip install -r results/oracle_experiment/requirements.txt
python results/oracle_experiment/scripts/validate_retained_data.py
python results/oracle_experiment/scripts/make_final_figure.py
```

Bulk-reference recovery figure:

```bash
python results/reference_recovery/scripts/generate_figure.py
```

## Result-module documentation

- [`results/hierarchy/README.md`](results/hierarchy/README.md)
- [`results/rare_celltype/README.md`](results/rare_celltype/README.md)
- [`results/oracle_experiment/README.md`](results/oracle_experiment/README.md)
- [`results/reference_recovery/README.md`](results/reference_recovery/README.md)

The release excludes caches, notebook checkpoints, local temporary files, and
author-specific paths. Raw inputs that cannot be redistributed are identified
in the relevant module README.
