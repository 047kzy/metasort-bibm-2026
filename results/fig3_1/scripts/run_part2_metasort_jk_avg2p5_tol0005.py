from __future__ import annotations

import csv
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path("results/fig3_1")
sys.path.insert(0, str(BASE_DIR))

from run_part2_metasort_jk_original_with_weights import (  # noqa: E402
    N_MIXTURES,
    TISSUES,
    WeightRecordingMetaSortSolver,
    make_config,
    normalize_nonnegative,
    preprocess,
    score_prediction,
    write_direct_weights,
)


METHOD = "MetaSort_jk_avg2p5_tol0005"


def make_test_config():
    config = make_config(hierarchical=False)
    config.averaging_old_weight = 2.5
    config.convergence_tol = 0.005
    return config


def run_case(tissue: str) -> dict[str, object]:
    case_dir = BASE_DIR / tissue
    signature_df = pd.read_csv(case_dir / "signature.txt", sep="\t", index_col=0)
    bulk_df = pd.read_csv(case_dir / "bulk.txt", sep="\t", index_col=0)
    truth_df = pd.read_csv(case_dir / "bulkRatio.txt", sep="\t", index_col=0)

    cell_types = signature_df.columns.astype(str).to_list()
    truth_df.index = truth_df.index.astype(str)
    missing_truth = sorted(set(cell_types) - set(truth_df.index))
    if missing_truth:
        raise ValueError(f"{tissue}: truth is missing cell types: {missing_truth}")

    mixtures = bulk_df.columns.astype(str).to_list()[:N_MIXTURES]
    solver = WeightRecordingMetaSortSolver(make_test_config())

    result = pd.DataFrame(index=mixtures, columns=cell_types, dtype=float)
    metric_rows: list[dict[str, object]] = []
    weight_columns: dict[str, np.ndarray] = {}
    weight_summary_rows: list[dict[str, object]] = []
    genes_expected: list[str] | None = None
    start = time.time()

    for idx, mixture in enumerate(mixtures, start=1):
        signature, bulk, genes = preprocess(signature_df, bulk_df, mixture)
        truth = normalize_nonnegative(truth_df.loc[cell_types, mixture].to_numpy(dtype=float))

        solver.last_meta_weights = None
        solved = solver.solve(signature, bulk, cell_types=cell_types)
        if solver.last_meta_weights is None:
            raise RuntimeError(f"{tissue}/{mixture}: no meta weights were recorded")

        if genes_expected is None:
            genes_expected = genes
        elif genes != genes_expected:
            raise ValueError(f"{tissue}: shared gene order changed between mixtures")

        pred = normalize_nonnegative(np.asarray(solved.proportions, dtype=float))
        result.loc[mixture, :] = pred
        weight_columns[mixture] = solver.last_meta_weights
        weight_summary_rows.append(
            {
                "tissue": tissue,
                "method": METHOD,
                "mixture": mixture,
                "mixture_index": idx,
                "n_genes": len(genes),
                "iterations": int(solved.iterations),
                "converged": bool(solved.converged),
                "meta_weight_min": float(solved.meta_weight_min),
                "meta_weight_mean": float(solved.meta_weight_mean),
                "meta_weight_max": float(solved.meta_weight_max),
                "total_weight_min": float(solved.total_weight_min),
                "total_weight_mean": float(solved.total_weight_mean),
                "total_weight_max": float(solved.total_weight_max),
                **(solved.meta_weight_metrics or {}),
            }
        )
        metric_rows.append(
            {
                "tissue": tissue,
                "method": METHOD,
                "mixture": mixture,
                "mixture_index": idx,
                "n_genes": len(genes),
                "iterations": int(solved.iterations),
                **score_prediction(pred, truth),
            }
        )
        print(f"{METHOD} {tissue} {idx}/{len(mixtures)} {mixture}", flush=True)

    result_dir = BASE_DIR / tissue / "result"
    result_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(result_dir / f"{METHOD}.txt", sep="\t")
    pd.DataFrame(metric_rows).to_csv(result_dir / f"{METHOD}_metrics.csv", index=False)
    if genes_expected is None:
        raise ValueError(f"{tissue}: no weights were recorded")
    write_direct_weights(
        result_dir=result_dir,
        method=METHOD,
        weight_columns=weight_columns,
        genes=genes_expected,
        rows=weight_summary_rows,
    )
    return {
        "tissue": tissue,
        "method": METHOD,
        "status": "OK",
        "rows": result.shape[0],
        "cols": result.shape[1],
        "seconds": round(time.time() - start, 3),
        "error": "",
    }


def main() -> None:
    rows: list[dict[str, object]] = []
    for tissue in TISSUES:
        try:
            rows.append(run_case(tissue))
        except Exception as exc:
            rows.append(
                {
                    "tissue": tissue,
                    "method": METHOD,
                    "status": "ERROR",
                    "rows": "",
                    "cols": "",
                    "seconds": "",
                    "error": repr(exc),
                }
            )
            print(f"{METHOD} {tissue} ERROR: {exc!r}", flush=True)

    log_path = BASE_DIR / f"{METHOD}_run_log.tsv"
    with log_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    config = make_test_config()
    config_row = asdict(config) | {
        "n_mixtures": N_MIXTURES,
        "confidence_weights": "off",
        "output_dir": str(BASE_DIR),
        "source": "code",
        "script": str(Path(__file__)),
    }
    pd.DataFrame([config_row]).to_csv(BASE_DIR / f"{METHOD}_config.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
