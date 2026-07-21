from __future__ import annotations

import csv
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = BASE_DIR / "data" / "inputs"
RESULTS_DIR = BASE_DIR / "results"
sys.path.insert(0, str(BASE_DIR))

from _metasort_common import (  # noqa: E402
    HIERARCHY_DIR,
    N_MIXTURES,
    TISSUES,
    WeightRecordingHierarchicalMetaSortSolver,
    load_hierarchy,
    make_config,
    normalize_nonnegative,
    preprocess,
    score_prediction,
    write_hierarchical_weights,
)


METHOD = "MetaSort_jk_hierarchy_avg2p5_tol0005"


def make_test_config():
    config = make_config(hierarchical=True)
    config.averaging_old_weight = 2.5
    config.convergence_tol = 0.005
    return config


def run_case(tissue: str) -> dict[str, object]:
    case_dir = INPUT_DIR / tissue
    signature_df = pd.read_csv(case_dir / "signature.tsv.gz", sep="\t", index_col=0)
    bulk_df = pd.read_csv(case_dir / "bulk.tsv.gz", sep="\t", index_col=0)
    truth_df = pd.read_csv(case_dir / "bulkRatio.tsv", sep="\t", index_col=0)

    cell_types = signature_df.columns.astype(str).to_list()
    truth_df.index = truth_df.index.astype(str)
    missing_truth = sorted(set(cell_types) - set(truth_df.index))
    if missing_truth:
        raise ValueError(f"{tissue}: truth is missing cell types: {missing_truth}")

    mixtures = bulk_df.columns.astype(str).to_list()[:N_MIXTURES]
    solver = WeightRecordingHierarchicalMetaSortSolver(make_test_config())
    hierarchy = load_hierarchy(HIERARCHY_DIR / f"{tissue}_hierarchy.json")

    result = pd.DataFrame(index=mixtures, columns=cell_types, dtype=float)
    metric_rows: list[dict[str, object]] = []
    weight_rows: list[dict[str, object]] = []
    weight_summary_rows: list[dict[str, object]] = []
    start = time.time()

    for idx, mixture in enumerate(mixtures, start=1):
        signature, bulk, genes = preprocess(signature_df, bulk_df, mixture)
        truth = normalize_nonnegative(truth_df.loc[cell_types, mixture].to_numpy(dtype=float))

        solver.weight_records = []
        solved = solver.solve(signature, bulk, cell_types=cell_types, hierarchy=hierarchy)
        stage_results = getattr(solved, "stage_results", [])
        iterations = sum(
            int(stage.result.iterations)
            for stage in stage_results
            if getattr(stage, "result", None) is not None
        )

        for record in solver.weight_records:
            selected = np.asarray(record["selected_genes"], dtype=int)
            meta_weights = np.asarray(record["meta_weights"], dtype=float)
            metrics = dict(record["meta_metrics"])
            weight_summary_rows.append(
                {
                    "tissue": tissue,
                    "method": METHOD,
                    "mixture": mixture,
                    "mixture_index": idx,
                    "solve_index": record["solve_index"],
                    "visible_nodes": record["visible_nodes"],
                    "n_genes": int(len(selected)),
                    "iterations": int(record["iterations"]),
                    "converged": bool(record["converged"]),
                    "meta_weight_min": float(np.min(meta_weights)),
                    "meta_weight_mean": float(np.mean(meta_weights)),
                    "meta_weight_max": float(np.max(meta_weights)),
                    **metrics,
                }
            )
            weight_rows.extend(
                {
                    "tissue": tissue,
                    "method": METHOD,
                    "mixture": mixture,
                    "mixture_index": idx,
                    "solve_index": record["solve_index"],
                    "visible_nodes": record["visible_nodes"],
                    "gene": genes[int(gene_idx)],
                    "meta_weight": float(weight),
                }
                for gene_idx, weight in zip(selected, meta_weights)
            )

        pred = normalize_nonnegative(np.asarray(solved.proportions, dtype=float))
        result.loc[mixture, :] = pred
        metric_rows.append(
            {
                "tissue": tissue,
                "method": METHOD,
                "mixture": mixture,
                "mixture_index": idx,
                "n_genes": len(genes),
                "iterations": iterations,
                **score_prediction(pred, truth),
            }
        )
        print(f"{METHOD} {tissue} {idx}/{len(mixtures)} {mixture}", flush=True)

    result_dir = RESULTS_DIR / tissue
    result_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(result_dir / f"{METHOD}.tsv", sep="\t")
    pd.DataFrame(metric_rows).to_csv(result_dir / f"{METHOD}_metrics.csv", index=False)
    write_hierarchical_weights(
        result_dir=result_dir,
        method=METHOD,
        rows=weight_rows,
        summary_rows=weight_summary_rows,
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

    log_path = RESULTS_DIR / f"{METHOD}_run_log.tsv"
    with log_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    config = make_test_config()
    config_row = asdict(config) | {
        "n_mixtures": N_MIXTURES,
        "confidence_weights": "off",
        "output_dir": str(BASE_DIR),
        "source": "METASORT_JK_PATH or installed metasort package",
        "script": str(Path(__file__)),
    }
    pd.DataFrame([config_row]).to_csv(RESULTS_DIR / f"{METHOD}_config.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
