from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from scipy.optimize import nnls

from metasort import MetaSortConfig, MetaSortSolver
from metasort.algorithm import (
    _joint_gene_zscore,
    _normalize_columns_to_sum_one,
    _normalize_vector_to_sum_one,
    solve_simplex_constrained_ls,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-output-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--mixture", required=True)
    parser.add_argument("--top-genes", type=int, required=True)
    parser.add_argument(
        "--gene-list-file",
        type=Path,
        default=None,
        help="Optional one-column gene list. If provided, --top-genes is ignored for selection.",
    )
    parser.add_argument(
        "--score-mode",
        choices=["meta_weight", "meta_weight_residual"],
        default="meta_weight_residual",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def normalize_nonnegative(values: np.ndarray) -> np.ndarray:
    values = np.clip(np.asarray(values, dtype=float), 0.0, None)
    total = float(np.sum(values))
    if total <= 0.0 or not np.isfinite(total):
        return np.full(values.shape, 1.0 / len(values), dtype=float)
    return values / total


def score_prediction(pred: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    pred = normalize_nonnegative(pred)
    truth = normalize_nonnegative(truth)
    diff = pred - truth
    l1 = float(np.sum(np.abs(diff)))
    pearson = float("nan")
    if pred.size > 1 and np.std(pred) > 0.0 and np.std(truth) > 0.0:
        pearson = float(np.corrcoef(pred, truth)[0, 1])
    return {
        "l1": l1,
        "accuracy": 1.0 - l1 / 2.0,
        "mae": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "pearson": pearson,
    }


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def preprocess(reference: pd.DataFrame, bulk: pd.Series, cell_types: list[str], genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    signature = _normalize_columns_to_sum_one(
        reference.loc[genes, cell_types].to_numpy(dtype=float),
        name="signature",
    )
    bulk_vec = _normalize_vector_to_sum_one(
        bulk.loc[genes].to_numpy(dtype=float),
        name="bulk",
    )
    return _joint_gene_zscore(signature, bulk_vec)


def raw_normalized(reference: pd.DataFrame, bulk: pd.Series, cell_types: list[str], genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    signature = _normalize_columns_to_sum_one(
        reference.loc[genes, cell_types].to_numpy(dtype=float),
        name="signature",
    )
    bulk_vec = _normalize_vector_to_sum_one(
        bulk.loc[genes].to_numpy(dtype=float),
        name="bulk",
    )
    return signature, bulk_vec


def solve_nnls_normalized(signature: np.ndarray, bulk: np.ndarray) -> np.ndarray:
    coef, _ = nnls(signature, bulk)
    return normalize_nonnegative(coef)


def solve_dwls_style(signature: np.ndarray, bulk: np.ndarray) -> np.ndarray:
    initial = solve_simplex_constrained_ls(signature, bulk)
    fitted = signature @ initial
    positive = fitted[(fitted > 0.0) & np.isfinite(fitted)]
    if positive.size == 0:
        return initial
    safe_fitted = np.where((fitted > 0.0) & np.isfinite(fitted), fitted, float(np.min(positive)))
    weights = 1.0 / np.square(safe_fitted)
    finite = weights[np.isfinite(weights) & (weights > 0.0)]
    if finite.size == 0:
        return initial
    weights = weights / float(np.min(finite))
    cap = float(np.quantile(weights[np.isfinite(weights)], 0.95))
    weights = np.minimum(weights, cap)
    weights[~np.isfinite(weights)] = cap
    return solve_simplex_constrained_ls(signature, bulk, weights=weights, initial=initial)


def evaluate_methods(
    reference: pd.DataFrame,
    bulk: pd.Series,
    truth: np.ndarray,
    cell_types: list[str],
    genes: list[str],
    solver: MetaSortSolver,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw_sig, raw_bulk = raw_normalized(reference, bulk, cell_types, genes)
    z_sig, z_bulk = _joint_gene_zscore(raw_sig, raw_bulk)
    rows: list[dict[str, Any]] = []
    pred_rows: list[dict[str, Any]] = []

    methods: dict[str, np.ndarray] = {}
    result = solver.solve(z_sig, z_bulk, cell_types=cell_types)
    methods["MetaSort"] = normalize_nonnegative(np.asarray(result.proportions, dtype=float))
    methods["NNLS_zscore"] = solve_nnls_normalized(z_sig, z_bulk)
    methods["NNLS_raw"] = solve_nnls_normalized(raw_sig, raw_bulk)
    methods["SimplexLS_raw"] = solve_simplex_constrained_ls(raw_sig, raw_bulk)
    methods["DWLS_style_raw"] = solve_dwls_style(raw_sig, raw_bulk)

    for method, pred in methods.items():
        rows.append(
            {
                "method": method,
                "n_genes": len(genes),
                "n_cell_types": len(cell_types),
                **score_prediction(pred, truth),
            }
        )
        for cell_type, truth_value, pred_value in zip(cell_types, truth, pred):
            pred_rows.append(
                {
                    "method": method,
                    "cell_type": cell_type,
                    "truth": float(truth_value),
                    "prediction": float(pred_value),
                    "abs_error": float(abs(pred_value - truth_value)),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(pred_rows)


def main() -> int:
    args = parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"{args.output_dir} exists and is not empty. Use --overwrite.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except Exception:
        pass

    source_case = args.source_output_dir / "sampled_cases" / args.dataset / args.tissue
    reference = read_matrix(source_case / "selected_cellTypeExpr.txt")
    bulk_df = read_matrix(source_case / "sampled_bulk.txt")
    truth_df = read_matrix(source_case / "sampled_bulkRatio.txt")
    if args.mixture not in bulk_df.columns:
        raise ValueError(f"{args.mixture} not found in sampled_bulk.txt")

    reference.index = reference.index.astype(str)
    reference.columns = reference.columns.astype(str)
    bulk_df.index = bulk_df.index.astype(str)
    truth_df.index = truth_df.index.astype(str)
    cell_types = [cell_type for cell_type in reference.columns if cell_type in truth_df.index]
    all_genes = reference.index.intersection(bulk_df.index).astype(str).tolist()
    bulk = bulk_df.loc[all_genes, args.mixture]
    truth = normalize_nonnegative(truth_df.loc[cell_types, args.mixture].to_numpy(dtype=float))

    cfg = MetaSortConfig()
    solver = MetaSortSolver(cfg)
    z_sig, z_bulk = preprocess(reference, bulk, cell_types, all_genes)
    full_result = solver.solve(z_sig, z_bulk, cell_types=cell_types)
    meta_weights = np.asarray(full_result.meta_weights, dtype=float)
    residual = np.abs(z_bulk - z_sig @ normalize_nonnegative(np.asarray(full_result.proportions, dtype=float)))
    if args.score_mode == "meta_weight":
        score = meta_weights
    else:
        score = meta_weights / (1.0 + residual)
    ranked = (
        pd.DataFrame(
            {
                "gene": all_genes,
                "meta_weight": meta_weights,
                "abs_residual": residual,
                "score": score,
            }
        )
        .sort_values(["score", "meta_weight"], ascending=False, kind="mergesort")
        .reset_index(drop=True)
    )
    if args.gene_list_file is not None:
        gene_list = [line.strip() for line in args.gene_list_file.read_text().splitlines()]
        selected_genes = [gene for gene in gene_list if gene and gene in all_genes]
        if not selected_genes:
            raise ValueError(f"No genes from {args.gene_list_file} overlap source genes")
    else:
        selected_genes = ranked["gene"].head(args.top_genes).tolist()
    metrics, predictions = evaluate_methods(reference, bulk, truth, cell_types, selected_genes, solver)

    case_out = args.output_dir / "sampled_cases" / args.dataset / args.tissue
    case_out.mkdir(parents=True, exist_ok=True)
    reference.loc[selected_genes, cell_types].to_csv(case_out / "selected_cellTypeExpr.txt", sep="\t")
    bulk_df.loc[selected_genes, [args.mixture]].to_csv(case_out / "sampled_bulk.txt", sep="\t")
    truth_df.loc[cell_types, [args.mixture]].to_csv(case_out / "sampled_bulkRatio.txt", sep="\t")
    pd.Series(selected_genes, name="gene").to_csv(case_out / "selected_genes.txt", sep="\t", index=False)
    ranked.to_csv(case_out / "metasort_gene_scores.csv", index=False)

    for frame in (metrics, predictions):
        frame.insert(0, "mixture_index", 1)
        frame.insert(0, "mixture", args.mixture)
        frame.insert(0, "case", f"{args.dataset}/{args.tissue}")
        frame.insert(0, "tissue", args.tissue)
        frame.insert(0, "dataset", args.dataset)

    metrics["runtime_sec"] = np.nan
    metrics.to_csv(args.output_dir / "metrics_long.csv", index=False)
    predictions.to_csv(args.output_dir / "predictions_long.csv", index=False)
    pd.DataFrame(
        [
            {
                "dataset": args.dataset,
                "tissue": args.tissue,
                "case": f"{args.dataset}/{args.tissue}",
                "mixture": args.mixture,
                "selected_genes": len(selected_genes),
                "source_genes": len(all_genes),
                "top_genes": args.top_genes,
                "gene_list_file": "" if args.gene_list_file is None else str(args.gene_list_file),
                "score_mode": args.score_mode,
                "source_output_dir": str(args.source_output_dir),
                "metasort_config": json.dumps(asdict(cfg), sort_keys=True),
            }
        ]
    ).to_csv(args.output_dir / "gene_selection_summary.csv", index=False)
    print(metrics.sort_values("l1").to_string(index=False), flush=True)
    print(f"output_dir: {args.output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
