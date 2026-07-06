from __future__ import annotations

import argparse
import hashlib
import json
import time
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


DEFAULT_DATA_ROOT = Path("data/benchmark_inputs")
DEFAULT_OUTPUT_DIR = Path("code/outputs/marker_nonmarker_sampled_bulk_search")
DEFAULT_DATASETS = ("CELLxGENE", "Tabula-Sapiens")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sample pseudo-bulk from single cells using each tissue bulkRatio, "
            "run MetaSort and lightweight baselines on marker + non-marker genes, "
            "and rank cases where MetaSort is more accurate."
        )
    )
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--n-mixtures", type=int, default=10)
    parser.add_argument("--total-cells", type=int, default=5000)
    parser.add_argument("--marker-genes", type=int, default=250)
    parser.add_argument("--nonmarker-genes", type=int, default=250)
    parser.add_argument(
        "--gene-selection",
        choices=["ranked", "random"],
        default="ranked",
        help="ranked keeps the historical marker-order/low-CV non-marker selection; random samples both pools.",
    )
    parser.add_argument("--min-genes", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260630)
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--tissues", nargs="+", default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--lambda-hessian", type=float, default=0.03)
    parser.add_argument("--lambda-residual", type=float, default=0.02)
    parser.add_argument("--convergence-tol", type=float, default=0.005)
    parser.add_argument("--averaging-old-weight", type=float, default=2.5)
    parser.add_argument("--meta-weight-baseline", type=float, default=10.0)
    parser.add_argument("--meta-weight-floor", type=float, default=1.0)
    parser.add_argument("--lambda3", type=float, default=0.0001)
    parser.add_argument("--lambda4", type=float, default=1e-5)
    parser.add_argument("--sqrt-sphere-hessian", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--normalize-meta-weight-mean", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def normalize_nonnegative(values: np.ndarray) -> np.ndarray:
    values = np.clip(np.asarray(values, dtype=float), 0.0, None)
    total = float(np.sum(values))
    if total <= 0.0 or not np.isfinite(total):
        return np.full(values.shape, 1.0 / len(values), dtype=float)
    return values / total


def pearson_safe(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or np.std(x) <= 0.0 or np.std(y) <= 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def score_prediction(pred: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    pred = normalize_nonnegative(pred)
    truth = normalize_nonnegative(truth)
    diff = pred - truth
    l1 = float(np.sum(np.abs(diff)))
    return {
        "l1": l1,
        "accuracy": 1.0 - l1 / 2.0,
        "mae": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "pearson": pearson_safe(pred, truth),
    }


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def find_case_dirs(data_root: Path, datasets: list[str], tissues: list[str] | None) -> list[Path]:
    case_dirs: list[Path] = []
    tissue_filter = None if tissues is None else set(tissues)
    required = {
        "bulkRatio.txt",
        "cellTypeExpr.txt",
        "singleCellExpr.txt",
        "singleCellLabels.txt",
    }
    for dataset in datasets:
        dataset_dir = data_root / dataset
        if not dataset_dir.exists():
            continue
        for case_dir in sorted(dataset_dir.iterdir()):
            if not case_dir.is_dir() or case_dir.name == "result":
                continue
            if tissue_filter is not None and case_dir.name not in tissue_filter:
                continue
            if all((case_dir / name).exists() for name in required):
                case_dirs.append(case_dir)
    return case_dirs


def load_marker_genes(case_dir: Path, cell_types: list[str]) -> list[str]:
    marker_path = case_dir / "markers.txt"
    if marker_path.exists():
        markers = pd.read_csv(marker_path, sep="\t")
        if {"cell_type", "gene"}.issubset(markers.columns):
            markers["cell_type"] = markers["cell_type"].astype(str)
            markers["gene"] = markers["gene"].astype(str)
            keep_types = set(cell_types)
            ordered = markers.loc[markers["cell_type"].isin(keep_types), "gene"].tolist()
            return list(dict.fromkeys(gene for gene in ordered if gene))

    sig_genes_path = case_dir / "sigGenes.txt"
    if sig_genes_path.exists():
        sig_genes = pd.read_csv(sig_genes_path, sep="\t")
        col = sig_genes.columns[0]
        return list(dict.fromkeys(sig_genes[col].dropna().astype(str).tolist()))

    signature_path = case_dir / "signature.txt"
    if signature_path.exists():
        return read_matrix(signature_path).index.astype(str).tolist()
    return []


def select_gene_subset(
    reference: pd.DataFrame,
    marker_genes: list[str],
    marker_gene_count: int,
    nonmarker_gene_count: int,
    mode: str = "ranked",
    rng: np.random.Generator | None = None,
) -> tuple[list[str], list[str], list[str]]:
    reference.index = reference.index.astype(str)
    marker_available = [gene for gene in marker_genes if gene in reference.index]
    if mode == "random":
        if rng is None:
            raise ValueError("rng is required when mode='random'")
        marker_n = min(marker_gene_count, len(marker_available))
        marker_selected = (
            rng.choice(marker_available, size=marker_n, replace=False).tolist()
            if marker_n > 0
            else []
        )
    else:
        marker_selected = marker_available[: min(marker_gene_count, len(marker_available))]

    marker_set = set(marker_available)
    nonmarker_candidates = reference.index.difference(pd.Index(marker_set), sort=False)
    nonmarker_ref = reference.loc[nonmarker_candidates]
    means = nonmarker_ref.mean(axis=1).astype(float)
    stds = nonmarker_ref.std(axis=1).astype(float)
    if mode == "random":
        expressed_nonmarkers = means.loc[means > 0.0].index.astype(str).tolist()
        if not expressed_nonmarkers:
            expressed_nonmarkers = nonmarker_candidates.astype(str).tolist()
        nonmarker_n = min(nonmarker_gene_count, len(expressed_nonmarkers))
        nonmarker_selected = (
            rng.choice(expressed_nonmarkers, size=nonmarker_n, replace=False).tolist()
            if nonmarker_n > 0
            else []
        )
    else:
        expressed = means > max(float(means.quantile(0.25)), 0.0)
        cv = stds / means.clip(lower=1e-12)
        ranked_nonmarkers = (
            pd.DataFrame({"mean": means, "cv": cv})
            .loc[expressed]
            .sort_values(["cv", "mean"], ascending=[True, False], kind="mergesort")
            .index.astype(str)
            .tolist()
        )
        nonmarker_selected = ranked_nonmarkers[: min(nonmarker_gene_count, len(ranked_nonmarkers))]
    selected = list(dict.fromkeys(marker_selected + nonmarker_selected))
    return selected, marker_selected, nonmarker_selected


def sample_counts_for_truth(truth: np.ndarray, total_cells: int, rng: np.random.Generator) -> np.ndarray:
    truth = normalize_nonnegative(truth)
    return rng.multinomial(total_cells, truth)


def sample_bulk_from_single_cells(
    sc_expr: pd.DataFrame,
    labels: pd.Series,
    cell_types: list[str],
    truth: np.ndarray,
    total_cells: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    counts = sample_counts_for_truth(truth, total_cells, rng)
    sampled_columns: list[str] = []
    replacement_by_type: dict[str, bool] = {}
    available_by_type: dict[str, int] = {}
    for cell_type, n_cells in zip(cell_types, counts):
        matching = labels.index[labels.astype(str) == cell_type].to_numpy()
        available_by_type[cell_type] = int(len(matching))
        if n_cells <= 0:
            replacement_by_type[cell_type] = False
            continue
        if len(matching) == 0:
            raise ValueError(f"no single cells available for cell type {cell_type!r}")
        replace = int(n_cells) > len(matching)
        chosen = rng.choice(matching, size=int(n_cells), replace=replace)
        sampled_columns.extend(chosen.tolist())
        replacement_by_type[cell_type] = bool(replace)

    if not sampled_columns:
        raise ValueError("sampled zero cells")

    sampled = sc_expr.loc[:, sampled_columns]
    bulk = sampled.sum(axis=1).to_numpy(dtype=float)
    sampled_truth = counts.astype(float) / float(np.sum(counts))
    info = {
        "sampled_cells": int(np.sum(counts)),
        "used_replacement": bool(any(replacement_by_type.values())),
        "replacement_by_type": replacement_by_type,
        "available_by_type": available_by_type,
    }
    return bulk, sampled_truth, info


def preprocess_metasort_inputs(
    reference: pd.DataFrame,
    bulk: pd.Series,
    cell_types: list[str],
    genes: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    common_genes = [gene for gene in genes if gene in reference.index and gene in bulk.index]
    signature = _normalize_columns_to_sum_one(
        reference.loc[common_genes, cell_types].to_numpy(dtype=float),
        name="signature",
    )
    bulk_vec = _normalize_vector_to_sum_one(
        bulk.loc[common_genes].to_numpy(dtype=float),
        name="bulk",
    )
    return _joint_gene_zscore(signature, bulk_vec)


def normalized_raw_inputs(
    reference: pd.DataFrame,
    bulk: pd.Series,
    cell_types: list[str],
    genes: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    common_genes = [gene for gene in genes if gene in reference.index and gene in bulk.index]
    signature = _normalize_columns_to_sum_one(
        reference.loc[common_genes, cell_types].to_numpy(dtype=float),
        name="signature",
    )
    bulk_vec = _normalize_vector_to_sum_one(
        bulk.loc[common_genes].to_numpy(dtype=float),
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


def run_methods(
    solver: MetaSortSolver,
    reference: pd.DataFrame,
    bulk: pd.Series,
    cell_types: list[str],
    genes: list[str],
) -> dict[str, tuple[np.ndarray, dict[str, Any]]]:
    raw_signature, raw_bulk = normalized_raw_inputs(reference, bulk, cell_types, genes)
    z_signature, z_bulk = _joint_gene_zscore(raw_signature, raw_bulk)

    results: dict[str, tuple[np.ndarray, dict[str, Any]]] = {}

    start = time.time()
    metasort_result = solver.solve(z_signature, z_bulk, cell_types=cell_types)
    results["MetaSort"] = (
        normalize_nonnegative(np.asarray(metasort_result.proportions, dtype=float)),
        {
            "runtime_sec": time.time() - start,
            "converged": bool(metasort_result.converged),
            "iterations": int(metasort_result.iterations),
            "error": "",
        },
    )

    start = time.time()
    results["NNLS_zscore"] = (
        solve_nnls_normalized(z_signature, z_bulk),
        {"runtime_sec": time.time() - start, "converged": True, "iterations": 1, "error": ""},
    )

    start = time.time()
    results["NNLS_raw"] = (
        solve_nnls_normalized(raw_signature, raw_bulk),
        {"runtime_sec": time.time() - start, "converged": True, "iterations": 1, "error": ""},
    )

    start = time.time()
    results["SimplexLS_raw"] = (
        solve_simplex_constrained_ls(raw_signature, raw_bulk),
        {"runtime_sec": time.time() - start, "converged": True, "iterations": 1, "error": ""},
    )

    start = time.time()
    results["DWLS_style_raw"] = (
        solve_dwls_style(raw_signature, raw_bulk),
        {"runtime_sec": time.time() - start, "converged": True, "iterations": 2, "error": ""},
    )
    return results


def process_case(
    case_dir: Path,
    output_dir: Path,
    solver: MetaSortSolver,
    n_mixtures: int,
    total_cells: int,
    marker_gene_count: int,
    nonmarker_gene_count: int,
    min_genes: int,
    seed: int,
    gene_selection: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    dataset = case_dir.parent.name
    tissue = case_dir.name
    case_name = f"{dataset}/{tissue}"

    reference = read_matrix(case_dir / "cellTypeExpr.txt")
    truth_df = read_matrix(case_dir / "bulkRatio.txt")
    sc_expr = read_matrix(case_dir / "singleCellExpr.txt")
    labels_df = pd.read_csv(case_dir / "singleCellLabels.txt", sep="\t")

    reference.columns = reference.columns.astype(str)
    reference.index = reference.index.astype(str)
    truth_df.index = truth_df.index.astype(str)
    sc_expr.index = sc_expr.index.astype(str)
    labels = pd.Series(labels_df.iloc[:, 0].astype(str).to_numpy(), index=sc_expr.columns.astype(str))

    cell_types = [
        cell_type
        for cell_type in reference.columns.astype(str).tolist()
        if cell_type in truth_df.index and cell_type in set(labels.to_numpy())
    ]
    if len(cell_types) < 2:
        raise ValueError(f"{case_name}: fewer than 2 usable cell types")

    case_seed = int(hashlib.sha256(case_name.encode("utf-8")).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed + case_seed)
    marker_genes = load_marker_genes(case_dir, cell_types)
    reference = reference.loc[reference.index.intersection(sc_expr.index)]
    selected_genes, marker_selected, nonmarker_selected = select_gene_subset(
        reference=reference,
        marker_genes=marker_genes,
        marker_gene_count=marker_gene_count,
        nonmarker_gene_count=nonmarker_gene_count,
        mode=gene_selection,
        rng=rng,
    )
    selected_genes = [gene for gene in selected_genes if gene in sc_expr.index]
    marker_selected = [gene for gene in marker_selected if gene in selected_genes]
    nonmarker_selected = [gene for gene in nonmarker_selected if gene in selected_genes]
    if len(selected_genes) < min_genes:
        raise ValueError(f"{case_name}: only {len(selected_genes)} selected genes")

    case_output = output_dir / "sampled_cases" / dataset / tissue
    case_output.mkdir(parents=True, exist_ok=True)
    pd.Series(selected_genes, name="gene").to_csv(case_output / "selected_genes.txt", sep="\t", index=False)
    pd.Series(marker_selected, name="gene").to_csv(case_output / "selected_marker_genes.txt", sep="\t", index=False)
    pd.Series(nonmarker_selected, name="gene").to_csv(
        case_output / "selected_nonmarker_genes.txt",
        sep="\t",
        index=False,
    )

    sc_expr = sc_expr.loc[selected_genes]
    reference = reference.loc[selected_genes]
    mixtures = truth_df.columns.astype(str).tolist()[:n_mixtures]

    metric_rows: list[dict[str, Any]] = []
    prediction_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    gene_rows = [
        {
            "dataset": dataset,
            "tissue": tissue,
            "case": case_name,
            "selected_genes": len(selected_genes),
            "selected_marker_genes": len(marker_selected),
            "selected_nonmarker_genes": len(nonmarker_selected),
            "available_marker_genes": len([gene for gene in marker_genes if gene in reference.index]),
            "reference_genes": int(reference.shape[0]),
            "cell_types": len(cell_types),
        }
    ]

    sampled_bulk_columns: dict[str, np.ndarray] = {}
    sampled_truth_columns: dict[str, np.ndarray] = {}
    for mixture_index, mixture in enumerate(mixtures, start=1):
        truth = normalize_nonnegative(truth_df.loc[cell_types, mixture].to_numpy(dtype=float))
        bulk_values, sampled_truth, sample_info = sample_bulk_from_single_cells(
            sc_expr=sc_expr,
            labels=labels,
            cell_types=cell_types,
            truth=truth,
            total_cells=total_cells,
            rng=rng,
        )
        bulk = pd.Series(bulk_values, index=selected_genes, name=mixture)
        sampled_bulk_columns[mixture] = bulk_values
        sampled_truth_columns[mixture] = sampled_truth

        sample_rows.append(
            {
                "dataset": dataset,
                "tissue": tissue,
                "case": case_name,
                "mixture": mixture,
                "mixture_index": mixture_index,
                "total_cells": sample_info["sampled_cells"],
                "used_replacement": sample_info["used_replacement"],
                "replacement_by_type": json.dumps(sample_info["replacement_by_type"], sort_keys=True),
                "available_by_type": json.dumps(sample_info["available_by_type"], sort_keys=True),
            }
        )

        method_results = run_methods(solver, reference, bulk, cell_types, selected_genes)
        for method, (prediction, extra) in method_results.items():
            metrics = score_prediction(prediction, sampled_truth)
            metric_rows.append(
                {
                    "dataset": dataset,
                    "tissue": tissue,
                    "case": case_name,
                    "mixture": mixture,
                    "mixture_index": mixture_index,
                    "method": method,
                    "n_genes": len(selected_genes),
                    "n_marker_genes": len(marker_selected),
                    "n_nonmarker_genes": len(nonmarker_selected),
                    "n_cell_types": len(cell_types),
                    **metrics,
                    **extra,
                }
            )
            for cell_type, truth_value, predicted_value in zip(cell_types, sampled_truth, prediction):
                prediction_rows.append(
                    {
                        "dataset": dataset,
                        "tissue": tissue,
                        "case": case_name,
                        "mixture": mixture,
                        "method": method,
                        "cell_type": cell_type,
                        "truth": float(truth_value),
                        "prediction": float(predicted_value),
                        "abs_error": float(abs(predicted_value - truth_value)),
                    }
                )

        best = min(
            (row for row in metric_rows if row["mixture"] == mixture),
            key=lambda row: row["l1"],
        )
        print(
            f"{case_name} {mixture_index}/{len(mixtures)} {mixture}: "
            f"best={best['method']} l1={best['l1']:.4f}",
            flush=True,
        )

    pd.DataFrame(sampled_bulk_columns, index=selected_genes).to_csv(case_output / "sampled_bulk.txt", sep="\t")
    pd.DataFrame(sampled_truth_columns, index=cell_types).to_csv(case_output / "sampled_bulkRatio.txt", sep="\t")
    reference.loc[selected_genes, cell_types].to_csv(case_output / "selected_cellTypeExpr.txt", sep="\t")

    return metric_rows, prediction_rows, sample_rows, gene_rows


def write_outputs(
    output_dir: Path,
    metrics: list[dict[str, Any]],
    predictions: list[dict[str, Any]],
    samples: list[dict[str, Any]],
    genes: list[dict[str, Any]],
) -> None:
    metrics_df = pd.DataFrame(metrics)
    predictions_df = pd.DataFrame(predictions)
    samples_df = pd.DataFrame(samples)
    genes_df = pd.DataFrame(genes)

    metrics_df.to_csv(output_dir / "metrics_long.csv", index=False)
    predictions_df.to_csv(output_dir / "predictions_long.csv", index=False)
    samples_df.to_csv(output_dir / "sampling_summary.csv", index=False)
    genes_df.to_csv(output_dir / "gene_selection_summary.csv", index=False)

    if metrics_df.empty:
        return

    case_summary = (
        metrics_df.groupby(["dataset", "tissue", "case", "method"], as_index=False)
        .agg(
            mixtures=("mixture", "count"),
            mean_l1=("l1", "mean"),
            median_l1=("l1", "median"),
            mean_accuracy=("accuracy", "mean"),
            mean_mae=("mae", "mean"),
            mean_rmse=("rmse", "mean"),
            mean_pearson=("pearson", "mean"),
            mean_runtime_sec=("runtime_sec", "mean"),
        )
        .sort_values(["case", "mean_l1", "mean_rmse"], kind="mergesort")
    )
    case_summary.to_csv(output_dir / "case_method_summary.csv", index=False)

    pivot = case_summary.pivot_table(
        index=["dataset", "tissue", "case"],
        columns="method",
        values="mean_l1",
        aggfunc="first",
    )
    baseline_cols = [col for col in pivot.columns if col != "MetaSort"]
    pivot["best_baseline_method"] = pivot[baseline_cols].idxmin(axis=1)
    pivot["best_baseline_l1"] = pivot[baseline_cols].min(axis=1)
    pivot["metasort_l1"] = pivot["MetaSort"]
    pivot["metasort_l1_advantage_vs_best_baseline"] = pivot["best_baseline_l1"] - pivot["metasort_l1"]
    pivot["metasort_wins_all_baselines"] = pivot["metasort_l1_advantage_vs_best_baseline"] > 0
    pivot = pivot.sort_values("metasort_l1_advantage_vs_best_baseline", ascending=False, kind="mergesort")
    pivot.reset_index().to_csv(output_dir / "case_l1_pivot.csv", index=False)
    pivot.loc[pivot["metasort_wins_all_baselines"]].reset_index().to_csv(
        output_dir / "top_metasort_winning_cases.csv",
        index=False,
    )

    mixture_pivot = metrics_df.pivot_table(
        index=["dataset", "tissue", "case", "mixture", "mixture_index"],
        columns="method",
        values="l1",
        aggfunc="first",
    )
    mixture_baseline_cols = [col for col in mixture_pivot.columns if col != "MetaSort"]
    mixture_pivot["best_baseline_method"] = mixture_pivot[mixture_baseline_cols].idxmin(axis=1)
    mixture_pivot["best_baseline_l1"] = mixture_pivot[mixture_baseline_cols].min(axis=1)
    mixture_pivot["metasort_l1"] = mixture_pivot["MetaSort"]
    mixture_pivot["metasort_l1_advantage_vs_best_baseline"] = (
        mixture_pivot["best_baseline_l1"] - mixture_pivot["metasort_l1"]
    )
    mixture_pivot = mixture_pivot.sort_values(
        "metasort_l1_advantage_vs_best_baseline",
        ascending=False,
        kind="mergesort",
    )
    mixture_pivot.reset_index().to_csv(output_dir / "mixture_l1_pivot.csv", index=False)


def main() -> int:
    args = parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"{args.output_dir} already exists and is not empty. Use --overwrite.")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except Exception:
        pass

    cfg = MetaSortConfig(
        convergence_tol=args.convergence_tol,
        averaging_old_weight=args.averaging_old_weight,
        lambda_hessian=args.lambda_hessian,
        lambda_residual=args.lambda_residual,
        use_sqrt_sphere_hessian=args.sqrt_sphere_hessian,
        lambda3=args.lambda3,
        lambda4=args.lambda4,
        meta_weight_baseline=args.meta_weight_baseline,
        meta_weight_floor=args.meta_weight_floor,
        normalize_meta_weight_mean=args.normalize_meta_weight_mean,
    )
    solver = MetaSortSolver(cfg)

    config_payload = vars(args).copy()
    config_payload["data_root"] = str(args.data_root)
    config_payload["output_dir"] = str(args.output_dir)
    config_payload["metasort_config"] = asdict(cfg)
    (args.output_dir / "config.json").write_text(json.dumps(config_payload, indent=2, sort_keys=True) + "\n")

    case_dirs = find_case_dirs(args.data_root, args.datasets, args.tissues)
    if args.max_cases is not None:
        case_dirs = case_dirs[: args.max_cases]
    print(f"cases: {len(case_dirs)}", flush=True)

    all_metrics: list[dict[str, Any]] = []
    all_predictions: list[dict[str, Any]] = []
    all_samples: list[dict[str, Any]] = []
    all_genes: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    for case_index, case_dir in enumerate(case_dirs, start=1):
        print(f"[{case_index}/{len(case_dirs)}] {case_dir.parent.name}/{case_dir.name}", flush=True)
        start = time.time()
        try:
            metrics, predictions, samples, genes = process_case(
                case_dir=case_dir,
                output_dir=args.output_dir,
                solver=solver,
                n_mixtures=args.n_mixtures,
                total_cells=args.total_cells,
                marker_gene_count=args.marker_genes,
                nonmarker_gene_count=args.nonmarker_genes,
                min_genes=args.min_genes,
                seed=args.seed,
                gene_selection=args.gene_selection,
            )
            all_metrics.extend(metrics)
            all_predictions.extend(predictions)
            all_samples.extend(samples)
            all_genes.extend(genes)
        except Exception as exc:
            failures.append(
                {
                    "dataset": case_dir.parent.name,
                    "tissue": case_dir.name,
                    "case": f"{case_dir.parent.name}/{case_dir.name}",
                    "error": repr(exc),
                    "elapsed_sec": time.time() - start,
                }
            )
            print(f"  FAILED: {exc!r}", flush=True)

        write_outputs(args.output_dir, all_metrics, all_predictions, all_samples, all_genes)
        pd.DataFrame(failures).to_csv(args.output_dir / "failures.csv", index=False)

    print(f"completed cases: {len(all_genes)}", flush=True)
    print(f"failed cases: {len(failures)}", flush=True)
    print(f"output_dir: {args.output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
