from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from metasort import MetaSortConfig, MetaSortSolver
from metasort.algorithm import (
    _joint_gene_zscore,
    _normalize_columns_to_sum_one,
    _normalize_vector_to_sum_one,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rerun MetaSort on a sampled bulk case and record L1 after every iteration."
    )
    parser.add_argument(
        "--case-dir",
        type=Path,
        default=Path(
            "outputs/search_mae_tongue_m10_try135/sampled_cases/Tabula-Sapiens/Tongue"
        ),
    )
    parser.add_argument("--dataset", default="Tabula-Sapiens")
    parser.add_argument("--tissue", default="Tongue")
    parser.add_argument("--mixture", default="Mixture10")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/fig4"),
    )
    parser.add_argument("--prefix", default="tongue_mixture10_try135_metasort")
    parser.add_argument(
        "--max-iter",
        type=int,
        default=None,
        help="Override MetaSortConfig.max_iter for this trace.",
    )
    parser.add_argument(
        "--ignore-convergence",
        action="store_true",
        help="Continue recording until --max-iter even after the convergence tolerance is met.",
    )
    parser.add_argument(
        "--initial-proportions",
        choices=["solver", "truth"],
        default="solver",
        help="Initial proportions for iteration 0. Use 'truth' to start directly from sampled_bulkRatio.",
    )
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
    return {
        "l1": l1,
        "mae": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "accuracy": 1.0 - l1 / 2.0,
    }


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def load_case(case_dir: Path, mixture: str) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    reference = read_matrix(case_dir / "selected_cellTypeExpr.txt")
    bulk_df = read_matrix(case_dir / "sampled_bulk.txt")
    truth_df = read_matrix(case_dir / "sampled_bulkRatio.txt")
    if mixture not in bulk_df.columns:
        raise ValueError(f"{mixture} not found in {case_dir / 'sampled_bulk.txt'}")
    if mixture not in truth_df.columns:
        raise ValueError(f"{mixture} not found in {case_dir / 'sampled_bulkRatio.txt'}")

    reference.index = reference.index.astype(str)
    reference.columns = reference.columns.astype(str)
    bulk_df.index = bulk_df.index.astype(str)
    truth_df.index = truth_df.index.astype(str)
    return reference, bulk_df[mixture], truth_df


def preprocess(
    reference: pd.DataFrame,
    bulk: pd.Series,
    cell_types: list[str],
    genes: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    signature = _normalize_columns_to_sum_one(
        reference.loc[genes, cell_types].to_numpy(dtype=float),
        name="signature",
    )
    bulk_vec = _normalize_vector_to_sum_one(
        bulk.loc[genes].to_numpy(dtype=float),
        name="bulk",
    )
    return _joint_gene_zscore(signature, bulk_vec)


def add_iteration_rows(
    rows: list[dict[str, Any]],
    proportion_rows: list[dict[str, Any]],
    *,
    dataset: str,
    tissue: str,
    mixture: str,
    iteration: int,
    phase: str,
    solution: np.ndarray,
    truth: np.ndarray,
    cell_types: list[str],
    change: float | None,
    meta_weights: np.ndarray,
    meta_metrics: dict[str, float] | None,
) -> None:
    proportions = normalize_nonnegative(solution)
    metrics = score_prediction(proportions, truth)
    row: dict[str, Any] = {
        "dataset": dataset,
        "tissue": tissue,
        "mixture": mixture,
        "iteration": iteration,
        "phase": phase,
        "change_l2": change,
        "n_cell_types": len(cell_types),
        "meta_weight_min": float(np.min(meta_weights)),
        "meta_weight_mean": float(np.mean(meta_weights)),
        "meta_weight_max": float(np.max(meta_weights)),
        **metrics,
    }
    if meta_metrics:
        for key, value in meta_metrics.items():
            row[f"meta_{key}"] = value
    rows.append(row)

    for cell_type, truth_value, pred_value in zip(cell_types, truth, proportions):
        proportion_rows.append(
            {
                "dataset": dataset,
                "tissue": tissue,
                "mixture": mixture,
                "iteration": iteration,
                "phase": phase,
                "cell_type": cell_type,
                "truth": float(truth_value),
                "prediction": float(pred_value),
                "abs_error": float(abs(pred_value - truth_value)),
            }
        )


def write_summary(
    path: Path,
    *,
    args: argparse.Namespace,
    cfg: MetaSortConfig,
    cell_types: list[str],
    genes: list[str],
    l1_df: pd.DataFrame,
    converged: bool,
) -> None:
    initial = l1_df.iloc[0]
    final = l1_df.iloc[-1]
    lines = [
        "# MetaSort iteration L1 trace",
        "",
        f"- Dataset: {args.dataset}",
        f"- Tissue: {args.tissue}",
        f"- Mixture: {args.mixture}",
        f"- Source case dir: `{args.case_dir}`",
        f"- Selected genes: {len(genes)}",
        f"- Cell types: {len(cell_types)} ({', '.join(cell_types)})",
        f"- Iterations recorded after initial: {int(final['iteration'])}",
        f"- Converged by change <= {cfg.convergence_tol}: {converged}",
        f"- Initial L1: {float(initial['l1']):.6f}",
        f"- Final L1: {float(final['l1']):.6f}",
        f"- Final per-celltype MAE: {float(final['mae']):.6f}",
        "",
        "## Files",
        "",
        f"- L1 trace: `{args.prefix}_iteration_l1.csv`",
        f"- Proportions: `{args.prefix}_iteration_proportions.csv`",
        f"- MetaSort config: `{args.prefix}_metasort_config.json`",
        f"- Genes: `{args.prefix}_genes.txt`",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass

    reference, bulk, truth_df = load_case(args.case_dir, args.mixture)
    cell_types = [cell_type for cell_type in reference.columns if cell_type in truth_df.index]
    if not cell_types:
        raise ValueError("No cell types overlap between reference and truth.")
    genes = reference.index.intersection(bulk.index).astype(str).tolist()
    if not genes:
        raise ValueError("No genes overlap between reference and bulk.")

    truth = normalize_nonnegative(truth_df.loc[cell_types, args.mixture].to_numpy(dtype=float))
    signature, bulk_vec = preprocess(reference, bulk, cell_types, genes)

    cfg = MetaSortConfig()
    if args.max_iter is not None:
        cfg.max_iter = args.max_iter
    solver = MetaSortSolver(cfg)
    confidence = np.ones(signature.shape[0], dtype=float)
    if args.initial_proportions == "truth":
        solution = truth.copy()
    else:
        solution = solver.solve_initial(signature, bulk_vec)
    meta_weights = np.full(signature.shape[0], cfg.meta_weight_baseline, dtype=float)

    rows: list[dict[str, Any]] = []
    proportion_rows: list[dict[str, Any]] = []
    add_iteration_rows(
        rows,
        proportion_rows,
        dataset=args.dataset,
        tissue=args.tissue,
        mixture=args.mixture,
        iteration=0,
        phase="initial_truth" if args.initial_proportions == "truth" else "initial",
        solution=solution,
        truth=truth,
        cell_types=cell_types,
        change=None,
        meta_weights=meta_weights,
        meta_metrics=None,
    )

    converged = False
    for iteration in range(1, cfg.max_iter + 1):
        new_solution, meta_step, meta_metrics = solver.solve_meta_weighted_step(
            signature=signature,
            bulk=bulk_vec,
            current_solution=solution,
            confidence_weights=confidence,
            prev_meta_weights=meta_weights,
        )
        solution_average = (new_solution + cfg.averaging_old_weight * solution) / float(
            cfg.averaging_old_weight + 1
        )
        change = float(np.linalg.norm(solution_average - solution))
        solution = solution_average
        meta_weights = meta_step
        add_iteration_rows(
            rows,
            proportion_rows,
            dataset=args.dataset,
            tissue=args.tissue,
            mixture=args.mixture,
            iteration=iteration,
            phase="iteration",
            solution=solution,
            truth=truth,
            cell_types=cell_types,
            change=change,
            meta_weights=meta_weights,
            meta_metrics=meta_metrics,
        )
        if change <= cfg.convergence_tol:
            converged = True
            if not args.ignore_convergence:
                break

    l1_df = pd.DataFrame(rows)
    proportions_df = pd.DataFrame(proportion_rows)

    l1_path = args.output_dir / f"{args.prefix}_iteration_l1.csv"
    proportions_path = args.output_dir / f"{args.prefix}_iteration_proportions.csv"
    config_path = args.output_dir / f"{args.prefix}_metasort_config.json"
    genes_path = args.output_dir / f"{args.prefix}_genes.txt"
    summary_path = args.output_dir / f"{args.prefix}_summary.md"

    l1_df.to_csv(l1_path, index=False)
    proportions_df.to_csv(proportions_path, index=False)
    config_path.write_text(json.dumps(asdict(cfg), indent=2, sort_keys=True) + "\n")
    genes_path.write_text("\n".join(genes) + "\n")
    write_summary(
        summary_path,
        args=args,
        cfg=cfg,
        cell_types=cell_types,
        genes=genes,
        l1_df=l1_df,
        converged=converged,
    )

    cols = ["iteration", "phase", "change_l2", "l1", "mae", "accuracy"]
    print(l1_df[cols].to_string(index=False), flush=True)
    print(f"output_dir: {args.output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
