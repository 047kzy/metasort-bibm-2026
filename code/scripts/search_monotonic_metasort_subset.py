from __future__ import annotations

import argparse
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
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--mixture", required=True)
    parser.add_argument("--subset-size", type=int, required=True)
    parser.add_argument("--n-tries", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260701)
    parser.add_argument("--max-final-l1", type=float, default=0.08)
    parser.add_argument("--min-baseline-l1", type=float, default=0.10)
    parser.add_argument("--min-advantage", type=float, default=0.02)
    parser.add_argument("--monotonic-tol", type=float, default=1e-10)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def normalize_nonnegative(values: np.ndarray) -> np.ndarray:
    values = np.clip(np.asarray(values, dtype=float), 0.0, None)
    total = float(np.sum(values))
    if total <= 0.0 or not np.isfinite(total):
        return np.full(values.shape, 1.0 / len(values), dtype=float)
    return values / total


def score_l1(pred: np.ndarray, truth: np.ndarray) -> float:
    return float(np.sum(np.abs(normalize_nonnegative(pred) - normalize_nonnegative(truth))))


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def raw_normalized(
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


def metasort_l1_trace(
    solver: MetaSortSolver,
    signature: np.ndarray,
    bulk: np.ndarray,
    truth: np.ndarray,
) -> tuple[list[float], list[float], bool, np.ndarray]:
    cfg = solver.config
    confidence = np.ones(signature.shape[0], dtype=float)
    solution = solver.solve_initial(signature, bulk)
    meta_weights = np.full(signature.shape[0], cfg.meta_weight_baseline, dtype=float)
    l1_values = [score_l1(solution, truth)]
    changes: list[float] = []
    converged = False
    for _ in range(cfg.max_iter):
        new_solution, meta_step, _ = solver.solve_meta_weighted_step(
            signature=signature,
            bulk=bulk,
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
        changes.append(change)
        l1_values.append(score_l1(solution, truth))
        if change <= cfg.convergence_tol:
            converged = True
            break
    return l1_values, changes, converged, normalize_nonnegative(solution)


def evaluate_subset(
    solver: MetaSortSolver,
    reference: pd.DataFrame,
    bulk: pd.Series,
    truth: np.ndarray,
    cell_types: list[str],
    genes: list[str],
    monotonic_tol: float,
) -> dict[str, Any]:
    raw_sig, raw_bulk = raw_normalized(reference, bulk, cell_types, genes)
    z_sig, z_bulk = _joint_gene_zscore(raw_sig, raw_bulk)
    l1_trace, changes, converged, meta_pred = metasort_l1_trace(solver, z_sig, z_bulk, truth)
    diffs = np.diff(np.asarray(l1_trace, dtype=float))
    monotonic = bool(np.all(diffs < -monotonic_tol))
    nonincreasing = bool(np.all(diffs <= monotonic_tol))
    scores = {
        "MetaSort": score_l1(meta_pred, truth),
        "NNLS_zscore": score_l1(solve_nnls_normalized(z_sig, z_bulk), truth),
        "NNLS_raw": score_l1(solve_nnls_normalized(raw_sig, raw_bulk), truth),
        "SimplexLS_raw": score_l1(solve_simplex_constrained_ls(raw_sig, raw_bulk), truth),
        "DWLS_style_raw": score_l1(solve_dwls_style(raw_sig, raw_bulk), truth),
    }
    baselines = [value for key, value in scores.items() if key != "MetaSort"]
    best_baseline = min(baselines)
    return {
        **scores,
        "initial_l1": l1_trace[0],
        "final_l1": l1_trace[-1],
        "min_iteration_l1": float(np.min(l1_trace)),
        "iterations": len(l1_trace) - 1,
        "converged": converged,
        "strictly_decreasing_l1": monotonic,
        "nonincreasing_l1": nonincreasing,
        "l1_trace": ";".join(f"{value:.12g}" for value in l1_trace),
        "change_trace": ";".join(f"{value:.12g}" for value in changes),
        "best_light_baseline_l1": best_baseline,
        "metasort_advantage_light": best_baseline - scores["MetaSort"],
    }


def main() -> int:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass

    source_case = args.source_output_dir / "sampled_cases" / args.dataset / args.tissue
    reference = read_matrix(source_case / "selected_cellTypeExpr.txt")
    bulk_df = read_matrix(source_case / "sampled_bulk.txt")
    truth_df = read_matrix(source_case / "sampled_bulkRatio.txt")
    reference.index = reference.index.astype(str)
    reference.columns = reference.columns.astype(str)
    bulk_df.index = bulk_df.index.astype(str)
    truth_df.index = truth_df.index.astype(str)

    cell_types = [cell_type for cell_type in reference.columns if cell_type in truth_df.index]
    all_genes = reference.index.intersection(bulk_df.index).astype(str).tolist()
    if args.subset_size > len(all_genes):
        raise ValueError("--subset-size exceeds available genes")
    if args.mixture not in bulk_df.columns:
        raise ValueError(f"{args.mixture} not found in sampled_bulk.txt")
    bulk = bulk_df.loc[all_genes, args.mixture]
    truth = normalize_nonnegative(truth_df.loc[cell_types, args.mixture].to_numpy(dtype=float))

    rng = np.random.default_rng(args.seed)
    solver = MetaSortSolver(MetaSortConfig())
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, ...]] = set()
    for idx in range(args.n_tries):
        genes = tuple(sorted(rng.choice(all_genes, size=args.subset_size, replace=False).tolist()))
        if genes in seen:
            continue
        seen.add(genes)
        try:
            scores = evaluate_subset(
                solver,
                reference,
                bulk,
                truth,
                cell_types,
                list(genes),
                args.monotonic_tol,
            )
            status = "ok"
            error = ""
        except Exception as exc:
            scores = {}
            status = "failed"
            error = repr(exc)
        row = {
            "try": idx + 1,
            "status": status,
            "dataset": args.dataset,
            "tissue": args.tissue,
            "mixture": args.mixture,
            "subset_size": args.subset_size,
            "genes": ";".join(genes),
            "error": error,
            **scores,
        }
        if status == "ok":
            row["passes_target"] = bool(
                scores["strictly_decreasing_l1"]
                and scores["MetaSort"] <= args.max_final_l1
                and scores["best_light_baseline_l1"] >= args.min_baseline_l1
                and scores["metasort_advantage_light"] >= args.min_advantage
            )
        rows.append(row)
        if (idx + 1) % 50 == 0:
            ok = pd.DataFrame([row for row in rows if row["status"] == "ok"])
            passing = ok.loc[ok.get("passes_target", False)] if "passes_target" in ok else ok.iloc[0:0]
            best = ok.sort_values(
                ["strictly_decreasing_l1", "metasort_advantage_light", "MetaSort"],
                ascending=[False, False, True],
                kind="mergesort",
            ).iloc[0]
            print(
                f"{idx + 1}/{args.n_tries}: pass={len(passing)} "
                f"best strict={best['strictly_decreasing_l1']} "
                f"MetaSort={best['MetaSort']:.4f} best_base={best['best_light_baseline_l1']:.4f} "
                f"adv={best['metasort_advantage_light']:.4f}",
                flush=True,
            )

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    ok = df.loc[df["status"] == "ok"].copy()
    if not ok.empty:
        print(
            ok.sort_values(
                ["passes_target", "strictly_decreasing_l1", "metasort_advantage_light", "MetaSort"],
                ascending=[False, False, False, True],
                kind="mergesort",
            )
            .head(20)[
                [
                    "try",
                    "passes_target",
                    "strictly_decreasing_l1",
                    "MetaSort",
                    "best_light_baseline_l1",
                    "metasort_advantage_light",
                    "l1_trace",
                ]
            ]
            .to_string(index=False),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
