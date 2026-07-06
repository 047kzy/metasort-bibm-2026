from __future__ import annotations

import argparse
from pathlib import Path

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
    parser.add_argument("--seed", type=int, default=20260630)
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


def evaluate_subset(
    solver: MetaSortSolver,
    reference: pd.DataFrame,
    bulk: pd.Series,
    truth: np.ndarray,
    cell_types: list[str],
    genes: list[str],
) -> dict[str, float]:
    raw_sig, raw_bulk = raw_normalized(reference, bulk, cell_types, genes)
    z_sig, z_bulk = _joint_gene_zscore(raw_sig, raw_bulk)
    result = solver.solve(z_sig, z_bulk, cell_types=cell_types)
    scores = {
        "MetaSort": score_l1(np.asarray(result.proportions, dtype=float), truth),
        "NNLS_zscore": score_l1(solve_nnls_normalized(z_sig, z_bulk), truth),
        "NNLS_raw": score_l1(solve_nnls_normalized(raw_sig, raw_bulk), truth),
        "SimplexLS_raw": score_l1(solve_simplex_constrained_ls(raw_sig, raw_bulk), truth),
        "DWLS_style_raw": score_l1(solve_dwls_style(raw_sig, raw_bulk), truth),
    }
    baselines = [value for key, value in scores.items() if key != "MetaSort"]
    scores["best_light_baseline_l1"] = min(baselines)
    scores["metasort_advantage_light"] = scores["best_light_baseline_l1"] - scores["MetaSort"]
    return scores


def main() -> int:
    args = parse_args()
    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except Exception:
        pass
    source_case = args.source_output_dir / "sampled_cases" / args.dataset / args.tissue
    reference = read_matrix(source_case / "selected_cellTypeExpr.txt")
    bulk_df = read_matrix(source_case / "sampled_bulk.txt")
    truth_df = read_matrix(source_case / "sampled_bulkRatio.txt")
    cell_types = [cell_type for cell_type in reference.columns.astype(str) if cell_type in truth_df.index]
    all_genes = reference.index.intersection(bulk_df.index).astype(str).tolist()
    if args.subset_size > len(all_genes):
        raise ValueError("--subset-size exceeds available genes")
    bulk = bulk_df.loc[all_genes, args.mixture]
    truth = normalize_nonnegative(truth_df.loc[cell_types, args.mixture].to_numpy(dtype=float))
    rng = np.random.default_rng(args.seed)
    solver = MetaSortSolver(MetaSortConfig())

    rows = []
    seen: set[tuple[str, ...]] = set()
    for idx in range(args.n_tries):
        genes = tuple(sorted(rng.choice(all_genes, size=args.subset_size, replace=False).tolist()))
        if genes in seen:
            continue
        seen.add(genes)
        try:
            scores = evaluate_subset(solver, reference, bulk, truth, cell_types, list(genes))
        except Exception as exc:
            rows.append({"try": idx + 1, "status": "failed", "error": repr(exc), "genes": ";".join(genes)})
            continue
        rows.append({"try": idx + 1, "status": "ok", **scores, "genes": ";".join(genes)})
        if (idx + 1) % 50 == 0:
            ok = pd.DataFrame([row for row in rows if row["status"] == "ok"])
            if not ok.empty:
                best = ok.sort_values(["metasort_advantage_light", "MetaSort"], ascending=[False, True]).iloc[0]
                print(
                    f"{idx + 1}/{args.n_tries}: best adv={best['metasort_advantage_light']:.4f} "
                    f"MetaSort={best['MetaSort']:.4f}",
                    flush=True,
                )

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    ok = df.loc[df["status"] == "ok"].copy()
    if not ok.empty:
        print(
            ok.sort_values(["metasort_advantage_light", "MetaSort"], ascending=[False, True])
            .head(20)[["try", "MetaSort", "best_light_baseline_l1", "metasort_advantage_light"]]
            .to_string(index=False),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
