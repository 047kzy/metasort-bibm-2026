from __future__ import annotations

import csv
import json
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

METASORT_SOURCE = os.environ.get("METASORT_JK_PATH")
if METASORT_SOURCE:
    sys.path.insert(0, METASORT_SOURCE)

from metasort import (  # noqa: E402
    HierarchicalMetaSortConfig,
    HierarchicalMetaSortSolver,
    HierarchyNode,
    MetaSortConfig,
    MetaSortSolver,
)
from metasort.algorithm import (  # noqa: E402
    _joint_gene_zscore,
    _normalize_columns_to_sum_one,
    _normalize_vector_to_sum_one,
)


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = BASE_DIR / "data" / "inputs"
RESULTS_DIR = BASE_DIR / "results"
HIERARCHY_DIR = BASE_DIR / "data" / "hierarchies"
TISSUES = ["Salivary_Gland", "Thymus"]
N_MIXTURES = 30


class WeightRecordingMetaSortSolver(MetaSortSolver):
    def __init__(self, config: MetaSortConfig | None = None) -> None:
        super().__init__(config)
        self.last_meta_weights: np.ndarray | None = None
        self.last_meta_metrics: dict[str, float] = {}

    def solve_meta_weighted_step(self, *args: Any, **kwargs: Any) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
        solution, meta_weights, meta_metrics = super().solve_meta_weighted_step(*args, **kwargs)
        self.last_meta_weights = np.asarray(meta_weights, dtype=float).copy()
        self.last_meta_metrics = dict(meta_metrics)
        return solution, meta_weights, meta_metrics


class WeightRecordingHierarchicalMetaSortSolver(HierarchicalMetaSortSolver):
    def __init__(self, config: HierarchicalMetaSortConfig | None = None) -> None:
        super().__init__(config)
        self.base_solver = WeightRecordingMetaSortSolver(self.config)
        self.weight_records: list[dict[str, object]] = []

    def _solve_stage(
        self,
        signature: np.ndarray,
        bulk: np.ndarray,
        nodes: list[HierarchyNode],
        cell_type_to_index: dict[str, int],
    ) -> tuple[np.ndarray, Any, dict[str, float]]:
        selected_genes = self._select_stage_genes(
            signature=signature,
            visible_nodes=nodes,
            cell_type_to_index=cell_type_to_index,
        )
        stage_signature = self._aggregate_signature(
            signature=signature[selected_genes, :],
            nodes=nodes,
            cell_type_to_index=cell_type_to_index,
        )
        stage_bulk = bulk[selected_genes]
        stage_names = [node.name for node in nodes]
        result = self.base_solver.solve(stage_signature, stage_bulk, cell_types=stage_names)
        meta_weights = self.base_solver.last_meta_weights
        if meta_weights is None:
            raise RuntimeError("MetaSort stage finished without recorded meta weights.")
        raw_props = {
            node.name: float(prop)
            for node, prop in zip(nodes, result.proportions)
        }
        self.weight_records.append(
            {
                "solve_index": len(self.weight_records) + 1,
                "visible_nodes": "|".join(stage_names),
                "selected_genes": np.asarray(selected_genes, dtype=int).copy(),
                "meta_weights": meta_weights.copy(),
                "meta_metrics": dict(self.base_solver.last_meta_metrics),
                "iterations": int(result.iterations),
                "converged": bool(result.converged),
            }
        )
        return selected_genes, result, raw_props


def normalize_nonnegative(values: np.ndarray) -> np.ndarray:
    values = np.clip(np.asarray(values, dtype=float), 0.0, None)
    total = float(np.sum(values))
    if total <= 0.0 or not np.isfinite(total):
        return np.full(values.shape, 1.0 / len(values), dtype=float)
    return values / total


def hierarchy_node_from_dict(data: Any) -> HierarchyNode:
    if not isinstance(data, dict):
        raise ValueError("Hierarchy JSON nodes must be objects.")
    name = data.get("name")
    if not isinstance(name, str) or not name:
        raise ValueError("Every hierarchy node must have a non-empty string name.")

    children_data = data.get("children") or []
    if not isinstance(children_data, list):
        raise ValueError(f"Hierarchy node {name} has non-list children.")
    children = [hierarchy_node_from_dict(child_data) for child_data in children_data]

    if children:
        inferred_cell_types = [
            cell_type
            for child in children
            for cell_type in child.cell_types
        ]
        provided_cell_types = data.get("cell_types")
        if provided_cell_types is not None and set(provided_cell_types) != set(inferred_cell_types):
            raise ValueError(f"Hierarchy node {name} cell_types do not match its children.")
        cell_types = inferred_cell_types
    else:
        provided_cell_types = data.get("cell_types")
        if provided_cell_types is None:
            cell_types = [name]
        elif isinstance(provided_cell_types, list) and provided_cell_types == [name]:
            cell_types = [name]
        else:
            raise ValueError(f"Leaf hierarchy node {name} has invalid cell_types.")

    return HierarchyNode(
        name=name,
        cell_types=cell_types,
        children=children,
        distance=float(data.get("distance", 0.0)),
    )


def load_hierarchy(path: Path) -> HierarchyNode:
    return hierarchy_node_from_dict(json.loads(path.read_text()))


def preprocess(
    signature_df: pd.DataFrame,
    bulk_df: pd.DataFrame,
    mixture_name: str,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    common_genes = signature_df.index.intersection(bulk_df.index)
    if len(common_genes) == 0:
        raise ValueError("No shared genes between signature and bulk.")
    signature = _normalize_columns_to_sum_one(
        signature_df.loc[common_genes].to_numpy(dtype=float),
        name="signature",
    )
    bulk = _normalize_vector_to_sum_one(
        bulk_df.loc[common_genes, mixture_name].to_numpy(dtype=float),
        name="bulk",
    )
    signature, bulk = _joint_gene_zscore(signature, bulk)
    return signature, bulk, common_genes.to_list()


def score_prediction(pred: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    pred = normalize_nonnegative(pred)
    truth = normalize_nonnegative(truth)
    diff = pred - truth
    l1 = float(np.sum(np.abs(diff)))
    return {
        "L1": l1,
        "Accuracy": 1.0 - l1 / 2.0,
        "MAE": float(np.mean(np.abs(diff))),
        "RMSE": float(np.sqrt(np.mean(diff**2))),
    }


def make_config(hierarchical: bool) -> MetaSortConfig:
    config_class = HierarchicalMetaSortConfig if hierarchical else MetaSortConfig
    return config_class(
        convergence_tol=0.02,
        averaging_old_weight=4,
        lambda_hessian=0.03,
        lambda_residual=0.02,
        use_sqrt_sphere_hessian=False,
        lambda3=0.0001,
        lambda4=0.00001,
        meta_weight_baseline=10.0,
        meta_weight_floor=1.0,
        normalize_meta_weight_mean=True,
    )


def write_direct_weights(
    result_dir: Path,
    method: str,
    weight_columns: dict[str, np.ndarray],
    genes: list[str],
    rows: list[dict[str, object]],
) -> None:
    weights = pd.DataFrame(weight_columns, index=genes)
    weights.index.name = "gene"
    weights.to_csv(result_dir / f"{method}_gene_weights.tsv.gz", sep="\t", compression="gzip")
    pd.DataFrame(rows).to_csv(result_dir / f"{method}_weight_summary.csv", index=False)


def write_hierarchical_weights(
    result_dir: Path,
    method: str,
    rows: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
) -> None:
    pd.DataFrame(rows).to_csv(
        result_dir / f"{method}_gene_weights.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(summary_rows).to_csv(result_dir / f"{method}_weight_summary.csv", index=False)


def run_case(tissue: str, hierarchical: bool) -> dict[str, object]:
    method = "MetaSort_jk_hierarchy" if hierarchical else "MetaSort_jk"
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
    config = make_config(hierarchical)
    solver = (
        WeightRecordingHierarchicalMetaSortSolver(config)
        if hierarchical
        else WeightRecordingMetaSortSolver(config)
    )
    hierarchy = load_hierarchy(HIERARCHY_DIR / f"{tissue}_hierarchy.json") if hierarchical else None

    result = pd.DataFrame(index=mixtures, columns=cell_types, dtype=float)
    metric_rows: list[dict[str, object]] = []
    direct_weight_columns: dict[str, np.ndarray] = {}
    direct_weight_summary_rows: list[dict[str, object]] = []
    hierarchical_weight_rows: list[dict[str, object]] = []
    hierarchical_weight_summary_rows: list[dict[str, object]] = []
    genes_expected: list[str] | None = None
    start = time.time()

    for idx, mixture in enumerate(mixtures, start=1):
        signature, bulk, genes = preprocess(signature_df, bulk_df, mixture)
        truth = normalize_nonnegative(truth_df.loc[cell_types, mixture].to_numpy(dtype=float))
        if hierarchical:
            assert isinstance(solver, WeightRecordingHierarchicalMetaSortSolver)
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
                summary = {
                    "tissue": tissue,
                    "method": method,
                    "mixture": mixture,
                    "mixture_index": idx,
                    "solve_index": record["solve_index"],
                    "visible_nodes": record["visible_nodes"],
                    "n_genes": int(len(selected)),
                    "iterations": record["iterations"],
                    "converged": record["converged"],
                    "meta_weight_min": float(np.min(meta_weights)),
                    "meta_weight_mean": float(np.mean(meta_weights)),
                    "meta_weight_max": float(np.max(meta_weights)),
                    **metrics,
                }
                hierarchical_weight_summary_rows.append(summary)
                hierarchical_weight_rows.extend(
                    {
                        "tissue": tissue,
                        "method": method,
                        "mixture": mixture,
                        "mixture_index": idx,
                        "solve_index": record["solve_index"],
                        "visible_nodes": record["visible_nodes"],
                        "gene": genes[int(gene_idx)],
                        "meta_weight": float(weight),
                    }
                    for gene_idx, weight in zip(selected, meta_weights)
                )
        else:
            assert isinstance(solver, WeightRecordingMetaSortSolver)
            solver.last_meta_weights = None
            solved = solver.solve(signature, bulk, cell_types=cell_types)
            iterations = int(solved.iterations)
            if genes_expected is None:
                genes_expected = genes
            elif genes != genes_expected:
                raise ValueError(f"{tissue}: shared gene order changed between mixtures")
            if solver.last_meta_weights is None:
                raise RuntimeError(f"{tissue}/{mixture}: no meta weights were recorded")
            direct_weight_columns[mixture] = solver.last_meta_weights
            direct_weight_summary_rows.append(
                {
                    "tissue": tissue,
                    "method": method,
                    "mixture": mixture,
                    "mixture_index": idx,
                    "n_genes": len(genes),
                    "iterations": iterations,
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

        pred = normalize_nonnegative(np.asarray(solved.proportions, dtype=float))
        result.loc[mixture, :] = pred
        metric_rows.append(
            {
                "tissue": tissue,
                "method": method,
                "mixture": mixture,
                "mixture_index": idx,
                "n_genes": len(genes),
                "iterations": iterations,
                **score_prediction(pred, truth),
            }
        )
        print(f"{method} {tissue} {idx}/{len(mixtures)} {mixture}", flush=True)

    result_dir = RESULTS_DIR / tissue
    result_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(result_dir / f"{method}.tsv", sep="\t")
    pd.DataFrame(metric_rows).to_csv(result_dir / f"{method}_metrics.csv", index=False)
    if hierarchical:
        write_hierarchical_weights(
            result_dir=result_dir,
            method=method,
            rows=hierarchical_weight_rows,
            summary_rows=hierarchical_weight_summary_rows,
        )
    else:
        if genes_expected is None:
            raise ValueError(f"{tissue}: no direct MetaSort weights were recorded")
        write_direct_weights(
            result_dir=result_dir,
            method=method,
            weight_columns=direct_weight_columns,
            genes=genes_expected,
            rows=direct_weight_summary_rows,
        )
    return {
        "tissue": tissue,
        "method": method,
        "status": "OK",
        "rows": result.shape[0],
        "cols": result.shape[1],
        "seconds": round(time.time() - start, 3),
        "error": "",
    }


def main() -> None:
    rows: list[dict[str, object]] = []
    for tissue in TISSUES:
        for hierarchical in (False, True):
            method = "MetaSort_jk_hierarchy" if hierarchical else "MetaSort_jk"
            try:
                rows.append(run_case(tissue, hierarchical))
            except Exception as exc:
                rows.append(
                    {
                        "tissue": tissue,
                        "method": method,
                        "status": "ERROR",
                        "rows": "",
                        "cols": "",
                        "seconds": "",
                        "error": repr(exc),
                    }
                )
                print(f"{method} {tissue} ERROR: {exc!r}", flush=True)

    log_path = RESULTS_DIR / "metasort_jk_original_with_weights_run_log.tsv"
    with log_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    config = make_config(hierarchical=True)
    config_row = asdict(config) | {
        "n_mixtures": N_MIXTURES,
        "confidence_weights": "off",
        "output_dir": str(BASE_DIR),
        "source": METASORT_SOURCE or "installed metasort package",
        "script": str(Path(__file__)),
    }
    pd.DataFrame([config_row]).to_csv(
        RESULTS_DIR / "metasort_jk_original_with_weights_config.csv", index=False
    )
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
