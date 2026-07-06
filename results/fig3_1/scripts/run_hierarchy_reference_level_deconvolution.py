from __future__ import annotations

import json
import re
import sys
import time
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

BASE_DIR = Path("results/fig3_1")
OUT_DIR = BASE_DIR / "hierarchy_reference_level_analysis"
VERSION2_CONFIG_PATH = BASE_DIR / "version_2" / "MetaSort_jk_avg2p5_tol0005_config.csv"

sys.path.insert(0, str(BASE_DIR))

from run_part2_metasort_jk_original_with_weights import (  # noqa: E402
    HIERARCHY_DIR,
    N_MIXTURES,
    TISSUES,
    MetaSortSolver,
    make_config,
    normalize_nonnegative,
    preprocess,
    score_prediction,
)


METHOD_PREFIX = "MetaSort_jk_ref_level_avg2p5_tol0005"
COARSE_LEVEL = "coarse"
LEAF_LEVEL = "leaf"
LEVELS = (COARSE_LEVEL, LEAF_LEVEL)


def clean_node_name(name: str) -> str:
    return re.sub(r"^L\d+:\s*", "", name).replace(" / ", "__")


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    return data


def leaf_names(node: dict[str, Any]) -> list[str]:
    children = node.get("children") or []
    if not children:
        return [str(node["name"])]
    leaves: list[str] = []
    for child in children:
        leaves.extend(leaf_names(child))
    return leaves


def nodes_at_depth(root: dict[str, Any], depth: int) -> list[dict[str, Any]]:
    if depth < 1:
        raise ValueError("depth must be >= 1")
    current = root.get("children") or []
    for _ in range(1, depth):
        next_nodes: list[dict[str, Any]] = []
        for node in current:
            children = node.get("children") or []
            next_nodes.extend(children if children else [node])
        current = next_nodes
    return current


def hierarchy_mappings(hierarchy: dict[str, Any], fine_cell_types: list[str]) -> dict[str, dict[str, list[str]]]:
    mappings: dict[str, dict[str, list[str]]] = {}

    mapping: dict[str, list[str]] = {}
    used_names: set[str] = set()
    for node in hierarchy.get("children") or []:
        name = clean_node_name(str(node["name"]))
        if name in used_names:
            name = f"{name}_{len(used_names) + 1}"
        used_names.add(name)
        leaves = [leaf for leaf in leaf_names(node) if leaf in fine_cell_types]
        if leaves:
            mapping[name] = leaves
    mappings[COARSE_LEVEL] = mapping

    mappings[LEAF_LEVEL] = {cell_type: [cell_type] for cell_type in fine_cell_types}
    return mappings


def leaf_weights_from_labels(case_dir: Path, fine_cell_types: list[str]) -> pd.Series:
    labels = pd.read_csv(case_dir / "singleCellLabels.txt", sep="\t")
    if labels.shape[1] != 1:
        return pd.Series(1.0, index=fine_cell_types)
    counts = labels.iloc[:, 0].astype(str).value_counts()
    weights = counts.reindex(fine_cell_types).fillna(0.0).astype(float)
    if float(weights.sum()) <= 0.0:
        weights[:] = 1.0
    return weights


def reference_from_single_cell(
    single_cell_df: pd.DataFrame,
    labels: pd.Series,
    mapping: dict[str, list[str]],
) -> pd.DataFrame:
    label_to_node = {
        leaf: node_name
        for node_name, leaves in mapping.items()
        for leaf in leaves
    }
    grouped_labels = labels.astype(str).map(label_to_node)
    keep = grouped_labels.notna()
    if not bool(keep.any()):
        raise ValueError("No single cells matched the requested hierarchy mapping.")

    # Grouping cells by the requested hierarchy node gives a reference profile
    # directly from the single-cell matrix at that resolution.
    reference = single_cell_df.loc[:, keep.to_numpy()].T.groupby(grouped_labels.loc[keep].to_numpy()).mean().T
    ordered_columns = [node_name for node_name in mapping if node_name in reference.columns]
    return reference.loc[:, ordered_columns]


def single_cell_counts(
    labels: pd.Series,
    mapping: dict[str, list[str]],
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for node_name, leaves in mapping.items():
        counts[node_name] = int(labels.astype(str).isin(leaves).sum())
    return counts


def aggregate_truth(truth_df: pd.DataFrame, mapping: dict[str, list[str]]) -> pd.DataFrame:
    columns: dict[str, pd.Series] = {}
    for node_name, leaves in mapping.items():
        columns[node_name] = truth_df.loc[leaves].sum(axis=0)
    return pd.DataFrame(columns).T


def partition_key(mapping: dict[str, list[str]]) -> tuple[tuple[str, ...], ...]:
    return tuple(tuple(leaves) for leaves in mapping.values())


def make_test_config():
    config = make_config(hierarchical=False)
    config_row = pd.read_csv(VERSION2_CONFIG_PATH).iloc[0].to_dict()
    config_fields = {field.name: field.type for field in fields(config)}
    bool_fields = {
        field.name
        for field in fields(config)
        if field.type is bool or str(field.type) == "bool"
    }
    int_fields = {
        field.name
        for field in fields(config)
        if field.type is int or str(field.type) == "int"
    }

    for key in config_fields:
        if key not in config_row or pd.isna(config_row[key]):
            continue
        value = config_row[key]
        if key in bool_fields:
            value = str(value).strip().lower() == "true"
        elif key in int_fields:
            value = int(value)
        else:
            value = float(value)
        setattr(config, key, value)
    return config


def run_level(
    tissue: str,
    level: str,
    signature_df: pd.DataFrame,
    bulk_df: pd.DataFrame,
    truth_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mixtures = bulk_df.columns.astype(str).to_list()[:N_MIXTURES]
    cell_types = signature_df.columns.astype(str).to_list()
    solver = MetaSortSolver(make_test_config())
    result = pd.DataFrame(index=mixtures, columns=cell_types, dtype=float)
    metric_rows: list[dict[str, object]] = []
    start = time.time()

    for idx, mixture in enumerate(mixtures, start=1):
        signature, bulk, genes = preprocess(signature_df, bulk_df, mixture)
        truth = normalize_nonnegative(truth_df.loc[cell_types, mixture].to_numpy(dtype=float))
        solved = solver.solve(signature, bulk, cell_types=cell_types)
        pred = normalize_nonnegative(np.asarray(solved.proportions, dtype=float))
        result.loc[mixture, :] = pred
        metric_rows.append(
            {
                "tissue": tissue,
                "level": level,
                "method": METHOD_PREFIX,
                "mixture": mixture,
                "mixture_index": idx,
                "n_genes": len(genes),
                "n_reference_celltypes": len(cell_types),
                "iterations": int(solved.iterations),
                **score_prediction(pred, truth),
            }
        )
        print(f"{METHOD_PREFIX} {tissue} {level} {idx}/{len(mixtures)} {mixture}", flush=True)

    seconds = round(time.time() - start, 3)
    print(f"{METHOD_PREFIX} {tissue} {level} finished in {seconds}s", flush=True)
    return result, pd.DataFrame(metric_rows)


def aggregate_prediction(pred_df: pd.DataFrame, mapping: dict[str, list[str]]) -> pd.DataFrame:
    columns: dict[str, pd.Series] = {}
    for node_name, leaves in mapping.items():
        columns[node_name] = pred_df.loc[:, leaves].sum(axis=1)
    return pd.DataFrame(columns, index=pred_df.index)


def consistency_rows(
    tissue: str,
    level: str,
    direct_pred: pd.DataFrame,
    leaf_pred: pd.DataFrame,
    mapping: dict[str, list[str]],
) -> list[dict[str, object]]:
    leaf_aggregated = aggregate_prediction(leaf_pred, mapping)
    rows: list[dict[str, object]] = []
    for idx, mixture in enumerate(direct_pred.index.astype(str), start=1):
        direct = normalize_nonnegative(direct_pred.loc[mixture].to_numpy(dtype=float))
        aggregated = normalize_nonnegative(leaf_aggregated.loc[mixture, direct_pred.columns].to_numpy(dtype=float))
        diff = direct - aggregated
        l1 = float(np.sum(np.abs(diff)))
        rows.append(
            {
                "tissue": tissue,
                "level": level,
                "comparison": "direct_ref_level_vs_leaf_aggregated",
                "mixture": mixture,
                "mixture_index": idx,
                "L1": l1,
                "Consistency": 1.0 - l1 / 2.0,
                "MAE": float(np.mean(np.abs(diff))),
                "RMSE": float(np.sqrt(np.mean(diff**2))),
            }
        )
    return rows


def relabel_prediction(
    pred: pd.DataFrame,
    source_mapping: dict[str, list[str]],
    target_mapping: dict[str, list[str]],
) -> pd.DataFrame:
    source_by_leaves = {
        tuple(leaves): node_name
        for node_name, leaves in source_mapping.items()
    }
    columns = {
        target_node: pred.loc[:, source_by_leaves[tuple(leaves)]]
        for target_node, leaves in target_mapping.items()
    }
    return pd.DataFrame(columns, index=pred.index)


def load_complete_outputs(
    pred_path: Path,
    metrics_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    if not (pred_path.exists() and metrics_path.exists()):
        return None
    pred = pd.read_csv(pred_path, sep="\t", index_col=0)
    metrics = pd.read_csv(metrics_path, sep="\t")
    if pred.shape[0] != N_MIXTURES or metrics.shape[0] != N_MIXTURES:
        return None
    return pred, metrics


def summarize(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    value_cols = ["Accuracy", "MAE", "RMSE", "L1"]
    available = [col for col in value_cols if col in df.columns]
    summary = df.groupby(group_cols, as_index=False)[available].agg(["mean", "median", "std"])
    summary.columns = [
        "_".join(col).rstrip("_") if isinstance(col, tuple) else col
        for col in summary.columns.to_flat_index()
    ]
    return summary


def run_tissue(tissue: str) -> tuple[list[pd.DataFrame], list[pd.DataFrame], list[dict[str, object]]]:
    case_dir = BASE_DIR / tissue
    tissue_out = OUT_DIR / tissue
    tissue_out.mkdir(parents=True, exist_ok=True)

    single_cell_df = pd.read_csv(case_dir / "singleCellExpr.txt", sep="\t", index_col=0)
    labels = pd.read_csv(case_dir / "singleCellLabels.txt", sep="\t").iloc[:, 0].astype(str)
    bulk_df = pd.read_csv(case_dir / "bulk.txt", sep="\t", index_col=0)
    truth_df = pd.read_csv(case_dir / "bulkRatio.txt", sep="\t", index_col=0)
    truth_df.index = truth_df.index.astype(str)
    fine_cell_types = truth_df.index.astype(str).to_list()
    if single_cell_df.shape[1] != len(labels):
        raise ValueError(
            f"{tissue}: singleCellExpr has {single_cell_df.shape[1]} cells but "
            f"singleCellLabels has {len(labels)} labels"
        )
    hierarchy = load_json(HIERARCHY_DIR / f"{tissue}_hierarchy.json")
    mappings = hierarchy_mappings(hierarchy, fine_cell_types)

    pred_by_level: dict[str, pd.DataFrame] = {}
    metric_tables: list[pd.DataFrame] = []
    consistency_tables: list[pd.DataFrame] = []
    mapping_rows: list[dict[str, object]] = []
    cached_by_partition: dict[
        tuple[tuple[str, ...], ...],
        tuple[str, dict[str, list[str]], pd.DataFrame, pd.DataFrame],
    ] = {}

    for level in LEVELS:
        mapping = mappings[level]
        level_signature = reference_from_single_cell(single_cell_df, labels, mapping)
        level_truth = aggregate_truth(truth_df, mapping)
        counts = single_cell_counts(labels, mapping)
        for node_name, leaves in mapping.items():
            mapping_rows.append(
                {
                    "tissue": tissue,
                    "level": level,
                    "node": node_name,
                    "n_leaf_celltypes": len(leaves),
                    "leaf_celltypes": "|".join(leaves),
                    "n_single_cells": counts[node_name],
                }
            )

        level_signature.to_csv(tissue_out / f"{tissue}_{level}_reference_signature.tsv", sep="\t")
        level_truth.to_csv(tissue_out / f"{tissue}_{level}_truth.tsv", sep="\t")
        pred_path = tissue_out / f"{tissue}_{level}_prediction.tsv"
        metrics_path = tissue_out / f"{tissue}_{level}_metrics.tsv"
        key = partition_key(mapping)
        complete_outputs = load_complete_outputs(pred_path, metrics_path)
        if complete_outputs is not None:
            pred, metrics = complete_outputs
            print(f"{METHOD_PREFIX} {tissue} {level} loaded existing outputs", flush=True)
        elif key in cached_by_partition:
            source_level, source_mapping, source_pred, source_metrics = cached_by_partition[key]
            pred = relabel_prediction(source_pred, source_mapping, mapping)
            metrics = source_metrics.copy()
            metrics["level"] = level
            metrics["n_reference_celltypes"] = len(mapping)
            print(
                f"{METHOD_PREFIX} {tissue} {level} reused {source_level} results "
                "because the hierarchy partition is identical",
                flush=True,
            )
        else:
            pred, metrics = run_level(tissue, level, level_signature, bulk_df, level_truth)
        pred_by_level[level] = pred
        pred.to_csv(pred_path, sep="\t")
        metrics.to_csv(metrics_path, sep="\t", index=False)
        metric_tables.append(metrics)
        cached_by_partition[key] = (level, mapping, pred, metrics)

    leaf_pred = pred_by_level[LEAF_LEVEL]
    rows = consistency_rows(
        tissue,
        COARSE_LEVEL,
        pred_by_level[COARSE_LEVEL],
        leaf_pred,
        mappings[COARSE_LEVEL],
    )
    table = pd.DataFrame(rows)
    table.to_csv(tissue_out / f"{tissue}_{COARSE_LEVEL}_consistency_vs_leaf.tsv", sep="\t", index=False)
    consistency_tables.append(table)

    return metric_tables, consistency_tables, mapping_rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_metrics: list[pd.DataFrame] = []
    all_consistency: list[pd.DataFrame] = []
    all_mapping_rows: list[dict[str, object]] = []

    for tissue in TISSUES:
        metrics, consistency, mapping_rows = run_tissue(tissue)
        all_metrics.extend(metrics)
        all_consistency.extend(consistency)
        all_mapping_rows.extend(mapping_rows)

    metrics_df = pd.concat(all_metrics, ignore_index=True)
    consistency_df = pd.concat(all_consistency, ignore_index=True)
    mapping_df = pd.DataFrame(all_mapping_rows)

    metrics_df.to_csv(OUT_DIR / "hierarchy_reference_level_metrics_long.tsv", sep="\t", index=False)
    consistency_df.to_csv(OUT_DIR / "hierarchy_reference_level_consistency_vs_leaf_long.tsv", sep="\t", index=False)
    mapping_df.to_csv(OUT_DIR / "hierarchy_reference_level_node_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            asdict(make_test_config())
            | {
                "n_mixtures": N_MIXTURES,
                "config_source": str(VERSION2_CONFIG_PATH),
                "reference_source": "singleCellExpr.txt grouped by hierarchy level and singleCellLabels.txt",
                "script": str(Path(__file__)),
            }
        ]
    ).to_csv(OUT_DIR / "hierarchy_reference_level_version2_config.csv", index=False)

    metric_summary = summarize(metrics_df, ["tissue", "level", "n_reference_celltypes"])
    consistency_summary = summarize(consistency_df.rename(columns={"Consistency": "Accuracy"}), ["tissue", "level"])
    metric_summary.to_csv(OUT_DIR / "hierarchy_reference_level_metric_summary.tsv", sep="\t", index=False)
    consistency_summary.to_csv(OUT_DIR / "hierarchy_reference_level_consistency_summary.tsv", sep="\t", index=False)

    print("\nAccuracy summary:", flush=True)
    print(metric_summary.to_string(index=False), flush=True)
    print("\nConsistency vs leaf-aggregated summary:", flush=True)
    print(consistency_summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
