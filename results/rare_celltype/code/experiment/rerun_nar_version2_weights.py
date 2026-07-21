from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from metasort import MetaSortConfig, MetaSortSolver
from run_nar30_first10_compare import load_case, normalize_nonnegative, preprocess_inputs


DEFAULT_RUN_DIR = Path("results/rare_celltype/work/metasort_benchmark")
DEFAULT_OUTPUT_DIR = Path("results/rare_celltype/work/metasort_weight_rerun")
DEFAULT_DATA_ROOT = Path("results/rare_celltype/data/inputs")
DEFAULT_TARGETS = (
    ("CELlxGENE", "Basal_Zone_Of_Heart"),
    ("Tabula-Sapiens", "Skin"),
)


class WeightRecordingSolver(MetaSortSolver):
    def __init__(self, config: MetaSortConfig) -> None:
        super().__init__(config)
        self.last_meta_weights: np.ndarray | None = None

    def solve_meta_weighted_step(self, *args, **kwargs):
        result = super().solve_meta_weighted_step(*args, **kwargs)
        self.last_meta_weights = np.asarray(result[1], dtype=float).copy()
        return result


def parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def load_config(config_file: Path) -> tuple[MetaSortConfig, dict[str, str]]:
    with config_file.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"{config_file} should contain exactly one config row")

    raw = rows[0]
    values: dict[str, Any] = {}
    defaults = MetaSortConfig()
    for field in fields(MetaSortConfig):
        if field.name not in raw or raw[field.name] == "":
            continue
        default_value = getattr(defaults, field.name)
        if isinstance(default_value, bool):
            values[field.name] = parse_bool(raw[field.name])
        elif isinstance(default_value, int):
            values[field.name] = int(float(raw[field.name]))
        elif isinstance(default_value, float):
            values[field.name] = float(raw[field.name])
        else:
            values[field.name] = raw[field.name]
    return MetaSortConfig(**values), raw


def score_vectors(pred: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    pred = normalize_nonnegative(pred)
    truth = normalize_nonnegative(truth)
    diff = pred - truth
    if len(pred) > 1 and np.std(pred) > 0.0 and np.std(truth) > 0.0:
        pearson = float(np.corrcoef(pred, truth)[0, 1])
    else:
        pearson = float("nan")
    return {
        "l1": float(np.sum(np.abs(diff))),
        "accuracy": float(1.0 - 0.5 * np.sum(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "mae": float(np.mean(np.abs(diff))),
        "pearson": pearson,
    }


def run_case(
    data_root: Path,
    output_dir: Path,
    saved_prediction_root: Path,
    dataset: str,
    tissue: str,
    cfg: MetaSortConfig,
    n_mixtures: int,
) -> dict[str, Any]:
    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except Exception:
        pass

    case_dir = data_root / dataset / tissue
    if not case_dir.exists() and dataset == "CELlxGENE":
        case_dir = data_root / "CELLxGENE" / tissue
    signature_df, bulk_df, truth_df = load_case(case_dir)
    cell_types = [cell_type for cell_type in signature_df.columns if cell_type in truth_df.index]
    if not cell_types:
        raise ValueError(f"{case_dir}: no common cell types between signature and truth")
    mixtures = bulk_df.columns.astype(str).to_list()[:n_mixtures]
    if not mixtures:
        raise ValueError(f"{case_dir}: no mixtures found")

    solver = WeightRecordingSolver(cfg)
    prediction_columns: dict[str, np.ndarray] = {}
    weight_columns: dict[str, np.ndarray] = {}
    metric_rows: list[dict[str, Any]] = []
    genes_expected: list[str] | None = None
    started = time.time()

    for mixture_index, mixture in enumerate(mixtures, start=1):
        signature, bulk, genes = preprocess_inputs(signature_df, bulk_df, mixture, cell_types)
        if genes_expected is None:
            genes_expected = genes
        elif genes != genes_expected:
            raise ValueError(f"{case_dir}: shared gene order changed between mixtures")

        solver.last_meta_weights = None
        mixture_started = time.time()
        result = solver.solve(signature, bulk, cell_types=cell_types)
        if solver.last_meta_weights is None:
            raise RuntimeError(f"{dataset}/{tissue}/{mixture}: no meta weights were recorded")

        props = normalize_nonnegative(np.asarray(result.proportions, dtype=float))
        prediction_columns[mixture] = props
        weight_columns[mixture] = solver.last_meta_weights
        truth = truth_df.loc[cell_types, mixture].to_numpy(dtype=float)
        metric_rows.append(
            {
                "dataset": dataset,
                "tissue": tissue,
                "case": f"{dataset}/{tissue}",
                "mixture": mixture,
                "mixture_index": mixture_index,
                "genes": len(genes),
                "cell_types": len(cell_types),
                "runtime_sec": time.time() - mixture_started,
                "converged": bool(result.converged),
                "iterations": int(result.iterations),
                "meta_weight_min": float(result.meta_weight_min),
                "meta_weight_mean": float(result.meta_weight_mean),
                "meta_weight_max": float(result.meta_weight_max),
                "total_weight_min": float(result.total_weight_min),
                "total_weight_mean": float(result.total_weight_mean),
                "total_weight_max": float(result.total_weight_max),
                **score_vectors(props, truth),
            }
        )
        if mixture_index % 20 == 0 or mixture_index == len(mixtures):
            print(f"{dataset}/{tissue}: {mixture_index}/{len(mixtures)} mixtures", flush=True)

    if genes_expected is None:
        raise ValueError(f"{case_dir}: no genes recorded")

    prediction = pd.DataFrame(prediction_columns, index=cell_types).T
    prediction_dir = output_dir / "predictions" / dataset / tissue
    prediction_dir.mkdir(parents=True, exist_ok=True)
    prediction_file = prediction_dir / "MetaSort.txt"
    prediction.to_csv(prediction_file, sep="\t")

    weights = pd.DataFrame(weight_columns, index=genes_expected)
    weights.index.name = "gene"
    weights_dir = output_dir / "gene_weight_matrices"
    weights_dir.mkdir(parents=True, exist_ok=True)
    weights_file = weights_dir / f"{dataset}__{tissue}__meta_weights.tsv.gz"
    weights.to_csv(weights_file, sep="\t", compression="gzip")

    metrics = pd.DataFrame(metric_rows)
    metrics_dir = output_dir / "mixture_metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    metrics_file = metrics_dir / f"{dataset}__{tissue}__mixture_metrics.csv"
    metrics.to_csv(metrics_file, index=False)

    max_saved_prediction_abs_difference = float("nan")
    saved_file = saved_prediction_root / dataset / tissue / "MetaSort.txt"
    if saved_file.exists():
        saved = pd.read_csv(saved_file, sep="\t", index_col=0)
        saved.index = saved.index.astype(str)
        saved.columns = saved.columns.astype(str)
        saved = saved.loc[prediction.index, prediction.columns]
        max_saved_prediction_abs_difference = float(
            np.max(np.abs(saved.to_numpy(dtype=float) - prediction.to_numpy(dtype=float)))
        )

    return {
        "dataset": dataset,
        "tissue": tissue,
        "case": f"{dataset}/{tissue}",
        "status": "ok",
        "mixtures": len(mixtures),
        "genes": len(genes_expected),
        "cell_types": len(cell_types),
        "elapsed_sec": time.time() - started,
        "prediction_file": str(prediction_file),
        "weights_file": str(weights_file),
        "metrics_file": str(metrics_file),
        "max_saved_prediction_abs_difference": max_saved_prediction_abs_difference,
    }


def parse_targets(values: list[str] | None) -> tuple[tuple[str, str], ...]:
    if not values:
        return DEFAULT_TARGETS
    targets: list[tuple[str, str]] = []
    for value in values:
        if "/" not in value:
            raise ValueError(f"Target must be DATASET/TISSUE: {value}")
        dataset, tissue = value.split("/", 1)
        targets.append((dataset, tissue))
    return tuple(dict.fromkeys(targets))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--targets", nargs="*", default=None, help="Optional DATASET/TISSUE list")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    cfg, raw_config = load_config(args.run_dir / "run_config.csv")
    n_mixtures = int(raw_config.get("n_mixtures") or 100)
    targets = parse_targets(args.targets)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        columns=[
            "dataset",
            "tissue",
            "case",
            "status",
            "mixtures",
            "genes",
            "cell_types",
            "elapsed_sec",
            "prediction_file",
            "weights_file",
            "metrics_file",
            "max_saved_prediction_abs_difference",
        ]
    ).to_csv(args.output_dir / "run_success.csv", index=False)
    pd.DataFrame(columns=["dataset", "tissue", "case", "status", "error"]).to_csv(
        args.output_dir / "run_failed.csv",
        index=False,
    )

    with (args.output_dir / "rerun_config.json").open("w") as handle:
        json.dump(
            {
                "data_root": str(args.data_root),
                "source_run_dir": str(args.run_dir),
                "output_dir": str(args.output_dir),
                "targets": [f"{dataset}/{tissue}" for dataset, tissue in targets],
                "n_mixtures": n_mixtures,
                "workers": args.workers,
                "metasort_config": asdict(cfg),
                "source_run_config": raw_config,
            },
            handle,
            indent=2,
        )
        handle.write("\n")

    success_rows: list[dict[str, Any]] = []
    failed_rows: list[dict[str, Any]] = []
    saved_prediction_root = args.run_dir / "predictions"

    def record_success(row: dict[str, Any]) -> None:
        success_rows.append(row)
        pd.DataFrame(success_rows).to_csv(args.output_dir / "run_success.csv", index=False)

    def record_failure(dataset: str, tissue: str, exc: BaseException) -> None:
        failed_rows.append(
            {
                "dataset": dataset,
                "tissue": tissue,
                "case": f"{dataset}/{tissue}",
                "status": "failed",
                "error": repr(exc),
            }
        )
        pd.DataFrame(failed_rows).to_csv(args.output_dir / "run_failed.csv", index=False)

    if args.workers <= 1:
        for index, (dataset, tissue) in enumerate(targets, start=1):
            print(f"[{index}/{len(targets)}] {dataset}/{tissue}", flush=True)
            try:
                record_success(
                    run_case(
                        args.data_root,
                        args.output_dir,
                        saved_prediction_root,
                        dataset,
                        tissue,
                        cfg,
                        n_mixtures,
                    )
                )
            except Exception as exc:
                record_failure(dataset, tissue, exc)
                print(f"FAILED {dataset}/{tissue}: {exc!r}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    run_case,
                    args.data_root,
                    args.output_dir,
                    saved_prediction_root,
                    dataset,
                    tissue,
                    cfg,
                    n_mixtures,
                ): (dataset, tissue)
                for dataset, tissue in targets
            }
            for completed, future in enumerate(as_completed(futures), start=1):
                dataset, tissue = futures[future]
                try:
                    row = future.result()
                    record_success(row)
                    print(f"[done {completed}/{len(targets)}] {dataset}/{tissue}", flush=True)
                except Exception as exc:
                    record_failure(dataset, tissue, exc)
                    print(f"FAILED {dataset}/{tissue}: {exc!r}", flush=True)

    metric_files = sorted((args.output_dir / "mixture_metrics").glob("*__mixture_metrics.csv"))
    if metric_files:
        pd.concat((pd.read_csv(path) for path in metric_files), ignore_index=True).to_csv(
            args.output_dir / "metasort_mixture_metrics.csv",
            index=False,
        )

    validation_columns = ["dataset", "tissue", "max_saved_prediction_abs_difference", "weights_file"]
    validation = pd.DataFrame(success_rows)
    if validation.empty:
        validation = pd.DataFrame(columns=validation_columns)
    else:
        validation = validation[validation_columns]
    validation.to_csv(args.output_dir / "rerun_prediction_validation.csv", index=False)
    print(f"Wrote rerun weights to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
