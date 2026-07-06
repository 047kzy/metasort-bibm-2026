from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


DEFAULT_DATA_ROOT = Path("data/benchmark_inputs")
DEFAULT_MUSIC_RSCRIPT = Path("<CONDA_ROOT>/bin/Rscript")
DEFAULT_LINDECONSEQ_RSCRIPT = Path("<CONDA_ROOT>/envs/r410_lindeconseq/bin/Rscript")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--script", type=Path, default=Path("scripts/run_sampled_bulk_r_baselines.R"))
    parser.add_argument("--music-rscript", type=Path, default=DEFAULT_MUSIC_RSCRIPT)
    parser.add_argument("--lindeconseq-rscript", type=Path, default=DEFAULT_LINDECONSEQ_RSCRIPT)
    parser.add_argument("--methods", nargs="+", default=["MuSiC", "LinDeconSeq"])
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def discover_cases(output_dir: Path) -> list[tuple[str, str, Path]]:
    sampled_root = output_dir / "sampled_cases"
    cases: list[tuple[str, str, Path]] = []
    for dataset_dir in sorted(sampled_root.iterdir()):
        if not dataset_dir.is_dir():
            continue
        for tissue_dir in sorted(dataset_dir.iterdir()):
            if not tissue_dir.is_dir():
                continue
            if (tissue_dir / "sampled_bulk.txt").exists():
                cases.append((dataset_dir.name, tissue_dir.name, tissue_dir))
    return cases


def rscript_for_method(args: argparse.Namespace, method: str) -> Path:
    if method == "MuSiC":
        return args.music_rscript
    if method == "LinDeconSeq":
        return args.lindeconseq_rscript
    raise ValueError(f"unsupported method: {method}")


def run_method(args: argparse.Namespace, dataset: str, tissue: str, case_dir: Path, method: str) -> None:
    out_file = case_dir / f"r_baseline_{method}_metrics_long.csv"
    if args.skip_existing and out_file.exists():
        print(f"[skip] {dataset}/{tissue} {method}", flush=True)
        return
    command = [
        str(rscript_for_method(args, method)),
        str(args.script),
        "--output-dir",
        str(args.output_dir),
        "--data-root",
        str(args.data_root),
        "--dataset",
        dataset,
        "--tissue",
        tissue,
        "--method",
        method,
    ]
    print(f"[run] {dataset}/{tissue} {method}", flush=True)
    subprocess.run(command, check=True)


def combine_outputs(output_dir: Path, methods: list[str]) -> None:
    base_metrics = pd.read_csv(output_dir / "metrics_long.csv")
    metric_frames = [base_metrics]
    summary_frames = []

    for case_dir in sorted((output_dir / "sampled_cases").glob("*/*")):
        for method in methods:
            metric_file = case_dir / f"r_baseline_{method}_metrics_long.csv"
            if metric_file.exists():
                metric_frames.append(pd.read_csv(metric_file))
            summary_file = case_dir / f"r_baseline_{method}_case_summary.csv"
            if summary_file.exists():
                summary_frames.append(pd.read_csv(summary_file))

    metrics = pd.concat(metric_frames, ignore_index=True)
    metrics.to_csv(output_dir / "metrics_long_with_music_lindeconseq.csv", index=False)

    case_summary = (
        metrics.groupby(["dataset", "tissue", "case", "method"], as_index=False)
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
    case_summary.to_csv(output_dir / "case_method_summary_with_music_lindeconseq.csv", index=False)

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
    pivot.reset_index().to_csv(output_dir / "case_l1_pivot_with_music_lindeconseq.csv", index=False)

    wins = pivot.loc[pivot["metasort_wins_all_baselines"]]
    wins.reset_index().to_csv(output_dir / "top_metasort_winning_cases_with_music_lindeconseq.csv", index=False)


def main() -> int:
    args = parse_args()
    cases = discover_cases(args.output_dir)
    if args.max_cases is not None:
        cases = cases[: args.max_cases]
    print(f"cases: {len(cases)}", flush=True)
    failures: list[dict[str, str]] = []
    for dataset, tissue, case_dir in cases:
        for method in args.methods:
            try:
                run_method(args, dataset, tissue, case_dir, method)
            except subprocess.CalledProcessError as exc:
                failures.append(
                    {
                        "dataset": dataset,
                        "tissue": tissue,
                        "method": method,
                        "error": repr(exc),
                    }
                )
                print(f"[failed] {dataset}/{tissue} {method}: {exc!r}", flush=True)
    if failures:
        pd.DataFrame(failures).to_csv(args.output_dir / "r_baseline_failures.csv", index=False)
    combine_outputs(args.output_dir, args.methods)
    print(f"failures: {len(failures)}", flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
