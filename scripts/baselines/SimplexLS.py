from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from metasort.algorithm import (
    _joint_gene_zscore,
    _normalize_columns_to_sum_one,
    _normalize_vector_to_sum_one,
    solve_simplex_constrained_ls,
)


REQUIRED_INPUTS = ("signature.txt", "bulk.txt")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run simplex-constrained least-squares deconvolution."
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        required=True,
        help="A tissue directory or a parent directory containing tissue inputs.",
    )
    parser.add_argument(
        "--mode",
        choices=("raw", "zscore"),
        default="raw",
        help="Use raw shared-gene matrices or MetaSort's normalized joint z-score preprocessing.",
    )
    parser.add_argument(
        "--output-name",
        default=None,
        help="Prediction filename under each tissue's result directory.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing prediction file instead of skipping it.",
    )
    return parser.parse_args()


def discover_cases(data_root: Path) -> list[Path]:
    data_root = data_root.resolve()
    if all((data_root / name).is_file() for name in REQUIRED_INPUTS):
        return [data_root]
    cases = {
        signature.parent
        for signature in data_root.rglob("signature.txt")
        if all((signature.parent / name).is_file() for name in REQUIRED_INPUTS)
    }
    return sorted(cases)


def read_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.apply(pd.to_numeric, errors="raise")


def prepare_inputs(
    signature: pd.DataFrame,
    bulk: pd.DataFrame,
    mixture: str,
    mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    genes = signature.index.intersection(bulk.index)
    if len(genes) == 0:
        raise ValueError("No genes are shared by signature.txt and bulk.txt.")
    signature_array = signature.loc[genes].to_numpy(dtype=float)
    bulk_array = bulk.loc[genes, mixture].to_numpy(dtype=float)
    if mode == "zscore":
        signature_array = _normalize_columns_to_sum_one(signature_array, name="signature")
        bulk_array = _normalize_vector_to_sum_one(bulk_array, name="bulk")
        signature_array, bulk_array = _joint_gene_zscore(signature_array, bulk_array)
    return signature_array, bulk_array


def run_case(case_dir: Path, mode: str, output_name: str, overwrite: bool) -> str:
    output_path = case_dir / "result" / output_name
    if output_path.exists() and not overwrite:
        return f"SKIP {case_dir}: {output_path.name} exists"

    signature = read_matrix(case_dir / "signature.txt")
    bulk = read_matrix(case_dir / "bulk.txt")
    prediction = pd.DataFrame(
        index=bulk.columns,
        columns=signature.columns,
        dtype=float,
    )
    for mixture in bulk.columns:
        signature_array, bulk_array = prepare_inputs(signature, bulk, mixture, mode)
        prediction.loc[mixture] = solve_simplex_constrained_ls(signature_array, bulk_array)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    prediction.to_csv(output_path, sep="\t")
    return f"OK {case_dir}: {prediction.shape[0]} mixtures x {prediction.shape[1]} cell types"


def main() -> int:
    args = parse_args()
    output_name = args.output_name or (
        "SimplexLS_zscore.txt" if args.mode == "zscore" else "SimplexLS_raw.txt"
    )
    cases = discover_cases(args.data_root)
    if not cases:
        raise FileNotFoundError(
            f"No directories containing {', '.join(REQUIRED_INPUTS)} under {args.data_root}"
        )
    failures = 0
    for case_dir in cases:
        try:
            print(run_case(case_dir, args.mode, output_name, args.overwrite), flush=True)
        except Exception as exc:
            failures += 1
            print(f"ERROR {case_dir}: {exc}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
