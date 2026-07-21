from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT_DIR / "data" / "figure_source"
CASE_DIR = ROOT_DIR / "data" / "selected_case"
CASE_METRICS_DIR = ROOT_DIR / "results" / "case_metrics"
METASORT_DIR = ROOT_DIR / "results" / "metasort_run"
TRUE_INIT_DIR = ROOT_DIR / "results" / "true_init_run"
METHOD_PREDICTIONS_DIR = ROOT_DIR / "results" / "method_predictions"
EXPECTED_METHODS = {
    "MetaSort",
    "MuSiC",
    "LinDeconSeq",
    "DWLS_style_raw",
    "NNLS_raw",
    "NNLS_zscore",
    "SimplexLS_raw",
}
EXPECTED_CASE_INPUTS = {
    "sampled_bulk.txt",
    "sampled_bulkRatio.txt",
    "selected_cellTypeExpr.txt",
    "selected_genes.txt",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def assert_close(actual: str, expected: str, label: str) -> None:
    if not math.isclose(float(actual), float(expected), rel_tol=1e-12, abs_tol=1e-12):
        raise AssertionError(f"{label}: {actual} != {expected}")


def validate_layout() -> None:
    case_inputs = {path.name for path in CASE_DIR.iterdir() if path.is_file()}
    case_subdirectories = [path.name for path in CASE_DIR.iterdir() if path.is_dir()]
    if case_inputs != EXPECTED_CASE_INPUTS or case_subdirectories:
        raise AssertionError(
            f"selected_case layout: files={sorted(case_inputs)}, dirs={sorted(case_subdirectories)}"
        )


def validate_panel_a() -> None:
    panel_rows = read_rows(SOURCE_DIR / "bibm_true_init_metasort_vs_simplex_zscore_l1_summary.csv")
    panel = {row["method"]: row["l1"] for row in panel_rows}

    true_init_rows = read_rows(TRUE_INIT_DIR / "tongue_mixture_true_init_metasort_iteration_l1.csv")
    one_step = next(row for row in true_init_rows if row["iteration"] == "1")
    assert_close(panel["Curvature step from true"], one_step["l1"], "Panel A true-init L1")

    metasort_rows = read_rows(METASORT_DIR / "tongue_mixture_metasort_iteration_l1.csv")
    initialization = next(row for row in metasort_rows if row["iteration"] == "0")
    assert_close(panel["Simplex LS (z-score)"], initialization["l1"], "Panel A simplex L1")


def validate_panel_b() -> None:
    panel_rows = read_rows(SOURCE_DIR / "bibm_metasort_tongue_case_combined_revised_values.csv")
    upstream_rows = read_rows(CASE_METRICS_DIR / "case_method_summary_with_music_lindeconseq.csv")
    upstream = {row["method"]: row for row in upstream_rows}
    for row in panel_rows:
        method = row["method"]
        assert_close(row["mean_l1"], upstream[method]["mean_l1"], f"Panel B {method} L1")
        assert_close(row["mean_mae"], upstream[method]["mean_mae"], f"Panel B {method} MAE")


def validate_panel_c() -> None:
    panel_rows = read_rows(SOURCE_DIR / "bibm_metasort_tongue_case_iteration_0_13.csv")
    run_rows = read_rows(METASORT_DIR / "tongue_mixture_metasort_iteration_l1.csv")
    run = {row["iteration"]: row for row in run_rows}
    expected_iterations = [str(index) for index in range(14)]
    observed_iterations = [row["iteration"] for row in panel_rows]
    if observed_iterations != expected_iterations:
        raise AssertionError(f"Panel C iterations: {observed_iterations} != {expected_iterations}")
    for row in panel_rows:
        assert_close(row["l1"], run[row["iteration"]]["l1"], f"Panel C iteration {row['iteration']}")


def validate_method_predictions() -> None:
    summary_rows = read_rows(CASE_METRICS_DIR / "case_method_summary_with_music_lindeconseq.csv")
    summary = {row["method"]: row for row in summary_rows}
    if set(summary) != EXPECTED_METHODS:
        raise AssertionError(f"Methods in summary: {sorted(summary)}")

    prediction_files = {path.stem: path for path in METHOD_PREDICTIONS_DIR.glob("*.csv")}
    if set(prediction_files) != EXPECTED_METHODS:
        raise AssertionError(f"Prediction files: {sorted(prediction_files)}")

    for method, path in prediction_files.items():
        rows = read_rows(path)
        if len(rows) != 5:
            raise AssertionError(f"{method}: expected 5 cell types, found {len(rows)}")
        l1 = sum(float(row["abs_error"]) for row in rows)
        assert_close(str(l1), summary[method]["mean_l1"], f"{method} prediction L1")


def main() -> int:
    validate_layout()
    validate_panel_a()
    validate_panel_b()
    validate_panel_c()
    validate_method_predictions()
    print("retained data validation: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
