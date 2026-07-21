from __future__ import annotations

import csv
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT_DIR / "data" / "selected_case"
CASE_METRICS_DIR = ROOT_DIR / "results" / "case_metrics"
OUT_DIR = ROOT_DIR / "results" / "method_predictions"
RAW_R_DIR = OUT_DIR / "raw_r_exports"
DATASET = "Tabula-Sapiens"
TISSUE = "Tongue"
MIXTURE = "Mixture10"

FIELDS = [
    "dataset",
    "tissue",
    "mixture",
    "method",
    "cell_type",
    "truth",
    "prediction",
    "abs_error",
]


def read_truth() -> dict[str, float]:
    path = CASE_DIR / "sampled_bulkRatio.txt"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return {row[0]: float(row[1]) for row in rows[1:]}


def read_long_predictions() -> dict[str, list[dict[str, str]]]:
    path = CASE_METRICS_DIR / "predictions_long.csv"
    grouped: dict[str, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            grouped.setdefault(row["method"], []).append(
                {
                    "dataset": row["dataset"],
                    "tissue": row["tissue"],
                    "mixture": row["mixture"],
                    "method": row["method"],
                    "cell_type": row["cell_type"],
                    "truth": row["truth"],
                    "prediction": row["prediction"],
                    "abs_error": row["abs_error"],
                }
            )
    return grouped


def read_r_prediction(method: str, truth: dict[str, float]) -> list[dict[str, str]]:
    path = RAW_R_DIR / f"{method}.txt"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    cell_types = rows[0][1:]
    predictions = rows[1][1:]
    return [
        {
            "dataset": DATASET,
            "tissue": TISSUE,
            "mixture": MIXTURE,
            "method": method,
            "cell_type": cell_type,
            "truth": str(truth[cell_type]),
            "prediction": prediction,
            "abs_error": str(abs(float(prediction) - truth[cell_type])),
        }
        for cell_type, prediction in zip(cell_types, predictions)
    ]


def write_method(method: str, rows: list[dict[str, str]]) -> None:
    path = OUT_DIR / f"{method}.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    truth = read_truth()
    grouped = read_long_predictions()
    grouped["MuSiC"] = read_r_prediction("MuSiC", truth)
    grouped["LinDeconSeq"] = read_r_prediction("LinDeconSeq", truth)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for method in sorted(grouped):
        write_method(method, grouped[method])
    print(f"exported {len(grouped)} method prediction files to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
