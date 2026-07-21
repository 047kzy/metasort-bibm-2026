from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metasort.algorithm import (  # noqa: E402
    _joint_gene_zscore,
    _normalize_columns_to_sum_one,
    _normalize_vector_to_sum_one,
    solve_simplex_constrained_ls,
)
from metasort import MetaSortConfig, MetaSortSolver  # noqa: E402


ROOT_DIR = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT_DIR / "data" / "selected_case"
TRUE_INIT_DIR = ROOT_DIR / "results" / "true_init_run"
FIGURE_SOURCE_DIR = ROOT_DIR / "data" / "figure_source"
OUT_DIR = ROOT_DIR / "results" / "true_init_diagnostics"
PREFIX = "bibm_true_init_metasort_vs_simplex_zscore"


def normalize_nonnegative(values: np.ndarray) -> np.ndarray:
    values = np.clip(np.asarray(values, dtype=float), 0.0, None)
    total = float(values.sum())
    if total <= 0.0 or not np.isfinite(total):
        return np.full(values.shape, 1.0 / len(values), dtype=float)
    return values / total


def score_l1(prediction: np.ndarray, truth: np.ndarray) -> float:
    return float(np.abs(normalize_nonnegative(prediction) - normalize_nonnegative(truth)).sum())


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def load_case_arrays() -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    reference = read_matrix(CASE_DIR / "selected_cellTypeExpr.txt")
    bulk_df = read_matrix(CASE_DIR / "sampled_bulk.txt")
    truth_df = read_matrix(CASE_DIR / "sampled_bulkRatio.txt")

    reference.index = reference.index.astype(str)
    reference.columns = reference.columns.astype(str)
    bulk_df.index = bulk_df.index.astype(str)
    truth_df.index = truth_df.index.astype(str)

    mixture = "Mixture10"
    cell_types = [cell_type for cell_type in reference.columns if cell_type in truth_df.index]
    genes = reference.index.intersection(bulk_df.index).astype(str).tolist()
    signature = _normalize_columns_to_sum_one(
        reference.loc[genes, cell_types].to_numpy(dtype=float),
        name="signature",
    )
    bulk = _normalize_vector_to_sum_one(
        bulk_df.loc[genes, mixture].to_numpy(dtype=float),
        name="bulk",
    )
    z_signature, z_bulk = _joint_gene_zscore(signature, bulk)
    prediction = solve_simplex_constrained_ls(z_signature, z_bulk)
    truth = normalize_nonnegative(truth_df.loc[cell_types, mixture].to_numpy(dtype=float))
    return cell_types, truth, normalize_nonnegative(prediction), z_signature, z_bulk


def load_true_init_metasort_prediction(cell_types: list[str]) -> np.ndarray:
    proportions = pd.read_csv(TRUE_INIT_DIR / "tongue_mixture_true_init_metasort_iteration_proportions.csv")
    one_step = proportions[proportions["iteration"].eq(1)].set_index("cell_type")
    return normalize_nonnegative(one_step.loc[cell_types, "prediction"].to_numpy(dtype=float))


def clean_cell_type(label: str) -> str:
    replacements = {
        "leukocyte": "Leukocyte",
        "fibroblast": "Fibroblast",
        "vein_endothelial_cell": "Vein\nendothelial",
        "pericyte_cell": "Pericyte",
        "capillary_endothelial_cell": "Capillary\nendothelial",
    }
    return replacements.get(label, label.replace("_", " ").title())


def build_long_table(cell_types: list[str], truth: np.ndarray, metasort: np.ndarray, simplex: np.ndarray) -> pd.DataFrame:
    rows = []
    for cell_type, true_value, metasort_value, simplex_value in zip(cell_types, truth, metasort, simplex):
        rows.extend(
            [
                {
                    "cell_type": cell_type,
                    "method": "True ratio",
                    "proportion": true_value,
                    "abs_error": 0.0,
                },
                {
                    "cell_type": cell_type,
                    "method": "Curvature step from true",
                    "proportion": metasort_value,
                    "abs_error": abs(metasort_value - true_value),
                },
                {
                    "cell_type": cell_type,
                    "method": "Simplex LS (z-score)",
                    "proportion": simplex_value,
                    "abs_error": abs(simplex_value - true_value),
                },
            ]
        )
    return pd.DataFrame(rows)


def sphere_tangent_basis(point: np.ndarray) -> np.ndarray:
    point = np.asarray(point, dtype=float)
    point = point / max(float(np.linalg.norm(point)), 1e-12)
    n_cell_types = point.shape[0]
    anchor = int(np.argmax(np.abs(point)))
    tangent_indices = [idx for idx in range(n_cell_types) if idx != anchor]
    seed = np.eye(n_cell_types, dtype=float)[:, tangent_indices]
    seed = seed - point[:, None] * (point @ seed)[None, :]
    basis, _ = np.linalg.qr(seed, mode="reduced")
    return basis


def curvature_eigenvalues(signature: np.ndarray, proportions: np.ndarray, weights: np.ndarray, epsilon: float) -> np.ndarray:
    proportions = normalize_nonnegative(proportions)
    sqrt_proportions = np.sqrt(np.clip(proportions, epsilon, None))
    sqrt_proportions = sqrt_proportions / max(float(np.linalg.norm(sqrt_proportions)), epsilon)
    basis = sphere_tangent_basis(sqrt_proportions)
    h_p = signature.T @ (np.asarray(weights, dtype=float)[:, None] * signature)
    h_p = 0.5 * (h_p + h_p.T)
    h_sphere = 4.0 * (sqrt_proportions[:, None] * h_p * sqrt_proportions[None, :])
    projected = basis.T @ h_sphere @ basis
    projected = 0.5 * (projected + projected.T)
    return np.clip(np.linalg.eigvalsh(projected), epsilon, None)


def load_curvature_tables(
    signature: np.ndarray,
    bulk: np.ndarray,
    truth: np.ndarray,
    simplex: np.ndarray,
    metasort: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = MetaSortConfig()
    solver = MetaSortSolver(cfg)
    baseline_weights = np.full(signature.shape[0], cfg.meta_weight_baseline, dtype=float)
    optimized_weights, _ = solver.optimize_meta_weights(
        signature=signature,
        bulk=bulk,
        proportions=truth,
        confidence_weights=np.ones(signature.shape[0], dtype=float),
        prev_meta_weights=baseline_weights,
    )
    summary_rows = []
    eig_rows = []
    for stage, proportions, weights in [
        ("Before optimization", simplex, baseline_weights),
        ("After optimization", metasort, optimized_weights),
    ]:
        eigvals = np.sort(curvature_eigenvalues(signature, proportions, weights, cfg.hessian_epsilon))
        summary_rows.append(
            {
                "stage": stage,
                "min_eig": float(np.min(eigvals)),
                "mean_eig": float(np.mean(eigvals)),
                "max_eig": float(np.max(eigvals)),
                "log10_min_eig": float(np.log10(np.min(eigvals))),
                "log10_mean_eig": float(np.log10(np.mean(eigvals))),
                "log10_max_eig": float(np.log10(np.max(eigvals))),
            }
        )
        for idx, eigval in enumerate(eigvals, start=1):
            eig_rows.append(
                {
                    "stage": stage,
                    "eigen_index": idx,
                    "eigenvalue": float(eigval),
                    "log10_eigenvalue": float(np.log10(eigval)),
                }
            )
    return pd.DataFrame(summary_rows), pd.DataFrame(eig_rows)


def plot(truth: np.ndarray, metasort: np.ndarray, simplex: np.ndarray, curvature_eigs: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    colors = {
        "Curvature step from true": "#0072B2",
        "Simplex LS (z-score)": "#D55E00",
        "Before optimization": "#8A8F98",
        "After optimization": "#0072B2",
    }
    l1_values = {
        "Curvature step from true": score_l1(metasort, truth),
        "Simplex LS (z-score)": score_l1(simplex, truth),
    }

    fig = plt.figure(figsize=(7.1, 3.3), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25])
    ax_l1 = fig.add_subplot(gs[0, 0])
    ax_radar = fig.add_subplot(gs[0, 1], projection="polar")

    l1_methods = ["Curvature step from true", "Simplex LS (z-score)"]
    l1_x = np.arange(len(l1_methods))
    bars = ax_l1.bar(
        l1_x,
        [l1_values[method] for method in l1_methods],
        width=0.58,
        color=[colors[method] for method in l1_methods],
        edgecolor="white",
        linewidth=0.6,
    )
    ax_l1.set_ylabel("L1 error")
    ax_l1.set_xticks(l1_x)
    ax_l1.set_xticklabels(["Curvature\nstep", "Simplex LS\n(z-score)"])
    ax_l1.set_ylim(0, max(l1_values.values()) * 1.22)
    ax_l1.grid(axis="y", color="#D8DEE9", linewidth=0.6, alpha=0.8)
    ax_l1.set_axisbelow(True)
    for bar, method in zip(bars, l1_methods):
        height = bar.get_height()
        ax_l1.text(
            bar.get_x() + bar.get_width() / 2,
            height + max(l1_values.values()) * 0.035,
            f"{height:.3f}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    ax_l1.text(
        0.0,
        1.04,
        "A",
        transform=ax_l1.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
    )

    axes = [f"λ{i}" for i in sorted(curvature_eigs["eigen_index"].unique())]
    angles = np.linspace(0, 2 * np.pi, len(axes), endpoint=False).tolist()
    closed_angles = angles + angles[:1]
    all_log_values = curvature_eigs["log10_eigenvalue"].to_numpy(dtype=float)
    radial_min = float(np.floor(np.min(all_log_values) * 2.0) / 2.0)
    radial_max = float(np.ceil(np.max(all_log_values) * 2.0) / 2.0)
    for stage, stage_df in curvature_eigs.groupby("stage", sort=False):
        values = stage_df.sort_values("eigen_index")["log10_eigenvalue"].to_list()
        closed_values = values + values[:1]
        ax_radar.plot(closed_angles, closed_values, color=colors[stage], linewidth=1.8, label=stage)
        ax_radar.fill(closed_angles, closed_values, color=colors[stage], alpha=0.13)
    ax_radar.set_xticks(angles)
    ax_radar.set_xticklabels(axes)
    ax_radar.set_ylim(radial_min, radial_max)
    ticks = np.linspace(radial_min, radial_max, 5)
    ax_radar.set_yticks(ticks)
    ax_radar.set_yticklabels([f"{tick:.1f}" for tick in ticks], fontsize=7)
    ax_radar.set_rlabel_position(25)
    ax_radar.set_title("log10 projected Hessian eigenvalues", pad=18, fontsize=8)
    ax_radar.grid(color="#D8DEE9", linewidth=0.7)
    ax_radar.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.35, 1.12), fontsize=7)
    ax_radar.text(
        -0.12,
        1.08,
        "B",
        transform=ax_radar.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
    )

    fig.suptitle(
        "Overall error and local curvature spectrum",
        fontsize=10,
        y=1.05,
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png", "svg", "tiff"):
        fig.savefig(OUT_DIR / f"{PREFIX}.{suffix}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    cell_types, truth, simplex, signature, bulk = load_case_arrays()
    metasort = load_true_init_metasort_prediction(cell_types)
    curvature_summary, curvature_eigs = load_curvature_tables(signature, bulk, truth, simplex, metasort)
    data = build_long_table(cell_types, truth, metasort, simplex)
    summary = pd.DataFrame(
        [
            {
                "method": "Curvature step from true",
                "l1": score_l1(metasort, truth),
                "mae": float(np.abs(metasort - truth).mean()),
            },
            {
                "method": "Simplex LS (z-score)",
                "l1": score_l1(simplex, truth),
                "mae": float(np.abs(simplex - truth).mean()),
            },
        ]
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUT_DIR / f"{PREFIX}_proportions.csv", index=False)
    summary.to_csv(OUT_DIR / f"{PREFIX}_l1_summary.csv", index=False)
    summary.to_csv(FIGURE_SOURCE_DIR / f"{PREFIX}_l1_summary.csv", index=False)
    curvature_summary.to_csv(OUT_DIR / f"{PREFIX}_curvature_log_eig_summary.csv", index=False)
    curvature_eigs.to_csv(OUT_DIR / f"{PREFIX}_curvature_log_eig_values.csv", index=False)
    plot(truth, metasort, simplex, curvature_eigs)
    print(summary.to_string(index=False), flush=True)
    print(curvature_summary.to_string(index=False), flush=True)
    print(curvature_eigs.to_string(index=False), flush=True)
    print(f"output_prefix: {OUT_DIR / PREFIX}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
