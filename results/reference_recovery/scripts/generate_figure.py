"""Regenerate the BIBM reference-recovery figure and its component panels."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
PAIR_NAME = "bulk_TSP14__reference_TSP10"
DATA_DIR = BASE_DIR / "data"
BIBM_DIR = BASE_DIR / "figures"

# Optional: point BIBM_FONT_DIR to a directory containing Times New Roman fonts.
# The script otherwise uses the first installed serif font in FONT_PREFERENCE.
font_dir = os.environ.get("BIBM_FONT_DIR")
if font_dir:
    for path in Path(font_dir).iterdir():
        if path.suffix.lower() in {".otf", ".ttf"}:
            font_manager.fontManager.addfont(str(path))

FONT_PREFERENCE = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]
available_fonts = {entry.name for entry in font_manager.fontManager.ttflist}
figure_font = next((name for name in FONT_PREFERENCE if name in available_fonts), "DejaVu Serif")

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": [figure_font],
        "font.size": 9.0,
        "axes.titlesize": 9.6,
        "axes.labelsize": 9.2,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.4,
        "legend.fontsize": 8.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "legend.frameon": False,
    }
)

COLORS = {
    "Input TSP10": "#F58518",
    "Input reference": "#4C78A8",
    "True TSP14": "#4C78A8",
    "Recovered": "#54A24B",
    "Recovered reference": "#54A24B",
}


def short_cell_label(label: str) -> str:
    mapping = {
        "macrophage": "Macrophage",
        "stromal_cell": "Stromal cell",
        "mast_cell": "Mast cell",
        "muscle_cell": "Muscle cell",
        "CD1c_positive_myeloid_dendritic_cell": "CD1c+ myeloid DC",
        "endothelial_cell": "Endothelial cell",
        "naive_thymus_derived_CD8_positive_alpha_beta_T_cell": "Naive CD8 T cell",
        "regulatory_T_cell": "Regulatory T cell",
    }
    if label in mapping:
        return mapping[label]
    label = label.replace("_", " ")
    label = label.replace(
        "naive thymus derived CD8 positive alpha beta T cell",
        "Naive CD8 T cell",
    )
    label = label.replace("CD1c positive myeloid dendritic cell", "CD1c+ myeloid DC")
    return label


def compact_cell_label(label: str) -> str:
    mapping = {
        "macrophage": "Mac",
        "stromal_cell": "Stromal",
        "mast_cell": "Mast",
        "muscle_cell": "Muscle",
        "CD1c_positive_myeloid_dendritic_cell": "CD1c DC",
        "endothelial_cell": "Endo",
        "naive_thymus_derived_CD8_positive_alpha_beta_T_cell": "CD8 T",
        "regulatory_T_cell": "Treg",
    }
    return mapping.get(label, short_cell_label(label))


def read_matrix(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / name, sep="\t", index_col=0)


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    BIBM_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(BIBM_DIR / f"{stem}_bibm.pdf")
    fig.savefig(BIBM_DIR / f"{stem}_bibm.png", dpi=600)
    plt.close(fig)


def draw_celltype_rmse(
    ax: mpl.axes.Axes,
    *,
    show_legend: bool = True,
    compact_labels: bool = False,
) -> None:
    metrics = pd.read_csv(DATA_DIR / "celltype_rmse.tsv", sep="\t")
    label_fn = compact_cell_label if compact_labels else short_cell_label
    labels = [label_fn(x) for x in metrics["cell_type"]]
    y = np.arange(len(metrics))
    bar_h = 0.34

    ax.barh(
        y - bar_h / 2,
        metrics["input_rmse"],
        height=bar_h,
        color=COLORS["Input reference"],
        label="Input reference",
    )
    ax.barh(
        y + bar_h / 2,
        metrics["recovered_rmse"],
        height=bar_h,
        color=COLORS["Recovered reference"],
        label="Recovered reference",
    )

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("RMSE to true TSP14 reference")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.35)
    ax.set_axisbelow(True)
    if show_legend:
        ax.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.04),
            ncol=2,
            handlelength=1.2,
            columnspacing=0.8,
            borderaxespad=0.0,
        )


def plot_celltype_rmse() -> None:
    fig, ax = plt.subplots(figsize=(3.35, 3.35))
    draw_celltype_rmse(ax, show_legend=False)
    fig.suptitle("TSP14 from TSP10: cell-type RMSE", fontsize=9.0, fontweight="bold", y=0.985)
    fig.legend(
        *ax.get_legend_handles_labels(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=2,
        handlelength=1.2,
        columnspacing=0.9,
        borderaxespad=0.0,
    )
    fig.subplots_adjust(left=0.39, right=0.985, bottom=0.13, top=0.80)
    save_figure(fig, f"{PAIR_NAME}_celltype_rmse")


def load_representative_gene_profile_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = pd.read_csv(DATA_DIR / "selected_representative_genes.tsv", sep="\t").head(4)
    true_ref = read_matrix("true_bulk_donor_reference.tsv")
    input_ref = read_matrix("input_reference_donor_reference.tsv")
    recovered_ref = read_matrix("recovered_bulk_domain_reference.tsv")
    return selected, true_ref, input_ref, recovered_ref


def draw_representative_gene_profiles(
    axes: np.ndarray,
    *,
    fig: mpl.figure.Figure,
    legend_anchor: tuple[float, float],
) -> None:
    selected, true_ref, input_ref, recovered_ref = load_representative_gene_profile_data()

    genes = list(selected["gene"])
    x = np.arange(true_ref.shape[1])
    cell_labels = [compact_cell_label(c) for c in true_ref.columns]
    series = [
        ("Input TSP10", input_ref, "s"),
        ("True TSP14", true_ref, "o"),
        ("Recovered", recovered_ref, "^"),
    ]

    axes_flat = axes.ravel()

    for idx, (ax, gene) in enumerate(zip(axes_flat, genes)):
        row = selected.loc[selected["gene"] == gene].iloc[0]
        for label, matrix, marker in series:
            ax.plot(
                x,
                matrix.loc[gene, true_ref.columns].to_numpy(dtype=float),
                marker=marker,
                markersize=2.15,
                markeredgewidth=0.25,
                linewidth=0.72,
                color=COLORS[label],
                label=label,
            )

        gene_name = row["feature_name"] if isinstance(row["feature_name"], str) and row["feature_name"] else gene
        ax.set_title(
            f"{gene_name}\nRMSE {row['reference_gene_rmse']:.2f} -> {row['recovered_gene_rmse']:.2f}",
            fontsize=8.4,
            fontweight="bold",
            pad=3,
        )
        ax.grid(axis="y", color="#E1E1E1", linewidth=0.35)
        ax.set_axisbelow(True)
        ax.margins(x=0.02)
        if idx % 2 == 0:
            ax.set_ylabel("log1p(CP10000)", labelpad=1.5)

    for ax in axes_flat:
        ax.set_xticks(x)
        ax.set_xticklabels(cell_labels, rotation=42, ha="right", rotation_mode="anchor")
    for ax in axes_flat[:2]:
        ax.tick_params(axis="x", labelbottom=False)

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=legend_anchor,
        ncol=3,
        handlelength=1.25,
        columnspacing=0.55,
        borderaxespad=0.0,
    )


def plot_representative_gene_profiles() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(3.35, 3.35), sharex=True)
    draw_representative_gene_profiles(axes, fig=fig, legend_anchor=(0.5, 0.895))
    fig.suptitle(
        "Representative gene profiles",
        fontsize=10.0,
        fontweight="bold",
        y=0.968,
    )
    fig.subplots_adjust(left=0.16, right=0.985, bottom=0.18, top=0.70, hspace=0.96, wspace=0.38)
    save_figure(fig, f"{PAIR_NAME}_representative_gene_profiles")


def plot_combined_celltype_rmse_and_gene_profiles() -> None:
    fig = plt.figure(figsize=(6.70, 3.35))
    outer = fig.add_gridspec(
        1,
        2,
        left=0.092,
        right=0.992,
        bottom=0.17,
        top=0.74,
        width_ratios=[1, 1],
        wspace=0.30,
    )

    ax_left = fig.add_subplot(outer[0, 0])
    draw_celltype_rmse(ax_left, show_legend=False, compact_labels=True)
    fig.text(
        0.264,
        0.962,
        "Cell-type RMSE",
        ha="center",
        va="center",
        fontsize=10.0,
        fontweight="bold",
    )

    right = outer[0, 1].subgridspec(2, 2, hspace=0.82, wspace=0.34)
    right_axes = np.array(
        [
            [fig.add_subplot(right[0, 0]), fig.add_subplot(right[0, 1])],
            [fig.add_subplot(right[1, 0]), fig.add_subplot(right[1, 1])],
        ]
    )
    draw_representative_gene_profiles(right_axes, fig=fig, legend_anchor=(0.735, 0.885))
    fig.text(
        0.735,
        0.962,
        "Representative gene profiles",
        ha="center",
        va="center",
        fontsize=10.0,
        fontweight="bold",
    )

    handles, labels = ax_left.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.264, 0.885),
        ncol=2,
        handlelength=1.2,
        columnspacing=0.8,
        borderaxespad=0.0,
    )

    fig.text(0.012, 0.965, "A", ha="left", va="center", fontsize=16, fontweight="bold")
    fig.text(0.505, 0.965, "B", ha="left", va="center", fontsize=16, fontweight="bold")
    save_figure(fig, f"{PAIR_NAME}_celltype_rmse_representative_gene_profiles")


def main() -> None:
    plot_celltype_rmse()
    plot_representative_gene_profiles()
    plot_combined_celltype_rmse_and_gene_profiles()


if __name__ == "__main__":
    main()
