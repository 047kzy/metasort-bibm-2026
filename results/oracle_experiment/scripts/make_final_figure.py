from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "data" / "figure_source"
OUT_DIR = ROOT_DIR / "results" / "final_figure"
OUT_PREFIX = OUT_DIR / "bibm_metasort_tongue_case_combined_revised_with_true_init_l1"
PANEL_FONT = ["Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif"]
TITLE_FONT_SIZE = 14
LABEL_FONT_SIZE = 12
TICK_FONT_SIZE = 11
VALUE_FONT_SIZE = 11
ANNOTATION_FONT_SIZE = 11
PANEL_LABEL_FONT_SIZE = 16


def load_method_values() -> pd.DataFrame:
    values = pd.read_csv(DATA_DIR / "bibm_metasort_tongue_case_combined_revised_values.csv")
    order = ["MetaSort", "MuSiC", "DWLS", "LinDeconSeq", "NNLS"]
    values["plot_method"] = pd.Categorical(values["plot_method"], categories=order, ordered=True)
    return values.sort_values("plot_method")


def load_iteration_trace() -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / "bibm_metasort_tongue_case_iteration_0_13.csv")


def load_true_init_values() -> pd.DataFrame:
    summary = pd.read_csv(DATA_DIR / "bibm_true_init_metasort_vs_simplex_zscore_l1_summary.csv")
    label_map = {
        "Curvature step from true": "MetaSort\ntrue init\n1 iter",
        "Simplex LS (z-score)": "Simplex LS\nz-score",
    }
    summary["label"] = summary["method"].map(label_map)
    order = ["MetaSort\ntrue init\n1 iter", "Simplex LS\nz-score"]
    summary["label"] = pd.Categorical(summary["label"], categories=order, ordered=True)
    return summary.sort_values("label")


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.22,
        1.25,
        label,
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_FONT_SIZE,
        fontweight="bold",
        fontfamily=PANEL_FONT,
        va="bottom",
        ha="left",
    )


def add_inset_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.18,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=7,
        fontweight="bold",
        fontfamily=PANEL_FONT,
        va="bottom",
        ha="left",
    )


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.9)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.7)
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE, width=1.0, length=4)
    ax.set_axisbelow(True)


def plot_methods(ax: plt.Axes, methods: pd.DataFrame, panel_label: str = "B") -> None:
    colors = ["#D62728" if method == "MetaSort" else "#7EA9C2" for method in methods["plot_method"]]
    bars = ax.bar(
        np.arange(len(methods)),
        methods["mean_l1"],
        color=colors,
        edgecolor="#333333",
        linewidth=0.7,
        width=0.72,
    )
    ax.set_title("Final deconvolution error", fontsize=TITLE_FONT_SIZE, pad=10)
    ax.set_ylabel("L1 error", fontsize=LABEL_FONT_SIZE)
    ax.set_xticks(np.arange(len(methods)))
    ax.set_xticklabels(methods["plot_method"], rotation=35, ha="right", fontsize=10.5)
    ax.set_ylim(0, 0.56)
    style_axis(ax)
    for bar, value in zip(bars, methods["mean_l1"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.012,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=VALUE_FONT_SIZE,
        )

    add_panel_label(ax, panel_label)


def plot_iteration(ax: plt.Axes, trace: pd.DataFrame, panel_label: str = "C") -> None:
    ax.plot(
        trace["iteration"],
        trace["l1"],
        color="#1F5EA8",
        linewidth=2.8,
        marker="o",
        markersize=6.0,
        markerfacecolor="white",
        markeredgewidth=1.8,
    )
    final = trace.iloc[-1]
    ax.scatter([final["iteration"]], [final["l1"]], s=86, color="#D62728", zorder=5)
    ax.annotate(
        f"iter. {int(final['iteration'])}\nL1 = {final['l1']:.3f}",
        xy=(final["iteration"], final["l1"]),
        xytext=(final["iteration"] - 2.8, final["l1"] + 0.065),
        arrowprops={"arrowstyle": "-", "color": "#6B6B6B", "lw": 1.2},
        fontsize=ANNOTATION_FONT_SIZE,
        ha="left",
        va="center",
        color="#2E2E2E",
    )
    ax.set_title("MetaSort convergence", fontsize=TITLE_FONT_SIZE, pad=10)
    ax.set_xlabel("MetaSort iteration", fontsize=LABEL_FONT_SIZE)
    ax.set_ylabel("L1 error", fontsize=LABEL_FONT_SIZE)
    ax.set_ylim(0, 0.35)
    ax.set_xlim(-0.35, float(trace["iteration"].max()) + 0.45)
    style_axis(ax)
    add_panel_label(ax, panel_label)


def plot_true_init_panel(ax: plt.Axes, values: pd.DataFrame, panel_label: str = "A") -> None:
    colors = ["#D62728", "#7EA9C2"]
    x = np.array([0.0, 0.60])
    bars = ax.bar(
        x,
        values["l1"],
        color=colors,
        edgecolor="#333333",
        linewidth=0.7,
        width=0.14,
    )
    ax.set_title("True-init L1", fontsize=TITLE_FONT_SIZE, pad=10)
    ax.set_ylabel("L1 error", fontsize=LABEL_FONT_SIZE)
    ax.set_xticks(x)
    ax.set_xticklabels(["MetaSort\ntrue init", "Simplex\nz-score"], fontsize=10)
    ax.set_xlim(-0.30, 0.90)
    ax.set_ylim(0, 0.36)
    style_axis(ax)
    for bar, value in zip(bars, values["l1"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.012,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=VALUE_FONT_SIZE,
        )
    add_panel_label(ax, panel_label)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Regenerate the retained Part 4 figure.")
    parser.add_argument(
        "--with-tiff",
        action="store_true",
        help="also export a 600-dpi TIFF (large; omitted by default)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": PANEL_FONT,
            "font.size": TICK_FONT_SIZE,
            "axes.linewidth": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    methods = load_method_values()
    trace = load_iteration_trace()
    true_init = load_true_init_values()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(7.25, 2.85))
    left = 0.070
    bottom = 0.34
    height = 0.51
    gap_ab = 0.078
    gap_bc = 0.095
    width_a = 0.170
    width_b = 0.235
    width_c = 0.220
    ax_a = fig.add_axes([left, bottom, width_a, height])
    ax_b = fig.add_axes([left + width_a + gap_ab, bottom, width_b, height])
    ax_c = fig.add_axes([left + width_a + gap_ab + width_b + gap_bc, bottom, width_c, height])

    plot_true_init_panel(ax_a, true_init, panel_label="A")
    plot_methods(ax_b, methods, panel_label="B")
    plot_iteration(ax_c, trace, panel_label="C")

    suffixes = ["pdf", "png", "svg"]
    if args.with_tiff:
        suffixes.append("tiff")
    for suffix in suffixes:
        fig.savefig(f"{OUT_PREFIX}.{suffix}", dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"output_prefix: {OUT_PREFIX}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
