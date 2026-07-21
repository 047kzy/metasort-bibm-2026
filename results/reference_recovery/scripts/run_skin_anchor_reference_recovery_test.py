from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from metasort import AnchorRecoveryConfig, AnchorRecoverySolver, MetaSortConfig, MetaSortSolver


DEFAULT_H5AD_PATH = Path("results/reference_recovery/external_inputs/Skin.h5ad")
DEFAULT_SIGNATURE_PATH = Path(
    "results/rare_celltype/data/inputs/Tabula-Sapiens/Skin/signature.txt"
)
DEFAULT_OUTPUT_DIR = Path("results/reference_recovery/rerun_outputs")

CELL_TYPE_MAP = {
    "macrophage": ["macrophage"],
    "stromal_cell": ["fibroblast", "stromal cell"],
    "mast_cell": ["mast cell"],
    "muscle_cell": ["smooth muscle cell", "muscle cell"],
    "CD1c_positive_myeloid_dendritic_cell": ["myeloid dendritic cell"],
    "endothelial_cell": ["endothelial cell"],
    "naive_thymus_derived_CD8_positive_alpha_beta_T_cell": [
        "CD8-positive, alpha-beta T cell"
    ],
    "regulatory_T_cell": ["regulatory T cell"],
}

TARGET_SUM = 1e4
MIN_DONOR_CELLTYPE_CELLS = 10
N_PSEUDOBULK = 16
DIRICHLET_ALPHA = 1.5
RANDOM_SEED = 7


def _log_normalize_counts(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.tocsr().astype(float)
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        factors = np.divide(
            TARGET_SUM,
            totals,
            out=np.zeros_like(totals, dtype=float),
            where=totals > 0.0,
        )
        return np.log1p(matrix.multiply(factors[:, None]).toarray())

    dense = np.asarray(matrix, dtype=float)
    totals = dense.sum(axis=1)
    factors = np.divide(
        TARGET_SUM,
        totals,
        out=np.zeros_like(totals, dtype=float),
        where=totals > 0.0,
    )
    return np.log1p(dense * factors[:, None])


def _map_cell_types(obs_cell_types: pd.Series, signature_cell_types: list[str]) -> pd.Series:
    mapped = pd.Series(index=obs_cell_types.index, dtype=object)
    for signature_cell_type in signature_cell_types:
        source_labels = CELL_TYPE_MAP[signature_cell_type]
        mapped.loc[obs_cell_types.isin(source_labels)] = signature_cell_type
    return mapped


def _donor_reference(
    expr: np.ndarray,
    obs: pd.DataFrame,
    genes: list[str],
    cell_types: list[str],
    donor_id: str,
) -> pd.DataFrame:
    donor_mask = obs["donor_id"].to_numpy() == donor_id
    values = {}
    for cell_type in cell_types:
        mask = donor_mask & (obs["signature_cell_type"].to_numpy() == cell_type)
        values[cell_type] = expr[mask].mean(axis=0)
    return pd.DataFrame(values, index=genes)


def _eligible_donors(counts: pd.DataFrame, cell_types: list[str]) -> list[str]:
    table = counts.pivot_table(
        index="donor_id",
        columns="signature_cell_type",
        values="n_cells",
        fill_value=0,
        aggfunc="sum",
    )
    table = table.reindex(columns=cell_types, fill_value=0)
    eligible = table.index[(table >= MIN_DONOR_CELLTYPE_CELLS).all(axis=1)].astype(str).tolist()
    return sorted(eligible)


def _make_pseudobulk(
    donor_signature: np.ndarray,
    rng: np.random.Generator,
    n_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_cell_types = donor_signature.shape[1]
    proportions = rng.dirichlet(
        np.full(n_cell_types, DIRICHLET_ALPHA, dtype=float),
        size=n_samples,
    ).T
    bulk = donor_signature @ proportions
    return bulk, proportions


def _nnls_proportions(signature: np.ndarray, bulk: np.ndarray) -> np.ndarray:
    n_cell_types = signature.shape[1]
    n_samples = bulk.shape[1]
    proportions = np.zeros((n_cell_types, n_samples), dtype=float)
    for sample_idx in range(n_samples):
        solution, _ = nnls(signature, bulk[:, sample_idx])
        total = float(solution.sum())
        if total <= 0.0:
            solution = np.full(n_cell_types, 1.0 / n_cell_types, dtype=float)
        else:
            solution = solution / total
        proportions[:, sample_idx] = solution
    return proportions


def _proportion_metrics(estimate: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    diff = estimate - truth
    l1_by_sample = np.sum(np.abs(diff), axis=0)
    return {
        "proportion_mae": float(np.mean(np.abs(diff))),
        "proportion_mean_l1": float(np.mean(l1_by_sample)),
        "proportion_accuracy_1_minus_l1_over_2": float(np.mean(1.0 - l1_by_sample / 2.0)),
    }


def _rmse(left: np.ndarray, right: np.ndarray, mask: np.ndarray | None = None) -> float:
    if mask is None:
        diff = left - right
    else:
        diff = left[mask, :] - right[mask, :]
    return float(np.sqrt(np.mean(diff ** 2)))


def _write_matrix(path: Path, matrix: np.ndarray, index: list[str], columns: list[str]) -> None:
    pd.DataFrame(matrix, index=index, columns=columns).to_csv(path, sep="\t")


def _run_pair(
    bulk_donor: str,
    reference_donor: str,
    donor_signatures: dict[str, pd.DataFrame],
    genes: list[str],
    cell_types: list[str],
    rng: np.random.Generator,
    output_dir: Path,
) -> dict[str, float | int | str]:
    pair_dir = output_dir / f"bulk_{bulk_donor}__reference_{reference_donor}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    true_bulk_signature = donor_signatures[bulk_donor].to_numpy(dtype=float)
    reference_signature = donor_signatures[reference_donor].to_numpy(dtype=float)
    bulk_matrix, true_proportions = _make_pseudobulk(
        donor_signature=true_bulk_signature,
        rng=rng,
        n_samples=N_PSEUDOBULK,
    )

    all_gene_nnls = _nnls_proportions(reference_signature, bulk_matrix)
    oracle_nnls = _nnls_proportions(true_bulk_signature, bulk_matrix)

    initial_solver = MetaSortSolver(MetaSortConfig(print_info=False))
    solver = AnchorRecoverySolver(
        AnchorRecoveryConfig(
            lambda_reference=0.1,
            anchor_fit_quantile=0.35,
            anchor_weight_quantile=0.65,
            min_anchor_genes=50,
            max_anchor_fraction=0.5,
            enforce_nonnegative_signature=True,
        ),
        initial_solver=initial_solver,
    )
    result = solver.solve(
        reference_signature,
        bulk_matrix,
        cell_types=cell_types,
        genes=genes,
    )
    initial_proportions = np.asarray(result.initial_proportions, dtype=float)
    anchor_proportions = np.asarray(result.proportions, dtype=float)
    recovered_signature = np.asarray(result.recovered_signature, dtype=float)
    anchor_mask = np.asarray(result.anchor_mask, dtype=bool)
    non_anchor_mask = ~anchor_mask

    _write_matrix(pair_dir / "true_bulk_donor_reference.tsv", true_bulk_signature, genes, cell_types)
    _write_matrix(pair_dir / "input_reference_donor_reference.tsv", reference_signature, genes, cell_types)
    _write_matrix(pair_dir / "recovered_bulk_domain_reference.tsv", recovered_signature, genes, cell_types)
    _write_matrix(pair_dir / "batch_shift_recovered_minus_reference.tsv", recovered_signature - reference_signature, genes, cell_types)
    _write_matrix(
        pair_dir / "pseudobulk.tsv",
        bulk_matrix,
        genes,
        [f"pseudobulk_{idx + 1}" for idx in range(N_PSEUDOBULK)],
    )
    _write_matrix(
        pair_dir / "true_proportions.tsv",
        true_proportions,
        cell_types,
        [f"pseudobulk_{idx + 1}" for idx in range(N_PSEUDOBULK)],
    )
    _write_matrix(
        pair_dir / "initial_metasort_proportions.tsv",
        initial_proportions,
        cell_types,
        [f"pseudobulk_{idx + 1}" for idx in range(N_PSEUDOBULK)],
    )
    _write_matrix(
        pair_dir / "anchor_simplex_proportions.tsv",
        anchor_proportions,
        cell_types,
        [f"pseudobulk_{idx + 1}" for idx in range(N_PSEUDOBULK)],
    )
    _write_matrix(
        pair_dir / "all_gene_nnls_proportions.tsv",
        all_gene_nnls,
        cell_types,
        [f"pseudobulk_{idx + 1}" for idx in range(N_PSEUDOBULK)],
    )
    _write_matrix(
        pair_dir / "oracle_true_reference_nnls_proportions.tsv",
        oracle_nnls,
        cell_types,
        [f"pseudobulk_{idx + 1}" for idx in range(N_PSEUDOBULK)],
    )

    pd.Series(result.anchor_genes, name="gene").to_csv(
        pair_dir / "anchor_genes.txt",
        sep="\t",
        index=False,
    )
    pd.DataFrame(result.gene_metrics, index=genes).assign(anchor=anchor_mask).to_csv(
        pair_dir / "gene_metrics.tsv",
        sep="\t",
    )
    reference_gene_rmse = np.sqrt(np.mean((reference_signature - true_bulk_signature) ** 2, axis=1))
    recovered_gene_rmse = np.sqrt(np.mean((recovered_signature - true_bulk_signature) ** 2, axis=1))
    gene_recovery = pd.DataFrame(
        {
            "gene": genes,
            "anchor": anchor_mask,
            "reference_gene_rmse": reference_gene_rmse,
            "recovered_gene_rmse": recovered_gene_rmse,
            "rmse_improvement": reference_gene_rmse - recovered_gene_rmse,
        }
    )
    gene_recovery["rmse_improvement_pct"] = np.divide(
        100.0 * gene_recovery["rmse_improvement"],
        gene_recovery["reference_gene_rmse"],
        out=np.full(gene_recovery.shape[0], np.nan, dtype=float),
        where=gene_recovery["reference_gene_rmse"].to_numpy() > 0.0,
    )
    gene_recovery.to_csv(pair_dir / "gene_recovery_metrics.tsv", sep="\t", index=False)

    metrics: dict[str, float | int | str] = {
        "bulk_donor": bulk_donor,
        "reference_donor": reference_donor,
        "n_genes": len(genes),
        "n_cell_types": len(cell_types),
        "n_pseudobulk": N_PSEUDOBULK,
        "n_anchor_genes": int(anchor_mask.sum()),
        "n_non_anchor_genes": int(non_anchor_mask.sum()),
        "reference_rmse_all": _rmse(reference_signature, true_bulk_signature),
        "recovered_rmse_all": _rmse(recovered_signature, true_bulk_signature),
        "reference_rmse_non_anchor": _rmse(reference_signature, true_bulk_signature, non_anchor_mask),
        "recovered_rmse_non_anchor": _rmse(recovered_signature, true_bulk_signature, non_anchor_mask),
        "reference_rmse_anchor": _rmse(reference_signature, true_bulk_signature, anchor_mask),
        "recovered_rmse_anchor": _rmse(recovered_signature, true_bulk_signature, anchor_mask),
        "bulk_reconstruction_rmse_input_reference_all_gene_nnls": _rmse(reference_signature @ all_gene_nnls, bulk_matrix),
        "bulk_reconstruction_rmse_recovered_anchor_nnls": _rmse(recovered_signature @ anchor_proportions, bulk_matrix),
        "bulk_reconstruction_rmse_oracle": _rmse(true_bulk_signature @ oracle_nnls, bulk_matrix),
    }
    metrics.update({f"all_gene_nnls_{k}": v for k, v in _proportion_metrics(all_gene_nnls, true_proportions).items()})
    metrics.update({f"initial_metasort_{k}": v for k, v in _proportion_metrics(initial_proportions, true_proportions).items()})
    metrics.update({f"anchor_nnls_{k}": v for k, v in _proportion_metrics(anchor_proportions, true_proportions).items()})
    metrics.update({f"oracle_{k}": v for k, v in _proportion_metrics(oracle_nnls, true_proportions).items()})

    with open(pair_dir / "metrics.json", "w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    return metrics


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the donor-mismatched Skin reference-recovery experiment used for Figure 4. "
            "All default paths are resolved from the repository root."
        )
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=DEFAULT_H5AD_PATH,
        help=f"Skin single-cell h5ad file (default: {DEFAULT_H5AD_PATH})",
    )
    parser.add_argument(
        "--signature",
        type=Path,
        default=DEFAULT_SIGNATURE_PATH,
        help=f"Skin marker signature TSV (default: {DEFAULT_SIGNATURE_PATH})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    return parser.parse_args()


def _resolve_repo_path(path: Path) -> Path:
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (REPO_ROOT / path).resolve()


def main() -> None:
    args = _parse_args()
    h5ad_path = _resolve_repo_path(args.h5ad)
    signature_path = _resolve_repo_path(args.signature)
    output_dir = _resolve_repo_path(args.output_dir)

    if not h5ad_path.is_file():
        raise FileNotFoundError(
            f"Skin h5ad file not found: {h5ad_path}. Pass its location with --h5ad."
        )
    if not signature_path.is_file():
        raise FileNotFoundError(
            f"Skin signature file not found: {signature_path}. Pass its location with --signature."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    signature = pd.read_csv(signature_path, sep="\t", index_col=0)
    marker_genes = signature.index.astype(str).tolist()
    cell_types = signature.columns.astype(str).tolist()

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        common_genes = [gene for gene in marker_genes if gene in set(adata.var_names)]
        mapped = _map_cell_types(adata.obs["cell_type"], cell_types)
        selected = mapped.notna() & adata.obs["donor_id"].notna()
        subset = adata[selected.to_numpy(), common_genes].to_memory()
        obs = subset.obs[["donor_id", "cell_type"]].copy()
        obs["signature_cell_type"] = mapped.loc[obs.index].astype(str)
        matrix = subset.layers["decontXcounts"] if "decontXcounts" in subset.layers else subset.X
        expr = _log_normalize_counts(matrix)
    finally:
        adata.file.close()

    counts = (
        obs.groupby(["donor_id", "signature_cell_type"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    counts.to_csv(output_dir / "donor_celltype_counts.tsv", sep="\t", index=False)
    eligible_donors = _eligible_donors(counts, cell_types)
    if len(eligible_donors) < 2:
        raise RuntimeError("Need at least two donors with all target cell types represented.")

    donor_signatures = {
        donor_id: _donor_reference(expr, obs, common_genes, cell_types, donor_id)
        for donor_id in eligible_donors
    }
    for donor_id, donor_signature in donor_signatures.items():
        donor_signature.to_csv(output_dir / f"donor_{donor_id}_reference_log1p_cp10k.tsv", sep="\t")

    rng = np.random.default_rng(RANDOM_SEED)
    metrics = []
    for bulk_donor in eligible_donors:
        for reference_donor in eligible_donors:
            if bulk_donor == reference_donor:
                continue
            metrics.append(
                _run_pair(
                    bulk_donor=bulk_donor,
                    reference_donor=reference_donor,
                    donor_signatures=donor_signatures,
                    genes=common_genes,
                    cell_types=cell_types,
                    rng=rng,
                    output_dir=output_dir,
                )
            )

    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(output_dir / "summary_metrics.tsv", sep="\t", index=False)

    summary = {
        "h5ad_path": str(h5ad_path),
        "signature_path": str(signature_path),
        "output_dir": str(output_dir),
        "normalization": f"log1p(CP{int(TARGET_SUM)})",
        "expression_source": "decontXcounts",
        "eligible_donors": eligible_donors,
        "n_marker_genes": len(marker_genes),
        "n_common_genes": len(common_genes),
        "n_pseudobulk": N_PSEUDOBULK,
        "dirichlet_alpha": DIRICHLET_ALPHA,
        "random_seed": RANDOM_SEED,
        "cell_type_map": CELL_TYPE_MAP,
        "method": "bulk donor pseudobulk generated from true donor-specific marker reference; anchor module uses a different donor reference to infer proportions and recover bulk-domain reference",
    }
    with open(output_dir / "experiment_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(metrics_df.to_string(index=False))


if __name__ == "__main__":
    main()
