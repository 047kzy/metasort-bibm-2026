# Figure 4: Donor-Mismatched Reference Recovery

This directory contains the `bulk_TSP14__reference_TSP10` experiment used for
Figure 4.

## Experiment

Pseudo-bulk mixtures were generated from donor TSP14. The initial deconvolution
reference was derived from donor TSP10. The retained files compare the recovered
bulk-domain reference against the true TSP14 reference and the original TSP10
reference.

## Key Files

- `input_reference_donor_reference.tsv`: initial donor-mismatched TSP10
  reference.
- `true_bulk_donor_reference.tsv`: true TSP14 bulk-domain reference used for
  evaluation.
- `recovered_bulk_domain_reference.tsv`: recovered reference profile.
- `gene_recovery_metrics.tsv`: gene-level recovery metrics.
- `gene_metrics.tsv`: retained gene-level diagnostic metrics.
- `recovered_reference_nnls_iteration_metrics.tsv`: iteration-level reference
  recovery metrics.
- `initial_metasort_proportions.tsv`: initial MetaSort proportion estimates.
- `anchor_genes.txt`: anchor genes used for recovery.
- `metrics.json`: compact experiment summary.
