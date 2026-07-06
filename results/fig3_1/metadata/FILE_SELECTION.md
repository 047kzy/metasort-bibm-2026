# File selection notes

## Included

- Version 2 result tables, metrics, run logs, configuration files, and plots.
- Per-tissue Version 2 MetaSort prediction tables, gene weights, metrics, and
  weight summaries.
- Lightweight processed inputs needed to interpret the result matrices:
  `bulk.txt`, `bulkRatio.txt`, `cellTypeExpr.txt`, `signature.txt`,
  `markers.txt`, `sigGenes.txt`, and `nCellTypes.txt`.
- Hierarchy JSON files.
- Selected supplementary summary tables and figure files that support the
  Figure 3A-B hierarchy analysis.
- Scripts needed to generate the retained Version 2 result summaries and plots.

## Excluded from this cleaned release

- `__pycache__/` and compiled Python cache files.
- Historical or superseded first30/original-weight runs that are not part of the
  final Version 2 result set.
- Large single-cell expression matrices (`singleCellExpr.txt`) and associated
  labels/subjects. These should be deposited separately if the publication
  requires raw or source-level processed matrices.
- Large all-gene intermediate caches from weight/DE exploratory analyses, such
  as long all-gene compressed tables. Retained files are summary tables and
  figure-source summaries.
- Per-level intermediate prediction/reference matrices from hierarchy reference
  analyses. Retained files are the consolidated long tables and summaries.

## Public deposition notes

- Add a final open data licence before deposit.
- Add repository DOI/accession and preferred dataset citation after upload.
- If the paper claims full reproducibility from raw data, deposit or reference
  the excluded single-cell matrices and any external source data.
