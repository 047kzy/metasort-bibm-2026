# MetaSort_jk low-abundance cell-type comparison

Rare threshold: `0 < truth < 0.05`; minimum rare observations per tissue/cell type: 5.
Every winner comparison uses the intersection of mixtures available for all seven methods within that tissue.

- Eligible tissue/cell-type pairs under the rare-observation definition: 133
- MetaSort_jk winners under the rare-observation definition: 16
- Cell types with mean truth below threshold: 77
- MetaSort_jk winners among low-mean-proportion cell types: 7

## Low-mean-proportion cell types won by MetaSort_jk

Here a cell type is rare when its mean true proportion across common mixtures is below the threshold; MAE uses all common mixtures.

| Dataset / tissue | Cell type | n | Mean truth | MetaSort MAE | Runner-up (MAE) | Improvement |
|---|---|---:|---:|---:|---|---:|
| Tabula-Sapiens / Pancreas | B_cell | 100 | 4.53% | 0.01355 | bayesprism (0.03614) | 62.5% |
| Tabula-Sapiens / Salivary_Gland | mature_NK_T_cell | 100 | 3.31% | 0.01509 | DWLS (0.02100) | 28.1% |
| Tabula-Sapiens / Skin | CD1c_positive_myeloid_dendritic_cell | 100 | 4.07% | 0.01310 | bayesprism (0.01813) | 27.7% |
| Tabula-Sapiens / Lymph_Node | intermediate_monocyte | 100 | 4.85% | 0.03633 | HiDecon (0.04019) | 9.6% |
| CELlxGENE / Primary_Somatosensory_Cortex | oligodendrocyte_precursor_cell | 100 | 3.37% | 0.01912 | bayesprism (0.02264) | 15.6% |
| CELlxGENE / Fimbria_Of_Uterine_Tube | ciliated_epithelial_cell | 100 | 3.88% | 0.02678 | DWLS (0.02977) | 10.0% |
| CELlxGENE / Basal_Zone_Of_Heart | vein_endothelial_cell | 100 | 3.04% | 0.00613 | DWLS (0.00752) | 18.5% |

## Rare-observation comparisons won by MetaSort_jk

Here MAE is restricted to mixtures where that cell type has a positive true proportion below the threshold.

| Dataset / tissue | Cell type | n | Mean truth | MetaSort MAE | Runner-up (MAE) | Improvement |
|---|---|---:|---:|---:|---|---:|
| Tabula-Sapiens / Pancreas | B_cell | 59 | 3.44% | 0.00892 | bayesprism (0.02998) | 70.2% |
| Tabula-Sapiens / Trachea | endothelial_cell | 13 | 4.10% | 0.00283 | DWLS (0.02226) | 87.3% |
| CELlxGENE / Small_Intestine | blood_vessel_endothelial_cell | 7 | 4.49% | 0.00600 | DWLS (0.02537) | 76.4% |
| Tabula-Sapiens / Lung | adventitial_cell | 17 | 4.23% | 0.01542 | MuSic (0.03345) | 53.9% |
| Tabula-Sapiens / Bladder | pericyte_cell | 14 | 4.31% | 0.00407 | DWLS (0.01747) | 76.7% |
| Tabula-Sapiens / Bladder | mature_NK_T_cell | 26 | 4.02% | 0.00375 | bayesprism (0.01417) | 73.5% |
| Tabula-Sapiens / Salivary_Gland | endothelial_cell | 30 | 3.79% | 0.01356 | bayesprism (0.02313) | 41.4% |
| Tabula-Sapiens / Salivary_Gland | mature_NK_T_cell | 76 | 3.31% | 0.01210 | DWLS (0.01998) | 39.4% |
| CELlxGENE / Primary_Visual_Cortex | lamp5_gabaergic_cortical_interneuron | 17 | 4.16% | 0.00271 | bayesprism (0.00690) | 60.6% |
| CELlxGENE / Primary_Somatosensory_Cortex | oligodendrocyte_precursor_cell | 74 | 3.26% | 0.01766 | bayesprism (0.02148) | 17.8% |
| CELlxGENE / Basal_Zone_Of_Heart | vein_endothelial_cell | 85 | 3.26% | 0.00655 | DWLS (0.00859) | 23.7% |
| CELlxGENE / Small_Intestine | neuron | 22 | 3.95% | 0.01526 | bayesprism (0.01728) | 11.7% |
| Tabula-Sapiens / Skin | CD1c_positive_myeloid_dendritic_cell | 66 | 3.51% | 0.01240 | bayesprism (0.01412) | 12.1% |
| Tabula-Sapiens / Lymph_Node | intermediate_monocyte | 37 | 3.75% | 0.02940 | bayesprism (0.03037) | 3.2% |
| CELlxGENE / Primary_Visual_Cortex | chandelier_pvalb_gabaergic_cortical_interneuron | 69 | 3.54% | 0.00726 | DWLS (0.00758) | 4.3% |
| CELlxGENE / Fimbria_Of_Uterine_Tube | endothelial_cell_of_lymphatic_vessel | 19 | 4.18% | 0.00657 | bayesprism (0.00659) | 0.4% |

## Overall MAE on rare observations

| Method | n | MAE | Median absolute error |
|---|---:|---:|---:|
| DWLS | 6542 | 0.01356 | 0.00825 |
| MuSic | 6542 | 0.02061 | 0.01438 |
| bayesprism | 6542 | 0.02184 | 0.01673 |
| LinDeconSeq | 6542 | 0.02384 | 0.01649 |
| MetaSort_jk | 6542 | 0.02511 | 0.01798 |
| NNLS | 6542 | 0.03789 | 0.02462 |
| HiDecon | 6542 | 0.05623 | 0.03820 |

The tissue/cell-type wins are local findings; MetaSort_jk is not the best method on rare observations overall.
Tabula-Sapiens/Thymus uses 10 common mixtures because its MuSic result contains only 10 rows; all other tissues use 100.

