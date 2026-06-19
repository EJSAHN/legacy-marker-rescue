# Source-data validation targets

`python scripts/validate_source_data.py --workbook source_data/Supplementary_Data_S1.xlsx` recomputes the following manuscript-level quantities from derived source-data tables.

| Area | Recomputed quantity | Expected value |
| --- | --- | --- |
| Mash calibration | RAPD/Mash pairwise comparisons | 253 |
| Mash calibration | RAPD-vs-Mash Pearson r | 0.756 |
| Mash calibration | RAPD-vs-Mash Spearman r | 0.636 |
| Mash calibration | Minimum leave-one-assembly-out Spearman r | 0.522 |
| Genome context | Observed gene-overlap fraction | 0.769 |
| Genome context | Matched-random gene-overlap mean | 0.658 |
| Genome context | Observed 10-kb contig-end fraction | 0.366 |
| Genome context | Observed mean GC fraction | 0.552 |
| Random-primer null | Accepted GC-matched random primer sets | 25 |
| Random-primer null | Historical predicted amplicons | 14,630 |
| Random-primer null | Historical non-duplicate Spearman r | 0.295 |
| Random-primer null | Historical fixed-core band fraction | 0.872 |
| Random-primer null | Historical polymorphic band fraction | 0.043 |
| Guthrie schematic rescue | Composite loci | 31 |
| Guthrie schematic rescue | Published-summary RMSE | 2.49 |
| Guthrie bridge | Legacy bins with any modern bridge | 27 of 31 |
| Guthrie bridge | Polymorphic bridge bins | 6 |
| Guthrie bridge | Core-only bridge bins | 21 |
| Guthrie bridge | No-modern-match bins | 4 |
| Denoyes-Rothan matrix audit | Total detected first-pass segments | 204 |
| Denoyes-Rothan matrix audit | OPC02 detected segments | 87 |
| Denoyes-Rothan matrix audit | OPA13 detected segments | 117 |
| Denoyes-Rothan matrix audit | Retained nonzero sample lanes | 33 |
| Denoyes-Rothan matrix audit | Manual-check sample lanes | 5 |

Values in the manuscript text are rounded for readability; the validator uses tighter numerical tolerances on the underlying source-data workbook.
