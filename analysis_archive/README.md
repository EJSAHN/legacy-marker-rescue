# Analysis provenance

This folder records the manuscript analysis stages at a public-facing level. It does not contain raw publisher PDFs, publisher-owned figure crops, manuscript display-item drawing code, environment test scripts, or exploratory run branches.

The public repository emphasizes two reproducibility layers:

1. reusable implementation code in `src/legacy_marker_rescue/`; and
2. source-data validation through `scripts/validate_source_data.py` and `data/source/Supplementary_Data_S1.xlsx`.

The manuscript was developed through multiple exploratory stages. The final reported numerical summaries are validated from the archived source-data workbook rather than from figure files. If exact run scripts are later deposited in a DOI archive, they should be placed in a separate `final_workflow_scripts/` directory after removing local paths, display-generation code, and exploratory troubleshooting branches.

## Final manuscript stages

| Manuscript stage | Public description | Public validation layer |
| --- | --- | --- |
| Genome inventory | Retain genomic FASTA entries and exclude CDS/non-genomic sequence files. | Workbook validation and genome-inventory metadata. |
| In-silico RAPD | Simulate one-mismatch 10-mer RAPD, retain inward-facing primer pairs, and bin predicted products. | RAPD/Mash distance and random-primer source tables. |
| Mash hardening | Run Mash with k = 21 and sketch size 10,000; join Mash distances to RAPD band-space distances. | `WGS distance pairs` and `WGS leave-one-out` workbook sheets. |
| Genome context | Compare predicted RAPD amplicons with size-matched random intervals from the same contigs. | `Genome context summary` workbook sheet. |
| Random-primer null | Simulate GC-matched random primer sets under the same RAPD condition and compare with the historical primer set. | `Random primer sets` workbook sheet. |
| Guthrie digitization | Convert schematic RAPD signatures into derived band-coordinate and composite matrix tables. | `Guthrie summary agreement` workbook sheet. |
| Guthrie bridge | Match digitized legacy band-size bins to modern in-silico RAPD band-size bins, including modern-polymorphic windows. | `Guthrie bridge tolerance` and `Guthrie bridge classes` workbook sheets. |
| Denoyes matrix audit | Convert representative OPC02/OPA13 gel panels into derived first-pass matrices and lane-level audit tables. | `Denoyes detected segments` and `Denoyes lane audit` workbook sheets. |

See `workflow_stage_manifest.tsv` for a compact provenance map.
