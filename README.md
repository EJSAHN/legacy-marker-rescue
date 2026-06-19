# legacy-marker-rescue

Reference implementation, source-data validation, and analysis provenance for recovering legacy random-polymorphism marker figures as auditable band-size fingerprint archives.

This repository accompanies a manuscript on rescuing published RAPD and related anonymous-marker figures. It is intentionally conservative: historical figures are treated as recoverable fingerprint archives, not as error-free genotype tables, locus-resolved whole-genome substitutes, or evidence that same-size bands correspond to the same genomic locus.

## What this repository provides

The public release separates three roles that are often mixed together.

1. **Reusable implementation.** `src/legacy_marker_rescue/` contains small, inspectable utilities for sequence handling, in-silico RAPD simulation, distance comparison, digitized segment tables, random-primer generation, and band-size bridging.
2. **Source-data validation.** `scripts/validate_source_data.py` recomputes headline manuscript results from the multisheet workbook in `data/source/Supplementary_Data_S1.xlsx`. This is the main public numerical check for the display items.
3. **Analysis provenance.** `analysis_archive/` records the final manuscript workflow stages and their derived-data targets. Raw publisher PDFs, publisher-owned figure crops, display-item drawing code, and exploratory troubleshooting scripts are not included.

The repository is **not** presented as a one-command reconstruction of every analysis from raw publisher PDFs. Instead, it provides clean implementation code, derived source data, a validation script, configuration templates, and a clear record of the manuscript workflow.

## Quick validation

Create the environment and run the source-data validator:

```bash
conda env create -f environment.yml
conda activate legacy-marker-rescue
python -m pip install -e .
python scripts/validate_source_data.py --workbook data/source/Supplementary_Data_S1.xlsx
```

The validator recomputes values from workbook rows and columns rather than copying summary cells. Expected headline checks include 253 RAPD/Mash pairwise comparisons, RAPD-vs-Mash Pearson r ≈ 0.76, Spearman r ≈ 0.64, 31 Guthrie composite loci, 27/31 legacy bands with a modern band-size bridge, 6/31 with a modern-polymorphic bridge, and 204 Denoyes-Rothan first-pass detected segments.

## Example implementation workflow

The example scripts illustrate reusable workflow components on user-supplied inputs. Edit `configs/guthrie_1992_cgraminicola.yaml` so that genome inventories, Mash distances, and source files point to local files, then run selected steps:

```bash
python scripts/01_prepare_genome_inventory.py --config configs/guthrie_1992_cgraminicola.yaml
python scripts/02_simulate_rapd.py --config configs/guthrie_1992_cgraminicola.yaml
python scripts/03_compare_bandspace_to_mash.py --config configs/guthrie_1992_cgraminicola.yaml
python scripts/06_digitize_legacy_figure.py --config configs/guthrie_1992_cgraminicola.yaml
```

Additional examples are in `scripts/examples/`. They are provided as building blocks and should not be confused with the manuscript source-data validation script.

Mash can be run outside Python using a FASTA list:

```bash
bash scripts/run_mash.sh genome_list.txt results/genomes_mash
```

## Data policy

Do not commit raw publisher PDFs, raw publisher figure crops, or images extracted from copyrighted sources. Commit only bibliographic metadata, crop templates, derived coordinate tables, derived band matrices, QC summaries, source-data workbooks, and scripts needed to reproduce or validate derived outputs. Public-domain material should still be handled with customary source credit.

## Citation and archive DOI

Temporary citation fields are provided in `CITATION.cff`. The software DOI and final repository URL should be filled after the archival release is deposited, for example through Zenodo. Until then, cite the associated manuscript and the marker-system and computational references listed in `docs/references.md`.
