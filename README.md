# legacy-marker-rescue

Reference implementation and source-data validation for recovering legacy RAPD and related random-polymorphism banding figures as auditable band-size fingerprint archives.

This repository accompanies a manuscript on genome-calibrated rescue of published anonymous-marker figures. The workflow is intentionally conservative: historical figures are treated as recoverable fingerprint records, not as error-free genotype tables, locus-resolved whole-genome substitutes, or evidence that same-size bands correspond to the same genomic locus.

## Repository organization

The public release separates three roles.

1. **Reusable code.** `src/legacy_marker_rescue/` contains small utilities for sequence handling, in-silico RAPD simulation, distance comparison, digitized segment tables, random-primer generation, and band-size bridging.
2. **Validation scripts.** `scripts/validate_source_data.py` recomputes headline manuscript values from the derived source-data workbook rather than copying summary cells. The `scripts/examples/` directory contains small illustrative scripts for local inputs.
3. **Derived source data and provenance.** `source_data/Supplementary_Data_S1.xlsx` contains the tables underlying the manuscript figures and tables. `docs/provenance/` maps manuscript analysis stages to the workbook sheets used for public validation.

The repository is **not** presented as a one-command reconstruction of every analysis from raw publisher PDFs. Raw publisher PDFs, figure crops, figure-generation code, and exploratory troubleshooting scripts are not included.

## Quick validation

Create the environment and run the source-data validator:

```bash
conda env create -f environment.yml
conda activate legacy-marker-rescue
python -m pip install -e .
python scripts/validate_source_data.py --workbook source_data/Supplementary_Data_S1.xlsx
```

The validator recomputes values from workbook rows and columns. Expected headline checks include 231 RAPD/Mash pairwise comparisons, RAPD-vs-Mash Pearson *r* approximately 0.76, Spearman *r* approximately 0.59, 31 Guthrie composite loci, 27/31 legacy bands with a modern band-size bridge, 6/31 with a modern-polymorphic bridge, and 204 Denoyes-Rothan first-pass detected segments.

## Example implementation workflow

The example scripts illustrate reusable workflow components on user-supplied inputs. Edit `configs/guthrie_1992_cgraminicola.yaml` so that genome inventories, Mash distances, and source files point to local files, then run selected steps:

```bash
python scripts/01_prepare_genome_inventory.py --config configs/guthrie_1992_cgraminicola.yaml
python scripts/02_simulate_rapd.py --config configs/guthrie_1992_cgraminicola.yaml
python scripts/03_compare_bandspace_to_mash.py --config configs/guthrie_1992_cgraminicola.yaml
python scripts/06_digitize_legacy_figure.py --config configs/guthrie_1992_cgraminicola.yaml
```

Additional examples are in `scripts/examples/`. They are building blocks for local analyses and should not be confused with the manuscript source-data validation script.

Mash can be run outside Python using a FASTA list:

```bash
bash scripts/run_mash.sh genome_list.txt results/genomes_mash
```

## Data policy

Do not commit raw publisher PDFs, raw publisher figure crops, or images extracted from copyrighted sources. Commit only bibliographic metadata, crop templates, derived coordinate tables, derived band matrices, QC summaries, source-data workbooks, and scripts needed to reproduce or validate derived outputs. Public-domain material should still be handled with customary source credit.

The copy of `Supplementary_Data_S1.xlsx` in `source_data/` is included so that the validator runs from a fresh clone. After journal publication, the journal-hosted supplementary workbook should be treated as the canonical copy.

## Citation and archive DOI

Temporary citation fields are provided in `CITATION.cff`. An archival DOI will be added after a versioned software release is deposited.
