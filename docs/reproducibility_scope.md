# Reproducibility scope

This repository is organized around two levels of reproducibility.

## Level 1: public source-data validation

The source-data workbook in `data/source/Supplementary_Data_S1.xlsx` contains the derived tables that underlie the manuscript display items. The script `scripts/validate_source_data.py` recomputes headline results directly from those tables and raises an error if the expected values are not recovered.

This validation layer is designed to run on any machine with the repository environment installed. It does not require local copies of genome FASTA files, GFF annotation files, raw publisher PDFs, or figure crops.

## Level 2: analysis provenance

The full analysis history used local public-genome downloads, local annotation files, and source PDFs. Some of those source materials are not redistributed here because of file size, source licensing, or publisher copyright. The `analysis_archive/` folder documents the final workflow stages and points to the derived outputs that are validated in Level 1.

Archived workflow scripts, if released in a DOI archive, should be treated as provenance records. They may require local FASTA, GFF, Mash output, or source-PDF inputs. They are not the primary public numerical validation layer; that role is served by `scripts/validate_source_data.py`.

## What is intentionally excluded

The public repository excludes raw publisher PDFs, publisher-owned figure crops, manuscript display-item drawing code, exploratory troubleshooting scripts, environment test scripts, and exploratory intermediate runs. This keeps the repository focused on implementation, derived data, and auditable numerical validation.
