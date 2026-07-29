# Workflow provenance

This directory is a provenance ledger for the manuscript workflow. It is not a raw-data rerun archive and does not contain publisher PDFs, figure crops, display-item drawing code, or exploratory troubleshooting scripts.

The public repository is organized around two reproducibility layers:

1. reusable workflow components in `src/legacy_marker_rescue/`; and
2. source-data validation through `scripts/validate_source_data.py` and `source_data/Supplementary_Data_S1.xlsx`.

The table `workflow_stage_manifest.tsv` maps each manuscript analysis stage to the derived workbook sheets used for public validation and notes which local inputs would be required for a full rerun. Exact exploratory run histories are intentionally not included in the GitHub companion repository.
